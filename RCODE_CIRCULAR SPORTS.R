# --------------------------------------------------------------------------------
# Circular Statistical Framework for a Simulated Network in R
# Based on: "A Novel Circular Statistical Framework for Adaptive Network Traffic Routing"
# --------------------------------------------------------------------------------

library(igraph)
library(circular)
library(ggplot2)
library(reshape2)

# 1. PRELIMINARIES & FOLDER
dir.create("results_simulated_network", showWarnings = FALSE)

# 2. SYNTHETIC NETWORK GENERATION
N <- 15        # number of nodes
p_edge <- 0.3  # probability of edge creation

set.seed(2023)
G <- sample_gnp(n = N, p = p_edge, directed = FALSE)
while (!is_connected(G)) {
  G <- sample_gnp(n = N, p = p_edge)
}

nodes <- 0:(N - 1)
A <- as.matrix(get.adjacency(G))

cat(sprintf("Generated a random connected graph with %d nodes and %d edges.\n", N, ecount(G)))

# 3. SIMULATE DATA (COVARIATE X_i AND PHASE theta_i)
logistic <- function(z) {
  return(1.0 / (1.0 + exp(-z)))
}

beta_true <- 1.7
kappa_true <- 4.0

set.seed(2023)
X <- rnorm(N, mean = 0.0, sd = 1.0)

eps <- as.numeric(rvonmises(n = N, mu = circular(0), kappa = kappa_true))
mu_true <- 2.0 * pi * logistic(X * beta_true)
theta_sim <- (mu_true + eps) %% (2.0 * pi)

df_nodes <- data.frame(
  node = nodes,
  X = X,
  theta_sim = theta_sim
)

cat("\nHead of the Node DataFrame:\n")
print(head(df_nodes))

# 4. FIT THE CIRCULAR REGRESSION (MLE)
neg_log_lik_circ_reg <- function(params, x, theta) {
  beta_ <- params[1]
  log_kappa <- params[2]
  kappa_ <- exp(log_kappa)
  mu_ <- 2.0 * pi * logistic(x * beta_)
  n <- length(theta)
  cos_term <- cos(theta - mu_)
  nll <- -kappa_ * sum(cos_term) + n * log(besselI(kappa_, 0) + 1e-12)
  return(nll)
}

init_params <- c(0.0, 0.0) # [beta, log(kappa)]
res <- optim(
  par = init_params,
  fn = neg_log_lik_circ_reg,
  x = df_nodes$X,
  theta = df_nodes$theta_sim,
  method = "BFGS"
)

beta_hat <- res$par[1]
kappa_hat <- exp(res$par[2])

cat("\n==========================\n")
cat("MLE RESULTS (CIRCULAR REGRESSION)\n")
cat("==========================\n")
cat(sprintf("Converged: %s, Value: %.4f\n", res$convergence == 0, res$value))
cat(sprintf("True beta:  %.4f,  Estimated beta:  %.4f\n", beta_true, beta_hat))
cat(sprintf("True kappa: %.4f, Estimated kappa: %.4f\n", kappa_true, kappa_hat))
cat("==========================\n\n")

df_mle <- data.frame(
  Parameter = c("beta", "kappa"),
  True_Value = c(beta_true, kappa_true),
  Estimate = c(beta_hat, kappa_hat)
)
cat("MLE Table:\n")
print(df_mle)
write.csv(df_mle, "results_simulated_network/mle_estimates.csv", row.names = FALSE)

# 5. COMPUTE CIRCULAR CENTRALITY
circular_centrality <- function(g, theta_dict) {
  cent <- numeric(vcount(g))
  v_list <- V(g)
  for (i in seq_along(v_list)) {
    v <- v_list[i]
    neighbors_v <- neighbors(g, v, mode = "all")
    if (length(neighbors_v) == 0) {
      cent[i] <- 0.0
      next
    }
    s <- 0.0
    for (u in neighbors_v) {
      # Match by index (1-based in R)
      u_idx <- which(v_list == u)
      s <- s + (1.0 - cos(theta_dict[i] - theta_dict[u_idx]))
    }
    cent[i] <- s / length(neighbors_v)
  }
  return(cent)
}

cent_vals <- circular_centrality(G, df_nodes$theta_sim)
df_nodes$circular_centrality <- cent_vals

cat("\nCircular Centrality (Head):\n")
print(head(df_nodes[, c("node", "theta_sim", "circular_centrality")]))
write.csv(df_nodes, "results_simulated_network/node_data_with_cent.csv", row.names = FALSE)

# 6. ADAPTIVE ROUTING DEMONSTRATION (TOY EXAMPLE)
T_period <- 24
arrival_times <- c(2, 8, 12, 15, 23, 29)

compute_phase_of_time <- function(t) {
  return(2.0 * pi * ((t %% T_period) / T_period))
}

predicted_mu <- 2.0 * pi * logistic(df_nodes$X * beta_hat)
df_nodes$mu_hat <- predicted_mu

route_decision <- function(t, df) {
  theta_t <- compute_phase_of_time(t)
  df$abs_diff <- abs(df$mu_hat - theta_t)
  best_idx <- which.min(df$abs_diff)
  return(list(chosen_node = df$node[best_idx], abs_diff = df$abs_diff[best_idx]))
}

routing_results <- list()
for (t in arrival_times) {
  res_route <- route_decision(t, df_nodes)
  routing_results[[length(routing_results) + 1]] <- data.frame(
    arrival_time = t,
    chosen_node = res_route$chosen_node,
    abs_phase_diff = res_route$abs_diff
  )
}

df_routing <- do.call(rbind, routing_results)
cat("\nToy Routing Decisions:\n")
print(df_routing)
write.csv(df_routing, "results_simulated_network/toy_routing_results.csv", row.names = FALSE)

# 7. PLOTTING RESULTS
# 1) Observed vs. Predicted
p1 <- ggplot(df_nodes, aes(x = theta_sim, y = mu_hat)) +
  geom_point(alpha = 0.6, color = "navy", size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Observed Phase (theta_sim)", y = "Predicted Phase (mu_hat)",
       title = "Observed vs. Predicted Node Phases") +
  theme_minimal()
ggsave("results_simulated_network/observed_vs_predicted.png", plot = p1, width = 6, height = 5, dpi = 300)

# 2) Heatmap of adjacency matrix
melted_adj <- melt(A)
colnames(melted_adj) <- c("Var1", "Var2", "value")
p2 <- ggplot(melted_adj, aes(x = factor(Var1 - 1), y = factor(Var2 - 1), fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "Node", y = "Node", title = "Adjacency Matrix of Synthetic Network") +
  theme_minimal()
ggsave("results_simulated_network/adjacency_matrix_heatmap.png", plot = p2, width = 6, height = 5, dpi = 300)

# 3) Bar plot of circular centrality
df_cent_sorted <- df_nodes[order(df_nodes$circular_centrality), ]
df_cent_sorted$node_f <- factor(df_cent_sorted$node, levels = df_cent_sorted$node)
p3 <- ggplot(df_cent_sorted, aes(x = node_f, y = circular_centrality)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(x = "Node", y = "C(v)", title = "Circular Centrality per Node (Synthetic Network)") +
  theme_minimal()
ggsave("results_simulated_network/circular_centrality_barplot.png", plot = p3, width = 8, height = 5, dpi = 300)

# 4) Toy routing demo plot
p4 <- ggplot(df_routing, aes(x = arrival_time, y = chosen_node)) +
  geom_point(color = "red", size = 4, alpha = 0.7) +
  geom_text(aes(label = sprintf("diff=%.2f", abs_phase_diff)), hjust = -0.1, vjust = -0.5, color = "blue", size = 3.5) +
  labs(x = "Arrival Time (Hours)", y = "Chosen Node", title = "Toy Adaptive Routing Decisions") +
  theme_minimal() +
  xlim(min(df_routing$arrival_time) - 1, max(df_routing$arrival_time) + 5)
ggsave("results_simulated_network/adaptive_routing_demo.png", plot = p4, width = 7, height = 4, dpi = 300)

# 8. PRINTING WRAP-UP
cat("\nAll figures saved in './results_simulated_network' folder:\n")
for (figname in c("observed_vs_predicted.png", "adjacency_matrix_heatmap.png", 
                  "circular_centrality_barplot.png", "adaptive_routing_demo.png")) {
  cat(sprintf("   --> %s\n", figname))
}

cat("\nTables saved in './results_simulated_network' folder:\n")
for (tabname in c("mle_estimates.csv", "node_data_with_cent.csv", "toy_routing_results.csv")) {
  cat(sprintf("   --> %s\n", tabname))
}

cat("\nScript completed successfully. The above approach is aligned with:\n")
cat("'A Novel Circular Statistical Framework for Adaptive Network Traffic Routing.'\n")

