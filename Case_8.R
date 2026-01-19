# ============================================================
# Case 8 (WORKING for L=3): Multiplex binary SBM (layer-wise collapse)
# Collapsed Gibbs + images + tables + summary
# Output: case8_outputs/
# ============================================================

# -----------------------------
# 0) Output directory + seed
# -----------------------------
out_dir <- "case8_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
set.seed(9)

# -----------------------------
# 1) Utilities (base R)
# -----------------------------
choose2 <- function(x) x * (x - 1) / 2

ari <- function(z_true, z_hat) {
  tab <- table(z_true, z_hat)
  n <- sum(tab)
  sum_choose2 <- function(x) sum(choose2(x))
  a <- sum_choose2(tab)
  b <- sum_choose2(rowSums(tab))
  c <- sum_choose2(colSums(tab))
  d <- choose2(n)
  expected <- b * c / d
  maxval <- (b + c) / 2
  if (maxval == expected) return(0)
  (a - expected) / (maxval - expected)
}

log_dirichlet_multinomial_prior <- function(nk, alpha) {
  K <- length(nk)
  n <- sum(nk)
  sum(lgamma(alpha / K + nk) - lgamma(alpha / K)) + lgamma(alpha) - lgamma(alpha + n)
}

# -----------------------------
# 2) Simulation: multiplex binary SBM
# -----------------------------
sim_sbm_binary <- function(z, P) {
  n <- length(z)
  A <- matrix(0L, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    p <- P[z[i], z[j]]
    aij <- rbinom(1, 1, p)
    A[i, j] <- aij
    A[j, i] <- aij
  }
  diag(A) <- 0L
  A
}

sim_multiplex_binary <- function(z, P_list) {
  lapply(P_list, function(P) sim_sbm_binary(z, P))
}

# -----------------------------
# 3) Diagnostics helpers
# -----------------------------
binary_block_densities <- function(A, z) {
  K <- max(z)
  within_edges <- 0
  within_dyads <- 0
  for (k in 1:K) {
    idx <- which(z == k)
    if (length(idx) >= 2) {
      within_edges <- within_edges + sum(A[idx, idx]) / 2
      within_dyads <- within_dyads + choose2(length(idx))
    }
  }
  
  between_edges <- 0
  between_dyads <- 0
  for (k in 1:(K - 1)) for (l in (k + 1):K) {
    ik <- which(z == k); il <- which(z == l)
    if (length(ik) && length(il)) {
      between_edges <- between_edges + sum(A[ik, il])
      between_dyads <- between_dyads + length(ik) * length(il)
    }
  }
  
  c(
    within_edges = within_edges,
    within_dyads = within_dyads,
    within_density = ifelse(within_dyads > 0, within_edges / within_dyads, NA),
    between_edges = between_edges,
    between_dyads = between_dyads,
    between_density = ifelse(between_dyads > 0, between_edges / between_dyads, NA)
  )
}

write_booktabs_table <- function(df, file, caption, label) {
  stopifnot(is.data.frame(df))
  cols <- ncol(df)
  align <- paste0(rep("l", cols), collapse = "")
  header <- paste(colnames(df), collapse = " & ")
  
  fmt_row <- function(row) paste(row, collapse = " & ")
  
  lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    sprintf("\\begin{tabular}{%s}", align),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule"
  )
  for (i in 1:nrow(df)) lines <- c(lines, paste0(fmt_row(df[i, , drop = TRUE]), " \\\\"))
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\end{table}"
  )
  writeLines(lines, con = file)
}

# -----------------------------
# 4) Spectral initialization on aggregated multiplex graph
# -----------------------------
spectral_init <- function(A_list, K) {
  A_sum <- Reduce("+", A_list)             # counts in {0,1,2,3}
  A_sum <- (A_sum + t(A_sum)) / 2
  diag(A_sum) <- 0
  
  # Top-K eigenvectors
  eig <- eigen(A_sum, symmetric = TRUE)
  X <- eig$vectors[, 1:K, drop = FALSE]
  
  # kmeans on embedding
  km <- kmeans(X, centers = K, nstart = 50)
  as.integer(km$cluster)
}

# -----------------------------
# 5) Collapsed Gibbs (multiplex) with trace + co-clustering
# -----------------------------
gibbs_multiplex_binary_trace <- function(A_list, K,
                                         a = 1, b = 1,
                                         alpha_part = 10,          # <- key change: discourage empty clusters
                                         n_iter = 6000, burn = 1500, thin = 10,
                                         z_init = NULL,
                                         z_true = NULL,
                                         store_cocluster = TRUE,
                                         verbose = TRUE) {
  L <- length(A_list)
  n <- nrow(A_list[[1]])
  stopifnot(all(sapply(A_list, function(A) is.matrix(A) && nrow(A) == n && ncol(A) == n)))
  
  # init z
  if (is.null(z_init)) z <- sample.int(K, n, replace = TRUE) else z <- as.integer(z_init)
  nk <- tabulate(z, K)
  
  # s1_list[[ell]]: KxK upper-tri counts of edges per block for each layer
  s1_list <- vector("list", L)
  for (ell in 1:L) {
    A <- A_list[[ell]]
    diag(A) <- 0
    s1 <- matrix(0, K, K)
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      r <- z[i]; s <- z[j]
      aa <- min(r, s); bb <- max(r, s)
      s1[aa, bb] <- s1[aa, bb] + A[i, j]
    }
    s1_list[[ell]] <- s1
  }
  
  n_dyad <- function(r, s, nk) if (r == s) choose2(nk[r]) else nk[r] * nk[s]
  log_m_bb <- function(m, n) lbeta(a + m, b + n - m) - lbeta(a, b)
  
  loglik_layers <- function(s1_list, nk) {
    ll <- numeric(L)
    for (ell in 1:L) {
      s1 <- s1_list[[ell]]
      tmp <- 0
      for (r in 1:K) for (s in r:K) tmp <- tmp + log_m_bb(s1[r, s], n_dyad(r, s, nk))
      ll[ell] <- tmp
    }
    ll
  }
  
  # storage
  trace <- data.frame(iter = integer(0), logprior = numeric(0),
                      loglik = numeric(0), logpost = numeric(0), ARI = numeric(0))
  
  cocluster <- if (store_cocluster) matrix(0, n, n) else NULL
  n_coclust <- 0
  
  best_lp <- -Inf
  z_map <- z
  
  for (it in 1:n_iter) {
    for (i in 1:n) {
      old <- z[i]
      
      # edge counts from node i to each community, for each layer
      e_list <- vector("list", L)
      for (ell in 1:L) {
        A <- A_list[[ell]]
        e <- numeric(K)
        for (l in 1:K) {
          idx <- which(z == l)
          if (l == old) idx <- idx[idx != i]
          if (length(idx)) e[l] <- sum(A[i, idx])
        }
        e_list[[ell]] <- e
      }
      
      # remove i: update sufficient statistics
      for (ell in 1:L) {
        e <- e_list[[ell]]
        s1 <- s1_list[[ell]]
        for (l in 1:K) {
          if (e[l] == 0) next
          aa <- min(old, l); bb <- max(old, l)
          if (old == l) s1[old, old] <- s1[old, old] - e[l]
          else s1[aa, bb] <- s1[aa, bb] - e[l]
        }
        s1_list[[ell]] <- s1
      }
      nk[old] <- nk[old] - 1
      z[i] <- 0
      
      # candidate log weights
      logw <- rep(-Inf, K)
      for (k in 1:K) {
        logprior_k <- log(nk[k] + alpha_part / K)
        delta <- 0
        for (l in 1:K) {
          for (ell in 1:L) {
            s1 <- s1_list[[ell]]
            e <- e_list[[ell]]
            if (l == k) {
              m_old <- s1[k, k]
              n_old <- n_dyad(k, k, nk)
              m_new <- m_old + e[l]
              n_new <- choose2(nk[k] + 1)
              delta <- delta + (log_m_bb(m_new, n_new) - log_m_bb(m_old, n_old))
            } else {
              aa <- min(k, l); bb <- max(k, l)
              m_old <- s1[aa, bb]
              n_old <- n_dyad(aa, bb, nk)
              m_new <- m_old + e[l]
              n_new <- (nk[k] + 1) * nk[l]
              delta <- delta + (log_m_bb(m_new, n_new) - log_m_bb(m_old, n_old))
            }
          }
        }
        logw[k] <- logprior_k + delta
      }
      
      # sample new label
      w <- exp(logw - max(logw))
      new <- sample.int(K, 1, prob = w)
      
      # add i to new: update sufficient statistics
      for (ell in 1:L) {
        e <- e_list[[ell]]
        s1 <- s1_list[[ell]]
        for (l in 1:K) {
          if (e[l] == 0) next
          aa <- min(new, l); bb <- max(new, l)
          if (new == l) s1[new, new] <- s1[new, new] + e[l]
          else s1[aa, bb] <- s1[aa, bb] + e[l]
        }
        s1_list[[ell]] <- s1
      }
      nk[new] <- nk[new] + 1
      z[i] <- new
    }
    
    # store
    if (it >= burn && ((it - burn) %% thin == 0)) {
      ll_layers <- loglik_layers(s1_list, nk)
      ll_total <- sum(ll_layers)
      lp <- log_dirichlet_multinomial_prior(nk, alpha_part)
      lpost <- lp + ll_total
      ari_val <- if (!is.null(z_true)) ari(z_true, z) else NA_real_
      
      trace <- rbind(trace, data.frame(
        iter = it, logprior = lp, loglik = ll_total, logpost = lpost, ARI = ari_val
      ))
      
      if (lpost > best_lp) { best_lp <- lpost; z_map <- z }
      
      if (store_cocluster) {
        cocluster <- cocluster + outer(z, z, "==")
        n_coclust <- n_coclust + 1
      }
    }
    
    if (verbose && (it %% 200 == 0)) message(sprintf("iter %d / %d", it, n_iter))
  }
  
  if (store_cocluster && n_coclust > 0) cocluster <- cocluster / n_coclust
  
  list(z_map = z_map, logpost_map = best_lp, trace = trace, cocluster = cocluster)
}

# -----------------------------
# 6) Plot helpers
# -----------------------------
plot_adj_heatmap <- function(A, ord, main = "") {
  M <- A[ord, ord]
  n <- nrow(M)
  image(1:n, 1:n, t(M)[, n:1],
        col = gray.colors(2, start = 1, end = 0),
        axes = FALSE, xlab = "", ylab = "", main = main)
  box()
}

plot_prob_heatmap <- function(P, ord, main = "") {
  M <- P[ord, ord]
  n <- nrow(M)
  image(1:n, 1:n, t(M)[, n:1],
        col = gray.colors(256, start = 1, end = 0),
        axes = FALSE, xlab = "", ylab = "", main = main)
  box()
}

# ============================================================
# 7) EXPERIMENT (L=3, K=3) -- FIXED so it recovers communities
# ============================================================
n <- 150
K_true <- 3
z_true <- rep(1:K_true, each = n / K_true)

# --- IMPORTANT: make the L=3 regime identifiable (p_in >> p_out) ---
# "Complementary layers": each layer emphasizes a different community a bit more
P_out <- 0.015

P1 <- matrix(P_out, K_true, K_true); diag(P1) <- c(0.130, 0.080, 0.080)
P2 <- matrix(P_out, K_true, K_true); diag(P2) <- c(0.080, 0.130, 0.080)
P3 <- matrix(P_out, K_true, K_true); diag(P3) <- c(0.080, 0.080, 0.130)

A_list <- sim_multiplex_binary(z_true, list(P1, P2, P3))
A_sum  <- Reduce("+", A_list)  # aggregated count in {0,1,2,3}

# Spectral init from multiplex aggregation
z0_multi  <- spectral_init(A_list, K_true)
z0_single <- spectral_init(list(A_list[[1]]), K_true)

# Hyperparameters
a <- 1; b <- 1
alpha_part <- 10
n_iter <- 6000
burn <- 1500
thin <- 10

# --- Single layer baseline ---
fit_single <- gibbs_multiplex_binary_trace(
  A_list = list(A_list[[1]]),
  K = K_true,
  a = a, b = b,
  alpha_part = alpha_part,
  n_iter = n_iter, burn = burn, thin = thin,
  z_init = z0_single,
  z_true = z_true,
  store_cocluster = TRUE,
  verbose = FALSE
)

# --- Multiplex inference (L=3) ---
fit_multi <- gibbs_multiplex_binary_trace(
  A_list = A_list,
  K = K_true,
  a = a, b = b,
  alpha_part = alpha_part,
  n_iter = n_iter, burn = burn, thin = thin,
  z_init = z0_multi,
  z_true = z_true,
  store_cocluster = TRUE,
  verbose = FALSE
)

z_hat_single <- fit_single$z_map
z_hat_multi  <- fit_multi$z_map

ari_single <- ari(z_true, z_hat_single)
ari_multi  <- ari(z_true, z_hat_multi)

cat("ARI single-layer:  ", round(ari_single, 3), "\n")
cat("ARI multiplex L=3: ", round(ari_multi, 3), "\n")
cat("MAP cluster sizes (multiplex): ", tabulate(z_hat_multi, K_true), "\n")

# -----------------------------
# 8) Save figures
# -----------------------------
ord_true <- order(z_true)

png(file.path(out_dir, "fig_case8_layers_trueorder.png"), width = 1400, height = 450)
par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))
plot_adj_heatmap(A_list[[1]], ord_true, main = "Layer 1 (ordered by true z*)")
plot_adj_heatmap(A_list[[2]], ord_true, main = "Layer 2 (ordered by true z*)")
plot_adj_heatmap(A_list[[3]], ord_true, main = "Layer 3 (ordered by true z*)")
dev.off()

png(file.path(out_dir, "fig_case8_trace.png"), width = 1400, height = 500)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(fit_single$trace$iter, fit_single$trace$logpost, type = "l",
     xlab = "Iteration (saved)", ylab = "Collapsed log-posterior",
     main = "Log-posterior trace")
lines(fit_multi$trace$iter, fit_multi$trace$logpost, type = "l", lty = 2)
legend("bottomright",
       legend = c("Single layer", "Multiplex (L=3)"),
       lty = c(1, 2), bty = "n")

plot(fit_single$trace$iter, fit_single$trace$ARI, type = "l",
     xlab = "Iteration (saved)", ylab = "ARI vs true z*",
     main = "Recovery over time",
     ylim = c(-0.1, 1.0))
lines(fit_multi$trace$iter, fit_multi$trace$ARI, type = "l", lty = 2)
legend("bottomright",
       legend = c("Single layer", "Multiplex (L=3)"),
       lty = c(1, 2), bty = "n")

dev.off()

ord_map <- order(z_hat_multi)
png(file.path(out_dir, "fig_case8_cocluster.png"), width = 900, height = 900)
par(mar = c(2, 2, 3, 1))
plot_prob_heatmap(fit_multi$cocluster, ord_map,
                  main = "Posterior co-clustering (multiplex, L=3)\nordered by MAP z")
dev.off()

# Optional igraph plot: show edges supported in >=2 layers (denoised view)
if (requireNamespace("igraph", quietly = TRUE)) {
  library(igraph)
  A_thr <- (A_sum >= 2) * 1
  diag(A_thr) <- 0
  g <- graph_from_adjacency_matrix(A_thr, mode = "undirected", diag = FALSE)
  
  pal <- c("#1b9e77", "#d95f02", "#7570b3")
  vcols <- pal[z_hat_multi]
  
  png(file.path(out_dir, "fig_case8_network.png"), width = 900, height = 900)
  plot(g,
       vertex.size = 6,
       vertex.label = NA,
       vertex.color = vcols,
       main = "Aggregated multiplex graph (edges in >=2 layers)\ncolored by MAP communities")
  dev.off()
}

# -----------------------------
# 9) Save summary tables + text summary
# -----------------------------
summary_comp <- data.frame(
  setting = c("Single layer (layer 1 only)", "Multiplex (L=3)"),
  n_layers = c(1, length(A_list)),
  ARI = c(round(ari_single, 3), round(ari_multi, 3)),
  logpost_MAP = c(round(fit_single$logpost_map, 1), round(fit_multi$logpost_map, 1))
)

write.csv(summary_comp, file = file.path(out_dir, "table_case8_summary.csv"),
          row.names = FALSE)

write_booktabs_table(
  df = summary_comp,
  file = file.path(out_dir, "table_case8_summary.tex"),
  caption = "Case 8 (multiplex, L=3): comparison of single-layer vs multiplex collapsed inference.",
  label = "tab:case8_multiplex_summary"
)

dens_rows <- list()
for (ell in 1:length(A_list)) {
  d_true <- binary_block_densities(A_list[[ell]], z_true)
  d_hat  <- binary_block_densities(A_list[[ell]], z_hat_multi)
  dens_rows[[ell]] <- data.frame(
    layer = ell,
    within_true = round(d_true["within_density"], 4),
    between_true = round(d_true["between_density"], 4),
    within_hat = round(d_hat["within_density"], 4),
    between_hat = round(d_hat["between_density"], 4)
  )
}
dens_df <- do.call(rbind, dens_rows)
write.csv(dens_df, file = file.path(out_dir, "table_case8_layer_densities.csv"), row.names = FALSE)

summary_txt <- c(
  "Case 8 (WORKING): Multiplex collapsed SBM (L=3)",
  "----------------------------------------------",
  sprintf("n = %d nodes, K_true = %d, L = %d layers", n, K_true, length(A_list)),
  "",
  "Key idea: conditional on z, layers are independent; collapsed log-evidence adds across layers.",
  "",
  sprintf("ARI single-layer:  %.3f", ari_single),
  sprintf("ARI multiplex L=3: %.3f", ari_multi),
  sprintf("MAP cluster sizes (multiplex): %s", paste(tabulate(z_hat_multi, K_true), collapse = ", ")),
  "",
  sprintf("Wrote outputs to: %s/", out_dir),
  "  - fig_case8_layers_trueorder.png",
  "  - fig_case8_trace.png",
  "  - fig_case8_cocluster.png",
  "  - fig_case8_network.png (optional; requires igraph)",
  "  - table_case8_summary.csv / table_case8_summary.tex",
  "  - table_case8_layer_densities.csv"
)

writeLines(summary_txt, con = file.path(out_dir, "summary_case8.txt"))
cat(paste(summary_txt, collapse = "\n"), "\n")
