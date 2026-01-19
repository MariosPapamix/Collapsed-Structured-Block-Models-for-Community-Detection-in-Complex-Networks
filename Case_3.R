# ============================================================
# Case 3: Gaussian weighted SBM (NIG collapse) — corrected + visuals
# ============================================================
# Run top-to-bottom. Outputs PNGs + CSVs in your working directory.

# ---- packages (auto-install) ----
need <- c("ggplot2", "dplyr", "tidyr")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

# ---- ARI helper (uses your ari() if you have it; otherwise mclust) ----
ARI <- function(a, b) {
  if (exists("ari", mode = "function")) return(ari(a, b))
  if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
  mclust::adjustedRandIndex(a, b)
}

# ---- vectorized simulator for Gaussian weighted SBM ----
sim_sbm_gaussian <- function(z, Mu, Sigma) {
  n <- length(z)
  Y <- matrix(0, n, n)
  idx <- which(upper.tri(Y), arr.ind = TRUE)
  a <- z[idx[, 1]]
  b <- z[idx[, 2]]
  mu <- Mu[cbind(a, b)]
  sd <- Sigma[cbind(a, b)]
  y <- rnorm(nrow(idx), mean = mu, sd = sd)
  Y[upper.tri(Y)] <- y
  Y <- Y + t(Y)
  diag(Y) <- 0
  Y
}

# ---- edge dataframe: within vs between for a given partition z ----
edge_df <- function(Y, z, label = "partition") {
  idx <- which(upper.tri(Y), arr.ind = TRUE)
  w <- Y[upper.tri(Y)]
  same <- z[idx[, 1]] == z[idx[, 2]]
  tibble(
    weight = w,
    type = if_else(same, "within", "between"),
    partition = label
  )
}

# ---- NIG posterior summary per block (conditioned on z) ----
# Prior: sigma^2 ~ Inv-Gamma(alpha0,beta0),  mu|sigma^2 ~ N(mu0, sigma^2/kappa0)
block_posterior_summary <- function(Y, z, K, mu0, kappa0, alpha0, beta0) {
  stopifnot(is.matrix(Y), nrow(Y) == ncol(Y), length(z) == nrow(Y))
  out <- vector("list", K * (K + 1) / 2)
  t <- 1
  
  for (a in 1:K) for (b in a:K) {
    Ia <- which(z == a); Ib <- which(z == b)
    y <- numeric(0)
    
    if (length(Ia) > 0 && length(Ib) > 0) {
      if (a == b) {
        if (length(Ia) >= 2) {
          sub <- Y[Ia, Ia, drop = FALSE]
          y <- sub[upper.tri(sub)]
        }
      } else {
        y <- as.vector(Y[Ia, Ib, drop = FALSE])
      }
    }
    
    m <- length(y)
    
    if (m == 0) {
      ybar <- NA_real_
      ysd  <- NA_real_
      kappa_n <- kappa0
      mu_n    <- mu0
      alpha_n <- alpha0
      beta_n  <- beta0
    } else {
      ybar <- mean(y)
      ysd  <- if (m > 1) sd(y) else NA_real_
      sse  <- sum((y - ybar)^2)
      
      kappa_n <- kappa0 + m
      mu_n    <- (kappa0 * mu0 + m * ybar) / kappa_n
      alpha_n <- alpha0 + m / 2
      beta_n  <- beta0 + 0.5 * sse + 0.5 * (kappa0 * m / kappa_n) * (ybar - mu0)^2
    }
    
    # Marginal posterior for mu_ab is Student-t with df = 2*alpha_n
    df <- 2 * alpha_n
    mu_scale <- sqrt(beta_n / (kappa_n * alpha_n))
    mu_ci <- mu_n + stats::qt(c(0.025, 0.975), df = df) * mu_scale
    
    # Posterior mean of sigma^2 exists for alpha_n > 1: E[sigma^2|y] = beta_n/(alpha_n-1)
    sigma2_mean <- if (alpha_n > 1) beta_n / (alpha_n - 1) else NA_real_
    sigma_mean  <- sqrt(sigma2_mean)
    
    out[[t]] <- tibble(
      a = a, b = b,
      m = m,
      y_mean = ybar,
      y_sd   = ysd,
      mu_post_mean = mu_n,
      mu_post_lo95 = mu_ci[1],
      mu_post_hi95 = mu_ci[2],
      sigma_post_mean = sigma_mean,
      kappa_n = kappa_n, alpha_n = alpha_n, beta_n = beta_n
    )
    t <- t + 1
  }
  
  bind_rows(out)
}

# ---- plotting helpers ----
plot_weight_heatmap <- function(Y, z, title) {
  n <- nrow(Y)
  ord <- order(z, seq_len(n))
  z_ord <- z[ord]
  Y_ord <- Y[ord, ord, drop = FALSE]
  
  df <- as.data.frame(as.table(Y_ord))
  names(df) <- c("i", "j", "w")
  df$i <- as.integer(df$i)
  df$j <- as.integer(df$j)
  
  b <- cumsum(as.integer(table(z_ord)))
  b <- b[-length(b)]
  
  p <- ggplot(df, aes(x = j, y = i, fill = w)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_gradient2(midpoint = 0) +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL, fill = "weight") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  if (length(b) > 0) {
    p <- p +
      geom_vline(xintercept = b + 0.5, linewidth = 0.25) +
      geom_hline(yintercept = b + 0.5, linewidth = 0.25)
  }
  
  p
}

plot_block_heatmap <- function(M, title) {
  K <- nrow(M)
  df <- expand_grid(a = 1:K, b = 1:K) %>%
    mutate(value = M[cbind(a, b)])
  
  ggplot(df, aes(x = b, y = a, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = title, x = "b", y = "a", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank())
}

draw_process_diagram <- function(file = "case3_process.png") {
  png(file, width = 1200, height = 500, res = 150)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  box <- function(x, y, w, h, txt) {
    rect(x - w/2, y - h/2, x + w/2, y + h/2, lwd = 2)
    text(x, y, txt, cex = 0.95)
  }
  arr <- function(x0, y0, x1, y1) arrows(x0, y0, x1, y1, lwd = 2, length = 0.10)
  
  box(0.12, 0.50, 0.22, 0.28, "Observed weights\nY = (y_ij)")
  box(0.36, 0.72, 0.24, 0.22, "Latent labels\nz_i in {1..K}")
  box(0.36, 0.28, 0.26, 0.22, "Block params\n(mu_ab, sigma^2_ab)")
  box(0.62, 0.50, 0.30, 0.28, "NIG prior\n+ conjugacy\n(integrate out params)")
  box(0.88, 0.50, 0.22, 0.28, "Collapsed Gibbs\nupdate z_i using\n(m_ab, sum, sumsq)")
  
  arr(0.23, 0.55, 0.26, 0.64)
  arr(0.23, 0.45, 0.26, 0.36)
  arr(0.48, 0.72, 0.52, 0.58)
  arr(0.49, 0.28, 0.52, 0.42)
  arr(0.77, 0.50, 0.80, 0.50)
  
  text(0.62, 0.10, "Outputs: MAP/consensus clustering + posterior block summaries", cex = 1.05)
  dev.off()
}

# ---- FIXED: sparse network plot (Fruchterman-Reingold needs positive weights) ----
plot_sparse_network <- function(Y, z, q = 0.99,
                                file = "case3_network_sparse.png",
                                seed = 1) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    message("Skipping sparse network plot (install.packages('igraph') to enable).")
    return(invisible(NULL))
  }
  set.seed(seed)
  
  n <- nrow(Y)
  idx <- which(upper.tri(Y), arr.ind = TRUE)
  w_signed_all <- Y[upper.tri(Y)]
  
  thr <- as.numeric(quantile(abs(w_signed_all), probs = q))
  keep <- which(abs(w_signed_all) >= thr)
  
  edges <- idx[keep, , drop = FALSE]
  w_signed <- w_signed_all[keep]
  w_abs <- abs(w_signed)
  
  # IMPORTANT:
  # - store positive weights in 'weight' (for layout)
  # - store signed weights separately for coloring
  eps <- 1e-12
  df_edges <- data.frame(
    from = edges[, 1],
    to   = edges[, 2],
    weight = w_abs + eps,
    signed_weight = w_signed
  )
  
  g <- igraph::graph_from_data_frame(
    df_edges,
    directed = FALSE,
    vertices = data.frame(name = 1:n)
  )
  
  K <- length(unique(z))
  pal <- if ("hcl.colors" %in% getNamespaceExports("grDevices")) {
    grDevices::hcl.colors(K, "Dark 3")
  } else {
    grDevices::rainbow(K)
  }
  
  igraph::V(g)$color <- pal[as.integer(factor(z))]
  igraph::V(g)$size  <- 4
  igraph::V(g)$label <- NA
  
  # edge widths use positive weights
  wpos <- igraph::E(g)$weight
  if (length(wpos) == 1) {
    wnorm <- 1
  } else {
    denom <- (max(wpos) - min(wpos))
    wnorm <- if (denom < 1e-12) rep(1, length(wpos)) else (wpos - min(wpos)) / denom
  }
  igraph::E(g)$width <- 0.5 + 2.5 * wnorm
  
  # color edges by sign (from separate attribute)
  wsgn <- igraph::E(g)$signed_weight
  igraph::E(g)$color <- ifelse(wsgn >= 0, "grey30", "grey70")
  
  # FIX: explicitly pass positive weights to FR layout
  lay <- igraph::layout_with_fr(g, weights = igraph::E(g)$weight)
  
  png(file, width = 900, height = 750, res = 150)
  plot(
    g,
    layout = lay,
    main = sprintf("Strongest |weights| (top %.1f%%), colored by inferred cluster", (1 - q) * 100)
  )
  dev.off()
  
  invisible(g)
}

# ---- OPTIONAL: if your sampler returns z samples, make “process” plots ----
extract_z_samples <- function(fit, n, K = NULL) {
  # Common names first
  for (nm in c("z_store", "z_samps", "z_samples", "z_trace", "z_chain", "Z_store")) {
    if (!is.null(fit[[nm]])) return(fit[[nm]])
  }
  
  # Safe fallback:
  # accept matrices that look like S x n or n x S (NOT n x n),
  # and (optionally) have integer-ish values in 1..K.
  mats <- Filter(function(x) is.matrix(x), fit)
  mats <- Filter(function(x) {
    dims <- dim(x)
    has_n <- (dims[2] == n && dims[1] != n) || (dims[1] == n && dims[2] != n)
    if (!has_n) return(FALSE)
    
    if (!is.null(K)) {
      vals <- as.numeric(x)
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(FALSE)
      intish <- mean(abs(vals - round(vals)) < 1e-8) > 0.95
      inrange <- (min(vals) >= 1) && (max(vals) <= K)
      return(intish && inrange)
    }
    TRUE
  }, mats)
  
  if (length(mats) > 0) return(mats[[1]])
  NULL
}

plot_ari_trace_and_coclustering <- function(z_samps, z_true, z_sort, file_trace, file_coclust) {
  n <- length(z_true)
  
  if (ncol(z_samps) != n && nrow(z_samps) == n) z_samps <- t(z_samps)
  if (ncol(z_samps) != n) return(invisible(NULL))
  
  S <- nrow(z_samps)
  
  # ARI trace
  ari_trace <- apply(z_samps, 1, function(z) ARI(z_true, z))
  df_trace <- tibble(iter = seq_len(S), ARI = ari_trace)
  
  p_trace <- ggplot(df_trace, aes(x = iter, y = ARI)) +
    geom_line() +
    labs(title = "ARI over saved Gibbs samples", x = "saved iteration", y = "ARI") +
    theme_minimal(base_size = 11)
  ggsave(file_trace, p_trace, width = 7, height = 3.5, dpi = 200)
  
  # Co-clustering matrix P_ij = Pr(z_i = z_j | data)
  P <- matrix(0, n, n)
  for (s in 1:S) {
    zs <- z_samps[s, ]
    P <- P + (outer(zs, zs, "==") * 1)
  }
  P <- P / S
  diag(P) <- 1
  
  ord <- order(z_sort, seq_len(n))
  P_ord <- P[ord, ord, drop = FALSE]
  
  dfP <- as.data.frame(as.table(P_ord))
  names(dfP) <- c("i", "j", "p")
  dfP$i <- as.integer(dfP$i)
  dfP$j <- as.integer(dfP$j)
  
  b <- cumsum(as.integer(table(z_sort[ord])))
  b <- b[-length(b)]
  
  pP <- ggplot(dfP, aes(x = j, y = i, fill = p)) +
    geom_raster() +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = "Posterior co-clustering probabilities", x = NULL, y = NULL, fill = "Pr(same)") +
    theme_minimal(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  
  if (length(b) > 0) {
    pP <- pP +
      geom_vline(xintercept = b + 0.5, linewidth = 0.25) +
      geom_hline(yintercept = b + 0.5, linewidth = 0.25)
  }
  
  ggsave(file_coclust, pP, width = 6.5, height = 5.5, dpi = 200)
  invisible(list(p_trace = p_trace, p_coclust = pP))
}

# ============================================================
# RUN: simulate -> fit -> report
# ============================================================

set.seed(3)
n <- 150
K_true <- 3
z_true <- rep(1:K_true, each = n / K_true)

Mu_true <- matrix(0, K_true, K_true); diag(Mu_true) <- 1.0
Sigma_true <- matrix(1.0, K_true, K_true)

Yw <- sim_sbm_gaussian(z_true, Mu_true, Sigma_true)

# Prior (match the collapse)
prior <- list(mu0 = 0, kappa0 = 0.1, alpha0 = 2, beta0 = 2)

# ---- uses YOUR existing functions ----
# make_log_m_nig(...)
# gibbs_sbm_undir_scalar(...)
K <- 3
log_m <- make_log_m_nig(
  mu0 = prior$mu0, kappa0 = prior$kappa0,
  alpha0 = prior$alpha0, beta0 = prior$beta0
)

fit <- gibbs_sbm_undir_scalar(
  Y = Yw, K = K, log_m_block = log_m,
  alpha = 1, n_iter = 5000, burn = 1000, thin = 10,
  verbose = FALSE
)

z_hat <- fit$z_map
ari_val <- ARI(z_true, z_hat)

# ---- concise text + tables ----
conf_mat <- table(true = z_true, inferred = z_hat)
sizes_hat <- as.data.frame(table(z_hat)) %>% dplyr::arrange(dplyr::desc(Freq))

df_true <- edge_df(Yw, z_true, "true")
df_hat  <- edge_df(Yw, z_hat,  "inferred")

within_between <- bind_rows(df_true, df_hat) %>%
  group_by(partition, type) %>%
  summarise(mean = mean(weight), sd = sd(weight), .groups = "drop")

cat("\n=== Case 3: Gaussian weighted SBM (NIG-collapsed) ===\n")
cat(sprintf("n = %d, K = %d\n", n, K))
cat(sprintf("ARI(true, MAP) = %.3f\n\n", ari_val))
cat("Cluster sizes (MAP):\n"); print(sizes_hat)
cat("\nConfusion matrix (true x inferred):\n"); print(conf_mat)
cat("\nWithin vs Between weights (mean ± sd):\n"); print(within_between)

# ---- block posterior summaries conditioned on z_hat ----
block_post <- block_posterior_summary(
  Y = Yw, z = z_hat, K = K,
  mu0 = prior$mu0, kappa0 = prior$kappa0,
  alpha0 = prior$alpha0, beta0 = prior$beta0
)

# Symmetric matrices for plotting
Mu_hat <- matrix(NA_real_, K, K)
Sigma_hat <- matrix(NA_real_, K, K)
for (r in seq_len(nrow(block_post))) {
  a <- block_post$a[r]; b <- block_post$b[r]
  Mu_hat[a, b] <- Mu_hat[b, a] <- block_post$mu_post_mean[r]
  Sigma_hat[a, b] <- Sigma_hat[b, a] <- block_post$sigma_post_mean[r]
}

# Save tables
write.csv(block_post, "case3_block_posterior_summary.csv", row.names = FALSE)
write.csv(as.data.frame(conf_mat), "case3_confusion_matrix.csv", row.names = FALSE)

# ---- figures ----
draw_process_diagram("case3_process.png")

p_true <- plot_weight_heatmap(Yw, z_true, "Edge weights sorted by TRUE communities")
p_hat  <- plot_weight_heatmap(Yw, z_hat,  "Edge weights sorted by INFERRED communities (MAP)")
ggsave("case3_heatmap_true.png", p_true, width = 6.5, height = 5.5, dpi = 200)
ggsave("case3_heatmap_inferred.png", p_hat, width = 6.5, height = 5.5, dpi = 200)

p_dens <- bind_rows(df_true, df_hat) %>%
  ggplot(aes(x = weight, fill = type)) +
  geom_density(alpha = 0.45) +
  facet_wrap(~ partition, nrow = 1) +
  labs(title = "Within vs between edge-weight distributions", x = "edge weight", y = "density", fill = NULL) +
  theme_minimal(base_size = 11)
ggsave("case3_within_between_density.png", p_dens, width = 10, height = 3.5, dpi = 200)

p_mu <- plot_block_heatmap(Mu_hat,    "Posterior E[mu_ab | z_hat] (NIG)")
p_si <- plot_block_heatmap(Sigma_hat, "Posterior E[sigma_ab | z_hat] (NIG)")
ggsave("case3_block_mu.png", p_mu, width = 5.5, height = 4.5, dpi = 200)
ggsave("case3_block_sigma.png", p_si, width = 5.5, height = 4.5, dpi = 200)

# FIXED call (no more FR negative-weight error)
plot_sparse_network(Yw, z_hat, q = 0.99, file = "case3_network_sparse.png")

# Optional “process” plots if z-samples are available in `fit`
z_samps <- extract_z_samples(fit, n, K = K)
if (!is.null(z_samps)) {
  plot_ari_trace_and_coclustering(
    z_samps = z_samps,
    z_true = z_true,
    z_sort = z_hat,
    file_trace = "case3_ari_trace.png",
    file_coclust = "case3_coclustering.png"
  )
}

cat("\nWrote:\n",
    "  case3_process.png\n",
    "  case3_heatmap_true.png\n",
    "  case3_heatmap_inferred.png\n",
    "  case3_within_between_density.png\n",
    "  case3_block_mu.png\n",
    "  case3_block_sigma.png\n",
    "  case3_network_sparse.png (if igraph installed)\n",
    "  case3_ari_trace.png + case3_coclustering.png (only if z-samples available)\n",
    "  case3_block_posterior_summary.csv\n",
    "  case3_confusion_matrix.csv\n", sep = "")
