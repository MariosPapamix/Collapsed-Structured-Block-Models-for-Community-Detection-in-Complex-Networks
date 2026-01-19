# ============================================================
# Case 2: Count SBM (Gamma–Poisson collapse)
#  - Simulate undirected count-SBM (Poisson)
#  - Fit collapsed Gibbs (your gibbs_sbm_undir_scalar + make_log_m_gamma_poisson)
#  - Produce: summary.txt, CSV tables, and PNG figures
#
# IMPORTANT:
#   This script assumes you have already sourced/loaded:
#     - make_log_m_gamma_poisson(a, b)
#     - gibbs_sbm_undir_scalar(Y, K, log_m_block, alpha, n_iter, burn, thin, verbose)
#   and optionally ari() (otherwise we use mclust::adjustedRandIndex).
# ============================================================

# ---- 0) Packages ----
load_pkgs <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
load_pkgs(c("ggplot2", "dplyr", "readr"))

# Optional: patchwork (side-by-side plots)
if (!requireNamespace("patchwork", quietly = TRUE)) {
  message("Optional package 'patchwork' not installed (side-by-side plot will be saved separately).")
}

# Optional: process diagram
have_diagram <- requireNamespace("DiagrammeR", quietly = TRUE) &&
  requireNamespace("DiagrammeRsvg", quietly = TRUE) &&
  requireNamespace("rsvg", quietly = TRUE)

# ---- 1) Metrics ----
ari_safe <- function(z_true, z_hat) {
  if (exists("ari", mode = "function")) return(ari(z_true, z_hat))
  if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
  mclust::adjustedRandIndex(z_true, z_hat)
}

# ---- 2) Simulation: undirected Poisson SBM counts ----
sim_sbm_poisson <- function(z, Lambda) {
  n <- length(z)
  Y <- matrix(0L, n, n)
  for (i in 1:(n - 1)) {
    zi <- z[i]
    for (j in (i + 1):n) {
      lam <- Lambda[zi, z[j]]
      yij <- rpois(1, lam)
      Y[i, j] <- yij
      Y[j, i] <- yij
    }
  }
  diag(Y) <- 0L
  Y
}

# ---- 3) Block posterior table ----
# This computes the conjugate posterior for each block intensity:
#   lambda_{kℓ} | Y,z  ~ Gamma(a0 + S_{kℓ}, b0 + N_{kℓ})
# where:
#   S_{kk} = sum_{i<j in k} Y_ij,    N_{kk} = choose(n_k,2)
#   S_{kℓ} = sum_{i in k, j in ℓ} Y_ij,  N_{kℓ} = n_k*n_ℓ,  for k<ℓ
#
# The mirroring (k,l)->(l,k) must NOT use `transmute(k=l, l=k)` because
# dplyr creates columns sequentially; we use temporary names to swap safely.
block_posterior_table <- function(Y, z, K, a0 = 1, b0 = 1) {
  n <- length(z)
  z <- as.integer(factor(z, levels = 1:K))
  idx <- split(seq_len(n), z)
  
  out <- list()
  for (k in 1:K) {
    ik <- idx[[k]]
    for (l in k:K) {
      il <- idx[[l]]
      
      if (k == l) {
        if (length(ik) < 2) {
          N <- 0
          S <- 0
        } else {
          sub <- Y[ik, ik, drop = FALSE]
          S <- sum(sub[upper.tri(sub)])   # each dyad once
          N <- choose(length(ik), 2)
        }
      } else {
        sub <- Y[ik, il, drop = FALSE]    # each dyad once (undirected)
        S <- sum(sub)
        N <- length(ik) * length(il)
      }
      
      shape_post <- a0 + S
      rate_post  <- b0 + N
      
      out[[length(out) + 1]] <- data.frame(
        k = k, l = l,
        dyads = N,
        sum_counts = S,
        mean_count_per_dyad = if (N > 0) S / N else NA_real_,
        post_shape = shape_post,
        post_rate  = rate_post,
        lambda_post_mean = shape_post / rate_post,
        lambda_post_q025 = qgamma(0.025, shape = shape_post, rate = rate_post),
        lambda_post_q975 = qgamma(0.975, shape = shape_post, rate = rate_post)
      )
    }
  }
  
  upper <- dplyr::bind_rows(out)
  
  # Safe swap via temp names 
  lower <- upper |>
    dplyr::filter(k != l) |>
    dplyr::rename(k_old = k, l_old = l) |>
    dplyr::transmute(
      k = l_old,
      l = k_old,
      dyads, sum_counts, mean_count_per_dyad,
      post_shape, post_rate,
      lambda_post_mean, lambda_post_q025, lambda_post_q975
    )
  
  dplyr::bind_rows(upper, lower) |>
    dplyr::arrange(k, l)
}

# ---- 4) Plots ----
plot_adjacency_heatmap <- function(Y, z_order, K = NULL, title = "") {
  n <- nrow(Y)
  
  if (is.null(K)) K <- length(unique(z_order))
  z_order <- as.integer(factor(z_order, levels = 1:K))
  ord <- order(z_order)
  Yord <- Y[ord, ord, drop = FALSE]
  
  df <- expand.grid(i = seq_len(n), j = seq_len(n))
  df$logy <- log1p(as.vector(Yord))
  
  # block boundaries (include empty clusters safely)
  sizes <- as.integer(table(factor(z_order[ord], levels = 1:K)))
  cuts <- cumsum(sizes)
  cuts <- cuts[cuts > 0]
  if (length(cuts) >= 2) cuts <- cuts[-length(cuts)] + 0.5 else cuts <- numeric(0)
  
  ggplot(df, aes(x = j, y = i, fill = logy)) +
    geom_raster() +
    geom_vline(xintercept = cuts, linewidth = 0.2) +
    geom_hline(yintercept = cuts, linewidth = 0.2) +
    coord_equal(expand = FALSE) +
    scale_y_reverse() +
    labs(title = title, fill = "log(1+count)") +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
}

plot_block_heatmap <- function(block_df, K, title = "Posterior mean block intensity") {
  block_df |>
    dplyr::mutate(
      k = factor(k, levels = 1:K),
      l = factor(l, levels = 1:K),
      label = sprintf("%.2f", lambda_post_mean)
    ) |>
    ggplot(aes(x = l, y = k, fill = lambda_post_mean)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 3) +
    coord_equal() +
    labs(x = "Block ℓ", y = "Block k", title = title, fill = "E[λ|Y,z]") +
    theme_minimal(base_size = 11)
}

save_process_diagram <- function(path_png) {
  if (!have_diagram) {
    message("Process diagram skipped (install DiagrammeR, DiagrammeRsvg, rsvg to enable).")
    return(invisible(FALSE))
  }
  dot <- "
  digraph sbm_process {
    graph [layout = dot, rankdir = LR]
    node [shape = box, style = rounded, fontsize = 12]

    A [label = 'Simulate counts Y\\n(Poisson SBM)']
    B [label = 'Gamma prior on block rates\\nλ_{kℓ} ~ Gamma(a,b)']
    C [label = 'Collapse λ\\n(integrate out)']
    D [label = 'Collapsed Gibbs sampler\\nupdate labels z']
    E [label = 'MAP labels z_hat\\n+ block posteriors']
    F [label = 'Summaries + figures\\nheatmaps, tables']

    A -> B -> C -> D -> E -> F
  }
  "
  g <- DiagrammeR::grViz(dot)
  svg_txt <- DiagrammeRsvg::export_svg(g)
  rsvg::rsvg_png(charToRaw(svg_txt), file = path_png, width = 1200)
  invisible(TRUE)
}

# ============================================================
# 5) RUN: simulate -> fit -> summarize -> save outputs
# ============================================================

set.seed(2)

out_dir <- "case2_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Ground truth ---
n <- 150
K_true <- 3
z_true <- rep(1:K_true, each = n / K_true)

Lambda <- matrix(0.15, K_true, K_true)
diag(Lambda) <- 1.5

Ycnt <- sim_sbm_poisson(z_true, Lambda)

# Scramble node order (makes heatmap less trivially blocky before sorting)
perm <- sample.int(n)
Ycnt <- Ycnt[perm, perm]
z_true <- z_true[perm]

# --- Fit collapsed SBM (your functions must exist) ---
K <- 3
a0 <- 1
b0 <- 1

stopifnot(exists("make_log_m_gamma_poisson", mode = "function"))
stopifnot(exists("gibbs_sbm_undir_scalar", mode = "function"))

log_m <- make_log_m_gamma_poisson(a = a0, b = b0)

fit <- gibbs_sbm_undir_scalar(
  Y = Ycnt, K = K, log_m_block = log_m,
  alpha = 1, n_iter = 5000, burn = 1000, thin = 10, verbose = FALSE
)

stopifnot(!is.null(fit$z_map))
z_hat <- fit$z_map

# --- Metrics + tables ---
ari_val <- ari_safe(z_true, z_hat)

conf_mat <- table(
  True = factor(z_true, levels = 1:K_true),
  Inferred = factor(z_hat, levels = 1:K)
)

block_tbl <- block_posterior_table(Ycnt, z_hat, K = K, a0 = a0, b0 = b0)

# Quick sanity: in this simulation, diagonals should be near 1.5, off-diagonals near 0.15
diag_est <- block_tbl |> dplyr::filter(k == l) |> dplyr::arrange(k) |> dplyr::pull(lambda_post_mean)
off_est  <- block_tbl |> dplyr::filter(k != l) |> dplyr::pull(lambda_post_mean)

# Edge-count summaries (true within vs between)
ut <- which(upper.tri(Ycnt), arr.ind = TRUE)
edges <- data.frame(
  y = Ycnt[ut],
  same_block = (z_true[ut[, 1]] == z_true[ut[, 2]])
)

edge_summary <- edges |>
  dplyr::group_by(same_block) |>
  dplyr::summarise(
    mean = mean(y),
    p_zero = mean(y == 0),
    q90 = as.numeric(stats::quantile(y, 0.90)),
    .groups = "drop"
  )

# --- Write outputs ---
readr::write_csv(as.data.frame(conf_mat), file.path(out_dir, "confusion_table.csv"))
readr::write_csv(block_tbl, file.path(out_dir, "block_posterior_table.csv"))
saveRDS(
  list(Y = Ycnt, z_true = z_true, z_hat = z_hat, fit = fit, Lambda = Lambda, perm = perm),
  file.path(out_dir, "case2_data_and_fit.rds")
)

summary_txt <- c(
  "=== Case 2: Count SBM (Gamma–Poisson collapsed) ===",
  sprintf("n = %d, K_true = %d, K_fit = %d", n, K_true, K),
  sprintf("Gamma prior: shape a = %g, rate b = %g", a0, b0),
  sprintf("ARI(true, MAP) = %.4f", ari_val),
  "",
  "Diagonal posterior means (E[λ_kk|Y,z_hat]) :",
  paste(sprintf("  k=%d: %.3f", 1:K, diag_est), collapse = "\n"),
  sprintf("Off-diagonal posterior means: mean=%.3f (min=%.3f, max=%.3f)",
          mean(off_est), min(off_est), max(off_est)),
  "",
  "Within vs between (TRUE blocks) edge-count summary over i<j:",
  paste(
    apply(edge_summary, 1, function(r) {
      sprintf("  same_block=%s: mean=%.3f, P(y=0)=%.3f, q90=%.0f",
              r[["same_block"]], r[["mean"]], r[["p_zero"]], r[["q90"]])
    }),
    collapse = "\n"
  ),
  "",
  "Wrote:",
  "  - confusion_table.csv",
  "  - block_posterior_table.csv",
  "  - case2_data_and_fit.rds"
)

writeLines(summary_txt, file.path(out_dir, "summary.txt"))

cat(paste(summary_txt, collapse = "\n"), "\n\n")
print(conf_mat)

# --- Figures ---
p_true <- plot_adjacency_heatmap(
  Ycnt, z_true, K = K_true,
  title = "Adjacency (counts) ordered by TRUE labels"
)
p_hat <- plot_adjacency_heatmap(
  Ycnt, z_hat, K = K,
  title = "Adjacency (counts) ordered by INFERRED (MAP) labels"
)

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_both <- p_true + p_hat + patchwork::plot_layout(ncol = 2)
  ggsave(file.path(out_dir, "adjacency_heatmaps_true_vs_hat.png"),
         p_both, width = 12, height = 5, dpi = 300)
} else {
  ggsave(file.path(out_dir, "adjacency_heatmap_true.png"),
         p_true, width = 6, height = 5, dpi = 300)
  ggsave(file.path(out_dir, "adjacency_heatmap_hat.png"),
         p_hat, width = 6, height = 5, dpi = 300)
}

p_block <- plot_block_heatmap(
  block_tbl, K = K,
  title = "Block intensities from collapsed Gamma–Poisson posterior "
)
ggsave(file.path(out_dir, "block_intensity_heatmap.png"),
       p_block, width = 6, height = 5, dpi = 300)

p_hist <- ggplot(edges, aes(x = y, fill = same_block)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.55) +
  labs(
    title = "Edge-count distribution (true within vs between blocks)",
    x = "count y_ij", y = "number of dyads", fill = "same block?"
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "edge_count_hist_within_between.png"),
       p_hist, width = 7, height = 4.5, dpi = 300)

save_process_diagram(file.path(out_dir, "process_diagram.png"))

cat("\nSaved outputs to: ", normalizePath(out_dir), "\n", sep = "")
