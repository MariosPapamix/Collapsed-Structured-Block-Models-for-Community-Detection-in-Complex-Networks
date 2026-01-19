# ============================================================
# Case 4a (enhanced): Gap-constrained Bernoulli SBM
#  - Runs the sampler
#  - Produces a narrative summary (console + summary.txt)
#  - Produces a block-level table (CSV + optional LaTeX)
#  - Produces multiple figures (PNG): adjacency heatmaps, block probs, network plot, flow diagram
# ============================================================

# ---- Packages (lightweight + common) ----
pkgs <- c("ggplot2", "igraph", "grid")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# ---- Output folder ----
out_dir <- "case4_gap_outputs"
dir.create(out_dir, showWarnings = FALSE)

# ---- Helpers: block counts for undirected adjacency ----
block_counts_undir <- function(A, z, K = max(z)) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  stopifnot(length(z) == nrow(A))
  if (!all(z %in% 1:K)) stop("z must take values in 1:K")
  
  # Ensure 0/1 and symmetric (defensive; keeps it robust)
  A <- (A > 0) * 1L
  A <- (A | t(A)) * 1L
  diag(A) <- 0L
  
  edges <- matrix(0, K, K)
  possible <- matrix(0, K, K)
  
  idx <- lapply(1:K, function(k) which(z == k))
  nk  <- vapply(idx, length, integer(1))
  
  for (k in 1:K) {
    # within-block
    if (nk[k] >= 2) {
      sub <- A[idx[[k]], idx[[k]], drop = FALSE]
      edges[k, k] <- sum(sub[upper.tri(sub)])
      possible[k, k] <- choose(nk[k], 2)
    } else {
      edges[k, k] <- 0
      possible[k, k] <- 0
    }
    
    # between-block
    if (k < K) {
      for (l in (k + 1):K) {
        if (nk[k] > 0 && nk[l] > 0) {
          sub <- A[idx[[k]], idx[[l]], drop = FALSE]
          edges[k, l] <- sum(sub)
          possible[k, l] <- nk[k] * nk[l]
          # symmetry
          edges[l, k] <- edges[k, l]
          possible[l, k] <- possible[k, l]
        } else {
          edges[k, l] <- 0
          possible[k, l] <- 0
          edges[l, k] <- 0
          possible[l, k] <- 0
        }
      }
    }
  }
  
  list(edges = edges, possible = possible)
}

# ---- Helpers: truncated-Beta posterior summaries for off-diagonal blocks ----
trunc_beta_mean <- function(a, b, upper) {
  # Mean of Beta(a,b) truncated to [0, upper]
  # E[p | p<=upper] = (a/(a+b)) * pbeta(upper, a+1, b) / pbeta(upper, a, b)
  # computed stably with log.p
  stopifnot(upper > 0, upper < 1)
  logF0 <- pbeta(upper, a, b, log.p = TRUE)
  if (is.infinite(logF0)) return(upper)  # extreme tail; mean ~ upper
  logF1 <- pbeta(upper, a + 1, b, log.p = TRUE)
  (a / (a + b)) * exp(logF1 - logF0)
}

trunc_beta_ci <- function(a, b, upper, probs = c(0.025, 0.975)) {
  # Quantiles of Beta(a,b) truncated to [0,upper]:
  # q_trunc(p) = qbeta( p * pbeta(upper,a,b), a, b )
  stopifnot(all(probs > 0 & probs < 1))
  F0 <- pbeta(upper, a, b)
  if (F0 == 0) return(c(NA_real_, NA_real_))
  qbeta(probs * F0, a, b)
}

posterior_theta_summary <- function(edges, possible, a_in, b_in, a_out, b_out, x_cap) {
  K <- nrow(edges)
  stopifnot(ncol(edges) == K, all(dim(possible) == c(K, K)))
  
  post_mean <- matrix(NA_real_, K, K)
  post_lo   <- matrix(NA_real_, K, K)
  post_hi   <- matrix(NA_real_, K, K)
  
  for (k in 1:K) {
    for (l in 1:K) {
      N <- possible[k, l]
      if (N <= 0) next
      s <- edges[k, l]
      f <- N - s
      
      if (k == l) {
        a <- a_in + s
        b <- b_in + f
        post_mean[k, l] <- a / (a + b)
        q <- qbeta(c(0.025, 0.975), a, b)
        post_lo[k, l] <- q[1]; post_hi[k, l] <- q[2]
      } else {
        a <- a_out + s
        b <- b_out + f
        post_mean[k, l] <- trunc_beta_mean(a, b, x_cap)
        q <- trunc_beta_ci(a, b, x_cap, probs = c(0.025, 0.975))
        post_lo[k, l] <- q[1]; post_hi[k, l] <- q[2]
      }
    }
  }
  
  list(mean = post_mean, lo = post_lo, hi = post_hi)
}

# ---- Helper: tidy block table ----
make_block_table <- function(edges, possible, post_mean, post_lo, post_hi, x_cap) {
  K <- nrow(edges)
  rows <- list()
  r <- 1
  for (k in 1:K) {
    for (l in k:K) {
      N <- possible[k, l]
      s <- edges[k, l]
      dens <- if (N > 0) s / N else NA_real_
      rows[[r]] <- data.frame(
        block_k = k,
        block_l = l,
        type    = if (k == l) "within" else "between",
        possible_edges = N,
        observed_edges = s,
        empirical_density = dens,
        post_mean = post_mean[k, l],
        post_lo   = post_lo[k, l],
        post_hi   = post_hi[k, l],
        cap_xmax  = if (k == l) NA_real_ else x_cap,
        stringsAsFactors = FALSE
      )
      r <- r + 1
    }
  }
  do.call(rbind, rows)
}

# ---- Plotting helpers ----
plot_adj_heatmap_base <- function(A, z, main = "") {
  A <- (A > 0) * 1L
  A <- (A | t(A)) * 1L
  diag(A) <- 0L
  
  ord <- order(z)
  A_ord <- A[ord, ord, drop = FALSE]
  n <- nrow(A_ord)
  
  # Heatmap-like display: 0=white, 1=black
  image(
    1:n, 1:n, t(A_ord[n:1, ]),
    col = gray.colors(2, start = 1, end = 0),
    axes = FALSE, xlab = "", ylab = "", main = main, useRaster = TRUE
  )
  
  # Block boundaries
  K <- max(z)
  sizes <- as.integer(table(factor(z[ord], levels = 1:K)))
  cuts <- cumsum(sizes)
  # vertical lines
  abline(v = cuts + 0.5, lwd = 1)
  # horizontal lines (flip because y is reversed)
  abline(h = n - cuts + 0.5, lwd = 1)
}

save_adj_heatmaps <- function(A, z_true, z_hat, file, res = 160) {
  png(file, width = 1800, height = 900, res = res)
  op <- par(mfrow = c(1, 2), mar = c(1, 1, 4, 1))
  plot_adj_heatmap_base(A, z_true, main = "Adjacency sorted by TRUE communities")
  plot_adj_heatmap_base(A, z_hat,  main = "Adjacency sorted by INFERRED communities")
  par(op)
  dev.off()
}

mat_to_df <- function(M) {
  K <- nrow(M)
  df <- expand.grid(from = 1:K, to = 1:K)
  df$val <- as.vector(M)
  df
}

save_block_heatmap <- function(M, file, title = "", lim = NULL) {
  df <- mat_to_df(M)
  p <- ggplot(df, aes(x = to, y = from, fill = val)) +
    geom_tile() +
    coord_equal() +
    scale_y_reverse(breaks = sort(unique(df$from))) +
    labs(x = "Block", y = "Block", fill = "prob", title = title) +
    theme_minimal()
  
  if (!is.null(lim)) {
    p <- p + scale_fill_continuous(limits = lim)
  }
  ggsave(file, p, width = 4.6, height = 4.1, dpi = 160)
}

save_network_plot <- function(A, z, file, title = "", seed = 1) {
  A <- (A > 0) * 1L
  A <- (A | t(A)) * 1L
  diag(A) <- 0L
  
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  V(g)$grp <- factor(z)
  
  set.seed(seed)
  lay <- igraph::layout_with_fr(g)
  
  png(file, width = 1200, height = 900, res = 160)
  plot(
    g, layout = lay,
    vertex.size = 4, vertex.label = NA,
    vertex.color = as.integer(V(g)$grp),
    edge.arrow.mode = 0
  )
  title(main = title)
  dev.off()
}

save_process_flow <- function(file, x_cap) {
  png(file, width = 1600, height = 420, res = 160)
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 3))
  rect_par <- list(col = "white")
  
  boxes <- data.frame(
    x0 = c(0.4, 2.4, 4.4, 6.4, 8.4),
    x1 = c(2.0, 4.0, 6.0, 8.0, 9.8),
    y0 = rep(1.0, 5),
    y1 = rep(2.6, 5),
    txt = c(
      "1) Simulate SBM\n(assortative)",
      sprintf("2) Priors\nwithin: Beta\nbetween: trunc Beta\n[0, %.3f]", x_cap),
      "3) Collapse Θ\n(Beta-Binomial)\n(trunc off-diag)",
      "4) Gibbs updates\nupdate z_i\nfrom counts",
      "5) Outputs\nARI, tables,\nheatmaps, network"
    ),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(boxes)) {
    rect(boxes$x0[i], boxes$y0[i], boxes$x1[i], boxes$y1[i], border = "black")
    text((boxes$x0[i] + boxes$x1[i]) / 2, (boxes$y0[i] + boxes$y1[i]) / 2, boxes$txt[i], cex = 0.95)
    if (i < nrow(boxes)) {
      arrows(boxes$x1[i], 1.8, boxes$x0[i + 1], 1.8, length = 0.10)
    }
  }
  title(main = "What is happening in Case 4a (gap-constrained SBM)")
  dev.off()
}

# ============================================================
# 1) Simulate data (same as your original case)
# ============================================================
set.seed(4)
n <- 150
K_true <- 3
z_true <- rep(1:K_true, each = n / K_true)

P <- matrix(0.005, K_true, K_true)  # very sparse between
diag(P) <- 0.12                     # moderate within

A <- sim_sbm_binary(z_true, P)

# ============================================================
# 2) Fit: diagonal Beta, off-diagonal truncated Beta on [0, x_cap]
# ============================================================
a_in <- 1; b_in <- 1
a_out <- 1; b_out <- 1
x_cap <- 0.02

log_m_in  <- make_log_m_beta_bernoulli(a_in, b_in)
log_m_out <- make_log_m_truncbeta(a_out, b_out, x_cap)

log_m_gap <- function(s1, s2, n, diag_flag) {
  if (diag_flag) log_m_in(s1, s2, n, TRUE) else log_m_out(s1, s2, n, FALSE)
}

fit <- gibbs_sbm_undir_scalar(
  Y = A, K = K_true, log_m_block = log_m_gap,
  alpha = 1, n_iter = 5000, burn = 1000, thin = 10, verbose = FALSE
)

z_hat <- fit$z_map

# ============================================================
# 3) Console summary + write summary.txt
# ============================================================
ari_val <- ari(z_true, z_hat)
size_tbl <- table(z_hat)

# overall density (upper triangle only)
overall_density <- mean(A[upper.tri(A)])

# block summaries under inferred labels
bc_hat <- block_counts_undir(A, z_hat, K_true)
theta_emp_hat <- with(bc_hat, edges / ifelse(possible > 0, possible, NA_real_))
post <- posterior_theta_summary(
  edges = bc_hat$edges, possible = bc_hat$possible,
  a_in = a_in, b_in = b_in,
  a_out = a_out, b_out = b_out,
  x_cap = x_cap
)

within_emp  <- sum(diag(bc_hat$edges)) / sum(diag(bc_hat$possible))
between_emp <- sum(bc_hat$edges[upper.tri(bc_hat$edges)]) / sum(bc_hat$possible[upper.tri(bc_hat$possible)])

within_post_min <- min(diag(post$mean), na.rm = TRUE)
between_post_max <- max(post$mean[upper.tri(post$mean)], na.rm = TRUE)
gap_post <- within_post_min - between_post_max

summary_lines <- c(
  "Gap-constrained Bernoulli SBM (Case 4a)",
  "--------------------------------------",
  sprintf("n = %d, K = %d", n, K_true),
  sprintf("True within prob (diag P) ≈ %.3f; true between prob ≈ %.3f", P[1,1], P[1,2]),
  sprintf("Between-block cap x_cap = %.3f", x_cap),
  "",
  sprintf("Overall observed edge density: %.4f", overall_density),
  sprintf("ARI(z_true, z_hat): %.4f", ari_val),
  sprintf("Inferred cluster sizes: %s",
          paste(names(size_tbl), as.integer(size_tbl), sep=":", collapse="  ")),
  "",
  sprintf("Empirical densities under z_hat: within = %.4f, between = %.4f", within_emp, between_emp),
  sprintf("Posterior separation (min within mean - max between mean): %.4f", gap_post),
  sprintf("Max between posterior mean: %.4f (must be <= x_cap = %.4f)", between_post_max, x_cap),
  "",
  "Tip: check the saved adjacency heatmaps and block-probability heatmaps to visually confirm the gap."
)

cat(paste0(summary_lines, collapse = "\n"), "\n")
writeLines(summary_lines, file.path(out_dir, "summary.txt"))

# A basic contingency table (label-invariant quality is ARI; this is just descriptive)
cat("\nContingency table (true vs inferred labels):\n")
print(table(true = z_true, inferred = z_hat))

# ============================================================
# 4) Block table (CSV + optional LaTeX)
# ============================================================
block_tbl <- make_block_table(
  edges = bc_hat$edges, possible = bc_hat$possible,
  post_mean = post$mean, post_lo = post$lo, post_hi = post$hi,
  x_cap = x_cap
)

# Save CSV
write.csv(block_tbl, file.path(out_dir, "block_summary.csv"), row.names = FALSE)

# Optional: save a LaTeX table if xtable is installed
if (requireNamespace("xtable", quietly = TRUE)) {
  xt <- xtable::xtable(block_tbl, digits = c(0,0,0,0,0,0,4,4,4,4,4))
  print(xt, include.rownames = FALSE, file = file.path(out_dir, "block_summary.tex"))
}

# ============================================================
# 5) Figures (PNG)
# ============================================================

# 5a) Adjacency heatmaps (side-by-side)
save_adj_heatmaps(
  A = A, z_true = z_true, z_hat = z_hat,
  file = file.path(out_dir, "adjacency_true_vs_hat.png")
)

# 5b) Block probability heatmaps (true P, empirical under z_hat, posterior mean under z_hat)
lim_all <- range(c(P, theta_emp_hat, post$mean), na.rm = TRUE)
save_block_heatmap(P,            file.path(out_dir, "P_true.png"),
                   title = "True block probabilities P", lim = lim_all)
save_block_heatmap(theta_emp_hat, file.path(out_dir, "P_empirical_under_zhat.png"),
                   title = "Empirical block densities under z_hat", lim = lim_all)
save_block_heatmap(post$mean,     file.path(out_dir, "P_posteriorMean_gapPrior.png"),
                   title = sprintf("Posterior mean (gap prior, cap=%.3f)", x_cap), lim = lim_all)

# 5c) Network plot colored by inferred communities
save_network_plot(
  A = A, z = z_hat,
  file = file.path(out_dir, "network_inferred.png"),
  title = "Network colored by inferred communities",
  seed = 4
)

# 5d) A simple “process diagram” image
save_process_flow(file.path(out_dir, "process_flow.png"), x_cap = x_cap)

cat("\nSaved outputs to: ", normalizePath(out_dir), "\n", sep = "")
cat(" - summary.txt\n - block_summary.csv (and block_summary.tex if xtable installed)\n")
cat(" - adjacency_true_vs_hat.png\n - P_true.png\n - P_empirical_under_zhat.png\n - P_posteriorMean_gapPrior.png\n - network_inferred.png\n - process_flow.png\n")
