# ============================================================
# Case 7: Zero-inflated Poisson SBM
#   - collapsed in (p, lambda)
#   - Z augmentation
#   - saves figures + a block-summary table + a short narrative summary
# ============================================================

# ---- Recommended packages for nicer plots/tables ----
# install.packages(c("ggplot2","dplyr","tidyr","patchwork","knitr"))
suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("dplyr",   quietly = TRUE) &&
    requireNamespace("tidyr",   quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE)
  if (!ok) stop("Please install ggplot2, dplyr, tidyr, patchwork (see install.packages line).")
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

# ----------------------------
# Helpers
# ----------------------------
choose2 <- function(n) n * (n - 1) / 2

ari <- function(x, y) {
  x <- as.integer(factor(x))
  y <- as.integer(factor(y))
  n <- length(x)
  tab <- table(x, y)
  a <- rowSums(tab)
  b <- colSums(tab)
  sum_ij <- sum(choose(tab, 2))
  sum_a  <- sum(choose(a, 2))
  sum_b  <- sum(choose(b, 2))
  expected <- (sum_a * sum_b) / choose(n, 2)
  max_index <- 0.5 * (sum_a + sum_b)
  (sum_ij - expected) / (max_index - expected)
}

log_dirichlet_multinomial_prior <- function(nk, alpha) {
  # symmetric Dirichlet with concentration alpha
  K <- length(nk)
  n <- sum(nk)
  a0 <- alpha / K
  lgamma(alpha) - lgamma(alpha + n) + sum(lgamma(a0 + nk) - lgamma(a0))
}

# ----------------------------
# ZIP-SBM simulation (vectorized)
# ----------------------------
sim_zip_sbm <- function(z,
                        p_active_in = 0.4, p_active_out = 0.1,
                        lam_in = 2.0, lam_out = 0.3) {
  n <- length(z)
  same <- outer(z, z, FUN = "==")
  pmat  <- ifelse(same, p_active_in, p_active_out)
  lammat <- ifelse(same, lam_in,      lam_out)
  
  ut <- upper.tri(pmat)
  p_ut <- pmat[ut]
  l_ut <- lammat[ut]
  
  active <- rbinom(length(p_ut), size = 1, prob = p_ut)
  y <- ifelse(active == 1, rpois(length(active), lambda = l_ut), 0L)
  
  A <- matrix(0L, n, n)
  A[ut] <- as.integer(y)
  A <- A + t(A)
  diag(A) <- 0L
  A
}

# ----------------------------
# Collapsed block log-marginal for ZIP
# stats: m_active (#active dyads), S_sum (sum of counts on active dyads), n_dyads
# ----------------------------
log_m_zip_block <- function(m_active, S_sum, n_dyads,
                            a_p = 1, b_p = 1, a_lam = 1, b_lam = 1) {
  # Beta-Bernoulli (Z) part
  lpZ <- lbeta(a_p + m_active, b_p + n_dyads - m_active) - lbeta(a_p, b_p)
  
  # Gamma-Poisson (A|Z=1) part (factorials omitted; constants wrt (z,Z) in this sampler)
  lpA <- a_lam * log(b_lam) - lgamma(a_lam) +
    lgamma(a_lam + S_sum) - (a_lam + S_sum) * log(b_lam + m_active)
  
  lpZ + lpA
}

zip_stats_from_state <- function(A, Z, z, K) {
  n <- nrow(A)
  m_active <- matrix(0, K, K)
  S_sum    <- matrix(0, K, K)
  
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    r <- z[i]; s <- z[j]
    a <- min(r, s); b <- max(r, s)
    m_active[a, b] <- m_active[a, b] + Z[i, j]
    S_sum[a, b]    <- S_sum[a, b] + Z[i, j] * A[i, j]
  }
  list(m_active = m_active, S_sum = S_sum)
}

# ----------------------------
# Collapsed Gibbs sampler (advanced: stores trace, MAP state incl. Z, saves stats)
# ----------------------------
gibbs_zip_sbm_adv <- function(A, K,
                              a_p = 1, b_p = 1, a_lam = 1, b_lam = 1,
                              alpha_part = 1,
                              n_iter = 3000, burn = 500, thin = 5,
                              z_true = NULL,
                              verbose = TRUE) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  n <- nrow(A)
  
  # dyad count helper
  n_dyad <- function(r, s, nk) if (r == s) choose2(nk[r]) else nk[r] * nk[s]
  
  # init z
  z  <- sample.int(K, n, replace = TRUE)
  nk <- tabulate(z, K)
  
  # init Z: if A>0 => active; if A==0 => random small probability
  Z <- matrix(0L, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    Z[i, j] <- if (A[i, j] > 0) 1L else rbinom(1, 1, 0.15)
    Z[j, i] <- Z[i, j]
  }
  diag(Z) <- 0L
  
  # block stats
  st <- zip_stats_from_state(A, Z, z, K)
  m_active <- st$m_active
  S_sum    <- st$S_sum
  
  loglik <- function(m_active, S_sum, nk) {
    out <- 0
    for (r in 1:K) for (s in r:K) {
      out <- out + log_m_zip_block(m_active[r, s], S_sum[r, s], n_dyad(r, s, nk),
                                   a_p, b_p, a_lam, b_lam)
    }
    out
  }
  
  # only need to update Z on zero edges
  zero_pairs <- which(A == 0 & upper.tri(A), arr.ind = TRUE)
  
  # trace storage (thinned, after burn)
  trace <- list(iter = integer(0), logpost = numeric(0), ari = numeric(0),
                n_active = integer(0), n_active_zeros = integer(0))
  
  best_lp <- -Inf
  z_map <- z
  Z_map <- Z
  nk_map <- nk
  m_map <- m_active
  S_map <- S_sum
  
  for (it in 1:n_iter) {
    
    # ---- Update Z for zero edges (A_ij == 0) ----
    if (nrow(zero_pairs) > 0) {
      for (k in 1:nrow(zero_pairs)) {
        i <- zero_pairs[k, 1]
        j <- zero_pairs[k, 2]
        
        r <- z[i]; s <- z[j]
        a <- min(r, s); b <- max(r, s)
        
        # remove current contribution if Z=1
        zcur <- Z[i, j]
        if (zcur == 1L) {
          m_active[a, b] <- m_active[a, b] - 1L
        }
        
        nrs <- n_dyad(a, b, nk)
        lp0 <- log_m_zip_block(m_active[a, b],     S_sum[a, b], nrs, a_p, b_p, a_lam, b_lam)
        lp1 <- log_m_zip_block(m_active[a, b] + 1, S_sum[a, b], nrs, a_p, b_p, a_lam, b_lam)
        
        p1 <- plogis(lp1 - lp0)
        Znew <- rbinom(1, 1, p1)
        
        if (Znew == 1L) m_active[a, b] <- m_active[a, b] + 1L
        Z[i, j] <- Znew
        Z[j, i] <- Znew
      }
    }
    
    # ---- Update z (collapsed in p, lambda given Z) ----
    for (i in 1:n) {
      old <- z[i]
      
      # per-group edge summaries from node i
      e_m <- numeric(K)  # sum Z_ij to group l
      e_S <- numeric(K)  # sum Z_ij * A_ij to group l
      for (l in 1:K) {
        idx <- which(z == l)
        if (l == old) idx <- idx[idx != i]
        if (length(idx) > 0) {
          e_m[l] <- sum(Z[i, idx])
          e_S[l] <- sum(Z[i, idx] * A[i, idx])
        }
      }
      
      # remove i from old
      for (l in 1:K) {
        if (e_m[l] == 0 && e_S[l] == 0) next
        a <- min(old, l); b <- max(old, l)
        m_active[a, b] <- m_active[a, b] - e_m[l]
        S_sum[a, b]    <- S_sum[a, b] - e_S[l]
      }
      nk[old] <- nk[old] - 1L
      z[i] <- 0L
      
      # candidate weights
      logw <- rep(-Inf, K)
      for (k in 1:K) {
        logprior_k <- log(nk[k] + alpha_part / K)
        delta <- 0
        
        for (l in 1:K) {
          if (l == k) {
            n_old <- n_dyad(k, k, nk)
            n_new <- choose2(nk[k] + 1L)
            lp_old <- log_m_zip_block(m_active[k, k],             S_sum[k, k],             n_old,
                                      a_p, b_p, a_lam, b_lam)
            lp_new <- log_m_zip_block(m_active[k, k] + e_m[l],    S_sum[k, k] + e_S[l],    n_new,
                                      a_p, b_p, a_lam, b_lam)
            delta <- delta + (lp_new - lp_old)
          } else {
            a <- min(k, l); b <- max(k, l)
            n_old <- n_dyad(a, b, nk)
            n_new <- (nk[k] + 1L) * nk[l]
            lp_old <- log_m_zip_block(m_active[a, b],          S_sum[a, b],          n_old,
                                      a_p, b_p, a_lam, b_lam)
            lp_new <- log_m_zip_block(m_active[a, b] + e_m[l], S_sum[a, b] + e_S[l], n_new,
                                      a_p, b_p, a_lam, b_lam)
            delta <- delta + (lp_new - lp_old)
          }
        }
        logw[k] <- logprior_k + delta
      }
      
      w <- exp(logw - max(logw))
      new <- sample.int(K, 1, prob = w)
      
      # add i to new
      for (l in 1:K) {
        if (e_m[l] == 0 && e_S[l] == 0) next
        a <- min(new, l); b <- max(new, l)
        m_active[a, b] <- m_active[a, b] + e_m[l]
        S_sum[a, b]    <- S_sum[a, b] + e_S[l]
      }
      nk[new] <- nk[new] + 1L
      z[i] <- new
    }
    
    # ---- Save trace + MAP (thinned) ----
    if (it >= burn && ((it - burn) %% thin == 0)) {
      lp <- log_dirichlet_multinomial_prior(nk, alpha_part) + loglik(m_active, S_sum, nk)
      
      trace$iter <- c(trace$iter, it)
      trace$logpost <- c(trace$logpost, lp)
      
      n_active <- sum(Z[upper.tri(Z)])
      n_active_zeros <- sum(Z[upper.tri(Z)] * (A[upper.tri(A)] == 0))
      trace$n_active <- c(trace$n_active, n_active)
      trace$n_active_zeros <- c(trace$n_active_zeros, n_active_zeros)
      
      if (!is.null(z_true)) trace$ari <- c(trace$ari, ari(z_true, z)) else trace$ari <- c(trace$ari, NA_real_)
      
      if (lp > best_lp) {
        best_lp <- lp
        z_map <- z
        Z_map <- Z
        nk_map <- nk
        m_map <- m_active
        S_map <- S_sum
      }
    }
    
    if (verbose && (it %% 100 == 0)) message(sprintf("iter %d / %d", it, n_iter))
  }
  
  list(
    z_map = z_map, Z_map = Z_map,
    nk_map = nk_map, m_active_map = m_map, S_sum_map = S_map,
    logpost_map = best_lp,
    z_last = z, Z_last = Z,
    trace = as.data.frame(trace)
  )
}

# ----------------------------
# Posterior block summaries + plotting helpers
# ----------------------------
block_posterior_summary <- function(m_active, S_sum, nk,
                                    a_p = 1, b_p = 1, a_lam = 1, b_lam = 1) {
  K <- length(nk)
  rows <- list()
  t <- 1
  for (r in 1:K) for (s in r:K) {
    n_dyads <- if (r == s) choose2(nk[r]) else nk[r] * nk[s]
    m <- m_active[r, s]
    S <- S_sum[r, s]
    
    ap <- a_p + m
    bp <- b_p + (n_dyads - m)
    p_mean <- ap / (ap + bp)
    p_ci <- qbeta(c(0.025, 0.975), ap, bp)
    
    shape <- a_lam + S
    rate  <- b_lam + m
    lam_mean <- shape / rate
    lam_ci <- qgamma(c(0.025, 0.975), shape = shape, rate = rate)
    
    rows[[t]] <- data.frame(
      r = r, s = s,
      n_dyads = n_dyads, m_active = m, S_sum = S,
      p_mean = p_mean, p_q025 = p_ci[1], p_q975 = p_ci[2],
      lambda_mean = lam_mean, lambda_q025 = lam_ci[1], lambda_q975 = lam_ci[2],
      mu_mean = p_mean * lam_mean
    )
    t <- t + 1
  }
  dplyr::bind_rows(rows)
}

block_heatmap_plot <- function(mat, title, fill_lab) {
  K <- nrow(mat)
  df <- expand.grid(r = 1:K, s = 1:K)
  df$value <- as.vector(mat)
  
  ggplot(df, aes(x = s, y = r, fill = value)) +
    geom_tile() +
    coord_equal() +
    scale_y_reverse(breaks = 1:K) +
    scale_x_continuous(breaks = 1:K) +
    labs(title = title, x = "block s", y = "block r", fill = fill_lab) +
    theme_minimal(base_size = 12)
}

adjacency_heatmap_plot <- function(A, z, title) {
  n <- nrow(A)
  ord <- order(z)
  Aord <- A[ord, ord, drop = FALSE]
  
  df <- as.data.frame(as.table(Aord))
  colnames(df) <- c("i", "j", "y")
  df$i <- as.integer(df$i)
  df$j <- as.integer(df$j)
  df$logy <- log1p(df$y)
  
  sizes <- as.numeric(table(z[ord]))
  bounds <- cumsum(sizes)
  if (length(bounds) >= 2) bounds <- bounds[-length(bounds)] + 0.5 else bounds <- numeric(0)
  
  ggplot(df, aes(i, j, fill = logy)) +
    geom_raster() +
    coord_equal() +
    scale_y_reverse() +
    labs(title = title, x = NULL, y = NULL, fill = "log(1 + A_ij)") +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    geom_vline(xintercept = bounds, linewidth = 0.3, alpha = 0.7) +
    geom_hline(yintercept = bounds, linewidth = 0.3, alpha = 0.7)
}

trace_plot <- function(trace_df) {
  ggplot(trace_df, aes(x = iter, y = logpost)) +
    geom_line() +
    labs(title = "Collapsed log-posterior trace (thinned)", x = "iteration", y = "log posterior") +
    theme_minimal(base_size = 12)
}

# ----------------------------
# One-command demo: run + summarize + save images/tables
# ----------------------------
run_case7_zip_demo <- function(seed = 8,
                               n = 120, K_true = 3, K_fit = 3,
                               p_active_in = 0.5, p_active_out = 0.1,
                               lam_in = 2.0, lam_out = 0.2,
                               a_p = 1, b_p = 1, a_lam = 1, b_lam = 1,
                               alpha_part = 1,
                               n_iter = 2000, burn = 500, thin = 5,
                               out_dir = "case7_outputs") {
  set.seed(seed)
  stopifnot(n %% K_true == 0)
  
  z_true <- rep(1:K_true, each = n / K_true)
  A <- sim_zip_sbm(z_true,
                   p_active_in = p_active_in, p_active_out = p_active_out,
                   lam_in = lam_in, lam_out = lam_out)
  
  fit <- gibbs_zip_sbm_adv(A, K = K_fit,
                           a_p = a_p, b_p = b_p, a_lam = a_lam, b_lam = b_lam,
                           alpha_part = alpha_part,
                           n_iter = n_iter, burn = burn, thin = thin,
                           z_true = z_true,
                           verbose = FALSE)
  
  z_hat <- fit$z_map
  
  # ---- summaries ----
  cat("\n================ Case 7: ZIP-SBM (collapsed) ================\n")
  cat(sprintf("n = %d, K_true = %d, K_fit = %d\n", n, K_true, K_fit))
  cat(sprintf("Edge sparsity (A_ij==0): %.1f%% of dyads\n", 100 * mean(A[upper.tri(A)] == 0)))
  cat(sprintf("ARI(z_true, z_map) = %.3f\n", ari(z_true, z_hat)))
  cat("Cluster sizes (MAP):\n"); print(table(z_hat))
  
  block_tbl <- block_posterior_summary(
    fit$m_active_map, fit$S_sum_map, fit$nk_map,
    a_p = a_p, b_p = b_p, a_lam = a_lam, b_lam = b_lam
  )
  
  # quick within vs between interpretation
  within <- dplyr::filter(block_tbl, r == s)
  between <- dplyr::filter(block_tbl, r < s)
  quick <- dplyr::bind_rows(
    data.frame(type = "within-block",  p = mean(within$p_mean),  lambda = mean(within$lambda_mean),  mu = mean(within$mu_mean)),
    data.frame(type = "between-block", p = mean(between$p_mean), lambda = mean(between$lambda_mean), mu = mean(between$mu_mean))
  )
  cat("\nPosterior means (averaged over blocks):\n")
  print(quick)
  
  # ---- make output dir ----
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- plots ----
  p_adj_true <- adjacency_heatmap_plot(A, z_true, "Adjacency heatmap (true ordering)")
  p_adj_map  <- adjacency_heatmap_plot(A, z_hat,  "Adjacency heatmap (MAP ordering)")
  p_trace    <- trace_plot(fit$trace)
  
  # build KxK matrices for heatmaps
  K <- K_fit
  Pmat <- matrix(NA_real_, K, K)
  Lmat <- matrix(NA_real_, K, K)
  Mumat <- matrix(NA_real_, K, K)
  for (ii in 1:nrow(block_tbl)) {
    r <- block_tbl$r[ii]; s <- block_tbl$s[ii]
    Pmat[r, s] <- block_tbl$p_mean[ii];      Pmat[s, r] <- Pmat[r, s]
    Lmat[r, s] <- block_tbl$lambda_mean[ii]; Lmat[s, r] <- Lmat[r, s]
    Mumat[r, s] <- block_tbl$mu_mean[ii];    Mumat[s, r] <- Mumat[r, s]
  }
  pP  <- block_heatmap_plot(Pmat,  "Posterior mean p_active (per block)", "p")
  pL  <- block_heatmap_plot(Lmat,  "Posterior mean lambda (per block)",   "lambda")
  pMu <- block_heatmap_plot(Mumat, "Posterior mean p*lambda (per dyad)",  "mu")
  
  # combine panels
  p_adj_panel <- p_adj_true | p_adj_map
  p_blk_panel <- pP | pL | pMu
  
  # ---- save images ----
  ggsave(file.path(out_dir, "adjacency_true.png"), p_adj_true, width = 5.5, height = 5.0, dpi = 300)
  ggsave(file.path(out_dir, "adjacency_map.png"),  p_adj_map,  width = 5.5, height = 5.0, dpi = 300)
  ggsave(file.path(out_dir, "trace_logpost.png"),  p_trace,    width = 6.5, height = 3.5, dpi = 300)
  ggsave(file.path(out_dir, "adjacency_panel.png"), p_adj_panel, width = 11.0, height = 5.0, dpi = 300)
  ggsave(file.path(out_dir, "block_heatmaps.png"),  p_blk_panel, width = 12.0, height = 4.0, dpi = 300)
  
  # ---- save tables ----
  write.csv(block_tbl, file.path(out_dir, "block_posterior_summary.csv"), row.names = FALSE)
  
  # Optional: write a LaTeX table file (nice if you're in RMarkdown/LaTeX workflow)
  if (requireNamespace("knitr", quietly = TRUE)) {
    tex <- knitr::kable(block_tbl, format = "latex", booktabs = TRUE, digits = 3,
                        caption = "ZIP-SBM block posterior summaries (MAP state).")
    writeLines(tex, con = file.path(out_dir, "block_posterior_summary.tex"))
  }
  
  # ---- “process diagram” (optional) ----
  # This shows the model graph (z -> p,lambda -> Z -> A). It will display in RStudio viewer.
  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    diagram <- DiagrammeR::grViz("
      digraph zip_sbm {
        rankdir=LR;
        node [shape=box, style=filled, color=lightgray];
        z[label='community labels z_i'];
        p[label='sparsity p_rs'];
        l[label='intensity lambda_rs'];
        Z[label='activity Z_ij'];
        A[label='counts A_ij'];
        z -> p; z -> l;
        p -> Z; l -> Z;
        Z -> A; l -> A;
      }
    ")
    # Print diagram to viewer
    print(diagram)
  }
  
  cat(sprintf("\nSaved outputs to: %s\n", normalizePath(out_dir)))
  cat("Key files:\n")
  cat("  - adjacency_map.png (block structure)\n")
  cat("  - trace_logpost.png (inference trace)\n")
  cat("  - block_heatmaps.png (p, lambda, p*lambda)\n")
  cat("  - block_posterior_summary.csv (+ .tex if knitr installed)\n")
  cat("=============================================================\n\n")
  
  invisible(list(A = A, z_true = z_true, fit = fit, block_tbl = block_tbl, quick = quick))
}

res <- run_case7_zip_demo()
