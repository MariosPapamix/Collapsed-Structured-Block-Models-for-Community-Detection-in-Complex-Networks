# ============================================================
# Case 1 (advanced): Binary SBM + collapsed Beta–Bernoulli Gibbs
# Produces:
#   - outputs_case1/figures/case1_overview.pdf  (multi-panel)
#   - outputs_case1/tables/case1_summary_table.tex
#   - outputs_case1/case1_summary.csv
# Also prints a short plain-language summary to console.
# ============================================================

# ----------------------------
# Utilities (base R only)
# ----------------------------
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
  # symmetric Dirichlet-multinomial prior integrated over pi
  K <- length(nk)
  n <- sum(nk)
  sum(lgamma(alpha / K + nk) - lgamma(alpha / K)) + lgamma(alpha) - lgamma(alpha + n)
}

# generate all permutations (K small; fine for K=3)
all_perms <- function(v) {
  if (length(v) == 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) {
    rest <- v[-i]
    subp <- all_perms(rest)
    for (p in subp) out[[length(out) + 1]] <- c(v[i], p)
  }
  out
}

best_perm_align <- function(z_hat, z_ref, K) {
  # find permutation of labels {1..K} that maximizes matches with z_ref
  perms <- all_perms(1:K)
  best_score <- -Inf
  best_p <- 1:K
  for (p in perms) {
    z_map <- p[z_hat]
    score <- sum(z_map == z_ref)
    if (score > best_score) { best_score <- score; best_p <- p }
  }
  best_p
}

order_by_labels <- function(z) order(z, seq_along(z))

block_counts <- function(A, z, K) {
  n <- nrow(A)
  nk <- tabulate(z, K)
  m <- matrix(0, K, K)  # store upper triangle counts
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    r <- z[i]; s <- z[j]
    a <- min(r, s); b <- max(r, s)
    m[a, b] <- m[a, b] + A[i, j]
  }
  # dyads per block
  nd <- matrix(0, K, K)
  for (r in 1:K) for (s in r:K) {
    nd[r, s] <- if (r == s) choose2(nk[r]) else nk[r] * nk[s]
  }
  list(m = m, nd = nd, nk = nk)
}

block_density_matrix <- function(m, nd, K) {
  D <- matrix(0, K, K)
  for (r in 1:K) for (s in 1:K) {
    a <- min(r, s); b <- max(r, s)
    D[r, s] <- if (nd[a, b] > 0) m[a, b] / nd[a, b] else NA_real_
  }
  D
}

plot_matrix_with_blocks <- function(M, ord, boundaries, main, col, zlim = NULL) {
  M2 <- M[ord, ord, drop = FALSE]
  n <- nrow(M2)
  
  # image uses lower-left origin; flip rows to look like a matrix
  if (is.null(zlim)) zlim <- range(M2, finite = TRUE)
  image(1:n, 1:n, t(M2[n:1, ]), col = col, axes = FALSE, xlab = "", ylab = "",
        main = main, zlim = zlim)
  
  # block separators
  if (length(boundaries) >= 2) {
    # boundaries are cumulative sizes, e.g., c(n1, n1+n2, n)
    for (b in boundaries[-length(boundaries)]) {
      abline(v = b + 0.5, col = "grey70", lwd = 1)
      abline(h = (n - b) + 0.5, col = "grey70", lwd = 1)
    }
  }
}

# ----------------------------
# Synthetic binary SBM
# ----------------------------
sim_sbm_binary <- function(z, P) {
  n <- length(z)
  A <- matrix(0, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    p <- P[z[i], z[j]]
    aij <- rbinom(1, 1, p)
    A[i, j] <- aij
    A[j, i] <- aij
  }
  diag(A) <- 0
  A
}

# ----------------------------
# Collapsed Beta–Bernoulli Gibbs with trace + samples
# ----------------------------
gibbs_bb_collapsed_trace <- function(A, K, a = 1, b = 1, alpha = 1,
                                     n_iter = 5000, burn = 1000, thin = 10,
                                     z_init = NULL, z_true = NULL,
                                     verbose = TRUE) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  n <- nrow(A)
  diag(A) <- 0
  
  # init labels
  if (is.null(z_init)) z <- sample.int(K, n, replace = TRUE) else z <- as.integer(z_init)
  nk <- tabulate(z, K)
  
  # init m counts (upper triangle only)
  m <- matrix(0, K, K)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    r <- z[i]; s <- z[j]
    aa <- min(r, s); bb <- max(r, s)
    m[aa, bb] <- m[aa, bb] + A[i, j]
  }
  
  n_dyad <- function(r, s, nk) if (r == s) choose2(nk[r]) else nk[r] * nk[s]
  log_m_block <- function(m_rs, n_rs) lbeta(a + m_rs, b + n_rs - m_rs) - lbeta(a, b)
  
  loglik_from_stats <- function(m, nk) {
    out <- 0
    for (r in 1:K) for (s in r:K) {
      out <- out + log_m_block(m[r, s], n_dyad(r, s, nk))
    }
    out
  }
  
  # storage
  keep <- floor((n_iter - burn) / thin)
  samples <- matrix(NA_integer_, n, keep)
  trace <- data.frame(iter = integer(keep), logpost = numeric(keep), ari = numeric(keep))
  
  best_lp <- -Inf
  z_map <- z
  
  kidx <- 0
  
  for (it in 1:n_iter) {
    for (i in 1:n) {
      old <- z[i]
      
      # edges from i to each community (excluding i if in old group)
      e <- numeric(K)
      for (l in 1:K) {
        idx <- which(z == l)
        if (l == old) idx <- idx[idx != i]
        if (length(idx)) e[l] <- sum(A[i, idx])
      }
      
      # remove i from old (update m)
      for (l in 1:K) {
        if (e[l] == 0) next
        aa <- min(old, l); bb <- max(old, l)
        m[aa, bb] <- m[aa, bb] - e[l]
      }
      nk[old] <- nk[old] - 1
      z[i] <- 0
      
      # candidate weights
      logw <- rep(-Inf, K)
      for (k in 1:K) {
        logprior_k <- log(nk[k] + alpha / K)
        delta <- 0
        
        # blocks involving k change when adding i to k
        for (l in 1:K) {
          if (l == k) {
            # diagonal block k,k
            n_old <- n_dyad(k, k, nk)
            n_new <- choose2(nk[k] + 1)
            m_old <- m[k, k]
            m_new <- m_old + e[l]
            delta <- delta + (log_m_block(m_new, n_new) - log_m_block(m_old, n_old))
          } else {
            aa <- min(k, l); bb <- max(k, l)
            n_old <- n_dyad(aa, bb, nk)
            n_new <- (nk[k] + 1) * nk[l]
            m_old <- m[aa, bb]
            m_new <- m_old + e[l]
            delta <- delta + (log_m_block(m_new, n_new) - log_m_block(m_old, n_old))
          }
        }
        logw[k] <- logprior_k + delta
      }
      
      w <- exp(logw - max(logw))
      new <- sample.int(K, 1, prob = w)
      
      # add i to new (update m)
      for (l in 1:K) {
        if (e[l] == 0) next
        aa <- min(new, l); bb <- max(new, l)
        m[aa, bb] <- m[aa, bb] + e[l]
      }
      nk[new] <- nk[new] + 1
      z[i] <- new
    }
    
    # record
    if (it > burn && ((it - burn) %% thin == 0)) {
      kidx <- kidx + 1
      lp <- log_dirichlet_multinomial_prior(nk, alpha) + loglik_from_stats(m, nk)
      samples[, kidx] <- z
      trace$iter[kidx] <- it
      trace$logpost[kidx] <- lp
      trace$ari[kidx] <- if (!is.null(z_true)) ari(z_true, z) else NA_real_
      
      if (lp > best_lp) {
        best_lp <- lp
        z_map <- z
      }
    }
    
    if (verbose && (it %% 250 == 0)) message(sprintf("iter %d / %d", it, n_iter))
  }
  
  list(z_map = z_map, logpost_map = best_lp, samples = samples, trace = trace)
}

# ----------------------------
# Posterior similarity matrix (PSM)
# ----------------------------
psm_from_samples <- function(samples) {
  n <- nrow(samples)
  S <- ncol(samples)
  P <- matrix(0, n, n)
  for (s in 1:S) {
    zs <- samples[, s]
    P <- P + outer(zs, zs, FUN = "==")
  }
  P / S
}

# ----------------------------
# Main run (you can tweak these)
# ----------------------------
set.seed(1)

# Simulation parameters
n <- 150
K_true <- 3
z_true <- rep(1:K_true, each = n / K_true)

p_in <- 0.15
p_out <- 0.02
P <- matrix(p_out, K_true, K_true); diag(P) <- p_in

A <- sim_sbm_binary(z_true, P)

# Prior + sampler settings
K <- 3
a <- 1; b <- 1               # Beta(a,b)
alpha <- 1                   # Dirichlet concentration for z prior
n_iter <- 5000; burn <- 1000; thin <- 10

# Run collapsed sampler
fit <- gibbs_bb_collapsed_trace(
  A = A, K = K, a = a, b = b, alpha = alpha,
  n_iter = n_iter, burn = burn, thin = thin,
  z_true = z_true, verbose = FALSE
)

z_map <- fit$z_map
trace <- fit$trace
samples <- fit$samples

# Align MAP labels to truth (only for prettier plots/tables; ARI itself is permutation-invariant)
perm <- best_perm_align(z_map, z_true, K)
z_map_aligned <- perm[z_map]

# Compute PSM
PSM <- psm_from_samples(samples)

# Block summaries (true vs MAP)
true_bc <- block_counts(A, z_true, K_true)
map_bc  <- block_counts(A, z_map_aligned, K)

D_true <- block_density_matrix(true_bc$m, true_bc$nd, K_true)
D_map  <- block_density_matrix(map_bc$m,  map_bc$nd,  K)

# Posterior mean of block probabilities under MAP:
# E[p_rs | A, z] = (a + m_rs) / (a + b + n_rs)
pmean <- matrix(NA_real_, K, K)
for (r in 1:K) for (s in 1:K) {
  aa <- min(r, s); bb <- max(r, s)
  nrs <- map_bc$nd[aa, bb]
  mrs <- map_bc$m[aa, bb]
  pmean[r, s] <- if (nrs > 0) (a + mrs) / (a + b + nrs) else NA_real_
}

# Global summaries
upper <- upper.tri(A)
density_global <- mean(A[upper])
avg_degree <- mean(rowSums(A))

ari_map <- ari(z_true, z_map_aligned)

within_true <- mean(diag(D_true), na.rm = TRUE)
between_true <- mean(D_true[upper.tri(D_true)], na.rm = TRUE)
within_map <- mean(diag(D_map), na.rm = TRUE)
between_map <- mean(D_map[upper.tri(D_map)], na.rm = TRUE)

# Confusion table
conf <- table(z_true, z_map_aligned)

# ----------------------------
# Output folders
# ----------------------------
out_dir <- "outputs_case1"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)

# ----------------------------
# Make multi-panel overview figure
# ----------------------------
ord_true <- order_by_labels(z_true)
ord_map  <- order_by_labels(z_map_aligned)

bound_true <- cumsum(tabulate(z_true, K_true))
bound_map  <- cumsum(tabulate(z_map_aligned, K))

pdf(file.path(fig_dir, "case1_overview.pdf"), width = 12, height = 7)
par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))

# (1) adjacency sorted by truth
plot_matrix_with_blocks(
  M = A, ord = ord_true, boundaries = bound_true,
  main = "Adjacency (sorted by truth)",
  col = c("white", "black"), zlim = c(0, 1)
)

# (2) adjacency sorted by MAP
plot_matrix_with_blocks(
  M = A, ord = ord_map, boundaries = bound_map,
  main = "Adjacency (sorted by MAP)",
  col = c("white", "black"), zlim = c(0, 1)
)

# (3) posterior mean block probs under MAP
image(1:K, 1:K, t(pmean[K:1, ]), axes = FALSE, xlab = "", ylab = "",
      main = "Posterior mean  E[p_rs | A, z_MAP]",
      col = gray(seq(1, 0, length.out = 256)),
      zlim = range(pmean, finite = TRUE))
axis(1, at = 1:K, labels = 1:K); axis(2, at = 1:K, labels = K:1)

# (4) PSM sorted by MAP
plot_matrix_with_blocks(
  M = PSM, ord = ord_map, boundaries = bound_map,
  main = "Posterior similarity matrix (PSM)",
  col = gray(seq(1, 0, length.out = 256)), zlim = c(0, 1)
)

# (5) log-posterior trace
plot(trace$iter, trace$logpost, type = "l",
     xlab = "Iteration", ylab = "Collapsed log-posterior",
     main = "Trace: collapsed log-posterior")

# (6) ARI trace (only meaningful because synthetic truth is known)
plot(trace$iter, trace$ari, type = "l", ylim = c(0, 1),
     xlab = "Iteration", ylab = "ARI(z_true, z_t)",
     main = "Trace: recovery (ARI)")

dev.off()

# ----------------------------
# Write LaTeX summary table (booktabs)
# ----------------------------
summary_df <- data.frame(
  Metric = c(
    "n", "K_true", "K_fit",
    "p_in (sim)", "p_out (sim)",
    "Global edge density", "Average degree",
    "Within density (true partition)", "Between density (true partition)",
    "Within density (MAP)", "Between density (MAP)",
    "ARI(z_true, z_MAP)"
  ),
  Value = c(
    n, K_true, K,
    p_in, p_out,
    density_global, avg_degree,
    within_true, between_true,
    within_map, between_map,
    ari_map
  )
)

# nice numeric formatting for LaTeX
fmt <- function(x) if (is.numeric(x)) sprintf("%.4f", x) else as.character(x)

tex_lines <- c(
  "\\begin{tabular}{lr}",
  "\\toprule",
  "Metric & Value\\\\",
  "\\midrule"
)
for (i in 1:nrow(summary_df)) {
  tex_lines <- c(tex_lines, paste0(summary_df$Metric[i], " & ", fmt(summary_df$Value[i]), "\\\\"))
}
tex_lines <- c(tex_lines, "\\bottomrule", "\\end{tabular}")

writeLines(tex_lines, file.path(tab_dir, "case1_summary_table.tex"))

# Also write CSV
write.csv(summary_df, file.path(out_dir, "case1_summary.csv"), row.names = FALSE)

# ----------------------------
# Console summary (plain-language)
# ----------------------------
cat("\n================ Case 1 summary ================\n")
cat(sprintf("Simulated SBM with K_true=%d, n=%d, p_in=%.3f, p_out=%.3f\n",
            K_true, n, p_in, p_out))
cat(sprintf("Collapsed inference (Beta–Bernoulli) with K=%d, a=%.2f, b=%.2f, alpha=%.2f\n",
            K, a, b, alpha))
cat(sprintf("MAP ARI: %.3f\n", ari_map))
cat(sprintf("True within/between densities: %.3f / %.3f\n", within_true, between_true))
cat(sprintf("MAP  within/between densities: %.3f / %.3f\n", within_map, between_map))
cat("\nConfusion table (truth rows x MAP columns; labels aligned for display):\n")
print(conf)

cat("\nWrote:\n")
cat(" - ", file.path(fig_dir, "case1_overview.pdf"), "\n", sep = "")
cat(" - ", file.path(tab_dir, "case1_summary_table.tex"), "\n", sep = "")
cat(" - ", file.path(out_dir, "case1_summary.csv"), "\n", sep = "")
cat("================================================\n\n")
