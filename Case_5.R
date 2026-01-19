# ============================================================
# Case 5+: Directed dyad-state SBM (collapsed Dirichlet–multinomial)
# + summaries, tables, and images
# ============================================================

# ----------------------------
# 0) Small utilities
# ----------------------------
choose2 <- function(x) x * (x - 1) / 2

ari <- function(z1, z2) {
  # Adjusted Rand Index (base R, no packages)
  stopifnot(length(z1) == length(z2))
  tab <- table(z1, z2)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  nij2 <- sum(choose2(tab))
  ai2  <- sum(choose2(rowSums(tab)))
  bj2  <- sum(choose2(colSums(tab)))
  n2   <- choose2(n)
  expected <- (ai2 * bj2) / n2
  max_idx  <- 0.5 * (ai2 + bj2)
  if (max_idx - expected == 0) return(0)
  (nij2 - expected) / (max_idx - expected)
}

log_multibeta <- function(alpha_vec) sum(lgamma(alpha_vec)) - lgamma(sum(alpha_vec))

log_m_dirichlet_multinomial <- function(counts_vec, alpha_vec) {
  # log p(counts | alpha) up to multinomial coefficient (cancels in ratios)
  log_multibeta(alpha_vec + counts_vec) - log_multibeta(alpha_vec)
}

log_dirichlet_multinomial_prior <- function(nk, alpha_part) {
  # Integrated-out mixing proportions: z ~ Dirichlet-multinomial(alpha_part/K)
  K <- length(nk)
  log_multibeta(nk + alpha_part / K) - log_multibeta(rep(alpha_part / K, K))
}

dyad_cat <- function(a_ij, a_ji) {
  # Map (A_ij, A_ji) to category index in {1,2,3,4} for (00,10,01,11)
  if (a_ij == 0 && a_ji == 0) return(1L)
  if (a_ij == 1 && a_ji == 0) return(2L)
  if (a_ij == 0 && a_ji == 1) return(3L)
  4L
}

dyad_counts_from_bits <- function(a1, a2) {
  # counts for categories (00,10,01,11) given two parallel bit-vectors
  c(
    sum(a1 == 0 & a2 == 0),  # 00
    sum(a1 == 1 & a2 == 0),  # 10
    sum(a1 == 0 & a2 == 1),  # 01
    sum(a1 == 1 & a2 == 1)   # 11
  )
}

within_counts_i_group <- function(i, idx, A) {
  # For within-group block counts we use a canonical node ordering (min index -> max index)
  idx_gt <- idx[idx > i]
  idx_lt <- idx[idx < i]
  c1 <- if (length(idx_gt)) dyad_counts_from_bits(A[i, idx_gt], A[idx_gt, i]) else c(0, 0, 0, 0)
  c2 <- if (length(idx_lt)) dyad_counts_from_bits(A[idx_lt, i], A[i, idx_lt]) else c(0, 0, 0, 0)
  c1 + c2
}

block_counts_unordered <- function(A, z, K) {
  # counts[a,b,] for a<=b, categories (00,10,01,11) relative to (a -> b, b -> a)
  n <- nrow(A)
  counts <- array(0L, dim = c(K, K, 4))
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    r <- z[i]; s <- z[j]
    a <- min(r, s); b <- max(r, s)
    if (r <= s) { a1 <- A[i, j]; a2 <- A[j, i] } else { a1 <- A[j, i]; a2 <- A[i, j] }
    cat <- dyad_cat(a1, a2)
    counts[a, b, cat] <- counts[a, b, cat] + 1L
  }
  counts
}

posterior_block_probs_unordered <- function(counts_unordered, alpha_dir) {
  K <- dim(counts_unordered)[1]
  Pi <- array(0, dim = c(K, K, 4))
  for (r in 1:K) for (s in r:K) {
    post <- alpha_dir + counts_unordered[r, s, ]
    Pi[r, s, ] <- post / sum(post)
  }
  Pi
}

expand_block_probs_ordered <- function(Pi_unordered) {
  # Build Pi_ordered[r,s,] meaning probabilities for (r->s, s->r) in that order.
  # For r>s we must swap categories 10 and 01.
  K <- dim(Pi_unordered)[1]
  Pi_ord <- array(0, dim = c(K, K, 4))
  for (r in 1:K) for (s in 1:K) {
    if (r <= s) {
      Pi_ord[r, s, ] <- Pi_unordered[r, s, ]
    } else {
      Pi_ord[r, s, ] <- Pi_unordered[s, r, c(1, 3, 2, 4)]  # swap 10 <-> 01
    }
  }
  Pi_ord
}

make_block_table <- function(Pi_unordered) {
  # One row per unordered block (r<=s)
  K <- dim(Pi_unordered)[1]
  rows <- list()
  for (r in 1:K) for (s in r:K) {
    p <- Pi_unordered[r, s, ]
    rows[[length(rows) + 1L]] <- data.frame(
      r = r, s = s,
      p00 = p[1], p10 = p[2], p01 = p[3], p11 = p[4],
      p_none   = p[1],
      p_oneway = p[2] + p[3],
      p_recip  = p[4],
      dir_bias = if (r < s) (p[2] - p[3]) else NA_real_,
      row.names = NULL
    )
  }
  do.call(rbind, rows)
}

# ----------------------------
# 1) Simulation (same model)
# ----------------------------
sim_directed_dyadstate <- function(z, Pi_unordered) {
  # Pi_unordered is K x K x 4 array for r<=s blocks with categories (00,10,01,11)
  # relative to (a -> b, b -> a), where a=min(r,s), b=max(r,s)
  n <- length(z)
  A <- matrix(0L, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    r <- z[i]; s <- z[j]
    a <- min(r, s); b <- max(r, s)
    probs <- Pi_unordered[a, b, ]
    state <- sample.int(4, 1, prob = probs) # 1:00, 2:10, 3:01, 4:11
    bits <- switch(state,
                   c(0L, 0L),
                   c(1L, 0L),
                   c(0L, 1L),
                   c(1L, 1L))
    if (r == a && s == b) {
      A[i, j] <- bits[1]; A[j, i] <- bits[2]
    } else {
      # swapped communities => swap interpretation of 10/01
      A[i, j] <- bits[2]; A[j, i] <- bits[1]
    }
  }
  diag(A) <- 0L
  A
}

# ----------------------------
# 2) Collapsed Gibbs sampler with trace
# ----------------------------
gibbs_dyadstate_sbm <- function(A, K,
                                alpha_dir = rep(1, 4),
                                alpha_part = 1,
                                n_iter = 5000, burn = 1000, thin = 10,
                                z_init = NULL, z_true = NULL,
                                verbose = TRUE, seed = NULL) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(A)
  diag(A) <- 0L
  
  # init labels
  z <- if (is.null(z_init)) sample.int(K, n, replace = TRUE) else as.integer(z_init)
  nk <- tabulate(z, K)
  
  # initialize unordered block counts
  counts <- block_counts_unordered(A, z, K)
  
  loglik_from_counts <- function(counts_arr) {
    out <- 0
    for (r in 1:K) for (s in r:K) {
      out <- out + log_m_dirichlet_multinomial(counts_arr[r, s, ], alpha_dir)
    }
    out
  }
  
  best_lp <- -Inf
  z_map <- z
  
  trace <- data.frame(iter = integer(0), logpost = numeric(0), ari = numeric(0))
  
  for (it in 1:n_iter) {
    
    for (i in 1:n) {
      old <- z[i]
      
      # remove i
      z[i] <- 0L
      nk[old] <- nk[old] - 1L
      
      idx_list <- lapply(1:K, function(l) which(z == l))
      
      # Counts between i and each group l, oriented as (i->j, j->i)
      cross_small <- lapply(1:K, function(l) {
        idx <- idx_list[[l]]
        if (!length(idx)) return(c(0, 0, 0, 0))
        dyad_counts_from_bits(A[i, idx], A[idx, i])
      })
      # If i's group index is larger than l, we need to swap 10/01 to match (min,max) orientation
      cross_large <- lapply(1:K, function(l) cross_small[[l]][c(1, 3, 2, 4)])
      
      # Within-group contribution if i were in group l
      within_cnt <- lapply(1:K, function(l) within_counts_i_group(i, idx_list[[l]], A))
      
      # subtract old contributions from counts
      for (l in 1:K) {
        if (l == old) {
          counts[old, old, ] <- counts[old, old, ] - within_cnt[[old]]
        } else if (old < l) {
          counts[old, l, ] <- counts[old, l, ] - cross_small[[l]]
        } else {
          counts[l, old, ] <- counts[l, old, ] - cross_large[[l]]
        }
      }
      
      # candidate weights for z[i]=k
      logw <- rep(-Inf, K)
      for (k in 1:K) {
        logprior_k <- log(nk[k] + alpha_part / K)
        delta <- 0
        
        for (l in 1:K) {
          if (l == k) {
            c_old <- counts[k, k, ]
            c_new <- c_old + within_cnt[[k]]
            delta <- delta + (log_m_dirichlet_multinomial(c_new, alpha_dir) -
                                log_m_dirichlet_multinomial(c_old, alpha_dir))
          } else if (k < l) {
            c_old <- counts[k, l, ]
            c_new <- c_old + cross_small[[l]]
            delta <- delta + (log_m_dirichlet_multinomial(c_new, alpha_dir) -
                                log_m_dirichlet_multinomial(c_old, alpha_dir))
          } else { # k > l
            c_old <- counts[l, k, ]
            c_new <- c_old + cross_large[[l]]
            delta <- delta + (log_m_dirichlet_multinomial(c_new, alpha_dir) -
                                log_m_dirichlet_multinomial(c_old, alpha_dir))
          }
        }
        logw[k] <- logprior_k + delta
      }
      
      w <- exp(logw - max(logw))
      new <- sample.int(K, 1, prob = w)
      
      # add i to new: update counts and nk
      for (l in 1:K) {
        if (l == new) {
          counts[new, new, ] <- counts[new, new, ] + within_cnt[[new]]
        } else if (new < l) {
          counts[new, l, ] <- counts[new, l, ] + cross_small[[l]]
        } else {
          counts[l, new, ] <- counts[l, new, ] + cross_large[[l]]
        }
      }
      
      z[i] <- new
      nk[new] <- nk[new] + 1L
    }
    
    # save trace after burn/thin
    if (it >= burn && ((it - burn) %% thin == 0)) {
      lp <- log_dirichlet_multinomial_prior(nk, alpha_part) + loglik_from_counts(counts)
      ari_it <- if (!is.null(z_true)) ari(z_true, z) else NA_real_
      trace <- rbind(trace, data.frame(iter = it, logpost = lp, ari = ari_it))
      if (lp > best_lp) { best_lp <- lp; z_map <- z }
    }
    
    if (verbose && (it %% 200 == 0)) message(sprintf("iter %d / %d", it, n_iter))
  }
  
  # Recompute counts & posterior means at MAP labeling
  counts_map <- block_counts_unordered(A, z_map, K)
  Pi_postmean_unordered <- posterior_block_probs_unordered(counts_map, alpha_dir)
  Pi_postmean_ordered <- expand_block_probs_ordered(Pi_postmean_unordered)
  
  list(
    z_map = z_map,
    logpost_map = best_lp,
    trace = trace,
    counts_map = counts_map,
    Pi_postmean_unordered = Pi_postmean_unordered,
    Pi_postmean_ordered = Pi_postmean_ordered
  )
}

# ----------------------------
# 3) Plotting helpers (all base R; saves PNGs)
# ----------------------------
plot_dyad_state_legend <- function(file = NULL, width = 7, height = 5, res = 160) {
  if (!is.null(file)) png(file, width = width, height = height, units = "in", res = res)
  op <- par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
  on.exit({ par(op); if (!is.null(file)) dev.off() }, add = TRUE)
  
  panels <- list(
    list(title = "00 : no ties", arrows = list()),
    list(title = "10 : i \u2192 j only", arrows = list(c(0.30, 0.50, 0.70, 0.50))),
    list(title = "01 : j \u2192 i only", arrows = list(c(0.70, 0.50, 0.30, 0.50))),
    list(title = "11 : mutual", arrows = list(c(0.30, 0.52, 0.70, 0.52),
                                              c(0.70, 0.48, 0.30, 0.48)))
  )
  
  for (p in panels) {
    plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
    points(c(0.3, 0.7), c(0.5, 0.5), pch = 21, bg = "white", cex = 3)
    text(c(0.3, 0.7), c(0.5, 0.5), labels = c("i", "j"), cex = 1.2)
    if (length(p$arrows)) {
      for (a in p$arrows) arrows(a[1], a[2], a[3], a[4], length = 0.12, lwd = 2)
    }
    title(p$title, cex.main = 1.1)
  }
}

plot_gibbs_flow <- function(file = NULL, width = 8, height = 4.3, res = 160) {
  if (!is.null(file)) png(file, width = width, height = height, units = "in", res = res)
  op <- par(mar = c(0.8, 0.8, 2.5, 0.8))
  on.exit({ par(op); if (!is.null(file)) dev.off() }, add = TRUE)
  
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  title("Collapsed Gibbs update (one node i)", cex.main = 1.15)
  
  boxes <- data.frame(
    x0 = rep(0.08, 4),
    y0 = c(0.78, 0.56, 0.34, 0.12),
    w  = rep(0.84, 4),
    h  = rep(0.14, 4),
    txt = c(
      "1) Remove i from its group\n   (update nk and block counts)",
      "2) Compute dyad-state counts\n   between i and each group",
      "3) For each candidate k:\n   prior + Dirichlet–multinomial count change",
      "4) Sample new group for i\n   then update counts"
    )
  )
  
  for (r in 1:nrow(boxes)) {
    rect(boxes$x0[r], boxes$y0[r], boxes$x0[r] + boxes$w[r], boxes$y0[r] + boxes$h[r], lwd = 1.6)
    text(boxes$x0[r] + boxes$w[r] / 2, boxes$y0[r] + boxes$h[r] / 2, boxes$txt[r], cex = 0.95)
    if (r < nrow(boxes)) {
      arrows(0.50, boxes$y0[r], 0.50, boxes$y0[r + 1] + boxes$h[r + 1], length = 0.10, lwd = 1.3)
    }
  }
}

plot_trace <- function(trace, file = NULL, width = 7.5, height = 5.5, res = 160) {
  if (!is.null(file)) png(file, width = width, height = height, units = "in", res = res)
  has_ari <- ("ari" %in% names(trace)) && any(!is.na(trace$ari))
  op <- par(mfrow = if (has_ari) c(2, 1) else c(1, 1), mar = c(4, 4, 2.5, 1))
  on.exit({ par(op); if (!is.null(file)) dev.off() }, add = TRUE)
  
  plot(trace$iter, trace$logpost, type = "l",
       xlab = "Iteration", ylab = "Log posterior", main = "Trace of collapsed log posterior")
  
  if (has_ari) {
    plot(trace$iter, trace$ari, type = "l", ylim = c(0, 1),
         xlab = "Iteration", ylab = "ARI", main = "Clustering accuracy trace (if z_true provided)")
  }
}

plot_adjacency_ordered <- function(A, z, file = NULL, width = 6.2, height = 6.2, res = 160) {
  if (!is.null(file)) png(file, width = width, height = height, units = "in", res = res)
  op <- par(mar = c(1, 1, 3, 1))
  on.exit({ par(op); if (!is.null(file)) dev.off() }, add = TRUE)
  
  ord <- order(z)
  Aord <- A[ord, ord]
  n <- nrow(Aord)
  
  image(1:n, 1:n, t(Aord[n:1, ]),
        col = c("white", "black"),
        axes = FALSE, xlab = "", ylab = "",
        main = "Adjacency matrix ordered by inferred communities")
  
  # draw block boundaries
  sizes <- tabulate(z[ord], max(z))
  cuts <- cumsum(sizes)
  for (c in cuts) {
    abline(v = c + 0.5, lwd = 1)
    abline(h = n - c + 0.5, lwd = 1)
  }
}

plot_block_heatmaps <- function(Pi_ordered, file = NULL, width = 8, height = 6.2, res = 160) {
  if (!is.null(file)) png(file, width = width, height = height, units = "in", res = res)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit({ par(op); if (!is.null(file)) dev.off() }, add = TRUE)
  
  K <- dim(Pi_ordered)[1]
  state_titles <- c("P(00): none", "P(10): r\u2192s only", "P(01): s\u2192r only", "P(11): mutual")
  pal <- gray.colors(200, start = 1, end = 0)  # white=0, black=1
  
  for (st in 1:4) {
    M <- Pi_ordered[, , st]
    image(1:K, 1:K, t(M[K:1, ]), zlim = c(0, 1), col = pal,
          axes = FALSE, xlab = "To community s", ylab = "From community r",
          main = state_titles[st])
    axis(1, at = 1:K, labels = 1:K)
    axis(2, at = 1:K, labels = K:1)
    
    if (K <= 6) {
      for (r in 1:K) for (s in 1:K) {
        text(s, K - r + 1, sprintf("%.2f", M[r, s]), cex = 0.95)
      }
    }
  }
}

plot_edgeprob_heatmap <- function(Pi_ordered, file = NULL, width = 6.8, height = 6, res = 160) {
  if (!is.null(file)) png(file, width = width, height = height, units = "in", res = res)
  op <- par(mar = c(4, 4, 3, 1))
  on.exit({ par(op); if (!is.null(file)) dev.off() }, add = TRUE)
  
  K <- dim(Pi_ordered)[1]
  # P(A_{ij}=1 | i in r, j in s) = P(10) + P(11) under ordered interpretation (r->s)
  Pedge <- Pi_ordered[, , 2] + Pi_ordered[, , 4]
  pal <- gray.colors(200, start = 1, end = 0)
  
  image(1:K, 1:K, t(Pedge[K:1, ]), zlim = c(0, 1), col = pal,
        axes = FALSE, xlab = "To community s", ylab = "From community r",
        main = "Directed edge probability: P(A_{ij}=1 | z_i=r, z_j=s)")
  axis(1, at = 1:K, labels = 1:K)
  axis(2, at = 1:K, labels = K:1)
  
  if (K <= 8) {
    for (r in 1:K) for (s in 1:K) {
      text(s, K - r + 1, sprintf("%.2f", Pedge[r, s]), cex = 0.9)
    }
  }
}

# ----------------------------
# 4) End-to-end runner: sim + fit + summary + images + CSV tables
# ----------------------------
demo_case5 <- function(output_dir = "case5_out",
                       seed = 6,
                       n = 120,
                       K_true = 2,
                       K_fit  = 2,
                       n_iter = 5000,
                       burn   = 1000,
                       thin   = 10) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(seed)
  
  # ---- simulate a directed dyad-state SBM ----
  z_true <- rep(1:K_true, each = n / K_true)
  
  Pi <- array(0, dim = c(K_true, K_true, 4))
  # within: more reciprocity
  Pi[1, 1, ] <- c(0.80, 0.06, 0.06, 0.08)
  Pi[2, 2, ] <- c(0.78, 0.07, 0.07, 0.08)
  # between: sparse and slightly asymmetric (1->2 more likely than 2->1)
  Pi[1, 2, ] <- c(0.93, 0.04, 0.02, 0.01)
  
  A <- sim_directed_dyadstate(z_true, Pi)
  
  # ---- fit ----
  fit <- gibbs_dyadstate_sbm(
    A, K = K_fit,
    alpha_dir = rep(1, 4),
    alpha_part = 1,
    n_iter = n_iter, burn = burn, thin = thin,
    z_true = z_true,
    verbose = FALSE,
    seed = seed
  )
  
  z_hat <- fit$z_map
  ari_val <- ari(z_true, z_hat)
  
  # ---- tables ----
  cluster_sizes <- data.frame(group = 1:K_fit, n = tabulate(z_hat, K_fit))
  block_table <- make_block_table(fit$Pi_postmean_unordered)
  conf <- table(true = z_true, inferred = z_hat)
  
  # save tables
  write.csv(cluster_sizes, file.path(output_dir, "cluster_sizes.csv"), row.names = FALSE)
  write.csv(block_table,   file.path(output_dir, "block_probs_unordered.csv"), row.names = FALSE)
  write.csv(as.data.frame(conf), file.path(output_dir, "confusion_matrix.csv"), row.names = FALSE)
  
  # ---- console summary ----
  cat("\n==================== Case 5+ Summary ====================\n")
  cat("Directed dyad-state SBM (collapsed Dirichlet–multinomial)\n")
  cat("n =", n, "  K_fit =", K_fit, "  seed =", seed, "\n")
  cat("ARI(z_true, z_hat) =", sprintf("%.3f", ari_val), "\n\n")
  
  cat("Inferred community sizes:\n")
  print(cluster_sizes)
  
  cat("\nPosterior mean dyad-state probabilities by unordered block (r<=s):\n")
  print(block_table)
  
  cat("\nConfusion matrix (true vs inferred):\n")
  print(conf)
  
  # ---- images ----
  plot_dyad_state_legend(file = file.path(output_dir, "dyad_states_legend.png"))
  plot_gibbs_flow(file = file.path(output_dir, "gibbs_flow.png"))
  plot_trace(fit$trace, file = file.path(output_dir, "trace.png"))
  plot_adjacency_ordered(A, z_hat, file = file.path(output_dir, "adjacency_ordered.png"))
  plot_block_heatmaps(fit$Pi_postmean_ordered, file = file.path(output_dir, "block_probs_states.png"))
  plot_edgeprob_heatmap(fit$Pi_postmean_ordered, file = file.path(output_dir, "block_edgeprob.png"))
  
  cat("\nSaved outputs to:", normalizePath(output_dir), "\n")
  cat("  - dyad_states_legend.png (what the 4 states mean)\n")
  cat("  - gibbs_flow.png (what one Gibbs update is doing)\n")
  cat("  - trace.png (log-posterior trace + ARI trace)\n")
  cat("  - adjacency_ordered.png (matrix ordered by inferred groups)\n")
  cat("  - block_probs_states.png (4 heatmaps: 00,10,01,11)\n")
  cat("  - block_edgeprob.png (heatmap of P(A_ij=1 | z_i=r,z_j=s))\n")
  cat("  - CSV tables: cluster_sizes.csv, block_probs_unordered.csv, confusion_matrix.csv\n")
  cat("==========================================================\n\n")
  
  invisible(list(
    A = A, z_true = z_true, z_hat = z_hat, fit = fit,
    tables = list(cluster_sizes = cluster_sizes, block_table = block_table, confusion = conf),
    output_dir = output_dir
  ))
}

# ============================================================
# What to run
# ============================================================
demo <- demo_case5(output_dir = "case5_out", seed = 6)
