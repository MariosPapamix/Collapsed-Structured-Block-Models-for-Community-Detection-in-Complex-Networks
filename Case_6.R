# ============================================================
# Case 6: Signed SBM (collapsed categorical Dirichlet–multinomial)
# + summary tables + plots (adj heatmap, block probs, trace, network, process)
# ============================================================

# ---- 0) Packages (auto-install if needed) ----
pkgs <- c("ggplot2", "dplyr", "tidyr", "igraph", "DiagrammeR", "mclust")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1) Math helpers ----
log_m_dirichlet_multinomial <- function(counts, alpha) {
  counts <- as.numeric(counts)
  alpha  <- as.numeric(alpha)
  stopifnot(length(counts) == length(alpha))
  lgamma(sum(alpha)) - lgamma(sum(alpha) + sum(counts)) +
    sum(lgamma(alpha + counts) - lgamma(alpha))
}

# Prior on labels with mixing proportions integrated out:
# z_i | pi ~ Categorical(pi), pi ~ Dirichlet(alpha_part/K,...)
log_dirichlet_multinomial_prior <- function(nk, alpha_part) {
  nk <- as.numeric(nk)
  K  <- length(nk)
  a  <- alpha_part / K
  lgamma(alpha_part) - lgamma(alpha_part + sum(nk)) +
    sum(lgamma(a + nk) - lgamma(a))
}

ari <- function(a, b) mclust::adjustedRandIndex(a, b)

edge_to_cat <- function(v) {           # map {-1,0,+1} -> {3,1,2}
  ifelse(v == 0, 1L, ifelse(v == 1, 2L, 3L))
}
cat_counts <- function(v) tabulate(edge_to_cat(v), nbins = 3) # (0, +, -)

# ---- 2) Simulator ----
sim_signed_sbm <- function(z, Prob_in, Prob_out) {
  # Prob_* are length-3 vectors for (0, +, -) summing to 1
  n <- length(z)
  Y <- matrix(0L, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    probs <- if (z[i] == z[j]) Prob_in else Prob_out
    state <- sample.int(3, 1, prob = probs)
    val   <- c(0L, 1L, -1L)[state]
    Y[i, j] <- val
    Y[j, i] <- val
  }
  diag(Y) <- 0L
  Y
}

# ---- 3) Block-pair counts + posterior means ----
counts_blockpairs <- function(Y, z, K) {
  n <- nrow(Y)
  ij <- which(upper.tri(Y), arr.ind = TRUE)
  r  <- z[ij[, 1]]
  s  <- z[ij[, 2]]
  a  <- pmin(r, s)
  b  <- pmax(r, s)
  cat <- edge_to_cat(Y[upper.tri(Y)])
  # R's linear indexing for array(dim=c(K,K,3)): a + (b-1)*K + (cat-1)*K*K
  idx <- a + (b - 1) * K + (cat - 1) * K * K
  tab <- tabulate(idx, nbins = K * K * 3)
  array(tab, dim = c(K, K, 3))
}

posterior_theta_mean <- function(counts, alpha_edge) {
  K <- dim(counts)[1]
  theta <- array(NA_real_, dim = dim(counts))
  for (r in 1:K) for (s in r:K) {
    a <- counts[r, s, ] + alpha_edge
    theta[r, s, ] <- a / sum(a)
  }
  theta
}

theta_df_from_upper <- function(theta_upper, mirror = TRUE) {
  K <- dim(theta_upper)[1]
  cats <- c("absent (0)", "positive (+1)", "negative (-1)")
  out <- list()
  t <- 1
  for (r in 1:K) for (s in r:K) {
    probs <- theta_upper[r, s, ]
    for (c in 1:3) {
      out[[t]] <- data.frame(r = r, s = s, category = cats[c], prob = probs[c])
      t <- t + 1
    }
    if (mirror && s != r) {
      for (c in 1:3) {
        out[[t]] <- data.frame(r = s, s = r, category = cats[c], prob = probs[c])
        t <- t + 1
      }
    }
  }
  dplyr::bind_rows(out)
}

# ---- 4) Collapsed Gibbs sampler (returns MAP + trace) ----
gibbs_signed_sbm_collapsed <- function(
    Y, K,
    alpha_edge = c(1, 1, 1),
    alpha_part = 1,
    n_iter = 5000, burn = 1000, thin = 10,
    z_init = NULL,
    store_chain = FALSE,
    verbose = TRUE
) {
  stopifnot(is.matrix(Y), nrow(Y) == ncol(Y))
  n <- nrow(Y)
  diag(Y) <- 0L
  
  z <- if (is.null(z_init)) sample.int(K, n, replace = TRUE) else as.integer(z_init)
  nk <- tabulate(z, nbins = K)
  
  # counts[r,s,cat] for r<=s, cat in (0,+,-) -> 1:3
  counts <- counts_blockpairs(Y, z, K)
  
  log_m_cat <- function(cvec) log_m_dirichlet_multinomial(cvec, alpha_edge)
  
  loglik_from_counts <- function(counts_arr) {
    out <- 0
    for (r in 1:K) for (s in r:K) out <- out + log_m_cat(counts_arr[r, s, ])
    out
  }
  
  # how many saved points?
  save_its <- seq.int(from = burn, to = n_iter, by = thin)
  save_its <- save_its[save_its >= burn]
  n_save <- length(save_its)
  
  logpost_trace <- rep(NA_real_, n_save)
  z_samples <- if (store_chain) matrix(NA_integer_, nrow = n_save, ncol = n) else NULL
  save_idx <- 0
  
  best_lp <- -Inf
  z_map <- z
  
  for (it in 1:n_iter) {
    
    for (i in 1:n) {
      old <- z[i]
      
      # counts of (0,+,-) between node i and each group l (excluding i itself)
      grp_counts <- matrix(0L, nrow = K, ncol = 3)
      for (l in 1:K) {
        idx <- which(z == l)
        if (l == old) idx <- idx[idx != i]
        if (length(idx)) grp_counts[l, ] <- cat_counts(Y[i, idx])
      }
      
      # remove i's contribution from counts
      for (l in 1:K) {
        if (sum(grp_counts[l, ]) == 0) next
        if (l == old) {
          counts[old, old, ] <- counts[old, old, ] - grp_counts[l, ]
        } else {
          a <- min(old, l); b <- max(old, l)
          counts[a, b, ] <- counts[a, b, ] - grp_counts[l, ]
        }
      }
      nk[old] <- nk[old] - 1L
      z[i] <- 0L
      
      # candidate weights for assigning i to k
      logw <- rep(-Inf, K)
      for (k in 1:K) {
        logprior_k <- log(nk[k] + alpha_part / K)
        
        delta <- 0
        for (l in 1:K) {
          if (sum(grp_counts[l, ]) == 0) next
          
          if (l == k) {
            c_old <- counts[k, k, ]
            c_new <- c_old + grp_counts[k, ]
            delta <- delta + (log_m_cat(c_new) - log_m_cat(c_old))
          } else {
            a <- min(k, l); b <- max(k, l)
            c_old <- counts[a, b, ]
            c_new <- c_old + grp_counts[l, ]
            delta <- delta + (log_m_cat(c_new) - log_m_cat(c_old))
          }
        }
        
        logw[k] <- logprior_k + delta
      }
      
      # sample new label
      w <- exp(logw - max(logw))
      new <- sample.int(K, 1, prob = w)
      
      # add i's contribution back under new label
      for (l in 1:K) {
        if (sum(grp_counts[l, ]) == 0) next
        if (l == new) {
          counts[new, new, ] <- counts[new, new, ] + grp_counts[l, ]
        } else {
          a <- min(new, l); b <- max(new, l)
          counts[a, b, ] <- counts[a, b, ] + grp_counts[l, ]
        }
      }
      nk[new] <- nk[new] + 1L
      z[i] <- new
    }
    
    # save trace / MAP on thinned iterations after burn
    if (it >= burn && ((it - burn) %% thin == 0)) {
      save_idx <- save_idx + 1L
      lp <- log_dirichlet_multinomial_prior(nk, alpha_part) + loglik_from_counts(counts)
      logpost_trace[save_idx] <- lp
      if (store_chain) z_samples[save_idx, ] <- z
      if (lp > best_lp) { best_lp <- lp; z_map <- z }
    }
    
    if (verbose && (it %% 200 == 0)) message(sprintf("iter %d / %d", it, n_iter))
  }
  
  list(
    z_map = z_map,
    z_last = z,
    logpost_map = best_lp,
    logpost_trace = logpost_trace,
    save_iters = save_its,
    z_samples = z_samples
  )
}

# ---- 5) Plot helpers ----
plot_signed_adj_heatmap <- function(Y, z, title = "Signed adjacency (ordered by z)") {
  n <- nrow(Y)
  ord <- order(z)
  Y_ord <- Y[ord, ord]
  
  df <- as.data.frame(as.table(Y_ord))
  colnames(df) <- c("row", "col", "value")
  df$row <- as.integer(df$row)
  df$col <- as.integer(df$col)
  df$value <- as.integer(as.character(df$value))
  
  z_ord <- z[ord]
  cuts <- which(diff(z_ord) != 0)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row,
                                        fill = factor(value, levels = c(-1, 0, 1)))) +
    ggplot2::geom_raster() +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = c(`-1` = "firebrick3", `0` = "grey95", `1` = "steelblue3"),
      name   = "edge",
      labels = c("-1 (negative)", "0 (absent)", "+1 (positive)")
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  
  if (length(cuts)) {
    for (c in cuts) {
      p <- p +
        ggplot2::geom_hline(yintercept = c + 0.5, linewidth = 0.3) +
        ggplot2::geom_vline(xintercept = c + 0.5, linewidth = 0.3)
    }
  }
  p
}

plot_block_prob_tiles <- function(theta_df, title = "Posterior mean block-pair edge-type probabilities") {
  ggplot2::ggplot(theta_df, ggplot2::aes(x = s, y = r, fill = prob)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~ category) +
    ggplot2::scale_x_continuous(breaks = unique(theta_df$s)) +
    ggplot2::scale_y_continuous(breaks = unique(theta_df$r)) +
    ggplot2::labs(title = title, x = "community s", y = "community r", fill = "prob") +
    ggplot2::theme_minimal(base_size = 11)
}

plot_logpost_trace <- function(save_iters, logpost_trace, title = "Collapsed Gibbs trace (log posterior)") {
  df <- data.frame(iter = save_iters, logpost = logpost_trace)
  ggplot2::ggplot(df, ggplot2::aes(x = iter, y = logpost)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::labs(title = title, x = "iteration", y = "log posterior") +
    ggplot2::theme_minimal(base_size = 11)
}

plot_signed_network <- function(Y, z, seed = 1, main = "Signed network (nonzero edges only)") {
  set.seed(seed)
  n <- nrow(Y)
  
  ij <- which(upper.tri(Y), arr.ind = TRUE)
  sign <- Y[upper.tri(Y)]
  keep <- sign != 0
  ij <- ij[keep, , drop = FALSE]
  sign <- sign[keep]
  
  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  if (nrow(ij) > 0) g <- igraph::add_edges(g, as.vector(t(ij)))
  igraph::E(g)$sign <- sign
  igraph::V(g)$comm <- z
  
  K <- length(unique(z))
  pal <- grDevices::hcl.colors(max(K, 3), "Set2")
  vcol <- pal[z]
  
  ecol <- ifelse(igraph::E(g)$sign == 1, "steelblue3", "firebrick3")
  
  lay <- igraph::layout_with_fr(g)
  plot(
    g, layout = lay,
    vertex.color = vcol,
    vertex.size = 6,
    vertex.label = NA,
    edge.color = ecol,
    edge.width = 1,
    main = main
  )
  invisible(g)
}

process_diagram <- function() {
  DiagrammeR::grViz("
    digraph signedSBM {
      graph [rankdir = LR]
      node [shape=box, style=filled, color=gray30, fillcolor=gray95, fontname=Helvetica]
      edge [color=gray40]

      PriorZ     [label='Dirichlet prior on π', shape=ellipse, fillcolor=white]
      PriorTheta [label='Dirichlet prior on θ_rs', shape=ellipse, fillcolor=white]
      Z     [label='Latent labels z_i (communities)']
      Theta [label='Block-pair tendencies θ_rs\\n(p0, p+, p−)']
      Y     [label='Observed edges Y_ij ∈ {0,+1,−1}']
      Gibbs [label='Collapsed Gibbs sampler\\n(update one z_i at a time)']
      Out   [label='Outputs\\nMAP z, posterior θ, diagnostics', shape=box3d, fillcolor=white]

      PriorZ -> Z
      PriorTheta -> Theta
      Z -> Y
      Theta -> Y
      Y -> Gibbs
      PriorZ -> Gibbs
      PriorTheta -> Gibbs
      Gibbs -> Out
    }
  ")
}

# ---- 6) One-shot runner: simulate + infer + summary + visuals ----
case6_signed_sbm_report <- function(
    n = 150,
    K_true = 3,
    Prob_in  = c(0.70, 0.25, 0.05),  # (0, +, -) within-block
    Prob_out = c(0.80, 0.05, 0.15),  # (0, +, -) between-block
    K_fit = 3,
    alpha_edge = c(1, 1, 1),
    alpha_part = 1,
    n_iter = 5000, burn = 1000, thin = 10,
    seed = 7,
    verbose = FALSE,
    save_dir = NULL
) {
  set.seed(seed)
  
  # --- truth
  stopifnot(n %% K_true == 0)
  z_true <- rep(1:K_true, each = n / K_true)
  
  # --- simulate
  Y <- sim_signed_sbm(z_true, Prob_in, Prob_out)
  
  # --- fit
  fit <- gibbs_signed_sbm_collapsed(
    Y, K = K_fit,
    alpha_edge = alpha_edge,
    alpha_part = alpha_part,
    n_iter = n_iter, burn = burn, thin = thin,
    verbose = verbose
  )
  z_hat <- fit$z_map
  ari_val <- ari(z_true, z_hat)
  
  # --- global edge composition (upper triangle only)
  vals <- Y[upper.tri(Y)]
  edge_comp <- table(factor(vals, levels = c(-1, 0, 1)))
  edge_comp_df <- data.frame(
    sign = c("negative (-1)", "absent (0)", "positive (+1)"),
    count = as.integer(edge_comp),
    prop = as.numeric(edge_comp) / sum(edge_comp)
  )
  
  # --- posterior mean block-pair probabilities under MAP z
  counts_hat <- counts_blockpairs(Y, z_hat, K_fit)
  theta_hat  <- posterior_theta_mean(counts_hat, alpha_edge)
  theta_df   <- theta_df_from_upper(theta_hat, mirror = TRUE)
  
  # --- key tables
  tab <- table(True = z_true, Estimated = z_hat)
  sizes_hat <- as.data.frame(table(z_hat))
  colnames(sizes_hat) <- c("community", "size")
  
  # --- narrative summary (printed)
  cat("\n================= CASE 6 SUMMARY =================\n")
  cat(sprintf("n = %d nodes | K_true = %d | K_fit = %d\n", n, K_true, K_fit))
  cat(sprintf("ARI (truth vs MAP labels) = %.3f\n\n", ari_val))
  
  cat("Edge-type composition (overall):\n")
  print(edge_comp_df, row.names = FALSE)
  
  cat("\nCommunity sizes (MAP):\n")
  print(sizes_hat, row.names = FALSE)
  
  cat("\nContingency table (True x Estimated):\n")
  print(tab)
  
  # Small intuitive statement comparing within vs between, if K_fit >= 2
  within_pos <- mean(sapply(1:K_fit, function(k) theta_hat[k, k, 2]))
  within_neg <- mean(sapply(1:K_fit, function(k) theta_hat[k, k, 3]))
  if (K_fit >= 2) {
    pairs <- which(upper.tri(matrix(1, K_fit, K_fit), diag = FALSE), arr.ind = TRUE)
    between_pos <- mean(apply(pairs, 1, function(rc) theta_hat[rc[1], rc[2], 2]))
    between_neg <- mean(apply(pairs, 1, function(rc) theta_hat[rc[1], rc[2], 3]))
    cat("\nInterpretation (posterior means under MAP z):\n")
    cat(sprintf("  Avg within-block:   P(+)=%.3f, P(-)=%.3f\n", within_pos, within_neg))
    cat(sprintf("  Avg between-block:  P(+)=%.3f, P(-)=%.3f\n", between_pos, between_neg))
  }
  
  # --- plots
  p_true <- plot_signed_adj_heatmap(Y, z_true, title = "Signed adjacency (ordered by TRUE blocks)")
  p_hat  <- plot_signed_adj_heatmap(Y, z_hat,  title = "Signed adjacency (ordered by MAP inferred blocks)")
  p_theta <- plot_block_prob_tiles(theta_df, title = "Posterior mean P(edge type | block pair)")
  p_trace <- plot_logpost_trace(fit$save_iters, fit$logpost_trace)
  
  # --- optional saving
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "adj_true.png"), p_true, width = 6, height = 5, dpi = 200)
    ggplot2::ggsave(file.path(save_dir, "adj_map.png"),  p_hat,  width = 6, height = 5, dpi = 200)
    ggplot2::ggsave(file.path(save_dir, "block_probs.png"), p_theta, width = 7, height = 4, dpi = 200)
    ggplot2::ggsave(file.path(save_dir, "logpost_trace.png"), p_trace, width = 7, height = 3.5, dpi = 200)
    
    png(file.path(save_dir, "network.png"), width = 900, height = 700)
    plot_signed_network(Y, z_hat, main = "Signed network (MAP communities)")
    dev.off()
  }
  
  list(
    Y = Y,
    z_true = z_true,
    fit = fit,
    z_hat = z_hat,
    ari = ari_val,
    tables = list(edge_comp = edge_comp_df, contingency = tab, sizes = sizes_hat, theta = theta_df),
    plots = list(adj_true = p_true, adj_map = p_hat, block_probs = p_theta, trace = p_trace),
    diagram = process_diagram()
  )
}

# ============================================================
# WHAT TO RUN (this produces summary + plots + tables)
# ============================================================
out <- case6_signed_sbm_report(
  n = 150, K_true = 3,
  Prob_in  = c(0.70, 0.25, 0.05),
  Prob_out = c(0.80, 0.05, 0.15),
  K_fit = 3,
  alpha_edge = c(1, 1, 1),
  alpha_part = 1,
  n_iter = 5000, burn = 1000, thin = 10,
  seed = 7,
  verbose = FALSE,
  save_dir = "case6_outputs"  # set NULL to not save PNGs
)

# View plots (RStudio / Rmd will display them)
out$plots$adj_true
out$plots$adj_map
out$plots$block_probs
out$plots$trace

# Signed network image (plotted via igraph)
plot_signed_network(out$Y, out$z_hat, main = "Signed network (MAP communities)")

# Process diagram (generative + inference pipeline)
out$diagram

# View key tables
out$tables$contingency
head(out$tables$theta)
