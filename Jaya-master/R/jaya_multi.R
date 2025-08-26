
#' Jaya Algorithm (Multi-Objective)
#'
#' Multi-objective variant with non-dominated sorting, crowding distance,
#' elitist archive, ε-dominance (optional), parallel evaluation (future.apply),
#' and tidy outputs.
#'
#' @param objectives List of objective functions f_k(x) returning scalars.
#' @param lower,upper Numeric vectors of bounds.
#' @param popSize Population size.
#' @param maxiter Maximum iterations.
#' @param n_var Number of variables (defaults to length(lower)).
#' @param seed Optional RNG seed.
#' @param suggestions Optional starting points (rows = candidates).
#' @param constraints Optional list of g_i(x) <= 0.
#' @param epsilon Optional numeric for ε-dominance archiving (per objective).
#' @param archive_size Max size of Pareto archive (default = popSize).
#' @param parallel Logical; future.apply backend if TRUE.
#' @param cores Integer worker count if parallel.
#' @param verbose Logical.
#'
#' @return An object of class "jaya_multi" with elements:
#'   \itemize{
#'     \item Pareto_Front data.frame of nondominated objective values.
#'     \item Solutions data.frame of decision variables and objectives.
#'     \item history Optional list with per-iteration fronts sizes.
#'     \item n_eval Total evaluations.
#'     \item runtime_sec Wall time seconds.
#'   }
#' @export
jaya_multi <- function(objectives, lower, upper,
                       popSize = 100L, maxiter = 200L,
                       n_var = length(lower), seed = NULL,
                       suggestions = NULL,
                       constraints = NULL,
                       epsilon = NULL, archive_size = popSize,
                       parallel = FALSE, cores = NULL,
                       verbose = FALSE) {

  if (!is.list(objectives) || length(objectives) < 2L)
    stop("'objectives' must be a list of at least 2 functions.")

  if (!is.null(seed)) {
    if (requireNamespace("withr", quietly = TRUE)) {
      withr::with_seed(seed, {})
    } else set.seed(seed)
  }

  D <- n_var; M <- length(objectives)
  lb <- as.numeric(lower); ub <- as.numeric(upper)
  clamp <- function(X) pmin(pmax(X, matrix(lb, nrow(X), D, byrow = TRUE)),
                            matrix(ub, nrow(X), D, byrow = TRUE))

  # Initial population
  n_sug <- if (!is.null(suggestions)) nrow(as.data.frame(suggestions)) else 0L
  if (n_sug > 0L) {
    X <- as.matrix(suggestions)
    if (ncol(X) != D) stop("suggestions must have D columns")
    if (nrow(X) > popSize) X <- X[seq_len(popSize), , drop = FALSE]
  } else {
    X <- matrix(runif(popSize * D, lb, ub), nrow = popSize, ncol = D)
  }
  if (nrow(X) < popSize) {
    X2 <- matrix(runif((popSize - nrow(X)) * D, lb, ub), ncol = D)
    X <- rbind(X, X2)
  }
  X <- clamp(X)

  # Constraints
  feasible <- function(x) {
    if (is.null(constraints)) return(TRUE)
    all(vapply(constraints, function(g) isTRUE(g(x) <= 0), TRUE))
  }
  feas_mask <- function(M) {
    if (is.null(constraints)) return(rep(TRUE, nrow(M)))
    apply(M, 1L, feasible)
  }

  # Prepare parallel plan once (if requested) ----
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Install 'future' and 'future.apply' for parallel evaluation.")
    }
    if (is.null(cores)) cores <- max(1L, future::availableCores() - 1L)
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
  }

  # Eval objectives matrix (N x M)
  eval_obj <- function(Mt) {
    if (!parallel) {
      t(apply(Mt, 1L, function(row) vapply(objectives, function(f) f(row), 0.0)))
    } else {
      t(future.apply::future_apply(Mt, 1L, function(row) {
        vapply(objectives, function(f) f(row), 0.0)
      }))
    }
  }

  # MO helpers ----
  dominates <- function(a, b) {
    all(a <= b) && any(a < b)
  }
  fast_nondom_sort <- function(F) {
    N <- nrow(F)
    S <- vector("list", N)
    n <- integer(N)
    rank <- integer(N)
    fronts <- list()
    F1 <- integer(0)
    for (p in seq_len(N)) {
      S[[p]] <- integer(0)
      n[p] <- 0L
      for (q in seq_len(N)) {
        if (p == q) next
        if (dominates(F[p, ], F[q, ])) S[[p]] <- c(S[[p]], q)
        else if (dominates(F[q, ], F[p, ])) n[p] <- n[p] + 1L
      }
      if (n[p] == 0L) {
        rank[p] <- 1L
        F1 <- c(F1, p)
      }
    }
    i <- 1L; fronts[[1]] <- F1
    while (length(fronts[[i]]) > 0L) {
      Q <- integer(0)
      for (p in fronts[[i]]) {
        for (q in S[[p]]) {
          n[q] <- n[q] - 1L
          if (n[q] == 0L) {
            rank[q] <- i + 1L
            Q <- c(Q, q)
          }
        }
      }
      i <- i + 1L
      fronts[[i]] <- Q
    }
    fronts[length(fronts)] <- NULL
    list(fronts = fronts, rank = rank)
  }
  crowding_distance <- function(F) {
    n <- nrow(F); m <- ncol(F)
    cd <- rep(0, n)
    for (k in seq_len(m)) {
      ord <- order(F[, k], decreasing = FALSE)
      cd[ord[c(1, n)]] <- Inf
      rng <- max(F[, k]) - min(F[, k])
      if (rng == 0) next
      for (i in 2:(n-1)) {
        cd[ord[i]] <- cd[ord[i]] + (F[ord[i+1], k] - F[ord[i-1], k]) / rng
      }
    }
    cd
  }
  epsilon_filter <- function(F, eps) {
    if (is.null(eps)) return(seq_len(nrow(F)))
    key <- apply(F, 2, function(col) floor(col / eps))
    idx <- !duplicated(as.data.frame(key))
    which(idx)
  }

  # Evaluate initial
  t0 <- proc.time()[["elapsed"]]
  F <- eval_obj(X); n_eval <- nrow(X)

  # Feasibility (reject infeasible by resampling until feasible or limit tries)
  if (length(constraints)) {
    ok <- feas_mask(X)
    tries <- 0L
    while (any(!ok) && tries < 10L) {
      X[!ok, ] <- matrix(runif(sum(!ok) * D, lb, ub), ncol = D)
      F[!ok, ] <- eval_obj(X[!ok, , drop = FALSE])
      n_eval <- n_eval + sum(!ok)
      ok <- feas_mask(X)
      tries <- tries + 1L
    }
  }

  # Archive (elitist non-dominated)
  fronts <- fast_nondom_sort(F)$fronts
  A <- X[fronts[[1]], , drop = FALSE]
  FA <- F[fronts[[1]], , drop = FALSE]

  # Reduce archive if needed
  trim_archive <- function(A, FA) {
    if (is.null(epsilon) && nrow(A) <= archive_size) return(list(A=A, FA=FA))
    idx <- seq_len(nrow(A))
    if (!is.null(epsilon)) {
      idx <- epsilon_filter(FA, epsilon)
    }
    if (length(idx) > archive_size) {
      # trim by crowding distance
      cd <- crowding_distance(FA[idx, , drop = FALSE])
      ord <- order(cd, decreasing = TRUE) # keep most diverse
      keep <- idx[ord][seq_len(archive_size)]
      list(A = A[keep, , drop = FALSE], FA = FA[keep, , drop = FALSE])
    } else {
      list(A = A[idx, , drop = FALSE], FA = FA[idx, , drop = FALSE])
    }
  }
  tmp <- trim_archive(A, FA); A <- tmp$A; FA <- tmp$FA

  # MO Jaya loop
  history <- integer(maxiter)
  for (it in seq_len(maxiter)) {

    # Normalize objectives to [0,1] for guidance selection
    mins <- apply(FA, 2L, min); maxs <- apply(FA, 2L, max)
    rng <- pmax(1e-12, maxs - mins)
    FAn <- sweep(FA, 2L, mins, "-"); FAn <- sweep(FAn, 2L, rng, "/")
    ideal <- rep(0, ncol(FA)); nadir <- rep(1, ncol(FA))
    # Select "best" (closest to ideal) and "worst" (closest to nadir)
    d_ideal <- sqrt(rowSums((FAn - matrix(ideal, nrow(FAn), ncol(FAn), TRUE))^2))
    d_nadir <- sqrt(rowSums((FAn - matrix(nadir, nrow(FAn), ncol(FAn), TRUE))^2))
    b_idx <- which.min(d_ideal); w_idx <- which.min(d_nadir)
    gbest <- A[b_idx, , drop = FALSE]
    gworst <- A[w_idx, , drop = FALSE]

    # Jaya update
    r1 <- matrix(runif(popSize * D), nrow = popSize)
    r2 <- matrix(runif(popSize * D), nrow = popSize)
    absX <- abs(X)
    Xnew <- X + r1 * (matrix(gbest, popSize, D, byrow = TRUE) - absX) -
      r2 * (matrix(gworst, popSize, D, byrow = TRUE) - absX)
    Xnew <- clamp(Xnew)

    # Evaluate new
    Fnew <- eval_obj(Xnew); n_eval <- n_eval + nrow(Xnew)

    # Greedy survival by Pareto dominance at pair level (keep better or equal)
    keep <- logical(popSize)
    for (i in seq_len(popSize)) {
      if (dominates(Fnew[i, ], F[i, ]) || all(Fnew[i, ] == F[i, ])) {
        keep[i] <- TRUE
      } else if (dominates(F[i, ], Fnew[i, ])) {
        keep[i] <- FALSE
      } else {
        # Neither dominates -> tie-break by distance to current archive ideal
        fn <- (Fnew[i, ] - mins) / rng
        f0 <- (F[i, ] - mins) / rng
        keep[i] <- (sum((fn - ideal)^2) <= sum((f0 - ideal)^2))
      }
    }
    X[keep, ] <- Xnew[keep, ]
    F[keep, ] <- Fnew[keep, ]

    # Update archive with current population
    AX <- rbind(A, X); AF <- rbind(FA, F)
    nd <- fast_nondom_sort(AF)
    first <- nd$fronts[[1]]
    A <- AX[first, , drop = FALSE]; FA <- AF[first, , drop = FALSE]
    tmp <- trim_archive(A, FA); A <- tmp$A; FA <- tmp$FA

    history[it] <- nrow(A)
    if (verbose) message(sprintf("Iter %03d | archive size = %d", it, nrow(A)))
  }

  runtime_sec <- proc.time()[["elapsed"]] - t0

  col_par <- paste0("x", seq_len(D))
  col_obj <- paste0("Obj", seq_len(M))
  Pareto_Front <- as.data.frame(FA); names(Pareto_Front) <- col_obj
  Solutions <- cbind(as.data.frame(A), as.data.frame(FA))
  names(Solutions) <- c(col_par, col_obj)

  out <- list(
    Pareto_Front = Pareto_Front,
    Solutions = Solutions,
    history = history[history != 0],
    n_eval = n_eval,
    runtime_sec = runtime_sec,
    call = match.call()
  )
  attr(out, "popSize") <- popSize
  attr(out, "maxiter") <- maxiter
  attr(out, "n_var") <- n_var
  attr(out, "lower") <- lower
  attr(out, "upper") <- upper
  attr(out, "objectives") <- objectives
  class(out) <- "jaya_multi"
  out
}
