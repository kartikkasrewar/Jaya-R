#' Jaya Algorithm (Single-Objective)
#'
#' Implements the Jaya optimization algorithm for single-objective problems with
#' optional parallel evaluation (future/furrr), constraint handling, early stopping,
#' and tidy outputs for plotting/integration.
#'
#' @param fun Function to optimize. Takes a numeric vector and returns a scalar.
#' @param lower,upper Numeric vectors of bounds (length D).
#' @param n_var Integer, number of decision variables (defaults to length(lower)).
#' @param popSize Integer population size.
#' @param maxiter Integer maximum iterations.
#' @param opt "minimize" (default) or "maximize".
#' @param seed Optional integer RNG seed.
#' @param suggestions Optional matrix/data.frame of starting points (rows = candidates).
#' @param constraints Optional list of functions g_i(x) <= 0 for feasibility.
#' @param constraint_handling One of "reject","repair","penalty".
#' @param penalty_fun Optional function p(x) returning penalty >= 0 if using "penalty".
#' @param early_stopping Logical, enable early stop on small improvement.
#' @param tolerance Numeric improvement tolerance for early stopping.
#' @param patience Integer number of stagnant iterations allowed.
#' @param parallel Logical; if TRUE, evaluate population with future/furrr.
#' @param cores Integer worker count for parallel; default uses availableCores()-1.
#' @param verbose Logical; print progress.
#'
#' @return An object of class "jaya" with elements:
#'   \itemize{
#'     \item par Best parameter vector.
#'     \item value Best objective value.
#'     \item history Numeric vector of best values by iteration.
#'     \item best_path Matrix of best parameters by iteration.
#'     \item n_eval Total objective evaluations.
#'     \item runtime_sec Wall time in seconds.
#'   }
#' @export
jaya <- function(fun, lower, upper,
                 n_var = length(lower),
                 popSize = 50L, maxiter = 200L,
                 opt = c("minimize","maximize"),
                 seed = NULL, suggestions = NULL,
                 constraints = NULL,
                 constraint_handling = c("reject","repair","penalty"),
                 penalty_fun = NULL,
                 early_stopping = TRUE, tolerance = 1e-8, patience = 25L,
                 parallel = FALSE, cores = NULL,
                 verbose = FALSE) {

  opt <- match.arg(opt)
  constraint_handling <- match.arg(constraint_handling)
  stopifnot(length(lower) == length(upper), length(lower) == n_var)

  # RNG scope
  if (!is.null(seed)) {
    if (requireNamespace("withr", quietly = TRUE)) {
      withr::with_seed(seed, {}) # set once
    } else {
      set.seed(seed)
    }
  }

  D <- n_var
  lb <- as.numeric(lower); ub <- as.numeric(upper)
  clamp <- function(X) pmin(pmax(X, matrix(lb, nrow(X), D, byrow = TRUE)),
                            matrix(ub, nrow(X), D, byrow = TRUE))

  # Build initial population
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

  # Helpers ----
  feasible <- function(x) {
    if (is.null(constraints)) return(TRUE)
    all(vapply(constraints, function(g) isTRUE(g(x) <= 0), TRUE))
  }
  feas_mask <- function(M) {
    if (is.null(constraints)) return(rep(TRUE, nrow(M)))
    apply(M, 1L, feasible)
  }

  penalized_eval <- function(M) {
    vals <- eval_pop(M)
    if (constraint_handling == "penalty" && length(constraints)) {
      if (is.null(penalty_fun)) {
        # default linear penalty on positive constraint violation
        pen_fun <- function(x) {
          sum(vapply(constraints, function(g) max(0, g(x)), 0.0))
        }
      } else pen_fun <- penalty_fun
      pen <- apply(M, 1L, pen_fun)
      vals <- vals + pen
    }
    vals
  }

  # Parallel evaluation ----
  eval_pop <- function(M) {
  #   if (!parallel) {
  #     apply(M, 1L, fun)
  #   } else {
  #     if (!requireNamespace("future", quietly = TRUE) ||
  #         !requireNamespace("furrr", quietly = TRUE)) {
  #       stop("Install 'future' and 'furrr' for parallel evaluation.")
  #     }
  #     if (is.null(cores)) cores <- max(1L, future::availableCores() - 1L)
  #     old_plan <- future::plan()
  #     on.exit(future::plan(old_plan), add = TRUE)
  #     future::plan(future::multisession, workers = cores)
  #     v <- furrr::future_map_dbl(seq_len(nrow(M)), function(i) fun(M[i, ]))
  #     unlist(v, use.names = FALSE)
  #   }
  # }

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
  }

  # Evaluation helpers ----
  eval_pop <- function(M) {
    if (!parallel) {
      apply(M, 1L, fun)
    } else {
      future.apply::future_apply(M, 1L, fun)
    }
  }

  # Evaluate initial population
  t0 <- proc.time()[["elapsed"]]
  vals <- penalized_eval(X)
  n_eval <- nrow(X)

  # Feasibility rejection
  if (constraint_handling == "reject" && length(constraints)) {
    ok <- feas_mask(X)
    # For infeasible rows, re-sample within bounds
    while (any(!ok)) {
      X[!ok, ] <- matrix(runif(sum(!ok) * D, lb, ub), ncol = D)
      vals[!ok] <- penalized_eval(X[!ok, , drop = FALSE])
      n_eval <- n_eval + sum(!ok)
      ok <- feas_mask(X)
    }
  }

  best_hist <- numeric(maxiter)
  best_path <- matrix(NA_real_, nrow = maxiter, ncol = D)
  improve_age <- 0L

  # Utility to choose best/worst
  which_best <- function(v) if (opt == "minimize") which.min(v) else which.max(v)
  which_worst <- function(v) if (opt == "minimize") which.max(v) else which.min(v)

  b_idx <- which_best(vals); w_idx <- which_worst(vals)
  gbest <- X[b_idx, , drop = FALSE]
  gworst <- X[w_idx, , drop = FALSE]
  gbest_val <- vals[b_idx]

  for (it in seq_len(maxiter)) {
    r1 <- matrix(runif(popSize * D), nrow = popSize)
    r2 <- matrix(runif(popSize * D), nrow = popSize)
    absX <- abs(X)
    Xnew <- X + r1 * (matrix(gbest, popSize, D, byrow = TRUE) - absX) -
                  r2 * (matrix(gworst, popSize, D, byrow = TRUE) - absX)
    # clamp
    Xnew <- clamp(Xnew)

    # Optional "repair": nothing beyond clamp; if constraints exist, try accept/reject at candidate level
    vnew <- penalized_eval(Xnew); n_eval <- n_eval + nrow(Xnew)

    if (length(constraints) && constraint_handling == "reject") {
      ok_new <- feas_mask(Xnew)
      # keep old where new infeasible
      Xnew[!ok_new, ] <- X[!ok_new, ]
      vnew[!ok_new] <- vals[!ok_new]
    }

    # Greedy survival
    better <- if (opt == "minimize") vnew <= vals else vnew >= vals
    X[better, ] <- Xnew[better, ]
    vals[better] <- vnew[better]

    # Global best & worst
    b_idx <- which_best(vals); w_idx <- which_worst(vals)
    gbest <- X[b_idx, , drop = FALSE]
    gworst <- X[w_idx, , drop = FALSE]
    gbest_val_new <- vals[b_idx]

    best_hist[it] <- gbest_val_new
    best_path[it, ] <- gbest

    if (verbose) {
      message(sprintf("Iter %03d | best = %.8g", it, gbest_val_new))
    }

    # early stopping
    if (early_stopping) {
      if (it == 1L || abs(gbest_val - gbest_val_new) > tolerance) {
        improve_age <- 0L
        gbest_val <- gbest_val_new
      } else {
        improve_age <- improve_age + 1L
        if (improve_age >= patience) break
      }
    } else {
      gbest_val <- gbest_val_new
    }
  }

  runtime_sec <- proc.time()[["elapsed"]] - t0
  out <- list(
    par = as.numeric(gbest),
    value = as.numeric(gbest_val),
    history = best_hist[best_hist != 0],
    best_path = best_path[!is.na(best_path[,1]), , drop = FALSE],
    n_eval = n_eval,
    runtime_sec = runtime_sec,
    call = match.call()
  )
  attr(out, "popSize") <- popSize
  attr(out, "maxiter") <- maxiter
  attr(out, "n_var") <- n_var
  attr(out, "opt") <- opt
  attr(out, "lower") <- lower
  attr(out, "upper") <- upper
  class(out) <- "jaya"
  out
}

#' @rdname jaya
#' @export
jaya_min <- function(fun, lower, upper, ...) {
  jaya(fun, lower, upper, opt = "minimize", ...)
}

#' @rdname jaya
#' @export
jaya_max <- function(fun, lower, upper, ...) {
  jaya(fun, lower, upper, opt = "maximize", ...)
}
