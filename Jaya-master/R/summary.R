#' Summary & Tidy Methods for Jaya
#'
#' @param object A "jaya" or "jaya_multi" object.
#' @param ... Unused.
#' @export
summary.jaya <- function(object, ...) {
  cat("Jaya (single-objective) summary\n")
  cat(sprintf("Variables: %d | popSize: %d | maxiter: %d\n",
              attr(object, "n_var"), attr(object, "popSize"), attr(object, "maxiter")))
  cat(sprintf("Best value: %.10g\n", object$value))
  cat("Best parameters:\n")
  par <- setNames(object$par, paste0("x", seq_along(object$par)))
  print(par)
  cat(sprintf("Evaluations: %d | Runtime: %.3f s\n", object$n_eval, object$runtime_sec))
  invisible(object)
}

#' @export
summary.jaya_multi <- function(object, ...) {
  cat("Jaya (multi-objective) summary\n")
  cat(sprintf("Variables: %d | popSize: %d | maxiter: %d\n",
              attr(object, "n_var"), attr(object, "popSize"), attr(object, "maxiter")))
  cat("Pareto front (first 6 rows):\n")
  print(utils::head(object$Pareto_Front))
  cat(sprintf("Archive size: %d | Evaluations: %d | Runtime: %.3f s\n",
              nrow(object$Pareto_Front), object$n_eval, object$runtime_sec))
  invisible(object)
}

#' @export
glance.jaya <- function(object, ...) {
  data.frame(best = object$value,
             iterations = length(object$history),
             n_eval = object$n_eval,
             runtime_sec = object$runtime_sec)
}

#' @export
glance.jaya_multi <- function(object, ...) {
  data.frame(archive = nrow(object$Pareto_Front),
             iterations = length(object$history),
             n_eval = object$n_eval,
             runtime_sec = object$runtime_sec)
}

#' @export
tidy.jaya <- function(object, ...) {
  data.frame(iter = seq_along(object$history),
             best = as.numeric(object$history))
}

#' @export
tidy.jaya_multi <- function(object, ...) {
  pf <- object$Pareto_Front
  pf$solution_id <- seq_len(nrow(pf))
  pf
}