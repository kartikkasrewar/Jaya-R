#' Plot & Autoplot Methods for Jaya Results
#'
#' Base and ggplot2 visualizations for single- and multi-objective Jaya runs.
#' ggplot2 will be used when available; otherwise base plots are used.
#'
#' @param x An object of class "jaya" or "jaya_multi".
#' @param ... Additional arguments.
#' @export
plot.jaya <- function(x, ...) {
  hist <- x$history
  if (length(hist) == 0) stop("No history to plot.")
  plot(seq_along(hist), hist, type = "o", pch = 19,
       xlab = "Iteration", ylab = "Best f(x)", main = "Jaya Convergence", ...)
}

#' @export
plot.jaya_multi <- function(x, ...) {
  PF <- x$Pareto_Front
  if (ncol(PF) < 2) stop("Need at least 2 objectives to plot.")
  if (ncol(PF) == 2) {
    plot(PF[,1], PF[,2], xlab = "Objective 1", ylab = "Objective 2",
         main = "Pareto Front", pch = 19, ...)
  } else {
    pairs(PF, main = "Pareto Pairs", pch = 19, ...)
  }
}

#' @export
autoplot.jaya <- function(object, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(plot.jaya(object, ...))
  }
  df <- data.frame(iter = seq_along(object$history),
                   best = as.numeric(object$history))
  ggplot2::ggplot(df, ggplot2::aes(iter, best)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Iteration", y = "Best f(x)", title = "Jaya Convergence")
}

#' @export
autoplot.jaya_multi <- function(object, type = c("scatter","pairs","parallel"), ...) {
  type <- match.arg(type)
  PF <- object$Pareto_Front
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(plot.jaya_multi(object, ...))
  }
  if (type == "pairs") {
    if (!requireNamespace("GGally", quietly = TRUE)) {
      warning("Install 'GGally' for pairs plot. Falling back to base pairs().")
      return(plot.jaya_multi(object, ...))
    }
    return(GGally::ggpairs(PF))
  }
  if (type == "parallel") {
    if (!requireNamespace("GGally", quietly = TRUE)) {
      warning("Install 'GGally' for parallel coord. Falling back to pairs().")
      return(plot.jaya_multi(object, ...))
    }
    return(GGally::ggparcoord(PF, columns = seq_len(ncol(PF)), scale = "globalminmax"))
  }
  # scatter: pick first two objectives
  ggplot2::ggplot(PF, ggplot2::aes(PF[[1]], PF[[2]])) +
    ggplot2::geom_point() +
    ggplot2::labs(x = names(PF)[1], y = names(PF)[2], title = "Pareto Front")
}

#' Save any base or ggplot figure
#' @param obj Plot object or a function that draws the plot.
#' @param file Output path (pdf/png).
#' @param width,height Numeric dimensions (inches for pdf).
#' @export
save_plot <- function(obj, file, width = 6, height = 4, ...) {
  ext <- tools::file_ext(file)
  if (ext %in% c("pdf", "png")) {
    if (ext == "pdf") grDevices::pdf(file, width = width, height = height)
    if (ext == "png") grDevices::png(file, width = width, height = height, units = "in", res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
    if (is.function(obj)) obj() else {
      if ("ggplot" %in% class(obj)) print(obj) else plot(obj, ...)
    }
  } else {
    stop("Unsupported file extension: ", ext)
  }
}