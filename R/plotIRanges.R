utils::globalVariables(c("start", "end", "ymin", "ymax"))

plotIRanges <- function(ir, xlim = NULL, main = "IRanges Plot", col = "blue", ...) {
  if (is.null(xlim)) {
    xlim <- range(IRanges::start(ir), IRanges::end(ir))
  }

  plot(0, 0, type = "n", xlim = xlim, ylim = c(0, length(ir) + 1),
       xlab = "Position", ylab = "Ranges", main = main, ...)

  for (i in seq_along(ir)) {
    graphics::arrows(IRanges::start(ir)[i], i, IRanges::end(ir)[i], i, code = 3, angle = 90, length = 0.1, col = col)
  }
}

plotRangesDF <- function(ir, xlim = NULL, main = "Ranges DF Plot", col = "blue", bar_height = .25, ...) {
  ir$ymin = rownames(ir) %>% as.integer() - (bar_height/2)
  ir$ymax = rownames(ir) %>% as.integer() + (bar_height/2)

  ggplot2::ggplot(data=ir %>% head()) + ggplot2::geom_rect( mapping=ggplot2::aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax))
}


