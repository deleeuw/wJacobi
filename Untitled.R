
jacobiPlot <- function(i, j, n) {
  y <- rep(1:n, n)
  x <- c()
  for (k in 1:n) {
    x <- c(x, rep(k, n))
  }
  plot(
    x,
    y,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    col = "RED",
    cex = 2,
    type = "n"
  )
  for (k in 1:(n * n)) {
    ni <- (n + 1) - i;
    nj <- (n + 1) - j
    if ((x[k] == i) || (y[k] == ni) || (x[k] == j) || (y[k] == nj)) {
      ppch = 15
      if (((x[k] == i) && (y[k] == ni)) ||
          ((x[k] == i) && (y[k] == nj)) ||
          ((x[k] == j) && (y[k] == ni)) ||
          ((x[k] == j) && (y[k] == nwj)))
      {
        pcex = 3
        pcol = "RED"
      } else {
        pcex = 2
        pcol = "BLUE"
      }
    } else {
      pcol = "BLACK"
      pcex = .5
      ppch = 0
    }
    points(x[k],
           y[k],
           col = pcol,
           cex = pcex,
           pch = ppch)
    axis(
      1,
      at = 1:n,
      side = 3,
      labels = as.character(1:n),
      tick = FALSE
    )
    axis(2,
         at = 1:n,
         labels = as.character(n:1),
         tick = FALSE)
  }
}
