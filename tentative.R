library(numDeriv)


mya <-
  matrix(c(9, -4, -2, -2, -4, 11, -3, -3, -2, -3, 6, 0, -2, -3, 0, 6), 4, 4)
myw <- 1 - diag(4)
myw[3, 4] <- myw[4, 3] <- 0
myx <- matrix(1:16 / 16, 4, 4)

fsum <- function(x = myx,
                 a = mya,
                 w = myw) {
  p <- ncol(x)
  f <- 0.0
  for (s in 1:p) {
    for (t in 1:p) {
      f <- f + w[s, t] * (sum(x[, s] * a %*% x[, t])) ^ 2
    }
  }
  return(f / 4)
}

value <- fsum(x = myx, a = mya, w = myw)

numgrad <- array(grad(fsum, x = myx, a = mya, w = myw), dim(myx))
numhess <- hessian(fsum, x = myx, a = mya, w = myw)

mygrad <- function(x = myx,
                   a = mya,
                   w = myw) {
  n <- nrow(a)
  p <- ncol(x)
  g <- matrix(0, n, p)
  for (s in 1:p) {
    for (t in 1:p) {
      g[, s] = g[, s] + w[s, t] * sum(x[, s] * (a %*% x[, t])) * x[, t]
    }
    g[, s] <- a %*% g[, s]
  }
  return(g)
}

myhess <- function(x = myx,
                   a = mya,
                   w = myw) {
  n <- nrow(a)
  p <- ncol(x)
  np = n * p
  h <- matrix(0, np, np)
  for (s in 1:p) {
    sind <- 1:n + (s - 1) * n
    as <- drop(a %*% x[, s])
    h[sind, sind] <-
      w[s, s] * (outer(as, as) + sum(x[, s] * (a %*% x[, s])) * a)
    for (v in 1:p) {
      av <- drop(a %*% x[, v])
      h[sind, sind] <- h[sind, sind] + w[s, v] * outer(av, av)
    }
    for (t in 1:p) {
      if (t == s) {
        next
      }
      tind <- 1:n + (t - 1) * n
      at <- drop(a %*% x[, t])
      h[sind, tind] <-
        w[s, t] * (outer(at, as) + sum(a * outer(x[, t], x[, s])) * a)
    }
  }
  return(h)
}


matgrad <- mygrad()
mathess <- myhess()
