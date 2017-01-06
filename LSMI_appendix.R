

# comparing 'true' density ratio with the estimation #

## true multivariate normal density ratio ##
true.mvt.dr <- function(x, y, m, sig2x2) {
  sig.x <- sig2x2[1,1]
  sig.y <- sig2x2[2,2]
  sig.xy <- sig2x2[2,1]
  m.x <- m[1]
  m.y <- m[2]

  rho <- sig.xy/(sig.x*sig.y)
  dr.val <- (1-rho^2)^(-1/2) * exp(rho*((x-m.x)/sig.x)*((y-m.y)/sig.y))
  dr.val # density ratio
}

# comparison:

## density ratio values:
estimated.dr <- rep(0, n)
true.dr <- rep(0, n)
for(i in 1:n) {
  estimated.dr[i] <- t(alpha.hat) %*% PHI(x.list[[i]], y.list[[i]], u.list, v.list, sig.kern)
  true.dr[i] <- true.mvt.dr(x.list[[i]], y.list[[i]], mean.mvt, sigma.mvt)
}

## 'MI' values
I.true.term <- function(x, y, m, sig2x2) {
  I.term.val <- (1 - true.mvt.dr(x, y, m, sig2x2))^2
  I.term.val
}

I.true <- mapply(I.true.term, x.list.pairs, y.list.pairs, MoreArgs = list(mean.mvt, sigma.mvt)) %>%
          sum %>%
          divide_by(n^2)

# I.dif[dif.i, dif.j] <- I.true - I.hat
#   }
#   dif.j <- 0
# }

# checking how much variance 'true' LSMI has #

exper.I.true <- rep(0, 100)
for(i in 1:100) {
  xy.list <- replicate(n, rmvnorm(1, mean.mvt, sigma.mvt), simplify = FALSE)
  x.list <- list()
  y.list <- list()
  for(z in 1:n) {
    x.list[[z]] <- xy.list[[z]][1]
    y.list[[z]] <- xy.list[[z]][2]
  }
  x.list.pairs <- rep(x.list, each = n)
# this vector contains all values of y repeated n times _as a whole_
  y.list.pairs <- rep(y.list, times = n)

  I.true <- mapply(I.true.term, x.list.pairs, y.list.pairs, MoreArgs = list(mean.mvt, sigma.mvt)) %>%
  sum %>%
  divide_by(n^2)
  exper.I.true[i] <- I.true
}




