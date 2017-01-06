#########################
# Implementaton of LSMI #
#########################

# very close to classical implementation except a few changes:
# 1. arbitrary multiplier c instead of 1/2
# 2. fixed 'mean' bug found in native Matlab implementation
# 3. slightly modified parameter space:
# excluded largest parameters lambda and sigma which are probably never actually used to reduce computational time
# 4. number of basis functions is 100 if sample size is greater than 100

lsmi.vanilla <- function(x, y, sigmas, lambdas, nbfuns, c = 1.5) {
  # x and y are LISTS for the sake of multi-dimensionality
  
  # dealing with the arguments #
  if(missing(sigmas)) sigmas <- 10^seq(-2, 1, length.out = 8)#10^seq(-2, 2, length.out = 9)
  if(missing(lambdas)) lambdas <- 10^seq(-3, 0, length.out = 8)
  if(missing(nbfuns)) nbfuns <- min(length(x), 100)
  
  ## some notations for convenience
  if(length(x) != length(y)) stop('Lengths of x and y are different!')
  n <- length(x)
  
  require(magrittr)
  ## basis functions - maybe will get expanded later ##
  ### we only take nbfuns out of total of n
  centroids <- 1:n
  if(nbfuns < n) centroids %<>% sample(.) %>% extract(1:nbfuns)
  
  ## squared distances are only computed once!
  dist2x <- dist(x) %>% as.matrix %>% extract(centroids, ) %>% .^2
  dist2y <- dist(y) %>% as.matrix %>% extract(centroids, ) %>% .^2
    
  if(length(sigmas) == 1 & length(lambdas) == 1) {
    cvScores <- -Inf
    # no CV
    sigma.chosen <- sigmas
    lambda.chosen <- lambdas
    phix.fin <- dist2x %>% divide_by(-2*sigma.chosen^2) %>% exp
    phiy.fin <- dist2y %>% divide_by(-2*sigma.chosen^2) %>% exp
    
  } else { # CV
    fold <- 5
    cvScores <- matrix(0, nrow = length(sigmas), ncol = length(lambdas))
    dimnames(cvScores) <- list(round(sigmas, 4), round(lambdas, 4))
    
    ### randomly splitting data into K disjoint subsets
    if(n %% fold == 0) {
      cvSubset <- rep(1:fold, n %/% fold)
    } else {
      cvSubset <- c(rep(1:fold, n %/% fold), c(1:fold)[1:(n %% fold)])
    }
    cvSubset %<>% sample 
    n_ <- rep(0, fold)
    for(k in 1:fold) n_[k] <- sum(cvSubset == k)
    ###
    ##
    phisx <- array(dist2x, dim = c(nbfuns, n, length(sigmas)))
    phisx %<>% sweep(3, -2*sigmas^2, '/') %>% exp
    
    phisy <- array(dist2y, dim = c(nbfuns, n, length(sigmas)))
    phisy %<>% sweep(3, -2*sigmas^2, '/') %>% exp
    ##
    H.cv.in <- array(0, dim = c(nbfuns, nbfuns, fold))
    h.cv.in <- array(0, dim = c(nbfuns, 1, fold))
    
    H.cv.out <- array(0, dim = c(nbfuns, nbfuns, fold))
    h.cv.out <- array(0, dim = c(nbfuns, 1, fold))
    
    ## cross-validation calculations ##
    ### !all centroids are used during CV, even those out of CV subset!
    for(sigma in sigmas) {
      sigma.ind <- which(sigmas == sigma)
      
      for(k in 1:fold) {
        kSlice <- cvSubset == k
        ### in-subset matrices
        h.cv.in[,,k] <- rowMeans(phisx[, kSlice, sigma.ind]*phisy[, kSlice, sigma.ind]) #%>% as.matrix
        H.cv.in[,,k] <- (phisx[, kSlice, sigma.ind] %*% t(phisx[, kSlice, sigma.ind])) *
                        (phisy[, kSlice, sigma.ind] %*% t(phisy[, kSlice, sigma.ind])) %>%
                        divide_by(n_[k]^2)
        ### out-of-subset matrices
        h.cv.out[,,k] <- rowMeans(phisx[, !kSlice, sigma.ind]*phisy[, !kSlice, sigma.ind]) #%>% as.matrix
        H.cv.out[,,k] <- (phisx[, !kSlice, sigma.ind] %*% t(phisx[, !kSlice, sigma.ind])) *
                         (phisy[, !kSlice, sigma.ind] %*% t(phisy[, !kSlice, sigma.ind])) %>%
                         divide_by((n-n_[k])^2)

        ## sum_i what^2(x_i,y_i) = alpha^t H alpha = alpha^t h h^t alpha
        ## this is simple, but not THAT trivial
        for(lambda in lambdas) {
          lambda.ind <- which(lambdas == lambda)
          alpha.cv <- solve(H.cv.out[,,k] + lambda*diag(nbfuns)) %*% h.cv.out[,,k]
          J.cv <- t(alpha.cv) %*% H.cv.in[,,k] %*% alpha.cv /2 - t(h.cv.in[,,k]) %*% alpha.cv
          
          cvScores[sigma.ind, lambda.ind] %<>% add(J.cv/fold)
        }
      }
    }
    # choosing the optimal parameters #
    chosen.pars.ind <- which(cvScores == min(cvScores), arr.ind = TRUE)
    sigma.chosen <- sigmas[chosen.pars.ind[1]]
    lambda.chosen <- lambdas[chosen.pars.ind[2]]

    phix.fin <- phisx[,,chosen.pars.ind[1]]
    phiy.fin <- phisy[,,chosen.pars.ind[1]]
  }
  #print('cvScores')
  #print(cvScores)
  
  #print(sigma.chosen)
  #print(lambda.chosen)
  
  phi.fin <- phix.fin*phiy.fin
  h.fin <- rowMeans(phi.fin)
  H.fin <- (phix.fin %*% t(phix.fin)) * (phiy.fin %*% t(phiy.fin)) %>% divide_by(n^2)
  alpha.fin <- solve(H.fin + lambda.chosen*diag(nbfuns)) %*% h.fin
  # final result #
  MI <- t(h.fin) %*% alpha.fin %>% as.numeric %>% subtract(1) %>% multiply_by(c)
  MI
}

