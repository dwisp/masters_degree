# Modified version of LSMI #
# Should include:
# -balanced classes base function selection
# -'feature space covering technique' base function selection
# -???

#########################
# Implementaton of LSMI #
#########################

# close to classical implementation, but more changes than in 'vanilla' script:
# 1. MI is calculated properly in bits / nats; introducing arbitrary multiplier c instead of 1/2
# 2. fixed 'mean' bug found in native Matlab implementation
# 3. slightly modified parameter space:
# excluded largest parameters lambda and sigma which are probably never actually used to reduce computational time
# it's actually a lot better to use adaptive sigmas w.r.t. distance distribution
# 4. number of kernel functions is 100 if sample size is greater than 100
# 5. kernel functions aren't permuted if max(nbfuns) > current sample size
# 6. added centroid selection options
# 7. it's possible to get log(w_hat(x,y)) as output

lsmi.extra <- function(x, y, method = c('bits', 'nats', 'suzuki'), y.discrete = FALSE, # main parameters
                       # custom base function selection & what retrieval
                       method.nbfuns = c('balanced classes', 'non-paired', 'uniform'), getWhat = FALSE, 
                       sigmas, lambdas, nbfuns, c = 1.5) { # rarely used parameters
  # x and y are LISTS for the sake of multi-dimensionality
  require(magrittr)
  
  # dealing with the arguments #
  if(missing(sigmas)) sigmas <- 10^seq(-2, 1, length.out = 8) # 10^seq(-2, 2, length.out = 9)
  if(missing(lambdas)) lambdas <- 10^seq(-3, 0, length.out = 8)
  if(missing(nbfuns)) nbfuns <- min(length(x), 100)
  MImethod <- match.arg(method, c('bits', 'nats', 'suzuki'))
  NBFmethod <- match.arg(method.nbfuns, c('balanced classes', 'non-paired', 'uniform'))
  ## some notations for convenience
  if(length(x) != length(y)) stop('Lengths of x and y are different!')
  ## support for one-dimensional numerics x and y
  ## if y is integer, it should be converted to factor manually
  if(is.numeric(x)) x %<>% as.list
  if(is.numeric(y) & not(y.discrete)) y %<>% as.list
  
  n <- length(x)
  if(nbfuns > n & NBFmethod != 'non-paired') {
    warning('Specified number of base functions is greater than sample size! Setting nbfuns <- n.')
    nbfuns <- n
  }
  ################################
  ## This section gets modified ##
  ################################
  ## kernel (a.k.a base) functions ##
  ### we only take nbfuns out of total of n
  
  #### discrete case handling; factor intuitively says y is discrete
  if(is.factor(y)) y.discrete <- TRUE
  if(y.discrete) {
    y.origLabs <- y
    ## vector labels are supported as well
    ## we already know that length(x) == length(y)
    if(class(y) != 'list') {
      y %<>% as.list
    } else {
      ## checking if all elements of y are of same length knowing it's a list
      if(y %>% sapply(length) %>% unique %>% length %>% equals(1) %>% not) {
        stop('Lengths of class vectors y should be the same for all x!')
      }
      ## checking if all elements of y are of same class knowing it's a list
      if(y %>% sapply(class) %>% unique %>% length %>% equals(1) %>% not ) {
        stop('All y should be of the same class!')
      }
    }
    ## transforming multi-dimensional class vectors to integers
    ## to later obtain phiY in a similar fashion as if length(y) was 1
    if(length(y[[1]]) != 1) y %<>% lapply(paste, collapse = '')
    ## if y is numeric, then checking for equality is performed faster
    ## so we transform y to numeric in any case
    ## we've just  checked if all elements of y are of same class so y[[1]] is enough
    
    if(!class(y[[1]]) %in% c('numeric', 'integer')) {
      y %<>% sapply(as.factor)
      levels(y) <- 1:length(levels(y))
      y %<>% as.numeric
    }
  }
  
  if(NBFmethod == 'uniform') {
    centroids <- 1:n
    ## if n is less than nbfuns, then no permutation is done on later stages
    if(nbfuns < n) centroids %<>% sample %>% extract(1:nbfuns)
  } else if(NBFmethod == 'balanced classes') {
    if(not(y.discrete)) stop('Balanced classes imply discrete y!')
    split.of.classes <- 
      y %>%
      #unlist %>%
      table %>% 
      multiply_by(nbfuns/n) %>%
      round %>%
      sort(decreasing = TRUE)
    
    # ad-hoc solution for now; sometimes 'fair' rounding can give more than nbfuns (or even n) samples
    if(sum(split.of.classes) > nbfuns) split.of.classes[1] %<>% subtract(sum(split.of.classes) - nbfuns)
    class.indices <- list()
    for(y.ind in seq_along(unique(y))) {
      subclass <- as.numeric(names(split.of.classes)[y.ind])
      subclass.nobjts <- split.of.classes[y.ind]
      
      class.indices[[y.ind]] <- sample(which(y == subclass), subclass.nobjts)
    }
    names(class.indices) <- names(split.of.classes)
    clas.ssbs.cardns <- sapply(class.indices, length)
    ## might be dealt with later
    if(0 %in% clas.ssbs.cardns) {
      warning(paste(sum(clas.ssbs.cardns == 0), 'Some class(es) do(es) not have kernel representation(s)!', collapse = ' '))
    }
    ## total selected centroids from all classes
    centroids <- unlist(class.indices)
  } else if(NBFmethod == 'non-paired') {
    x.centroids <- sample(1:n, nbfuns, TRUE)
    y.centroids <- sample(1:n, nbfuns, TRUE)
    
    ## a rare exception handling:
    ## how much base functions are actually duplicated
    duplBases <- duplicated(cbind(x.centroids, y.centroids))
    nduplBases <- sum(duplBases)
    ## removing the duplicates
    x.centroids %<>% extract(!duplBases)
    y.centroids %<>% extract(!duplBases)
    ## new centroid holders
    x.chain <- integer(0)
    y.chain <- integer(0)
    
    if(nduplBases > 0) {
      ## we need all duplicated centroids sorted out
    while(not(length(c(x.centroids, x.chain)) == nbfuns & 
              sum(duplicated(cbind(c(x.centroids, x.chain), c(y.centroids, y.chain)))) == 0)) {
      x.chain <- sample(1:n, nduplBases, TRUE)
      y.chain <- sample(1:n, nduplBases, TRUE)
      }
    }
    x.centroids %<>% c(x.chain)
    y.centroids %<>% c(y.chain)
  }
  
  if(NBFmethod != 'non-paired') {
    x.centroids <- centroids
    y.centroids <- centroids
  }
  
  if(y.discrete) {
    phiy.discr <- 
      sapply(as.list(y), function(arg) lapply(as.list(y), '==', arg)) %>% 
      extract(y.centroids, ) %>%
      as.matrix
    mode(phiy.discr) <- 'numeric'}
  
  ## squared distances are only computed once!
  require(fields) # faster euclidean distance matrix calculation
  if(length(x[[1]]) != 1) x %<>% Reduce(rbind, .)
  
  dist2x <- 
    rdist(x) %>%
    extract(x.centroids, ) %>%
    {.^2}
  
  if(!y.discrete) {
    if(length(y[[1]] != 1)) y %<>% Reduce(rbind, .) ## test!
    dist2y <-
      rdist(y) %>%
      extract(y.centroids, ) %>%
      {.^2}
  }
  
  if(length(sigmas) == 1 & length(lambdas) == 1) {
    cvScores <- -Inf
    # no CV
    sigma.chosen <- sigmas
    lambda.chosen <- lambdas
    phix.fin <- 
      dist2x %>%
      divide_by(-2*sigma.chosen^2) %>%
      exp
    if(!y.discrete) {
      phiy.fin <-
        dist2y %>%
        divide_by(-2*sigma.chosen^2) %>% 
        exp
    } else {
      phiy.fin <- phiy.discr
    }
    
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
    phisx %<>% 
      sweep(3, -2*sigmas^2, '/') %>%
      exp
    
    if(!y.discrete) {
      phisy <- array(dist2y, dim = c(nbfuns, n, length(sigmas)))
      phisy %<>% 
        sweep(3, -2*sigmas^2, '/') %>%
        exp 
    } else {
      ### this is an array so I don't have to rewrite array slicing below which corresponds to
      ### trying out various sigmas in continious y case; hope that doesn't take so much time and memory
      phisy <- array(phiy.discr, dim = c(dim(phiy.discr), length(sigmas)))
    } 
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
        h.cv.in[,,k] <- rowMeans(phisx[, kSlice, sigma.ind]*phisy[, kSlice, sigma.ind])
        H.cv.in[,,k] <- (phisx[, kSlice, sigma.ind] %*% t(phisx[, kSlice, sigma.ind])) *
          (phisy[, kSlice, sigma.ind] %*% t(phisy[, kSlice, sigma.ind])) %>%
          divide_by(n_[k]^2)
        ### out-of-subset matrices
        h.cv.out[,,k] <- rowMeans(phisx[, !kSlice, sigma.ind]*phisy[, !kSlice, sigma.ind])
        H.cv.out[,,k] <- (phisx[, !kSlice, sigma.ind] %*% t(phisx[, !kSlice, sigma.ind])) *
          (phisy[, !kSlice, sigma.ind] %*% t(phisy[, !kSlice, sigma.ind])) %>%
          divide_by((n-n_[k])^2)
        
        ## sum_i what^2(x_i,y_i) = alpha^t H alpha = alpha^t h h^t alpha
        ## this is simple, but not THAT trivial
        for(lambda in lambdas) {
          lambda.ind <- which(lambdas == lambda)
          alpha.cv <- solve(H.cv.out[,,k] + lambda*diag(nbfuns)) %*% h.cv.out[,,k]
          J.cv <- t(alpha.cv) %*% H.cv.in[,,k] %*% alpha.cv / 2 - t(h.cv.in[,,k]) %*% alpha.cv
          
          cvScores[sigma.ind, lambda.ind] %<>% add(J.cv/fold)
        }
      }
    }
    # choosing optimal parameters #
    chosen.pars.ind <- which(cvScores == min(cvScores), arr.ind = TRUE)
    sigma.chosen <- sigmas[chosen.pars.ind[1]]
    lambda.chosen <- lambdas[chosen.pars.ind[2]]
    
    phix.fin <- phisx[,,chosen.pars.ind[1]]
    phiy.fin <- phisy[,,chosen.pars.ind[1]]
  }
  
  phi.fin <- phix.fin*phiy.fin
  h.fin <- rowMeans(phi.fin)
  H.fin <- (phix.fin %*% t(phix.fin)) * (phiy.fin %*% t(phiy.fin)) %>% divide_by(n^2)
  alpha.fin <- solve(H.fin + lambda.chosen*diag(nbfuns)) %*% h.fin
  
  if(getWhat) {
    xc <- x[x.centroids]
    if(y.discrete) {
      ## keeping original labels to match strings rather than integers
      yc <- y.origLabs[y.centroids] 
    } else {
      yc <- y[y.centroids]
      }
    
    ## returning a function to calculate density ratio in a feasible region ##
    What <- function(x0, y0, yd = FALSE, units = c('log2', 'log10', 'I', 'alphas', 'sigma', 'centroids')) {
      require(fields)
      require(magrittr)
      ## slight automation for convenience
      if(inherits(y0, 'factor') | inherits(y0, 'character')) yd <- TRUE
      ## distance from point x0 to all the centroids
      #af <- alpha.fin
      #sf <- sigma.chosen
      units <- match.arg(units, c('log2', 'log10', 'I', 'alphas', 'sigma', 'centroids'))
      
      phi.x0xc <- 
        do.call(rbind, c(list(x0), xc)) %>%
        rdist() %>%
        extract(-1, 1) %>%
        {.^2} %>%
        divide_by(-2*sigma.chosen^2) %>%
        exp()
      
      if(!yd) {
        ## distance from point y0 to all the centroids assuming it's continuous
        phi.y0yc <-
          do.call(rbind, c(list(y0), yc)) %>%
          rdist() %>%
          extract(-1, 1) %>%
          {.^2} %>%
          divide_by(-2*sigma.chosen^2) %>%
          exp()
      } else {
          phi.y0yc <- as.numeric(yc == y0)
        }
      
      phi.x0y0 <- phi.x0xc*phi.y0yc
      ## w_hat at (x0, y0)
      what.x0y0 <- sum(phi.x0y0*alpha.fin)
      switch(units, log2 = log(what.x0y0, 2), log10 = log(what.x0y0, 10), I = what.x0y0, 
             alphas = alpha.fin, sigma = sigma.chosen, centroids = cbind(x = unlist(xc), y = unlist(yc)))
    }
    #return(Vectorize(What, c('x0', 'y0'))) # prevents a correct return of centroids
    return(What)
  }
  # final result as authors do#
  MI <- switch(MImethod,
               bits = 
                 t(alpha.fin) %*% phi.fin %>%
                 inset(. <= 0, 0.001) %>% # if this is removed, function will rarely produce NaNs
                 log(base = 2) %>%
                 mean,
               nats = 
                 t(alpha.fin) %*% phi.fin %>%
                 inset(. <= 0, 0.001) %>% # if this is removed, function will rarely produce NaNs
                 log(base = exp(1)) %>%
                 mean,
               suzuki =
                 t(h.fin) %*% alpha.fin %>%
                 as.numeric %>%
                 subtract(1) %>%
                 mean()
  )
  MI
}