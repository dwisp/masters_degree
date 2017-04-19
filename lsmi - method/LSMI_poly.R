#####################################################################
## implementation of LSMI with polynomial kernels of various order ##
#####################################################################


lsmi.poly <- function(x, y, pow.from = 1, pow.to = 3, pow.by = 1, MImethod = c('bits', 'nats', 'suzuki')) {
  library(magrittr)
  #! various technical operation and
  # foolproof code
  # settings etc
  n <- length(x)
  dotProd <- function(x, y) sum(x*y) # dot product utility function
  
  # calculating all necessary dot products
  #! assuming x and y are lists of vectors of dim 'd'
  allPairwiseDPs  <- matrix(mapply(dotProd, rep(x), rep(y, each = n)), 
                            nrow = n, byrow = TRUE) # <x_i, y_j>; y in the rows, x in the columns
  allPairwiseDPs_inArow <- as.numeric(t(allPairwiseDPs)) # a format more suitable for processing later on for creation of H_cv's
  
  # <x_i, y_i>, for all i, i.e. paired observations
  pairwiseDPs <- diag(allPairwiseDPs)
  #! here be parameters set by the user (on the later stages of development
  
  ## sequence of powers for the base funtion
  ker.pow.space <- seq(pow.from, pow.to, pow.by)
  nbfun <- length(ker.pow.space)
  ## sequence of addition parameter to search through; depends on the range of dot products
  ker.add.scale <- ceiling(median(pairwiseDPs) + mad(pairwiseDPs))
  ker.add.step <- ceiling(ceiling(mad(pairwiseDPs)) / ker.add.scale)
  ker.add.space <- seq(-ker.add.scale, ker.add.scale, by = ker.add.step)
  
  lambda.space <- 10 ^ seq(-3, 0, length.out = 8)
  
  # Cross-Validation
  fold <- 5
  ### randomly splitting data into K disjoint subsets
  if(n %% fold == 0) {
    cvSubset <- rep(1:fold, n %/% fold)
  } else {
    cvSubset <- c(rep(1:fold, n %/% fold), c(1:fold)[1:(n %% fold)])
  }
  cvSubset %<>% sample 
  # n_ <- rep(0, fold)
  # for(k in 1:fold) n_[k] <- sum(cvSubset == k)
  
  # utility function which is used repeatedly to calculate matrices during cross validation
  ## lesser utility function; defined in this envirorment to be initialized only once
  revPow <- function(x, y) y ^ x 
  doPolyKerPHI <- function(dotprods, const2add, powers) {
    polyKerMatrix <- matrix(nrow = length(powers), ncol = length(dotprods))
    firstRow <- dotprods + const2add
    polyKerMatrix <- sapply(as.list(powers), revPow, firstRow)
    t(polyKerMatrix)
  }
  
  # these are 'address' vectors for easier subsetting & calcultaion of H (to obtain all pairwise phi_ij)
  toSubsetX <- rep(1:n, n)
  toSubsetY <- rep(1:n, each = n)
  
  ## this utility function uses few variables in lsmi.poly's environment hence only one input argument
  ## it is used for creation of H inside of cross-validation
  doCV.H <- function(samplesInSubset, ker.const) {
    x2take <- toSubsetX %in% samplesInSubset
    y2take <- toSubsetY %in% samplesInSubset
    n_cv <- length(samplesInSubset)
    
    dp.relevant <- 
      allPairwiseDPs_inArow %>%
      extract(x2take & y2take)
    
    pariwisePHI <- doPolyKerPHI(dp.relevant, ker.const, ker.pow.space)
    cv.H <- 
      vapply(split(pariwisePHI, rep(1:ncol(pariwisePHI), each = nrow(pariwisePHI))), # splitting matrix into column-lists to feed in to apply
             function(x) x %*% t(x), # formula taken from the paper
             matrix(0, nbfun, nbfun)) %>% # format of the output: produces nbfun X nbfun matrix
      apply(c(1,2), sum) %>% # summation of all resulting matrices;
      divide_by(n_cv^2)
    ### vapply produces an array which is a lot easier to handle than lapply's list or sapply's stangely concantenated matrix
    cv.H
  }
  
  # vapply(split(allPairwiseDPs, rep(1:(n^2), each = nbfun)), function(x) x %*% t(x), matrix(0, nrow = nbfun, ncol = nbfun))
  
  # creating an array of PHI matrices for CV
  cv.KerPHIs <- array(sapply(as.list(ker.add.space), FUN = doPolyKerPHI, pairwiseDPs, ker.pow.space), 
                      dim = c(length(ker.pow.space), n, length(ker.add.space)))
  
  ## addition constants in rows, lambdas in columns
  cv.Scores <- matrix(0, nrow = length(ker.add.space), ncol = length(lambda.space),
                      dimnames = list(round(ker.add.space, 3), round(lambda.space, 3)))
  
  #for(bfuns in seq(2, length(ker.pow.space), by = 1)) { # will add nbfuns cv later
    for(ker.add in ker.add.space) {
      #print(ker.add)
      ker.add.ind <- which(ker.add.space == ker.add)
      ker.add.val <- ker.add.space[ker.add.ind]
      
      for(k in 1:fold) {
        kSlice <- cvSubset == k
        kSliceIndices <- which(kSlice)
        kSliceOutIndices <- which(!kSlice)
        ### PHIs
        cv.PHI.in <- cv.KerPHIs[, kSlice, ker.add.ind]
        cv.PHI.out <- cv.KerPHIs[, !kSlice, ker.add.ind]
        ### in-subset matrices
        cv.h.in <- rowMeans(cv.PHI.in)
        cv.H.in <- doCV.H(kSliceIndices, ker.add.val)
        ### out-of-subset matrices
        cv.h.out <- rowMeans(cv.PHI.out)
        cv.H.out <- doCV.H(kSliceOutIndices, ker.add.val)
        
        
        for(lambda in lambda.space) {
         #print(lambda)
          lambda.ind <- which(lambda.space == lambda)
          
          cv.alpha <- solve(cv.H.out + lambda*diag(nbfun)) %*% cv.h.out
          J.cv <- t(cv.alpha) %*% cv.H.in %*% cv.alpha / 2 - dotProd(cv.h.in, cv.alpha)
          
          cv.Scores[ker.add.ind, lambda.ind] %<>% add(J.cv/fold)
        }
      }
    }
  fin.par.ind <- 
    cv.Scores %>%
    {. == min(.)} %>%
    which(arr.ind = TRUE)
  #}
  fin.add <- ker.add.space[fin.par.ind[1]]
  fin.lambda <- lambda.space[fin.par.ind[2]]
  fin.PHI <- cv.KerPHIs[,,fin.par.ind[1]]
  fin.H <- doCV.H(1:n, fin.add)
  fin.h <- rowMeans(fin.PHI)
  fin.alha <- solve(fin.H + fin.lambda*diag(nbfun)) %*% fin.h
  
  MImethod <- match.arg(MImethod)
  # print('alpha = ')
  # print(fin.alha)
  # print('alpha^t * PHI = ')
  # print(t(fin.alha) %*% fin.PHI)
  # print('n of negative instances')
  # print(sum((t(fin.alha) %*% fin.PHI) <= 0))
  # print('final PHI')
  # print(fin.PHI)
  # print('final H')
  # print(fin.H)
  # print(rankMatrix(fin.H))
  MI <- switch(MImethod,
               bits = 
                 t(fin.alha) %*% fin.PHI %>%
                 inset(. <= 0, 0.001) %>% # if this is removed, function will rarely produce NaNs
                 log(base = 2) %>%
                 mean,
               nats = 
                 t(fin.alha) %*% fin.PHI %>%
                 inset(. <= 0, 0.001) %>% # if this is removed, function will rarely produce NaNs
                 log(base = exp(1)) %>%
                 mean,
               suzuki =
                 t(fin.h) %*% fin.alha %>%
                 as.numeric %>%
                 subtract(1)
  )
  MI
}