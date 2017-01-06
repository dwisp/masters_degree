######################################################################
# this script contains impelentation of LSMI method by Suzuki et al. #
######################################################################

# generating input data #
## list with iid vectors of 'gene expressions' (or discrete markers)##

n <- 50 # sample size
# dx <- 5 # input data dimensions
# dy <- 2 #
sig.kern <- 1 # sigma value for gaussian kernel
lambda.reg <- 0.01 # regularization parameter

# grid validation
# dif.i <- 0
# dif.j <- 0
# 
# for(sig.kern in seq(0.2, 5, by = 0.2)) {
#   dif.i <- dif.i + 1
#   for(lambda.reg in seq(0, 5, by = 0.2)) {
#   dif.j <- dif.j + 1
#     
# simplest multivariate normal case


# I commented this out, because atm I have to use pre-generated data

sigma.mvt <- matrix(c(1, 0.5,
                      0.5, 1), nrow = 2)
mean.mvt <- c(0, 0)

MI.real <- rep(0, 100)
MI.alt <- rep(0, 100)


#for(zzz in 1:100) {
  
library(mvtnorm)
xy.list <- replicate(n, rmvnorm(1, mean.mvt, sigma.mvt), simplify = FALSE)
x.list <- list()
y.list <- list()

for(i in 1:n) {
  x.list[[i]] <- xy.list[[i]][1]
  y.list[[i]] <- xy.list[[i]][2]
}

# instead I'll try to use simple linear dependence example #
# x.list <- x.list.01 #-24:25
# y.list <- y.list.01#2*x.list + 3

# x.list %<>% as.list
# y.list %<>% as.list

### currently we try this out on small sample sizes, so no randomization for kernel centroids is made
u.list <- x.list
v.list <- y.list
  
## implemeting kernel functions ##
### this function is a component of a vector of kernel function values ###
gauss.kernel <- function(x, y, u, v, sig.k) {
  gk <- exp(-(sum((x - u)^2) + sum((y - v)^2))/(2*sig.k^2) )
  gk
}

## gaussian kernel used in the paper ##
gauss.kernel.discr.y <- function(x, y, u, v, sig.k) {
  gk <- exp(-sum((x - u)^2)/(2*sig.k^2) ) * as.numeric(all(y == v))
  gk
}

### less computationally expensive ###
gauss.kernel.single <- function(x, u, sig.k) {
  gk <- exp(-(sum((x - u)^2))/(2*sig.k^2))
  gk
}

## this function is a vector of basis functions ##
## u.s and v.s should be lists of length b (=n for now) ###
PHI <- function(x, y, u.s, v.s, sig.k, kern2use = gauss.kernel) {
  phi.value <- rep(0, length(u.s)) # which equals n
  phi.value <- mapply(kern2use, u = u.s, v = v.s, MoreArgs = list(x, y, sig.k))
  phi.value
}

PHIsingle <- function(x, u, sig.k, kern2use = gauss.kernel.single) {
  phi.value <- vapply(u, kern2use, FUN.VALUE = numeric(1), x, sig.k)
  phi.value
}

## estimating matrices ##

### h.hat
library(magrittr)
h.hat <- mapply(PHI, x.list, y.list, MoreArgs = list(u.list, v.list, sig.kern)) %>% 
  rowSums %>%
  divide_by(n) %>%
  as.matrix

### H.hat

## creating all pairwise PHI(x_i, y_j) ##
# this vector contains all values of x repeated n times in _each_ _sequentially_
x.list.pairs <- rep(x.list, each = n) 
# this vector contains all values of y repeated n times _as a whole_
y.list.pairs <- rep(y.list, times = n)

# each column in this matrix represents PHI(x_i, y_j) for all i,j
H.hat.base <- mapply(PHI, x.list.pairs, y.list.pairs, MoreArgs = list(u.list, v.list, sig.kern))

H.hat <- matrix(0, nrow = n, ncol = n)
for(index.phi in 1:(n^2)) {
  current.PHI <- as.matrix(H.hat.base[,index.phi])
  H.hat <- H.hat + (current.PHI %*% t(current.PHI))
}
rm(current.PHI, H.hat.base) # removing large objects to free memory and workspace
H.hat %<>% divide_by(n^2)

## estimating alpha ##
alpha.hat <- solve(H.hat + lambda.reg*diag(n)) %*% h.hat

# estimated density ratio:

# I.hat.term equals (w.hat.alpha[i,j] - 1)^2
I.hat.term <- function(x, y, alpha, PHI_fun, ...) {
  w.hat.a.value <- t(alpha.hat) %*% PHI_fun(x, y, ...)
  term <- (w.hat.a.value - 1)^2
  term
}

# estimating LS-mutual information #
I.hat <- mapply(I.hat.term, x.list.pairs, y.list.pairs, 
                MoreArgs = list(alpha.hat, PHI, u.list, v.list, sig.kern)) %>%
         sum %>%
         divide_by(n^2)

I.hat.alt <- t(h.hat) %*% alpha.hat/2 - 1/2

# MI.alt[zzz] <- I.hat.alt
# MI.real[zzz] <- I.hat
# }
# asessing quality by estimating cost J #

## it seems this works differently in the original implementation

# 
# K <- 5#-fold CV
# # randomly splitting data into K disjoint subsets
# if(n%%K == 0) {
#   data.indexes <- rep(1:K, n%/%K)
# } else {
#   data.indexes <- c(rep(1:K, n%/%K), c(1:K)[1:(n%%K)])
# }
# 
# # perturbing indexes
# Kdivision <- sample(data.indexes)
# 
# # Z is a set of CV samples
# Z <- list()
# for(k in 1:K) Z[[k]] <- which(Kdivision == k)
# 
# # cost function for a particular k in 1:K #
# ### x and y are lists of ALL samples, they're subsetted later
# J.k <- function(x, y, Z, k, alpha, PHI_fun, ...) {
#   x.prime <- x[ Z[[k]] ]
#   y.prime <- y[ Z[[k]] ]
# 
#   n.k <- length(Z[[k]])
# 
#   x.p.pairs <- rep(x.prime, each = n.k)
#   y.p.pairs <- rep(y.prime, times = n.k)
# 
#   first.sum <- mapply(PHI_fun, x.p.pairs, y.p.pairs, MoreArgs = list(...)) %>%
#                {.^2} %>% # squaring it
#                sum %>% # sum of squares
#                divide_by(2*n.k^2)
#   print(first.sum)
#   second.sum <- mapply(PHI_fun, x.prime, y.prime, MoreArgs = list(...)) %>%
#                 sum %>%
#                 divide_by(n.k)
#   print(second.sum)
#   J.val <- first.sum - second.sum
#   J.val
# }
# 
# ## obtaining average cost over all k ##
# J_K_CV <- rep(0, K)
# for(k in 1:K) J_K_CV[k] <- J.k(x.list, y.list, Z, k, alpha.hat, PHI, u.list, v.list, sig.kern)
# J_K_CV %<>% mean













