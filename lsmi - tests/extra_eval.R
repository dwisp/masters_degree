###
# testing out new base functions selection systems
###

##
# non-paired centroids vs. discrete Y
##
nbfun.lsmi.ydiscr <- 
  data.frame(base.funs = rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps))

set.seed(451)
iris.cropS <- 
  iris$Species %>%
  equals('setosa') %>%
  which %>%
  sample(size = 25)

irisSpecies.cropS <- 
  iris$Species %>%
  extract(-iris.cropS)

iris4lsmi.cropS <-
  iris4lsmi %>%
  extract(-iris.cropS)

###

set.seed(451)
iris.cropVi <- 
  iris$Species %>%
  equals('virginica') %>%
  which %>%
  sample(size = 25)

irisSpecies.cropVi <- 
  iris$Species %>%
  extract(-iris.cropVi)

iris4lsmi.cropVi <-
  iris4lsmi %>%
  extract(-iris.cropVi)

###

indices.vals <- 1:nbfun.prec.reps
for(i in seq(15, 125, by = nbfun.prec.step)) {
  ### cropped data ###
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.S.np'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris4lsmi.cropS, irisSpecies.cropS, nbfuns = i, method.nbfuns = 'non-paired'))
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.S.bl'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris4lsmi.cropS, irisSpecies.cropS, nbfuns = i, method.nbfuns = 'balanced'))
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.S.Orig'] <-
    replicate(nbfun.prec.reps, lsmi.vanilla(iris4lsmi.cropS, irisSpecies.cropS, nbfuns = i))
  ##
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.Vi.np'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris4lsmi.cropVi, irisSpecies.cropVi, nbfuns = i, method.nbfuns = 'non-paired'))
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.Vi.bl'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris4lsmi.cropVi, irisSpecies.cropVi, nbfuns = i, method.nbfuns = 'balanced'))
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.Vi.Orig'] <-
    replicate(nbfun.prec.reps, lsmi.vanilla(iris4lsmi.cropVi, irisSpecies.cropVi, nbfuns = i))
  ### original data ###
  
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.F.np'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris4lsmi, iris$Species, nbfuns = i, method.nbfuns = 'non-paired'))
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.F.bl'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris4lsmi, iris$Species, nbfuns = i, method.nbfuns = 'balanced'))
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.F.Orig'] <-
    replicate(nbfun.prec.reps, lsmi.vanilla(iris4lsmi, iris$Species, nbfuns = i))

  indices.vals %<>% add(nbfun.prec.reps)
}
rm(indices.vals)

# mean squarred error compared to the median LSMI computed using all possible base functions
## splitting results to feed them as lists to mse
nbfun.mse.ydiscr <- data_frame(base.funs = seq(15, 125, by = nbfun.prec.step))
for(var in colnames(nbfun.lsmi.ydiscr)[-1]) {
  nbfun.mse.ydiscr[, var] <- sapply(
    split(nbfun.lsmi.ydiscr[, var], rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = median(nbfun.lsmi.ydiscr[, var] %>% tail(nbfun.prec.reps)))
}
library(reshape2)
nbfun.mse.ydiscr %<>%
  melt(id.vars = 'base.funs')

nbfun.mse.ydiscr$variable %<>% as.character()
###
# plotting RMSE
###
library(ggplot2)

### Full iris data ###
ggplot(filter(nbfun.mse.ydiscr, str_detect(nbfun.mse.ydiscr$variable, 'F')), 
       aes(x = base.funs, y = value, color = variable)) + 
  stat_smooth(method = 'loess', se = FALSE, size = 0.75, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.8) + 
  scale_x_continuous(breaks = seq(15, 125, by = 10)) + 
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.1), limits = c(0, 1.3)) +
  guides(color = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold')),
         shape = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold'))) +
  theme(legend.background = element_rect(colour = 'darkgrey', fill = 'lightgrey'), legend.position = c(0.75, 0.8), legend.key.height = unit(1.5, 'lines'), 
        title = element_text(size = 8, angle = 0)) +
  scale_color_discrete(labels = c('Balanced', 'Non-paired', 'Uniform')) +
  labs(x = 'Number of centroids considered',
       y = 'Root of mean squared error',
       color = 'Centroid selection method') +
  ggtitle(
'Comparison of non-paired and uniform centroid selection methods on Fisher iris data
X is a set of 4-dimensional vectors and Y are labels denoting species'
  ) +
  ggsave('nbfuns_mse_full_ydiscr.png', width = 8, height = 4)

### Virginica class sparsed ###
ggplot(filter(nbfun.mse.ydiscr, str_detect(nbfun.mse.ydiscr$variable, 'Vi')), 
       aes(x = base.funs, y = value, color = variable)) + 
  stat_smooth(method = 'loess', se = FALSE, size = 0.75, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.8) + 
  scale_x_continuous(breaks = seq(15, 125, by = 10)) + 
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.1), limits = c(0, 1.3)) +
  guides(color = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold')),
         shape = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold'))) +
  theme(legend.background = element_rect(colour = 'darkgrey', fill = 'lightgrey'), legend.position = c(0.75, 0.8), legend.key.height = unit(1.5, 'lines'), 
        title = element_text(size = 8, angle = 0)) +
  scale_shape_discrete(labels = c('Balanced', 'Non-paired', 'Uniform')) +
  scale_color_discrete(labels = c('Balanced', 'Non-paired', 'Uniform')) +
  labs(x = 'Number of centroids considered',
       y = 'Root of mean squared error',
       color = 'Centroid selection method',
       shape = 'Centroid selection method') +
  ggtitle(
'Comparison of non-paired and uniform centroid selection methods on Fisher iris data
Virginica class is half-sparsed. X is a set of 4-dimensional vectors and Y are labels denoting species'
  ) +
  ggsave('nbfuns_mse_virginica_ydiscr.png', width = 8, height = 4)

### Setosa class sparsed ###
ggplot(filter(nbfun.mse.ydiscr, str_detect(nbfun.mse.ydiscr$variable, 'S')), 
       aes(x = base.funs, y = value, color = variable)) + 
  stat_smooth(method = 'loess', se = FALSE, size = 0.75, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.8) + 
  scale_x_continuous(breaks = seq(15, 125, by = 10)) + 
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.1), limits = c(0, 1.3)) +
  guides(color = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold')),
         shape = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold'))) +
  theme(legend.background = element_rect(colour = 'darkgrey', fill = 'lightgrey'), legend.position = c(0.75, 0.8), 
        legend.key.height = unit(1.5, 'lines'), title = element_text(size = 8, angle = 0)) +
  scale_shape_discrete(labels = c('Balanced', 'Non-paired', 'Uniform')) +
  scale_color_discrete(labels = c('Balanced', 'Non-paired', 'Uniform')) +
  labs(x = 'Number of centroids considered',
       y = 'Root of mean squared error',
       color = 'Centroid selection method',
       shape = 'Centroid selection method') +
  ggtitle(
'Comparison of non-paired and uniform centroid selection methods on Fisher iris data
Setosa class is half-sparsed. X is a set of 4-dimensional vectors and Y are labels denoting species'
  ) +
  ggsave('nbfuns_mse_setosa_ydiscr.png', width = 8, height = 4)



##
# non-paired centroids vs. continuous Y
##
nbfun.lsmi.ycont <- data_frame(base.funs = rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps))

indices.vals <- 1:nbfun.prec.reps
for(i in seq(15, 125, by = nbfun.prec.step)) {
  nbfun.lsmi.ycont[indices.vals, 'lsmi.Full'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris[, 1], iris[, 3], nbfuns = i, method.nbfuns = 'non-paired'))
  nbfun.lsmi.ycont[indices.vals, 'lsmi.Full.Orig'] <-
    replicate(nbfun.prec.reps, lsmi.vanilla(iris[, 1], iris[, 3], nbfuns = i))
  
  nbfun.lsmi.ycont[indices.vals, 'lsmi.S'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris[-iris.cropS, 1], iris[-iris.cropS, 3], nbfuns = i, method.nbfuns = 'non-paired'))
  nbfun.lsmi.ycont[indices.vals, 'lsmi.S.Orig'] <-
    replicate(nbfun.prec.reps, lsmi.vanilla(iris[-iris.cropS, 1], iris[-iris.cropS, 3], nbfuns = i))

  nbfun.lsmi.ycont[indices.vals, 'lsmi.Vi'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris[-iris.cropVi, 1], iris[-iris.cropVi, 3], nbfuns = i, method.nbfuns = 'non-paired'))
  nbfun.lsmi.ycont[indices.vals, 'lsmi.Vi.Orig'] <-
    replicate(nbfun.prec.reps, lsmi.vanilla(iris[-iris.cropVi, 1], iris[-iris.cropVi, 3], nbfuns = i))
  
  
  indices.vals %<>% add(nbfun.prec.reps)
}
rm(indices.vals)
nbfun.lsmi.ycont %<>% as.data.frame()

nbfun.mse.ycont <- data_frame(base.funs = seq(15, 125, by = nbfun.prec.step))
for(var in colnames(nbfun.lsmi.ycont)[-1]) {
  nbfun.mse.ycont[, var] <- sapply(
    split(nbfun.lsmi.ycont[, var], rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = nbfun.lsmi.ycont[, var] %>% tail(nbfun.prec.reps) %>% median)
}
library(reshape2)
nbfun.mse.ycont %<>% melt(id.vars = 'base.funs')

ggplot(nbfun.mse.ycont, aes(x = base.funs, y = value, color = variable, shape = variable)) +
  stat_smooth(method = 'loess', se = FALSE, alpha = 0.8, size = 0.7) +
  geom_point(size = 2.5, alpha = 0.8) + 
  guides(color = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold'))) +
  scale_x_continuous(breaks = seq(15, 125, by = 10)) + 
  scale_y_continuous(breaks = seq(0, 9, by = 0.1), limits = c(0, 0.9)) +
  scale_color_brewer(palette = 'Paired', 
                     labels = str_c(c('Non-paired\n', 'Uniform\n'), rep(c('Full data', 'Setosa sparsed', 'Virginica sparsed'), each = 2))) + 
  scale_shape_manual(values = rep(c(17, 16), 3), 
                     labels = str_c(c('Non-paired\n', 'Uniform\n'), rep(c('Full data', 'Setosa sparsed', 'Virginica sparsed'), each = 2))) +
  labs(x = 'Number of centroids considered', y = 'Root of mean squared error', color = 'Centroid selection method', shape = 'Centroid selection method') +
  theme(legend.background = element_rect(colour = 'darkgrey', fill = 'lightgrey'), 
        legend.position = c(0.75, 0.8), legend.key.height = unit(1.5, 'lines'), 
        title = element_text(size = 8, angle = 0)) +
  guides(color = guide_legend(ncol = 3, title.theme = element_text(size = 12, angle = 0, face = 'bold')),
         shape = guide_legend(ncol = 3), title.theme = element_text(size = 12, angle = 0, face = 'bold')) +
  ggtitle(
'Comparison of non-paired and uniform centroid selection methods on Fisher iris data
Mutual information between "Sepal Length" & "Petal Length"'
  ) +
  ggsave('nbfuns_mse_nonpaired_ycont.png', width = 8, height = 4)


##
# comparing log what with log w
##

mvtsigma <- matrix(c(1, 0, 0, 1), ncol = 2)
## easier off-diag cor setting
mvtCor <- 0.3
mvtsigma %<>% inset(!diag(ncol(.)), mvtCor)

library(mvtnorm)
mvtnormXY <- rmvnorm(200, c(0, 0), mvtsigma)
mvtnormX <- mvtnormXY[, 1]
mvtnormY <- mvtnormXY[, 2]

mvtMethod <- 'uniform'
# mvtMethod <- 'non-paired'

mvtWhat <- lsmi.extra(as.list(mvtnormX), as.list(mvtnormY), method.nbfuns = mvtMethod, nbfuns = 100, getWhat = TRUE)

## x, y and z values
mvtWhatVals <- 
  expand.grid(x = seq(min(mvtnormX), max(mvtnormX), length.out = 100),
              y = seq(min(mvtnormY), max(mvtnormY), length.out = 100),
              z = 0) %>%
  as_data_frame()

mvtWhatVec <- Vectorize(mvtWhat, c('x0', 'y0'))
mvtWhatVals$z <- mvtWhatVec(mvtWhatVals$x, mvtWhatVals$y, units = 'I')
mvtWhatVals$logz <- log(mvtWhatVals$z, 2)

mvtWhatVals$z %<>% inset(. <= 0.125, 0.125)
mvtWhatVals$z %<>% log(2)

mvtWhatAlphas <- as.numeric(mvtWhat(1, 1, units = 'alphas'))
mvtWhatCentroids <- mvtWhat(1, 1, units = 'centroids')
mvtWhatSigma <- mvtWhat(1, 1, units = 'sigma')

mvtWhatAux <- data_frame(xcentroid = mvtWhatCentroids[, 'x'], ycentroid = mvtWhatCentroids[, 'y'], 
                         alphas = mvtWhatAlphas, signAlpha = ifelse(alphas <= 0, 'neg', 'pos'))
library(ggplot2)
library(RColorBrewer)
library(stringr)
ggplot() +
  geom_tile(data = mvtWhatVals, aes(x = x, y = y, fill = z)) +
  stat_contour(data = mvtWhatVals, bins = 6, aes(x = x, y = y, z = z), color = 'lightslategrey', size = 0.6, alpha = 0.4) +
  scale_fill_gradientn(colors = brewer.pal(6, 'RdYlGn'),
                       breaks = seq(-3, 3, length.out = 7),
                       labels = c(-3:-1, str_c(' ', 0:3)),
                       limits = c(-3, 3)) +
  geom_point(data = mvtWhatAux, aes(x = xcentroid, y = ycentroid, size = abs(alphas), color = factor(signAlpha)), alpha = 0.9) +
  scale_color_manual(values = c('orangered2', 'chartreuse4'), labels = c('Negative', 'Positive')) +
  scale_size_area(limits = c(0, 5)) + 
  guides(color = guide_legend(title = 'Sign of\nweight alpha', override.aes = list(size = 3.5), order = 1), 
         size = guide_legend(title = 'Absolute value\nof weight alpha', ncol = 2, override.aes = list(color = 'darkgrey', alpha = 1), order = 2),
         fill = guide_colorbar(title = 'Log2(Density ratio)', barheight = 4, barwidth = 3, order = 3)) +
  ggtitle(str_c('LSMI density ratio estimate for bivariate normal distribution, cor = ', mvtsigma[1, 2], 
                '\nKernel Sigma = ', round(mvtWhatSigma, 3), ', ', mvtMethod, ' centroid selection (', nrow(mvtWhatCentroids), ' centroids).')) +
  ggsave(str_c('mvtDensityRatio', 'Cor', mvtsigma[1,2], mvtMethod, '.png'), width = 10, height = 7)

mvtTrueDr <- function(x, y, m = c(0, 0), sig2x2) {
  sig.x <- sig2x2[1,1]
  sig.y <- sig2x2[2,2]
  sig.xy <- sig2x2[2,1]
  m.x <- m[1]
  m.y <- m[2]
    
  rho <- sig.xy/(sig.x*sig.y)
  dr.val <- (1-rho^2)^(-1/2) * exp(rho*((x-m.x)/sig.x)*((y-m.y)/sig.y))
  dr.val # density ratio
}
mvtTrueDr %<>% Vectorize(., c('x', 'y'))

mvtWtrueVals <- 
  expand.grid(x = seq(min(mvtnormX), max(mvtnormX), length.out = 100),
              y = seq(min(mvtnormY), max(mvtnormY), length.out = 100),
              z = 0) %>%
  as_data_frame()

mvtWtrueVals$z <- mvtTrueDr(mvtWtrueVals$x, mvtWtrueVals$y, sig2x2 = mvtsigma)
mvtWtrueVals %<>% mutate(logz = log(z, 2))

library(RColorBrewer)
ggplot(mvtWtrueVals, aes(x = x, y = y, z = logz)) + 
  geom_tile(data = mvtWtrueVals, aes(x = x, y = y, fill = z)) +
  stat_contour(data = mvtWtrueVals, bins = 24, aes(x = x, y = y, z = z), color = 'lightslategrey', size = 0.6, alpha = 0.4) +
  scale_fill_gradientn(colors = brewer.pal(6, 'RdYlGn'),
                       breaks = seq(-3, 3, length.out = 7),
                       labels = c(-3:-1, str_c(' ', 0:3)),
                       limits = c(-3, 3)) +
  guides(fill = guide_colorbar(title = 'Log2(Density ratio)', barheight = 4, barwidth = 3, order = 3)) +
  ggtitle(str_c('True Log2(Density ratio) for bivariate normal distribution, cor = ', mvtsigma[1, 2])) +
  ggsave(str_c('mvtDensityRatioTrue', 'Cor', mvtsigma[1,2], '.png'), width = 10, height = 7)