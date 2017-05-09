###
# testing out new base functions selection systems
###

##
# non-paired centroids vs. discrete Y
##
nbfun.lsmi.ydiscr <- 
  data.frame(base.funs = rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps))

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
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.S'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris4lsmi.cropS, irisSpecies.cropS, nbfuns = i, method.nbfuns = 'non-paired'))
  
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.S.Orig'] <-
    replicate(nbfun.prec.reps, lsmi.vanilla(iris4lsmi.cropS, irisSpecies.cropS, nbfuns = i))
  
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.Vi'] <-
    replicate(nbfun.prec.reps, lsmi.extra(iris4lsmi.cropVi, irisSpecies.cropVi, nbfuns = i, method.nbfuns = 'non-paired'))
  
  nbfun.lsmi.ydiscr[indices.vals, 'lsmi.Vi.Orig'] <- 
    replicate(nbfun.prec.reps, lsmi.vanilla(iris4lsmi.cropVi, irisSpecies.cropVi, nbfuns = i))
  
  indices.vals %<>% add(nbfun.prec.reps)
}
rm(indices.vals)

nbfun.mse.ydiscr <- 
  data_frame(base.funs = seq(15, 125, by = nbfun.prec.step))

# mean squarred error compared to the median LSMI computed using all possible base functions
## splitting results to feed them as lists to mse
nbfun.mse.ydiscr$mse.S <- 
  sapply(
    split(nbfun.lsmi.ydiscr$lsmi.S, rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = median(nbfun.lsmi.ydiscr$lsmi.S %>% tail(nbfun.prec.reps)))

nbfun.mse.ydiscr$mse.S.Orig <- 
  sapply(
    split(nbfun.lsmi.ydiscr$lsmi.S.Orig, rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = median(nbfun.lsmi.ydiscr$lsmi.S.Orig %>% tail(nbfun.prec.reps)))

nbfun.mse.ydiscr$mse.Vi <- 
  sapply(
    split(nbfun.lsmi.ydiscr$lsmi.Vi, rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = median(nbfun.lsmi.ydiscr$lsmi.Vi %>% tail(nbfun.prec.reps)))

nbfun.mse.ydiscr$mse.Vi.Orig <- 
  sapply(
    split(nbfun.lsmi.ydiscr$lsmi.Vi.Orig, rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = median(nbfun.lsmi.ydiscr$lsmi.Vi.Orig %>% tail(nbfun.prec.reps)))

library(reshape2)
nbfun.mse.ydiscr %<>% melt(id.vars = 'base.funs')
# nbfun.mse.ydiscr$var <- 
#   sapply(
#     split(nbfun.lsmi.ydiscr$lsmi.S, rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
#     sd)

## plotting mse
library(ggplot2)
ggplot(nbfun.mse.ydiscr, 
       aes(base.funs, value, 
           color = variable, 
           shape = variable)) + 
  stat_smooth(method = 'loess', se = FALSE) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_segment(aes(x = 38, xend = 15, y = 1.35, yend = 1.35), size = 1.2, color = 'olivedrab4') +
  annotate('text', x = 28, y = 1.4, label = '1.35 - optimal MI value', color = 'black') +
  guides(color = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold')),
         shape = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold'))) +
  scale_x_continuous(breaks = seq(15, 125, by = 10)) + 
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.25)) +
  labs(x = 'Number of centroids considered',
       y = 'Root of mean squared error',
       color = 'Centroid selection method',
       shape = 'Centroid selection method') +
  theme(legend.background = element_rect(colour = 'darkgrey', fill = 'lightgrey'), legend.position = c(0.75, 0.8), legend.key.height = unit(1.5, 'lines'), 
        title = element_text(size = 8, angle = 0)) +
  scale_shape_manual(values = c(18, 18, 17, 17), labels = c('Setosa sparsed\nNon-paired', 'Virginica sparsed\nNon-paired', 'Setosa sparsed\nUniform', 'Virginica sparsed\nUniform')) +
  scale_color_discrete(labels = c('Setosa sparsed\nNon-paired', 'Virginica sparsed\nNon-paired', 'Setosa sparsed\nUniform', 'Virginica sparsed\nUniform')) +
  ggtitle(
    'Comparison of non-paired and uniform centroid selection methods on Fisher iris data
    Setosa / Virginica classes are half-sparsed. X is a set of 4-dimensional vectors and Y are labels denoting species'
  ) +
  ggsave('nbfuns_mse_nonpaired_ydiscr.png')

##
# non-paired centroids vs. continuous Y
##
nbfun.lsmi.ycont <- data_frame(base.funs = rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps))

indices.vals <- 1:nbfun.prec.reps
for(i in seq(15, 125, by = nbfun.prec.step)) {
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

nbfun.mse.ycont <- data_frame(base.funs = seq(15, 125, by = nbfun.prec.step))
nbfun.mse.ycont$mse.S <- 
  sapply(
    split(nbfun.lsmi.ycont$lsmi.S, rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = median(nbfun.lsmi.ycont$lsmi.S %>% tail(nbfun.prec.reps)))

nbfun.mse.ycont$mse.S.Orig <- 
  sapply(
    split(nbfun.lsmi.ycont$lsmi.S.Orig, rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = median(nbfun.lsmi.ycont$lsmi.S.Orig %>% tail(nbfun.prec.reps)))

nbfun.mse.ycont$mse.Vi <- 
  sapply(
    split(nbfun.lsmi.ycont$lsmi.Vi, rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = median(nbfun.lsmi.ycont$lsmi.Vi %>% tail(nbfun.prec.reps)))

nbfun.mse.ycont$mse.Vi.Orig <- 
  sapply(
    split(nbfun.lsmi.ycont$lsmi.Vi.Orig, rep(seq(15, 125, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    rmse, eta = median(nbfun.lsmi.ycont$lsmi.Vi.Orig %>% tail(nbfun.prec.reps)))

nbfun.mse.ycont %<>% melt(id.vars = 'base.funs')

ggplot(nbfun.mse.ycont, aes(x = base.funs, y = value, color = variable)) +
  stat_smooth(method = 'loess', se = FALSE) +
  geom_point(size = 3, alpha = 0.8) + 
  guides(color = guide_legend(ncol = 2, title.theme = element_text(size = 12, angle = 0, face = 'bold'))) +
  scale_x_continuous(breaks = seq(10, 125, by = 10)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.15)) +
  labs(x = 'Number of centroids considered',
       y = 'Root of mean squared error',
       color = 'Centroid selection method') +
  theme(legend.background = element_rect(colour = 'darkgrey', fill = 'lightgrey'), 
        legend.position = c(0.75, 0.8), legend.key.height = unit(1.5, 'lines'), 
        title = element_text(size = 8, angle = 0)) +
  ggtitle(
'Comparison of non-paired and uniform centroid selection methods on Fisher iris data
Mutual information between "Sepal Length" & "Petal Length"'
  ) +
  ggsave('nbfuns_mse_nonpaired_ycont.png')


##
# comparing log what with log w
##

library(mvtnorm)
mvtsigma <-  matrix(c(1, 0, 0, 1), ncol = 2)
## easier off-diag cor setting
mvtCor <- 0
mvtsigma %<>% inset(!diag(ncol(.)), mvtCor)

mvtnormXY <- rmvnorm(100, c(0, 0), mvtsigma)
mvtnormX <- mvtnormXY[, 1]
mvtnormY <- mvtnormXY[, 2]

mvtWhat <- lsmi.extra(as.list(mvtnormX), as.list(mvtnormY), method.nbfuns = 'non-paired', nbfuns = 100, getWhat = TRUE)

## x, y and z values
mvtWhatVals <- 
  expand.grid(x = seq(min(mvtnormX), max(mvtnormX), length.out = 100),
              y = seq(min(mvtnormY), max(mvtnormY), length.out = 100),
              z = 0) %>%
  as_data_frame()

mvtWhatVals$z <- mvtWhat(mvtWhatVals$x, mvtWhatVals$y, units = 'I')

# mvtWhatVals$z %<>% inset(. <= 0, 0.001) # removing numerical artifacts
# mvtWhatVals %<>% mutate(logz = log(z, 10)) # adding alternate z column
mvtWhatAlphas <- as.numeric(mvtWhat(1, 1, units = 'alphas'))
mvtWhatCentroids <- mvtWhat(1, 1, units = 'centroids')

mvtWhatAux <- data_frame(xcentroid = mvtWhatCentroids[, 'x'], ycentroid = mvtWhatCentroids[, 'y'], 
                         alphas = mvtWhatAlphas, signAlpha = ifelse(alphas <= 0, 'neg', 'pos'))


library(ggplot2)
library(RColorBrewer)
library(stringr)
ggplot() +
  geom_tile(data = mvtWhatVals, aes(x = x, y = y, fill = z)) +
  stat_contour(data = mvtWhatVals, bins = 6, aes(x = x, y = y, z = z), color = 'lightslategrey', size = 0.6, alpha = 0.4) +
  scale_fill_gradientn(colors = brewer.pal(6, 'YlOrRd'), 
                       breaks = seq(min(0, min(mvtWhatVals$z)), max(mvtWhatVals$z), length.out = 6), 
                       labels = round(seq(min(0, min(mvtWhatVals$z)), max(mvtWhatVals$z), length.out = 6), 2),
                       limits = c(min(0, min(mvtWhatVals$z)), max(mvtWhatVals$z))) +
  geom_point(data = mvtWhatAux, aes(x = xcentroid, y = ycentroid, size = abs(alphas), color = factor(signAlpha)), 
             alpha = 0.7) +
  scale_color_manual(values = c('orangered2', 'olivedrab'), labels = c('Negative', 'Positive')) +
  scale_size_area(breaks = seq(0, max(abs(mvtWhatAux$alphas)), by = max(mvtWhatAux$alphas)/4),
                  labels = round(seq(0, max(abs(mvtWhatAux$alphas)), by = max(mvtWhatAux$alphas)/4), 2)) + 
  guides(color = guide_legend(title = 'Sign of\ncoefficient alpha', override.aes = list(size = 3.5)), 
         size = guide_legend(title = 'Absolute value\nof coefficient alpha', ncol = 2, override.aes = list(color = 'darkgrey', alpha = 1)),
         fill = guide_legend(title = 'Density ratio', ncol = 2)) +
  ggtitle(str_c('LSMI density ratio estimate for bivariate normal distribution, cor = ', mvtsigma[1,2]))
  ggsave(str_c('mvtDensityRatio', 'Cor', mvtsigma[1,2], '.png'))






