#######################################################
## tests performed on classification subtype of LSMI ##
## experiments on the famous iris dataset            ##
#######################################################

library(magrittr)
# two normal distributions
# what's the typical dispersion value for lsmi in discrete case?
# will it be able to differentiate original case from  one with permuted labels?

tst.disy.x1 <- c(rnorm(50, 2, 2), rnorm(50, -2, 2))
tst.disy.y1 <- c(rep('a', 50), rep('b', 50))
tst.disy.y1_wrong <- 
  tst.disy.y1 %>%
  unlist %>%
  sample

tst.disy.lsmi1 <- replicate(50, lsmi.vanilla(tst.disy.x1, tst.disy.y1, y.discrete = T))
tst.disy.lsmi1_wrong <- replicate(50, lsmi.vanilla(tst.disy.x1, tst.disy.y1_wrong, y.discrete = T))

df.disy.lsmi_comparison <- data.frame(original_labels = tst.disy.lsmi1, 
                                      permuted_labels = tst.disy.lsmi1_wrong)
library(reshape2)
library(dplyr)
df.disy.lsmi_comparison %<>%
  melt %>%
  rename(case = variable)

library(ggplot2)
ggplot(df.disy.lsmi_comparison, aes(case, value)) +
  geom_boxplot(aes(color = case)) +
  labs(y = 'lsmi estimate') + 
  ggtitle('Boxplot of discrete-case LSMI values for normal distributions N(-2, 2) / N(2, 2)
        right boxplot shows the result of calculating LSMI using permuted labels') +
  ggsave('discrete_lsmi_normal_distrs.png')

## testing discrete lsmi on iris dataset to validate clustering results
## core idea is to check if we'll be able to get noticeably higher lsmi values when
## providing 'right' clustering labels to the method

library(dendextend)
hc.iris <- iris[, -5] %>%
  as.matrix %>%
  dist %>%
  hclust(method = 'ward.D2') 

iris.cols <- as.factor(iris$Species)
levels(iris.cols) <- 1:length(levels(iris.cols))
iris.cols %<>% as.numeric

png('hclust_iris.png', width = 1920, height = 1080)
hc.iris %>% 
  as.dendrogram %>%
  set('branches_k_color', k = 3, c('black', 'forestgreen', 'firebrick2')) %>%
  set('leaves_pch', 19) %>%
  set('leaves_col', iris.cols[hc.iris$order]) %>%
  set('leaves_cex', 1) %>%
  plot 

hc.iris %>% 
  as.dendrogram %>%
  rect.dendrogram(k = 3, lty = 5, lwd = 1, xpd = T, lower_rect = -0.25, upper_rect = 1.6,
                  border = c('black', 'forestgreen', 'firebrick2'))
dev.off()

hc.iris.labels <- cutree(hc.iris, k = 3)
hc.iris.labels_wrong <- 
  hc.iris.labels %>%
  sample

## slicing the data by rows to produce list of 4-dimensional vectors
iris4lsmi <- split(as.matrix(iris[, -5]), row(iris[,-5]))
iris.2classes <- 
  iris[5] %>%
  unlist %>%
  as.character %>%
  inset(. == 'versicolor', 'virginica')

iris.2classes_wrong <- 
  iris[5] %>%
  unlist %>%
  as.character %>%
  inset(. == 'setosa', 'virginica')

iris.lsmi <- replicate(25, lsmi.vanilla(iris4lsmi, hc.iris.labels, y.discrete = TRUE, sigmas = 10^seq(-2, 5, length.out = 8)))
iris.lsmi_true <- replicate(25, lsmi.vanilla(iris4lsmi, unlist(iris[5]), y.discrete = TRUE))
iris.lsmi_twoclasses <- replicate(25, lsmi.vanilla(iris4lsmi, iris.2classes, y.discrete = TRUE))
iris.lsmi_twoclasses_wrong <- replicate(25, lsmi.vanilla(iris4lsmi, iris.2classes_wrong, y.discrete = TRUE))
iris.lsmi_wrong <- replicate(25, lsmi.vanilla(iris4lsmi, hc.iris.labels_wrong, y.discrete = TRUE))


iris.lsmi.results <- data.frame(hclust_labels = iris.lsmi,
                                true_labels = iris.lsmi_true,
                                merged_hybrids = iris.lsmi_twoclasses,
                                merged_wrong_classes = iris.lsmi_twoclasses_wrong,
                                permuted_labels = iris.lsmi_wrong) %>% melt %>% rename(case = variable)

qplot(factor(case), value, data = iris.lsmi.results, 
      geom = c('boxplot'), color = case, xlab = 'case', ylab = 'lsmi estimate',
      main = 'Boxplot of 25 discrete-case LSMI values for iris data using labels obtained from 
      a) hierarchical clustering b) original data and c) permuted clustering labels') +
  ggsave('iris_lsmi_hclust_validation.png')

###
# looking how lsmi changes with number of clusters in iris data #
###

lsmi.vs.clustnum <- 
  data_frame(lsmi = rep(0, 25*5), k = rep(2:6, each = 25))
  
subs.data <- 1:25
for(clustnum in 2:6) {
  loc.clustLabs <- cutree(hc.iris, k = clustnum)
  lsmi.vs.clustnum[subs.data, 'lsmi'] <- replicate(25, lsmi.vanilla(iris4lsmi, loc.clustLabs, nbfuns = 150, y.discrete = TRUE))
  subs.data %<>% add(25)
}
rm(subs.data)
rm(loc.clustLabs)

lsmi.vs.clustnum %<>% mutate(lsmi_adj = lsmi / sqrt(k), k = as.factor(k))

qplot(k, lsmi_adj, data = lsmi.vs.clustnum, geom = 'boxplot', color = k,
      xlab = 'number of clusters used in validation', ylab = 'lsmi estimate / sqrt(n. of clusters)', 
      main = 'Boxplot of discrete-case lsmi distributions vs. number of iris data clusters supplied 
      clusters are obtained by cutting hierarchical clustering tree') + 
  labs(color = 'n. of clusters') +
  ggsave('iris_lsmi_vs_clustnum.png')

###
# testing how features correlate / lsmi-ate #
###

library(GGally)

png('corrplot_iris.png', width = 900, height = 900, pointsize = 3)
ggpairs(iris[, 1:4])
dev.off()

# comparing how lsmi can detect dependencies between the features analogous to correlation
## lower triangular matrix of between-features lsmi 
iris.lsmi.matrix <- matrix(0, nrow = 4, ncol = 4)
dimnames(iris.lsmi.matrix) <- list(colnames(iris)[1:4], colnames(iris)[1:4])
for(i in 1:4) {
  for(j in 1:i) {
    iris.lsmi.matrix[i, j] <- 
      replicate(10, lsmi.vanilla(iris[,i], iris[,j])) %>%
      median
  }
}

## plotting result
iris.lsmi.vs.cor <- data.frame(lsmi = iris.lsmi.matrix %>% extract(lower.tri(.)),
                               pcor = cor(iris[,1:4]) %>% extract(lower.tri(.)))
tmp.cor_pcor_lsmi <- cor.test(iris.lsmi.vs.cor[,1], iris.lsmi.vs.cor[,2])

ggplot(iris.lsmi.vs.cor, aes(lsmi, pcor)) +
  geom_point(color = 'cornflowerblue', size = 2) +
  stat_smooth(method = 'lm', col = 'aquamarine3', se = FALSE) +
  ggtitle(paste('Comparison of pearson correlation vs LSMI in identifying dependent features\n',
                'Fisher iris dataset;',
                'cor(pearson_rho, lsmi) = ', tmp.cor_pcor_lsmi$estimate %>% round(3), 
                '; p-value = ', tmp.cor_pcor_lsmi$p.value %>% round(3))) + 
  ggsave('iris_pcor_vs_lsmi.png')

###
# measuring accuracy depending on the (relative?) number of base functions
###

nbfun.prec.reps <- 25
nbfun.prec.step <- 10

nbfun.lsmi.values <- data.frame(lsmi = rep(0, nbfun.prec.reps*(150 %/% nbfun.prec.step)),
                                base.funs = rep(seq(10, 150, by = nbfun.prec.step), each = nbfun.prec.reps))

indices.vals <- 1:nbfun.prec.reps
for(i in seq(10, 150, by = nbfun.prec.step)) {
  nbfun.lsmi.values[indices.vals, 'lsmi'] <- 
    replicate(nbfun.prec.reps, lsmi.vanilla(iris4lsmi, iris$Species, nbfuns = i))
  indices.vals %<>% add(nbfun.prec.reps)
}
rm(indices.vals)

## mean squared error auxilary function
mse <- function(obs, eta) {
  require(magrittr)
  obs %>%
    subtract(eta) %>%
    {.^2} %>%
    divide_by(length(.)) %>%
    sum
}
rmse <- function(obs, eta) sqrt(mse(obs, eta))

nbfun.mse.values <- data_frame(base.funs = seq(10, 150, by = nbfun.prec.step))

# mean squarred error compared to the median LSMI computed using all possible base functions
## splitting results to feed them as lists to mse
nbfun.mse.values$mse <- 
  lapply(
    split(nbfun.lsmi.values$lsmi, rep(seq(10, 150, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    mse, eta = median(nbfun.lsmi.values$lsmi %>% tail(nbfun.prec.reps))) %>% as.numeric

nbfun.mse.values$var <- 
  lapply(
    split(nbfun.lsmi.values$lsmi, rep(seq(10, 150, by = nbfun.prec.step), each = nbfun.prec.reps)), 
    sd) %>%
  as.numeric #%>% .^2

## plotting mse
ggplot(nbfun.mse.values, aes(base.funs, mse)) +
  geom_point(color = 'firebrick2') + 
  labs(x = 'number of base functions considered',
       y = 'mean squared error') +
  ggtitle('Mean squared error w.r.t. median LSMI calculated using all 150 base functions
          Fisher iris data; y are labels denoting iris species') +
  ggsave('nbfuns_mse.png')

## plotting variance
ggplot(nbfun.mse.values, aes(base.funs, var)) +
  geom_point(color = 'darkslategray') + 
  labs(x = 'number of base functions considered',
       y = 'sd') +
  ggtitle('LSMI standard deviation depending on number of base functions considered
          Fisher iris data; y are labels denoting iris species') +
  ggsave('nbfuns_sd.png')

# ggplot(df.mse, aes(nbfuns, mse)) +
#   geom_point(color = 'firebrick2') + 
#   labs(x = 'number of base functions considered',
#        y = 'mean squared error') +
#   ggtitle(plotTitle)