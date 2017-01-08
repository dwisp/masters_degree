#######################################################
## tests performed on classification subtype of LSMI ##
#######################################################

library(magrittr)
##########################
# two normal distributions
# what's the typical dispersion value for lsmi in discrete case?
# will it be able to differentiate original case from  one with permuted labels?

tst.disy.x1 <- as.list(c(rnorm(25, 2, 2), rnorm(25, -2, 2)))
tst.disy.y1 <- as.list(c(rep('a', 25), rep('b', 25)))
tst.disy.y1_wrong <- tst.disy.y1 %>% unlist %>% sample %>% as.list

tst.disy.lsmi1 <- replicate(50, lsmi.vanilla(tst.disy.x1, tst.disy.y1, y.discrete = T))
tst.disy.lsmi1_wrong <- replicate(50, lsmi.vanilla(tst.disy.x1, tst.disy.y1_wrong, y.discrete = T))

df.disy.lsmi_comparison <- data.frame(original_labels = tst.disy.lsmi1, 
                                      permuted_labels = tst.disy.lsmi1_wrong)
library(reshape2)
library(dplyr)
df.disy.lsmi_comparison %<>% melt %>% rename(case = variable)

library(ggplot2)
ggplot(df.disy.lsmi_comparison, aes(case, value)) +
  geom_boxplot(aes(color = case)) +
  ggtitle('Boxplot of discrete-case LSMI values for two normal distributions N(-2, 2) / N(2, 2)
        right boxplot shows the result of calculating LSMI using permuted labels') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggsave('discrete_lsmi_normal_distrs.png')

##########################
## testing discrete lsmi on iris dataset to validate clustering results
## core idea is to check if we'll be able to get noticeably higher lsmi values when
## providing 'right' clustering labels to the method


library(dendextend)
hc.iris <- iris[, -5] %>% as.matrix %>% dist %>% hclust(method = 'ward.D2') 
iris.cols <- iris$Species %>% as.factor
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

hc.iris.labels <- cutree(hc.iris, k = 3) %>% as.list
hc.iris.labels_wrong <- hc.iris.labels %>% unlist %>% sample %>% as.list

iris4lsmi <- split(as.matrix(iris[, -5]), row(iris[,-5]))

iris.lsmi <- replicate(25, lsmi.vanilla(iris4lsmi, hc.iris.labels, y.discrete = T),)
iris.lsmi_true <- replicate(25, lsmi.vanilla(iris4lsmi, iris[5] %>% unlist %>% as.character %>% as.list, y.discrete = T))
iris.lsmi_wrong <- replicate(25, lsmi.vanilla(iris4lsmi, hc.iris.labels_wrong, y.discrete = T))

iris.lsmi.results <- data.frame(hclust_labels = iris.lsmi,
                                true_labels = iris.lsmi_true,
                                permuted_labels = iris.lsmi_wrong) %>% melt %>% rename(case = variable)

qplot(factor(case), value, data = iris.lsmi.results, 
      geom = c('boxplot'), color = case, xlab = 'case', 
      main = 'Boxplot of 25 discrete-case LSMI values for iris data using labels obtained from 
      a) hierarchical clustering b) original data and c) permuted clustering labels') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggsave('iris_lsmi_hclust_validation.png')

## testing how features correlate / lsmi-ate
plotmatrix(iris[, 1:4])







