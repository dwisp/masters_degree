##
# trying to build a wrapper method using LSMI for wine data
##

library(readr)
library(stringr)
library(dplyr)
wqRed <- 
  read_csv('C:/Users/Ivan/Documents/R/experiments/datasets/food-wine-quality/winequality-red.csv') %>%
  na.omit() %>%
  mutate(id = str_c('red', 1:nrow(.)))
wqWhite <- 
  read_csv('C:/Users/Ivan/Documents/R/experiments/datasets/food-wine-quality/winequality-white.csv') %>%
  mutate(id = str_c('white', 1:nrow(.)))

wqFull <- 
  full_join(wqRed %>% mutate(wine = 'red'), 
            wqWhite %>% mutate(wine = 'white'))

set.seed(451)
wqRedSbs <- 
  sample_n(wqRed, 500) %>%
  mutate(wine = 'red')
wqWhiteSbs <- 
  sample_n(wqWhite, 1500) %>%
  mutate(wine = 'white')

wqSbs <- full_join(wqRedSbs, wqWhiteSbs)

# wine.MI <- numeric(12)
# wine.cor <- numeric(12)
# 
# names(wine.MI) <- colnames(wqSbs)[-13]
# names(wine.cor) <- colnames(wqSbs)[-13]
# 
# for(i in 1:12) {
#   wine.MI[i] <- lsmi.extra(wqSbs %>% select(i) %>% unlist(),
#                            factor(wqSbs$wine), nbfuns = 200,
#                            method.nbfuns = 'balanced', method = 'suzuki')
#   wine.cor[i] <- cor(wqSbs %>% select(i) %>% unlist(), 
#                      wqSbs$wine %>% factor() %>% as.numeric())
# }
# 
# # doesn't work! yay
# # library(FNN)
# # mutinfo(wqSbs %>% select(5) %>% unlist(), wqSbs$wine %>% factor() %>% as.numeric())
# 
# ggplot(data = data_frame(MI = wine.MI, feature = str_replace_all(names(wine.MI), ' ', '/n')), aes(feature, MI)) +
#   geom_col(fill = 'lightskyblue', alpha = 0.7) +
#   ggsave('MI_feature_importance.png', width = 8, height =  4)

wqSbsScaled <-
  wqSbs %>%
  select(-wine, -id) %>%
  mutate_all(funs((. - min(.))/(max(.) - min(.)))) %>%
  cbind(wine = wqSbs$wine, id = wqSbs$id) %>%
  as_data_frame()

source('./lsmi - method/LSMI_feature_selection.R')
wqFeaSelBW <- miFeaSelect(wqSbsScaled[, -c(13, 14)], wqSbsScaled$wine, selection = 'backward')
wqFeaSelFW <- miFeaSelect(wqSbsScaled[, -c(13, 14)], wqSbsScaled$wine, selection = 'forward')

### Backward selection
ggplot(wqFeaSelBW, aes(x = iteration, y = miIter, size = importance, label = feaIter)) +
  geom_line(color = 'lightcyan3', alpha = 0.8) +
  geom_point(shape = 18) +
  geom_text(color = 'palevioletred4', size = 3, vjust = 2) +
  ylab('Mutual Information') +
  scale_x_continuous(breaks = 0:11, limits = c(0, 11)) +
  scale_y_continuous(breaks = seq(round(min(wqFeaSelBW$miIter) - 0.05, 2), round(max(wqFeaSelBW$miIter) + 0.05, 2), by = 0.02)) +
  scale_size_area(limits = c(0, 0.8), breaks = seq(0, 0.6, by = 0.2)) +
  theme(legend.position = c(0.85, 0.2)) +
  guides(size = guide_legend(title = 'Importance\n(marginal MI)', override.aes = list(lty = NULL), ncol = 2)) +
  ggsave('MI_backward_selection.png')

### Forward selection
ggplot(wqFeaSelFW[-13, ], aes(x = iteration, y = miIter, size = importance, label = feaIter)) +
  geom_line(color = 'lightcyan3', alpha = 0.8) +
  geom_point(shape = 18) +
  geom_text(color = 'palevioletred4', size = 3, vjust = 2) +
  ylab('Mutual Information') +
  scale_x_continuous(breaks = 0:12, limits = c(0, 12)) +
  scale_y_continuous(breaks = seq(round(min(wqFeaSelFW$miIter) - 0.05, 2), round(max(wqFeaSelFW$miIter) + 0.05, 2), by = 0.02)) +
  scale_size_area(limits = c(0, 0.8), breaks = seq(0, 0.6, by = 0.2)) +
  theme(legend.position = c(0.85, 0.2)) +
  guides(size = guide_legend(title = 'Importance\n(marginal MI)', override.aes = list(lty = NULL), ncol = 2)) +
  ggsave('MI_forward_selection.png')

# irisFeaSelBW <- miFeaSelect(as_data_frame(iris[, -5]), iris$Species, selection = 'backward')
# irisFeaSelFW <- miFeaSelect(iris[, -5], iris$Species, selection = 'forward')

###
# fitting logistic regression ##
###

wineLogReg_full <- glm(data = select(wqSbsScaled, -id), formula = factor(wine) ~ ., family = 'binomial')
wineLogReg_bw <- glm(data = select(wqSbsScaled, -id), 
                      formula = factor(wine) ~ `volatile acidity` + `residual sugar` + density + `total sulfur dioxide`, 
                      family = 'binomial')
wineLogReg_fw <- glm(data = select(wqSbsScaled, -id),
                          formula = factor(wine) ~ chlorides + `total sulfur dioxide`, 
                          family = 'binomial')

wineGrid <- 
  expand.grid(chlorides = seq(0, 1, by = 0.05), 
              `total sulfur dioxide` = seq(0, 1, by = 0.05)) %>%
  as_data_frame()
wineGrid$z <- predict(wineLogReg_bw, newdata = wineGrid, type = 'response')

set.seed(1984)
wineTest <- sample_n(filter(wqFull, !id %in% wqSbs$id), 1200)
wineTestScaled <-
  wineTest %>%
  select(-wine, -id) %>%
  mutate_all(funs((. - min(.))/(max(.) - min(.)))) %>%
  cbind(wine = wineTest$wine) %>%
  as_data_frame()

wineLogRTest <- predict(wineLogReg_bw, newdata = wineTestScaled, type = 'response')
wineLogRTest01 <- as.numeric(wineLogRTest > 0.5)

## red = 0, white = 1
wineLogRerrors <- (as.numeric(factor(wineTest$wine)) - 1) != wineLogRTest01
wineConfMat <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c('pred Red', 'pred White'), c('Red', 'White')))

wineConfMat[1, 1] <- sum(wineLogRTest01 == 0 & wineTest$wine == 'red')
wineConfMat[2, 2] <- sum(wineLogRTest01 == 1 & wineTest$wine == 'white')
wineConfMat[1, 2] <- sum(wineLogRTest01 == 0 & wineTest$wine == 'white')
wineConfMat[2, 1] <- sum(wineLogRTest01 == 1 & wineTest$wine == 'red')

print(str_c('Accuracy = ', round(sum(diag(wineConfMat)) / sum(wineConfMat), 3)))
print(str_c('Red Sensitivity = ', round(sum(wineConfMat[1, 1]) / sum(wineConfMat[, 1]), 3)))
print(str_c('White Sensitivity = ', round(sum(wineConfMat[2, 2]) / sum(wineConfMat[, 2]), 3)))

library(ggthemes)
library(RColorBrewer)
library(forcats)
ggplot() +
  geom_tile(data = wineGrid, aes(x = `total sulfur dioxide`, y = chlorides, fill = z), alpha = 0.5) +
  scale_fill_gradientn(colors = brewer.pal(6, 'RdYlBu'), breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  geom_point(data = wqSbsScaled, aes(x = `total sulfur dioxide`, y = chlorides, color = fct_rev(wine)), alpha = 0.3, size = 2.5) +
  theme_solarized() +
  scale_color_solarized() +
  guides(color = guide_legend(title = 'Wine type', override.aes = list(size = 4, alpha = 1), order = 1),
         fill = guide_colorbar(title = 'Logit score\nP{wine is white}\n', barwidth = 2.4, order = 2)) +
  ggsave('wineLogRegR2d.png')








