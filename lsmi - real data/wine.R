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
miFlowPlot(wqFeaSelBW, colnames(wqSbsScaled[, -c(13, 14)]), nfeat = 2) + 
  ggtitle('Backward feature elimination for wine classification using LSMI') +
  ggsave('MI_backward_selection.png', width = 10, height = 4)

### Forward selection
miFlowPlot(wqFeaSelFW, nfeat = 5) + 
  ggtitle('Forward feature selection for wine classification using LSMI') +
  ggsave('MI_forward_selection.png', width = 10, height = 4)

###
# fitting logistic regression ##
###

wineLogReg_full <- glm(data = select(wqSbsScaled, -id), formula = factor(wine) ~ ., family = 'binomial')
wineLogReg_bw <- glm(data = select(wqSbsScaled, -id), 
                     formula = factor(wine) ~ `residual sugar` + density + `total sulfur dioxide`, family = 'binomial')
wineLogReg_fw <- glm(data = select(wqSbsScaled, -id),
                     formula = factor(wine) ~ `total sulfur dioxide` + chlorides + `free sulfur dioxide` + density + `residual sugar`, family = 'binomial')

wineGrid <- 
  expand.grid(chlorides = seq(0, 1, by = 0.05), `total sulfur dioxide` = seq(0, 1, by = 0.05)) %>%
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

wineLogRTest_full <- predict(wineLogReg_full, wineTestScaled, type = 'response')
wineLogRTest_bw <- predict(wineLogReg_bw, wineTestScaled, type = 'response')
wineLogRTest_fw <- predict(wineLogReg_fw, wineTestScaled, type = 'response')

###
# ggROC
###

wq.df.roc <- 
  data_frame(D = as.numeric(wineTestScaled$wine) - 1,
             D.str = as.character(wineTestScaled$wine),
             Full = wineLogRTest_full,
             `Backward feat. elimination` = wineLogRTest_bw,
             `Forward feat. selection` = wineLogRTest_fw) #%>%
library(reshape2)
wq.df.roc.melt <- melt(wq.df.roc, id.vars = c('D', 'D.str'), variable.name = 'Model', value.name = 'M')
library(caTools)
wq.auc <- colAUC(wq.df.roc[,3:5], wq.df.roc$D.str)

ggplot(wq.df.roc.melt, aes(d = D, m = M, color = Model)) + 
  geom_roc(linealpha = 0.8, labelround = 3) +
  style_roc() +
  labs(caption = str_c('AUC values are ', str_c(round(wq.auc, 4), collapse = ', '), ', respectively.')) +
  ggsave('roc_for_wine.png', width = 8, height = 6)

####
## Confusion matrices for these three classifiers ##
####
library(caret)
confusionMatrix(ifelse(wineLogRTest_full > 0.5, 'white', 'red'), wineTestScaled$wine)
confusionMatrix(ifelse(wineLogRTest_bw > 0.5, 'white', 'red'), wineTestScaled$wine)
confusionMatrix(ifelse(wineLogRTest_fw > 0.5, 'white', 'red'), wineTestScaled$wine)

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
