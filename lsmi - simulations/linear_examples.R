#####################################
## experimental validation of LSMI ##
#####################################
library(magrittr)
library(ggplot2)
#library(tictoc)

#################################
### independent r.v.'s example ##
#################################

ind.ex1.x <- rnorm(50, 0, 4)
ind.ex1.y <- rnorm(50, 5, 1)

ind.ex1.cor <- cor(ind.ex1.x, ind.ex1.y) %>% round(3)
ind.ex1.lsmi <- lsmi.vanilla(as.list(ind.ex1.x), as.list(ind.ex1.y)) %>% round(3)

ggplot(data.frame(x = ind.ex1.x, y = ind.ex1.y), aes(x, y)) + 
  geom_point(color = 'orange', size = 2) + 
  ggtitle(paste('Independent r.v.s example:\n x ~ N(0, 4), y ~ N(5, 1) \n', 
                'Cor = ', ind.ex1.cor, '; ', 'LSMI = ', ind.ex1.lsmi, sep = '')) + 
  theme(plot.background = element_rect(fill = "darkseagreen"), plot.title = element_text(hjust = 0.5)) +
  ggsave('indep_example_normal.png')

replicate(50, lsmi.vanilla(as.list(ind.ex1.x), as.list(ind.ex1.y))) %>% table

#######################################
### how lsmi depends on correlation ###
#######################################

library(mvtnorm)
mvtsigma.list <- list()

## various parameters for this section ##
lin.step <- 21
lin.dots <- 50
lin.mean <- c(0, 0)
lin.rep.lsmi <- 10
##

for(i in 1:lin.step) mvtsigma.list[[i]] <- matrix(c(1, (i-1)/(lin.step - 1), 
                                                    (i-1)/(lin.step - 1), 1), 
                                                  nrow = 2)

# preparing data frame for filling with data
cor.vs.lsmi.linear <- data.frame(x = rep(0, lin.dots*lin.step),
                                 y = rep(0, lin.dots*lin.step), 
                                 cor = rep(seq(0, 1, length.out = lin.step), each = lin.dots)
                                 )

# generating r.v.s
indices.data <- c(1, lin.dots)

for(i in 1:lin.step) {
  subs.data <- indices.data[1]:indices.data[2] 
  cor.vs.lsmi.linear[subs.data, 1:2] <- rmvnorm(lin.dots, lin.mean, mvtsigma.list[[i]])
  indices.data %<>% add(lin.dots)
}
rm(indices.data)
rm(subs.data)

lin.lsmi.values <- data.frame(lsmi = rep(0, lin.rep.lsmi*lin.step), 
                              cor = rep(seq(0, 1, length.out = lin.step), each = lin.rep.lsmi),
                              # correlation estimate
                              corhat = rep(0, lin.rep.lsmi*lin.step),
                              corhat.cimin = rep(0, lin.rep.lsmi*lin.step),
                              corhat.cimax = rep(0, lin.rep.lsmi*lin.step)
                              )

# calculating LSMI values
indices.data <- c(1, lin.dots)
indices.vals <- c(1, lin.rep.lsmi)

for(i in 1:lin.step) {
  subs.data <- indices.data[1]:indices.data[2]
  subs.vals <- indices.vals[1]:indices.vals[2]
  
  lin.lsmi.values[subs.vals,'lsmi'] <- replicate(lin.rep.lsmi,
                                        
                                        lsmi.vanilla(as.list(cor.vs.lsmi.linear$x[subs.data]),
                                             as.list(cor.vs.lsmi.linear$y[subs.data])))
  
  current.cortest <- cor.test(cor.vs.lsmi.linear$x[subs.data],
                              cor.vs.lsmi.linear$y[subs.data])
  
  lin.lsmi.values[subs.vals, 'corhat'] <- current.cortest$estimate
  lin.lsmi.values[subs.vals, c('corhat.cimin', 'corhat.cimax')] <- current.cortest$conf.int %>% 
                                                                   rep(., lin.rep.lsmi) %>%
                                                                   matrix(., ncol = 2, byrow = TRUE) %>%
                                                                   as.data.frame

  indices.data %<>% add(lin.dots)
  indices.vals %<>% add(lin.rep.lsmi)
}
rm(indices.data)
rm(indices.vals)
rm(subs.data)
rm(current.cortest)

# plotting results

## dependence of LSMI / Cor on true cor
library(reshape2)
lin.lsmi.values %<>% melt(., id = c('cor', 'corhat.cimin', 'corhat.cimax'))

library(ggplot2)
ggplot(data = lin.lsmi.values, aes(cor, value)) + 
  geom_point(aes(color = factor(variable)), size = 3, alpha = 0.3) + 
  geom_errorbar(aes(ymin = corhat.cimin, ymax = corhat.cimax, colour = factor(variable)), width = 0.025) + 
  scale_color_hue(labels = c("LSMI estimate", "Cor estimate")) +
  
  theme(plot.background = element_rect(fill = "aliceblue"), legend.title = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, 16), alpha = 1))) + 
  ggtitle(paste('LSMI / Cor estimates (95% CI) vs. true Cor\n bivariate normal distribution, n = ', lin.dots, ';', sep = '')) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  
  ggsave('corhat_vs_lsmi.png')

# zoomed version
ggplot(data = lin.lsmi.values, aes(cor, value)) + 
  coord_cartesian(ylim = c(0, 1)) +
  geom_point(aes(color = factor(variable)), size = 3, alpha = 0.3) + 
  geom_errorbar(aes(ymin = corhat.cimin, ymax = corhat.cimax, colour = factor(variable)), width = 0.025) + 
  geom_abline(slope = 1, intercept = 0, color = 'blue', size = 1) + 
  
  theme(plot.background = element_rect(fill = "aliceblue"), legend.title = element_blank()) + 
  scale_color_hue(labels = c("LSMI estimate", "Cor estimate")) +
  guides(colour = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, 16), alpha = 1))) + 
  ggtitle(paste('LSMI / Cor estimates (95% CI) vs. true Cor\n bivariate normal distribution, n = ', lin.dots, ';', sep = '')) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  
  ggsave('corhat_vs_lsmi_zoom.png')

# Mean squared error vs. dependence strength  #

mse <- function(obs, eta) {
  require(magrittr)
  obs %>%
    subtract(eta) %>%
    {.^2} %>%
    mean
}

## Function for calculaing MSE vs. nbfuns
MSEvsNbfuns <- function(x, y, nrep = 25, nstep = 10, ydisc = FALSE) {
  require(magrittr)
  require(dplyr)
  if(!is.list(x)) x %<>% as.list
  if(!is.list(y)) y %<>% as.list
  
  maxnbfuns <- length(x)
  nbfun.seq <- seq(10, maxnbfuns, by = nstep)
  
  df.vals <- data_frame(nbfuns = nbfun.seq %>% rep(each = nrep), 
                        lsmi = numeric(length(nbfun.seq)*nrep))
  
  df.mse <- data_frame(nbfuns = nbfun.seq, 
                       mse = numeric(length(nbfun.seq)))
  
  indices.vals <- 1:nrep
  for(i in seq(10, maxnbfuns, by = nstep)) {
    df.vals$lsmi[indices.vals] <- replicate(nrep, lsmi.vanilla(x, y, nbfuns = i, y.discrete = ydisc))
    indices.vals %<>% add(nrep)
  }
  
  df.mse$mse <- lapply(split(df.vals$lsmi, rep(seq(10, maxnbfuns, by = nstep), each = nrep)), 
                       mse, 
                       eta = median(df.vals$lsmi %>% tail(nrep))) %>% as.numeric
  df.mse
}

# Generating multivariate normal data with different correlation strengths
library(mvtnorm)
library(stringr)
listXY <- lapply(mvtsigma.list[seq(1, 21, by = 5)], rmvnorm, n = 150, mean = c(0, 0))

for(i in 1:5) {
  colnames(listXY[[i]]) <- 
    mvtsigma.list[seq(1, 21, by = 5)] %>%
    sapply(extract, 1, 2) %>%
    extract(i) %>%
    rep(2)
}

MSEvsDepStr <- function(xymat) {
  require(magrittr)
  require(stringr)
  x <- xymat[, 1]
  y <- xymat[, 2]
  sigma <- 
    colnames(xymat) %>%
    extract(1)
  
  df.lsmi <- MSEvsNbfuns(x = x, y = y)
  df.lsmi %<>% mutate(depstr = rep(sigma, nrow(df.lsmi)))
  df.lsmi
}

#listByDepStr <- lapply(listXY, MSEvsDepStr)
df.depStr <- 
  do.call(rbind, listByDepStr) %>%
  filter(depstr != 0) %>%
  mutate(depstr = as.factor(depstr), logMse = log(mse, 10), rmse = sqrt(mse), logRMse = log(rmse, 10))

ggplot(df.depStr, aes(nbfuns, logRMse)) +
  geom_point(aes(color = depstr), size = 3) + 
  labs(x = 'number of base functions considered',
       y = 'log10(RMSE)') +
  geom_smooth(method = 'loess', se = F, aes(color = depstr)) +
  scale_x_continuous(breaks = seq(0, 150, by = 10)) +
  scale_y_continuous(breaks = unique(c(seq(0, floor(min(df.depStr$logRMse)), by = -0.2), seq(0, ceiling(max(df.depStr$logRMse)), by = 0.2)))) + 
  scale_colour_discrete(name = 'Pearson\ncorrelation') +
  ggtitle('RMSE vs. number of base functions considered for RBF MI') +
  ggsave('nbfuns_depstr_log.png')

df.depStr %<>% filter(depstr != 1)

ggplot(df.depStr, aes(nbfuns, rmse)) +
  geom_point(aes(color = depstr), size = 3) + 
  labs(x = 'number of base functions considered',
       y = 'RMSE') +
  geom_smooth(method = 'loess', se = F, aes(color = depstr)) +
  scale_x_continuous(breaks = seq(0, 150, by = 10)) +
  scale_y_continuous(breaks = seq(0, ceiling(max(df.depStr$rmse)), by = 0.02)) + 
  scale_colour_discrete(name = 'Pearson\ncorrelation') +
  ggtitle('RMSE vs. number of base functions considered for RBF MI') +
  ggsave('nbfuns_depstr.png')


## parametric graph (LSMI, corhat) as a function of true cor
## is most likely obsolette

# library(ggplot2)
# ggplot(data = lin.lsmi.values, aes(corhat, lsmi)) +
#   geom_point(color = 'blueviolet', size = 3, alpha = 0.3) +
#   geom_errorbarh(aes(xmin = corhat.cimin, xmax = corhat.cimax), height = 0.02, colour = 'blueviolet') +
# #  geom_jitter(position = position_jitter(width = 0.04), colour = 'blueviolet', size = 3, alpha = 0.3) +
#   theme(plot.background = element_rect(fill = "aliceblue")) +
# #  stat_smooth(se = FALSE) +
#   ggtitle(paste('Plot of Pearson correlation estimates (95% CI) vs. LSMI
#           bivariate normal distribution, n = ', lin.dots, ' for each unique cor', sep = '')) +
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   
#   ggsave('cor_vs_lsmi.png')
# 
# ggplot(data = lin.lsmi.values, aes(corhat, lsmi)) +
#   coord_cartesian(ylim = c(0, 1)) +
#   geom_point(colour = 'blueviolet', size = 3, alpha = 0.3) +
#   geom_errorbarh(aes(xmin = corhat.cimin, xmax = corhat.cimax), height = 0.02, colour = 'blueviolet') +
#   #  stat_smooth(se = FALSE) +
#   # geom_jitter(position = position_jitter(width = 0.04), colour = 'blueviolet', size = 3, alpha = 0.3) +
#   
#   theme(plot.background = element_rect(fill = "aliceblue")) +
#   ggtitle(paste('Plot of Pearson correlation estimates (95% CI) vs. LSMI
#           bivariate normal distribution, n = ', lin.dots, ' for each unique cor\n zoomed in at LSMI in [0;1]', sep = '')) +
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   
#   ggsave('cor_vs_lsmi_zoom.png')

########################################
###  simplest linear dependence case ###
########################################

slin.step <- 21
slin.dots <- 50
slin.rep.lsmi <- 10

slin.slope <- 3
slin.intercept <- 1

slin.data2plot <- data.frame(x = runif(slin.dots, -1, 1), y = rep(0, slin.dots))
slin.data2plot$y <- slin.slope*slin.data2plot$x + slin.intercept + rnorm(slin.dots, 0, 2)
slin.cor2plot <- cor(slin.data2plot$x, slin.data2plot$y) %>% round(3)
slin.lsmi2plot <- lsmi.vanilla(as.list(slin.data2plot$x), as.list(slin.data2plot$y)) %>% round(3)


ggplot(slin.data2plot, aes(x, y)) + 
  geom_point(color = 'orange', size = 2) + 
  ggtitle(paste('Simplest linear case
                x ~ U[-1; 1]; y = ', slin.slope, '*x + ', slin.intercept, ' + N(0, sigma = 2) ; n = ', slin.dots, ';\n',
                'Cor = ', slin.cor2plot, '; LSMI = ', slin.lsmi2plot, sep = '')) + 
  theme(plot.background = element_rect(fill = "darkseagreen"), plot.title = element_text(hjust = 0.5)) +
  ggsave('linear_example.png')


# generating r.v.s
slin.lsmi.values <- data.frame(lsmi = rep(0, slin.rep.lsmi*slin.step),
                               noise = rep(seq(0, 5, length.out = slin.step), each = slin.rep.lsmi),
                               # correlation estimate
                               corhat = rep(0, slin.rep.lsmi*slin.step),
                               corhat.cimin = rep(0, slin.rep.lsmi*slin.step),
                               corhat.cimax = rep(0, slin.rep.lsmi*slin.step)
)

slin.lsmi <- data.frame(x = rep(x = runif(slin.dots, -1, 1), slin.step),
                               y = rep(0, slin.dots*slin.step),
                               noise = rep(seq(0, 5, length.out = slin.step), each = slin.dots)
)

slin.lsmi$y <- slin.slope*slin.lsmi$x + slin.intercept + sapply(as.list(slin.lsmi$noise), rnorm, n = 1, mean = 0)

# calculating LSMI
indices.vals <- c(1, slin.rep.lsmi)
indices.data <- c(1, slin.dots)

for(i in 1:slin.step) {
  subs.vals <- indices.vals[1]:indices.vals[2]
  subs.data <- indices.data[1]:indices.data[2]
  
  slin.lsmi.values[subs.vals, 'lsmi'] <- replicate(slin.rep.lsmi,
                                                   lsmi.vanilla(as.list(slin.lsmi[subs.data, 'x']), 
                                                        as.list(slin.lsmi[subs.data, 'y']))
  )
  
  current.cortest <- cor.test(slin.lsmi$x[subs.data],
                              slin.lsmi$y[subs.data])
  
  slin.lsmi.values[subs.vals, 'corhat'] <- current.cortest$estimate
  slin.lsmi.values[subs.vals, c('corhat.cimin', 'corhat.cimax')] <- current.cortest$conf.int %>% 
    rep(., slin.rep.lsmi) %>%
    matrix(., ncol = 2, byrow = TRUE) %>%
    as.data.frame
  
  indices.vals %<>% add(slin.rep.lsmi)
  indices.data %<>% add(slin.dots)
}
rm(indices.vals)
rm(indices.data)
rm(subs.vals)
rm(current.cortest)

## plotting results

library(reshape2)
slin.lsmi.values %<>% melt(., id = c('noise', 'corhat.cimin', 'corhat.cimax'))

library(ggplot2)
#library(ggforce)
ggplot(data = slin.lsmi.values, aes(noise, value)) + 
  geom_hline(yintercept = 1, color = 'darkgray', size = 1) +
  geom_point(aes(color = factor(variable)), size = 3, alpha = 0.3) + 
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 3)) +
  geom_errorbar(aes(ymin = corhat.cimin, ymax = corhat.cimax, colour = factor(variable)), width = 0.042) + 

  scale_color_hue(labels = c("LSMI estimate", "Cor estimate")) +
  theme(plot.background = element_rect(fill = "aliceblue"), legend.title = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, 16), alpha = 1))) + 
  ggtitle(paste('LSMI / Cor estimates (95% CI) vs. noise\n x ~ U[-1; 1]; y =', 
                slin.slope, '*x + ', slin.intercept, ' + N(0, sigma = noise)', '; n = ', slin.dots, ';', sep = '')) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  #facet_zoom(y = value < 1 & value > 0, x = noise > 0.5 & noise < 1.5, horizontal = FALSE) +
  
  ggsave('simplest_linear.png')


