
## miscellanous tests and experiments related to LSMI ##

##########################################################
### how lsmi computational time depends on sample size ###
##########################################################

# generating data

compt.rep.lsmi <- 10
compt.step.lsmi <- 10

compt.x <- rnorm(100, 0, 5)
compt.y <- compt.x^2 + rnorm(100, 0, 2)
compt.data <- data.frame(x = compt.x,
                         y = compt.y)
rm(compt.x)
rm(compt.y)

compt.lsmi.values <- data.frame(cmpt.time = rep(0, compt.rep.lsmi*compt.step.lsmi),
                                sampl.size = rep(seq(10, 100, by = compt.step.lsmi), each = compt.rep.lsmi))

indices.vals <- c(1, compt.rep.lsmi)
subs.data <- compt.step.lsmi

for(i in 1:compt.step.lsmi) {
  subs.vals <- indices.vals[1]:indices.vals[2]
  
  compt.lsmi.values[subs.vals, 'cmpt.time'] <- replicate(compt.rep.lsmi,
                                                         system.time(
                                                           lsmi.vanilla(as.list(
                                                             sample(compt.data$x, subs.data)
                                                           ),
                                                           as.list(
                                                             sample(compt.data$y, subs.data)
                                                           )
                                                           )
                                                         )[[3]]
  )
  
  indices.vals %<>% add(compt.rep.lsmi)
  subs.data %<>% add(compt.step.lsmi)
}
rm(subs.data)
rm(indices.vals)

library(ggplot2)
ggplot(data = compt.lsmi.values, aes(sampl.size, cmpt.time)) + 
  geom_point(color = 'red', size = 3, alpha = 0.3) +
  geom_smooth(method = 'lm',formula = y ~ poly(x, 2), colour = 'darkorange') + 
  geom_jitter(position = position_jitter(width = 1.5), colour = 'red', size = 3, alpha = 0.3) + 
  
  theme(plot.background = element_rect(fill = "lightgoldenrod1")) +
  ggtitle(paste('Plot of computation time dependence on sample size
                fitted quadratic model; repetitions = ', compt.rep.lsmi, ';', sep = '')) + 
  labs(x = 'sample size', y = 'computation time, sec') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  
  ggsave('time_vs_samplesize.png')

# how lsmi computational time depends on a number of basis functions #
## this can provide insight on what's the t(nbfuns) dependence type is and
## provide info which will help to decide how much basis functions we're going to set
## to save computational time;

## note: 'step' in code above (how much iterations to go through) means different thing
## than here (single step in a number of basis functions)
iris4lsmi <- split(as.matrix(iris[, -5]), row(iris[,-5])) ## in case ydiscr.R has not been sourced

compt.rep.bf.lsmi <- 25
compt.step.bf.lsmi <- 10

compt.bf.lsmi.values <- data.frame(cmpt.time = rep(0, compt.rep.bf.lsmi*(150 %/% compt.step.bf.lsmi)),
                                   basis.funs = rep(seq(10, 150, by = compt.step.bf.lsmi), each = compt.rep.bf.lsmi)
                                   )

indices.vals <- 1:compt.rep.bf.lsmi
for(i in seq(10, 150, by = compt.step.bf.lsmi)) {
  
  compt.bf.lsmi.values[indices.vals, 'cmpt.time'] <- replicate(compt.rep.bf.lsmi,
                                                                system.time(lsmi.vanilla(x = iris4lsmi, y = iris[5] %>% unlist,
                                                                                         y.discrete = TRUE, nbfuns = i)
                                                                            )[[3]]
                                                               )
  indices.vals %<>% add(compt.rep.bf.lsmi)
}
rm(indices.vals)

ggplot(data = compt.bf.lsmi.values, aes(basis.funs, cmpt.time)) + 
  geom_point(color = 'firebrick2') +
  stat_smooth(method = 'lm', formula = y ~ poly(x, 2), color = 'forestgreen', se = FALSE) +
  labs(x = 'number of base functions',
       y = 'computation time, sec') +
  ggtitle('LSMI computational time vs. number of base functions used
          fitted quadratic model; iris data, n = 150; dim(x) = 4; y are labels denoting iris species') +
  # geom_jitter(height = 0, width = 2) +
  ggsave('time_vs_nbfuns.png')


