
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
                                                           LSMI(as.list(
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

