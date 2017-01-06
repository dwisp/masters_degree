
###################################
## sine function, matlab example ##
###################################

nl.ex1.x <- runif(50, -1, 1)
nl.ex1.y <- sin(nl.ex1.x*4*pi) + rnorm(50, 0, 0.25)

nl.ex1.cor <- cor(nl.ex1.x, nl.ex1.y) %>% round(3)
nl.ex1.lsmi <- lsmi.vanilla(as.list(nl.ex1.x), as.list(nl.ex1.y)) %>% round(3)

# plotting dependence #
library(ggplot2)
ggplot(data.frame(x = nl.ex1.x, y = nl.ex1.y), aes(x, y)) + 
  geom_point(color = 'orange', size = 2) + 
  ggtitle(paste('Nonlinear dependence example 1: high-frequency sine function\n x ~ U[-1;1]; y = sin 4*pi*x + N(0, 0.25); n = 50\n',
                'Cor = ', nl.ex1.cor, '; ', 'LSMI = ', nl.ex1.lsmi, sep = '')) +
  theme(plot.background = element_rect(fill = "darkseagreen"), plot.title = element_text(hjust = 0.5)) +
  ggsave('sine_example.png')

# lsmi.reps.nl.ex1 <- replicate(50, lsmi.vanilla(as.list(nl.ex1.x), as.list(nl.ex1.y)), TRUE)
# table(lsmi.reps.nl.ex1)

#################################
###        a parabola         ###
#################################

nl.ex2.x <- runif(50, -1, 1)
nl.ex2.y <- nl.ex2.x^2 + rnorm(50, 0, 0.25)

nl.ex2.cor <- cor(nl.ex2.x, nl.ex2.y) %>% round(3)
nl.ex2.lsmi <- lsmi.vanilla(as.list(nl.ex2.x), as.list(nl.ex2.y)) %>% round(3)

ggplot(data.frame(x = nl.ex2.x, y = nl.ex2.y), aes(x, y)) + 
  geom_point(color = 'orange', size = 2) + 
  ggtitle(paste('Nonlinear dependence example 2: parabola\n x ~ U[-1;1]; y = x^2 + N(0, 0.25); n = 50\n',
                'Cor = ', nl.ex2.cor, '; ', 'LSMI = ', nl.ex2.lsmi, sep = '')) +
  theme(plot.background = element_rect(fill = "darkseagreen"), plot.title = element_text(hjust = 0.5)) +
  ggsave('parabola_example.png')

# lsmi.reps.nl.ex2 <- replicate(50, lsmi.vanilla(as.list(nl.ex2.x), as.list(nl.ex2.y)), TRUE)
# table(lsmi.reps.nl.ex2)

#################################
###         a circle          ###
#################################

nl.ex3.x <- runif(50, -1, 1)
nl.ex3.y <- sample(c(-1,1), 50, T, c(1/2,1/2))*sqrt(1 - nl.ex3.x^2) + rnorm(50, 0, 0.25)

nl.ex3.cor <- cor(nl.ex3.x, nl.ex3.y) %>% round(3)
nl.ex3.lsmi <- lsmi.vanilla(as.list(nl.ex3.x), as.list(nl.ex3.y)) %>% round(3)

ggplot(data.frame(x = nl.ex3.x, y = nl.ex3.y), aes(x, y)) + 
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1.1, 1.1)) +
  geom_point(color = 'orange', size = 2) +
  ggtitle(paste('Nonlinear dependence example 3: unit circle\n x ~ U[-1;1]; y = random(-1;1)*sqrt(1 - x^2) + N(0, 0.25), n = 50\n',
                'Cor = ', nl.ex3.cor, '; ', 'LSMI = ', nl.ex3.lsmi, sep = '')) +
  theme(plot.background = element_rect(fill = "darkseagreen"), plot.title = element_text(hjust = 0.5)) +
  ggsave('circle_example.png')#, width = 9, height = 9)

# lsmi.reps.nl.ex3 <- replicate(50, lsmi.vanilla(as.list(nl.ex3.x), as.list(nl.ex3.y)), TRUE)
# table(lsmi.reps.nl.ex3)


################################################################
### how lsmi indicates  n/l dependencies of various strength ###
################################################################

nlin.step <- 21
nlin.dots <- 50
nlin.rep.lsmi <- 10

# generating r.v.s
nlin.lsmi.values <- data.frame(lsmi = rep(0, nlin.rep.lsmi*nlin.step),
                               noise = rep(seq(0, 0.5, length.out = nlin.step), each = nlin.rep.lsmi),
                               # correlation estimate
                               corhat = rep(0, nlin.rep.lsmi*nlin.step),
                               corhat.cimin = rep(0, nlin.rep.lsmi*nlin.step),
                               corhat.cimax = rep(0, nlin.rep.lsmi*nlin.step)
)

noise.vs.lsmi.nlin <- data.frame(x = rep(x = runif(nlin.dots, -1, 1), nlin.step),
                                 y = rep(0, nlin.dots*nlin.step),
                                 noise = rep(seq(0, 0.5, length.out = nlin.step), each = nlin.dots)
)

# calculating LSMI
indices.vals <- c(1, nlin.rep.lsmi)
indices.data <- c(1, nlin.dots)

for(i in 1:nlin.step) {
  subs.vals <- indices.vals[1]:indices.vals[2]
  subs.data <- indices.data[1]:indices.data[2]
  
  ##########################################
  ## parabola or sine nonlinear functions ##
  
  #noise.vs.lsmi.nlin[subs.data, 'y'] <- noise.vs.lsmi.nlin[subs.data, 'x']^2 + rnorm(nlin.dots, 0, noise.vs.lsmi.nlin[subs.data[1], 'noise'])
  #noise.vs.lsmi.nlin[subs.data, 'y'] <- sin(4*pi*noise.vs.lsmi.nlin[subs.data, 'x']) + rnorm(nlin.dots, 0, noise.vs.lsmi.nlin[subs.data[1], 'noise'])
  noise.vs.lsmi.nlin[subs.data, 'y'] <- sample(c(-1,1), nlin.dots, TRUE)*sqrt(1 - noise.vs.lsmi.nlin[subs.data, 'x']^2) + rnorm(nlin.dots, 0, noise.vs.lsmi.nlin[subs.data[1], 'noise'])
  
  nlin.lsmi.values[subs.vals, 'lsmi'] <- replicate(nlin.rep.lsmi,
                                                   lsmi.vanilla(as.list(noise.vs.lsmi.nlin[subs.data, 'x']), 
                                                        as.list(noise.vs.lsmi.nlin[subs.data, 'y']))
  )
  
  current.cortest <- cor.test(noise.vs.lsmi.nlin$x[subs.data],
                              noise.vs.lsmi.nlin$y[subs.data])
  
  nlin.lsmi.values[subs.vals, 'corhat'] <- current.cortest$estimate
  nlin.lsmi.values[subs.vals, c('corhat.cimin', 'corhat.cimax')] <- current.cortest$conf.int %>% 
    rep(., nlin.rep.lsmi) %>%
    matrix(., ncol = 2, byrow = TRUE) %>%
    as.data.frame
  
  indices.vals %<>% add(nlin.rep.lsmi)
  indices.data %<>% add(nlin.dots)
}
rm(indices.vals)
rm(indices.data)
rm(subs.vals)
rm(current.cortest)

# plotting results
library(reshape2)
nlin.lsmi.values %<>% melt(., id = c('noise', 'corhat.cimin', 'corhat.cimax'))

ggplot(data = nlin.lsmi.values, aes(noise, value)) + 
  geom_point(aes(color = factor(variable)), size = 3, alpha = 0.3) + 
  geom_errorbar(aes(ymin = corhat.cimin, ymax = corhat.cimax, colour = factor(variable)), width = 0.005) +
  
  scale_color_hue(labels = c("LSMI estimate", "Cor estimate")) +
  theme(plot.background = element_rect(fill = "aliceblue"), legend.title = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(linetype=c(0, 1), shape=c(16, 16), alpha = 1))) + 
  #ggtitle(paste('LSMI / Cor estimates (95% CI) vs. noise strength\n x ~ U[-1;1]; Y = X^2 + N(0, sigma = noise), n = ', nlin.dots, ';', sep = '')) +
  # ggtitle(paste('LSMI / Cor estimates (95% CI) vs. noise strength\n x ~ U[-1;1]; y = Sin 4*pi*x + N(0, sigma = noise), n = ', nlin.dots, ';', sep = '')) + 
  ggtitle(paste('LSMI / Cor estimates (95% CI) vs. noise strength\n x ~ U[-1;1]; y = random(-1;1)*sqrt(1 - x^2) + N(0, sigma = noise), n = ', nlin.dots, ';', sep = '')) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggsave('nonlin_detection_noise_circle.png')
