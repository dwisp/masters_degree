
###
## sine function, matlab example
###
library(dplyr)
nl.ex1 <- data_frame(x = runif(50, -1, 1), y = sin(x*4*pi) + rnorm(50, 0, 0.25))

nl.ex1.cor <- cor(nl.ex1$x, nl.ex1$y) %>% round(3)
nl.ex1.lsmi <- lsmi.vanilla(nl.ex1$x, nl.ex1$y) %>% round(3)

# plotting dependence #
library(ggplot2)
library(ggthemes)
library(stringr)
ggplot(nl.ex1,  aes(x = x, y = y)) + 
  geom_point(color = 'orange', size = 2, alpha = 0.8) + 
  ggtitle(str_c('Nonlinear dependence example 1: high-frequency sine function\n x ~ U[-1;1]; y = sin 4*pi*x + N(0, 0.25); n = 50\n',
                'Cor = ', nl.ex1.cor, '; ', 'LSMI = ', nl.ex1.lsmi)) +
  theme_solarized() +
  #theme(plot.background = element_rect(fill = "darkseagreen"), plot.title = element_text(hjust = 0.5)) +
  ggsave('sine_example.png')

###
## parabola
###

nl.ex2 <- data_frame(x = runif(50, -1, 1), y = x^2 + rnorm(50, 0, 0.1))

nl.ex2.cor <- cor(nl.ex2$x, nl.ex2$y) %>% round(3)
nl.ex2.lsmi <- lsmi.vanilla(nl.ex2$x, nl.ex2$y) %>% round(3)

ggplot(nl.ex2, aes(x = x, y = y)) + 
  geom_point(color = 'orange', size = 2, alpha = 0.8) + 
  ggtitle(str_c('Nonlinear dependence example 2: parabola\n x ~ U[-1;1]; y = x^2 + N(0, 0.1); n = 50\n',
                'Cor = ', nl.ex2.cor, '; ', 'LSMI = ', nl.ex2.lsmi)) +
  #theme(plot.background = element_rect(fill = "darkseagreen"), plot.title = element_text(hjust = 0.5)) +
  theme_solarized() + 
  ggsave('parabola_example.png')

###
## circle
###

nl.ex3 <- data_frame(x = runif(50, -1, 1), y = sample(c(-1,1), 50, T, c(1,1)/2)*sqrt(1 - x^2) + rnorm(50, 0, 0.1))

nl.ex3.cor <- cor(nl.ex3$x, nl.ex3$y) %>% round(3)
nl.ex3.lsmi <- lsmi.vanilla(nl.ex3$x, nl.ex3$y) %>% round(3)

ggplot(nl.ex3, aes(x = x, y = y)) + 
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1.1, 1.1)) +
  geom_point(color = 'orange', size = 2, alpha = 0.8) +
  scale_x_continuous(limits = c(-1.3, 1.3)) +
  scale_y_continuous(limits = c(-1.3, 1.3)) +
  ggtitle(str_c('Nonlinear dependence example 3: unit circle\n x ~ U[-1;1]; y = random(-1;1)*sqrt(1 - x^2) + N(0, 0.1), n = 50\n',
                'Cor = ', nl.ex3.cor, '; ', 'LSMI = ', nl.ex3.lsmi)) +
  #theme(plot.background = element_rect(fill = "darkseagreen"), plot.title = element_text(hjust = 0.5)) +
  theme_solarized() +
  ggsave('circle_example.png', width = 8, height = 8)

###
## how lsmi indicates  n/l dependencies of various strength
###

nlin.step <- 21
nlin.dots <- 50
nlin.rep.lsmi <- 10

# generating r.v.s
nlin.lsmi.values <- 
  data.frame(lsmi_par = rep(0, nlin.rep.lsmi*nlin.step),
             lsmi_sin = rep(0, nlin.rep.lsmi*nlin.step),
                               lsmi_cir = rep(0, nlin.rep.lsmi*nlin.step),
             noise = rep(seq(0, 0.5, length.out = nlin.step), each = nlin.rep.lsmi),
             # correlation estimate
             ## parabola
             corhat_par = rep(0, nlin.rep.lsmi*nlin.step),
             corhat.cimin_par = rep(0, nlin.rep.lsmi*nlin.step),
             corhat.cimax_par = rep(0, nlin.rep.lsmi*nlin.step),
             ## sine
             corhat_sin = rep(0, nlin.rep.lsmi*nlin.step),
             corhat.cimin_sin = rep(0, nlin.rep.lsmi*nlin.step),
             corhat.cimax_sin = rep(0, nlin.rep.lsmi*nlin.step),
             ## circle
             corhat_cir = rep(0, nlin.rep.lsmi*nlin.step),
             corhat.cimin_cir = rep(0, nlin.rep.lsmi*nlin.step),
             corhat.cimax_cir = rep(0, nlin.rep.lsmi*nlin.step))

noise.vs.lsmi.nlin <- 
  data.frame(x = rep(x = runif(nlin.dots, -1, 1), nlin.step),
             y_par = rep(0, nlin.dots*nlin.step),
             y_sin = rep(0, nlin.dots*nlin.step),
             y_cir = rep(0, nlin.dots*nlin.step),
             noise = rep(seq(0, 0.5, length.out = nlin.step), each = nlin.dots))

# calculating LSMI
indices.vals <- c(1, nlin.rep.lsmi)
indices.data <- c(1, nlin.dots)

for(i in 1:nlin.step) {
  subs.vals <- indices.vals[1]:indices.vals[2]
  subs.data <- indices.data[1]:indices.data[2]
  
  ## parabola, sine and circle
  noise.vs.lsmi.nlin[subs.data, 'y_par'] <- 
    noise.vs.lsmi.nlin[subs.data, 'x']^2 + rnorm(nlin.dots, 0, noise.vs.lsmi.nlin[subs.data[1], 'noise'])
  
  noise.vs.lsmi.nlin[subs.data, 'y_sin'] <- 
    sin(4*pi*noise.vs.lsmi.nlin[subs.data, 'x']) + rnorm(nlin.dots, 0, noise.vs.lsmi.nlin[subs.data[1], 'noise'])
  
  noise.vs.lsmi.nlin[subs.data, 'y_cir'] <- 
    sample(c(-1,1), nlin.dots, TRUE)*sqrt(1 - noise.vs.lsmi.nlin[subs.data, 'x']^2) + rnorm(nlin.dots, 0, noise.vs.lsmi.nlin[subs.data[1], 'noise'])
  
  
  nlin.lsmi.values[subs.vals, 'lsmi_par'] <- 
    replicate(nlin.rep.lsmi, 
              lsmi.vanilla(noise.vs.lsmi.nlin[subs.data, 'x'], 
                           noise.vs.lsmi.nlin[subs.data, 'y_par']))
  
  nlin.lsmi.values[subs.vals, 'lsmi_sin'] <- 
    replicate(nlin.rep.lsmi, 
              lsmi.vanilla(noise.vs.lsmi.nlin[subs.data, 'x'], 
                           noise.vs.lsmi.nlin[subs.data, 'y_sin']))
  
  nlin.lsmi.values[subs.vals, 'lsmi_cir'] <- 
    replicate(nlin.rep.lsmi, 
              lsmi.vanilla(noise.vs.lsmi.nlin[subs.data, 'x'], 
                           noise.vs.lsmi.nlin[subs.data, 'y_cir']))
  
  
  current.cortest_par <- cor.test(noise.vs.lsmi.nlin$x[subs.data],
                                  noise.vs.lsmi.nlin$y_par[subs.data])
  current.cortest_sin <- cor.test(noise.vs.lsmi.nlin$x[subs.data],
                                  noise.vs.lsmi.nlin$y_sin[subs.data])
  current.cortest_cir <- cor.test(noise.vs.lsmi.nlin$x[subs.data],
                                  noise.vs.lsmi.nlin$y_cir[subs.data])
  
  nlin.lsmi.values[subs.vals, 'corhat_par'] <- current.cortest_par$estimate
  nlin.lsmi.values[subs.vals, c('corhat.cimin_par', 'corhat.cimax_par')] <- 
    current.cortest_par$conf.int %>% 
    rep(nlin.rep.lsmi) %>%
    matrix(ncol = 2, byrow = TRUE) %>%
    as.data.frame
  
  nlin.lsmi.values[subs.vals, 'corhat_sin'] <- current.cortest_sin$estimate
  nlin.lsmi.values[subs.vals, c('corhat.cimin_sin', 'corhat.cimax_sin')] <- 
    current.cortest_sin$conf.int %>% 
    rep(nlin.rep.lsmi) %>%
    matrix(ncol = 2, byrow = TRUE) %>%
    as.data.frame
  
  nlin.lsmi.values[subs.vals, 'corhat_cir'] <- current.cortest_cir$estimate
  nlin.lsmi.values[subs.vals, c('corhat.cimin_cir', 'corhat.cimax_cir')] <- 
    current.cortest_cir$conf.int %>% 
    rep(nlin.rep.lsmi) %>%
    matrix(ncol = 2, byrow = TRUE) %>%
    as.data.frame
  
  indices.vals %<>% add(nlin.rep.lsmi)
  indices.data %<>% add(nlin.dots)
}
rm(list = ls(pattern = 'indices'))
rm(subs.vals)
rm(list = ls(pattern = 'current.cortest'))

# plotting results
library(reshape2)
nlin.lsmi.values %<>% melt(id = c('corhat.cimin_par', 'corhat.cimax_par', 
                                  'corhat.cimin_sin', 'corhat.cimax_sin', 
                                  'corhat.cimin_cir', 'corhat.cimax_cir', 
                                  'noise'))

## parabola ##
library(dplyr)
nlin.lsmi.values_par <-
  nlin.lsmi.values %>%
  as_data_frame %>%
  select(value, variable, corhat.cimin_par, corhat.cimax_par, noise) %>%
  filter(variable %in% c('lsmi_par', 'corhat_par'))

library(ggplot2)
library(stringr)
ggplot(data = nlin.lsmi.values_par, aes(noise, value)) + 
  geom_point(aes(color = variable), size = 3, alpha = 0.26) + 
  # stat_smooth(aes(color = variable), se = FALSE, method = 'loess') +
  geom_errorbar(aes(ymin = corhat.cimin_par, ymax = corhat.cimax_par, color = variable), width = 0.005) +
  #scale_color_hue(labels = c("LSMI estimate", "Cor estimate")) +
  #theme(plot.background = element_rect(fill = 'aliceblue'), legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(linetype=c(0, 1), shape = c(16, 16), alpha = 1))) + 
  ggtitle(str_c('LSMI / Cor estimates (95% CI) vs. noise strength\n x ~ U[-1;1]; Y = X^2 + N(0, sigma = noise), n = ', nlin.dots, ';')) + 
  theme_solarized() +
  scale_colour_solarized('blue', labels = c('LSMI estimate', 'Cor estimate')) +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  ggsave('nonlin_detection_noise_parabola.png')

## sine ##
nlin.lsmi.values_sin <-
  nlin.lsmi.values %>%
  as_data_frame %>%
  select(value, variable, corhat.cimin_sin, corhat.cimax_sin, noise) %>%
  filter(variable %in% c('lsmi_sin', 'corhat_sin'))

ggplot(data = nlin.lsmi.values_sin, aes(noise, value)) + 
  geom_point(aes(color = variable), size = 3, alpha = 0.26) +
  # stat_smooth(aes(color = variable), se = FALSE, method = 'loess') +
  geom_errorbar(aes(ymin = corhat.cimin_sin, ymax = corhat.cimax_sin, color = variable), width = 0.005) +
  scale_colour_solarized('blue', labels = c('LSMI estimate', 'Cor estimate')) +
  # theme(plot.background = element_rect(fill = 'aliceblue'), legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, 16), alpha = 1))) + 
  ggtitle(str_c('LSMI / Cor estimates (95% CI) vs. noise strength\n x ~ U[-1;1]; y = Sin 4*pi*x + N(0, sigma = noise), n = ', nlin.dots, ';')) + 
  theme_solarized() +
  # theme(plot.title = element_text(hjust = 0.5)) + 
  ggsave('nonlin_detection_noise_sine.png')


## circle ##
nlin.lsmi.values_cir <-
  nlin.lsmi.values %>%
  as_data_frame %>%
  select(value, variable, corhat.cimin_cir, corhat.cimax_cir, noise) %>%
  filter(variable %in% c('lsmi_cir', 'corhat_cir'))

ggplot(data = nlin.lsmi.values_cir, aes(noise, value)) + 
  geom_point(aes(color = variable), size = 3, alpha = 0.26) + 
  # stat_smooth(aes(color = variable), se = FALSE, method = 'loess') +
  geom_errorbar(aes(ymin = corhat.cimin_cir, ymax = corhat.cimax_cir, color = variable), width = 0.005) +
  scale_colour_solarized('blue', labels = c('LSMI estimate', 'Cor estimate')) +
  # theme(plot.background = element_rect(fill = 'aliceblue'), legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(linetype=c(0, 1), shape = c(16, 16), alpha = 1))) + 
  theme_solarized() +
  ggtitle(str_c('LSMI / Cor estimates (95% CI) vs. noise strength\n x ~ U[-1;1]; y = random(-1;1)*sqrt(1 - x^2) + N(0, sigma = noise), n = ', nlin.dots, ';')) + 
  # theme(plot.title = element_text(hjust = 0.5)) + 
  ggsave('nonlin_detection_noise_circle.png')




