# ##############################################################
# ## testing and evaluating polynomial kernel MI performance  ##
# ##############################################################
# 
# ###################################
# ## sine function, matlab example ##
# ###################################
# 
# nlp.ex1 <- data_frame(x = runif(50, -1, 1), y = sin(x*4*pi) + rnorm(50, 0, 0.25))
# 
# nlp.ex1.cor <- cor(nlp.ex1$x, nlp.ex1$y) %>% round(3)
# nlp.ex1.lsmi <- lsmi.poly(nlp.ex1$x, nlp.ex1$y) %>% round(3)
# 
# # plotting dependence #
# library(ggplot2)
# ggplot(nlp.ex1, aes(x, y)) + 
#   geom_point(color = 'orange', size = 2) + 
#   ggtitle(paste('Nonlinear dependence example 1: high-frequency sine function\n x ~ U[-1;1]; y = sin 4*pi*x + N(0, 0.25); n = 50\n',
#                 'Cor = ', nlp.ex1.cor, '; ', 'polyLSMI = ', nlp.ex1.lsmi, sep = '')) +
#   ggsave('sine_example_poly.png')
# 
# #################################
# ###        a parabola         ###
# #################################
# 
# nlp.ex2 <- data_frame(x = runif(50, -1, 1), y = x^2 + rnorm(50, 0, 0.1))
# 
# nlp.ex2.cor <- cor(nlp.ex2$x, nlp.ex2$y) %>% round(3)
# nlp.ex2.lsmi <- lsmi.poly(nlp.ex2$x, nlp.ex2$y) %>% round(3)
# 
# ggplot(nlp.ex2, aes(x, y)) + 
#   geom_point(color = 'orange', size = 2) + 
#   ggtitle(paste('Nonlinear dependence example 2: parabola\n x ~ U[-1;1]; y = x^2 + N(0, 0.1); n = 50\n',
#                 'Cor = ', nlp.ex2.cor, '; ', 'polyLSMI = ', nlp.ex2.lsmi, sep = '')) +
#   ggsave('parabola_example_poly.png')
# 
# #################################
# ###         a circle          ###
# #################################
# 
# nlp.ex3 <- data_frame(x = runif(50, -1, 1),y = sample(c(-1,1), 50, T, c(1/2,1/2))*sqrt(1 - x^2) + rnorm(50, 0, 0.1))
# 
# nlp.ex3.cor <- cor(nlp.ex3$x, nlp.ex3$y) %>% round(3)
# nlp.ex3.lsmi <- lsmi.poly(nlp.ex3$x, nlp.ex3$y) %>% round(3)
# 
# ggplot(nlp.ex3, aes(x, y)) + 
#   coord_cartesian(xlim = c(-1, 1), ylim = c(-1.1, 1.1)) +
#   geom_point(color = 'orange', size = 2) +
#   ggtitle(paste('Nonlinear dependence example 3: unit circle\n x ~ U[-1;1]; y = random(-1;1)*sqrt(1 - x^2) + N(0, 0.1), n = 50\n',
#                 'Cor = ', nlp.ex3.cor, '; ', 'polyLSMI = ', nlp.ex3.lsmi, sep = '')) +
#   ggsave('circle_example_poly.png', width = 9, height = 9)

###
# how poly lsmi indicates  n/l dependencies of various strength
###

nlin.step <- 21
nlin.dots <- 50
nlin.rep.lsmi <- 10

# generating r.v.s
nlinp.lsmi.values <- 
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

# calculating polyLSMI
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
  
  
  nlinp.lsmi.values[subs.vals, 'lsmi_par'] <- replicate(nlin.rep.lsmi, 
                                                       lsmi.poly(noise.vs.lsmi.nlin[subs.data, 'x'], 
                                                                 noise.vs.lsmi.nlin[subs.data, 'y_par']))
  
  nlinp.lsmi.values[subs.vals, 'lsmi_sin'] <- replicate(nlin.rep.lsmi, 
                                                       lsmi.poly(noise.vs.lsmi.nlin[subs.data, 'x'], 
                                                                 noise.vs.lsmi.nlin[subs.data, 'y_sin']))
  
  nlinp.lsmi.values[subs.vals, 'lsmi_cir'] <- replicate(nlin.rep.lsmi, 
                                                       lsmi.poly(noise.vs.lsmi.nlin[subs.data, 'x'], 
                                                                 noise.vs.lsmi.nlin[subs.data, 'y_cir']))
  
  current.cortest_par <- cor.test(noise.vs.lsmi.nlin$x[subs.data],
                                  noise.vs.lsmi.nlin$y_par[subs.data])
  current.cortest_sin <- cor.test(noise.vs.lsmi.nlin$x[subs.data],
                                  noise.vs.lsmi.nlin$y_sin[subs.data])
  current.cortest_cir <- cor.test(noise.vs.lsmi.nlin$x[subs.data],
                                  noise.vs.lsmi.nlin$y_cir[subs.data])
  
  nlinp.lsmi.values[subs.vals, 'corhat_par'] <- current.cortest_par$estimate
  nlinp.lsmi.values[subs.vals, c('corhat.cimin_par', 'corhat.cimax_par')] <- 
    current.cortest_par$conf.int %>% 
    rep(nlin.rep.lsmi) %>%
    matrix(ncol = 2, byrow = TRUE) %>%
    as.data.frame
  
  nlinp.lsmi.values[subs.vals, 'corhat_sin'] <- current.cortest_sin$estimate
  nlinp.lsmi.values[subs.vals, c('corhat.cimin_sin', 'corhat.cimax_sin')] <- 
    current.cortest_sin$conf.int %>% 
    rep(nlin.rep.lsmi) %>%
    matrix(ncol = 2, byrow = TRUE) %>%
    as.data.frame
  
  nlinp.lsmi.values[subs.vals, 'corhat_cir'] <- current.cortest_cir$estimate
  nlinp.lsmi.values[subs.vals, c('corhat.cimin_cir', 'corhat.cimax_cir')] <- 
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
nlinp.lsmi.values %<>% melt(id = c('corhat.cimin_par', 'corhat.cimax_par', 
                                  'corhat.cimin_sin', 'corhat.cimax_sin', 
                                  'corhat.cimin_cir', 'corhat.cimax_cir', 
                                  'noise'))
##
# parabola
##

library(dplyr)
nlinp.lsmi.values_par <-
  nlinp.lsmi.values %>%
  as_data_frame %>%
  select(value, variable, corhat.cimin_par, corhat.cimax_par, noise) %>%
  filter(variable %in% c('lsmi_par', 'corhat_par')) %>%
  filter(value > -1.5) %>%
  rename(corhat.cimin = corhat.cimin_par, corhat.cimax = corhat.cimax_par)


## a shortcut for plotting polynomial LSMI
plotPolyLSMI <- function(dataset) {
  require(ggplot2)
  require(ggthemes)
  require(stringr)
  ggplot(data = dataset, aes(noise, value)) + 
    geom_point(aes(color = variable), size = 3, alpha = 0.26) + 
    geom_errorbar(aes(ymin = corhat.cimin, ymax = corhat.cimax, color = variable), width = 0.005) +
    geom_smooth(method = 'loess', se = FALSE, aes(color = variable, linetype = variable)) +
    theme_solarized() +
    scale_color_solarized(labels = c('Polynomial\nLSMI estimate', 'Cor estimate')) +
    scale_linetype_manual(values = c(1, 0), guide = FALSE) +
    guides(color = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, 16), alpha = 1)))
}

 plotPolyLSMI(nlinp.lsmi.values_par) + 
  ggtitle(str_c('Polynomial LSMI and Cor estimates (95% CI) vs. noise strength\n x ~ U[-1;1]; Y = X^2 + N(0, sigma = noise), n = ', nlin.dots, ';')) + 
  ggsave('nonlin_detection_noise_parabola_poly.png')

##
# sine
##
 
nlinp.lsmi.values_sin <-
  nlinp.lsmi.values %>%
  as_data_frame %>%
  select(value, variable, corhat.cimin_sin, corhat.cimax_sin, noise) %>%
  filter(variable %in% c('lsmi_sin', 'corhat_sin')) %>%
  rename(corhat.cimin = corhat.cimin_sin, corhat.cimax = corhat.cimax_sin)

plotPolyLSMI(nlinp.lsmi.values_sin) + 
  ggtitle(str_c('Polynomial MI and Cor estimates (95% CI) vs. noise strength\n x ~ U[-1;1]; y = Sin 4*pi*x + N(0, sigma = noise), n = ', nlin.dots, ';')) + 
  ggsave('nonlin_detection_noise_sine_poly.png')

##
# circle
##

nlinp.lsmi.values_cir <-
  nlinp.lsmi.values %>%
  as_data_frame %>%
  select(value, variable, corhat.cimin_cir, corhat.cimax_cir, noise) %>%
  filter(variable %in% c('lsmi_cir', 'corhat_cir')) %>%
  rename(corhat.cimin = corhat.cimin_cir, corhat.cimax = corhat.cimax_cir)

plotPolyLSMI(nlinp.lsmi.values_cir) + 
  ggtitle(str_c('Polynomial MI and Cor estimates (95% CI) vs. noise strength\n x ~ U[-1;1]; y = random(-1;1)*sqrt(1 - x^2) + N(0, sigma = noise), n = ', nlin.dots, ';')) + 
  ggsave('nonlin_detection_noise_circle_poly.png')
