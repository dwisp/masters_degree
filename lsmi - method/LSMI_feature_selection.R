###
## implementation of LSMI-based feature selection
###

###
## Input args:
###
# dataset: data_frame with cases in rows and variables in columns
# clases: a factor which gives labels of the data
# nbfuns: how many centroids to select
# selection: 
#  'forward' starts from nothing and gaters the mose useful features;
#  'backward' starts with all the variables and excludes the least relevant
# verbose: display processing messages?
# method: method of LSMI calculation
# method.nbfuns: method of centroid selection; 'balanced' or 'uniform'
# max.steps: how much steps to allow? The procedure won't terminate until
#            it reaches max.steps or processes all the features
# 
###
## Output values:
###
# returns a dplyr's data_frame with the following columns:
#   
# miStep: mutual information on each step 
# feaStep: names of the features which were included / excluded on each step
# step: step numbers
# importance: marginal importances of all included features
# 
# This result can be plotted with the miFlowPlot() function
# 

miFeaSelect <- function(dataset, classes, nbfuns, selection = c('backward', 'forward'), verbose = TRUE,
                        method = 'suzuki', method.nbfuns = 'balanced', max.steps) {
  require(magrittr)
  require(dplyr)
  require(stringr)
  
  selection %<>% match.arg(c('backward', 'forward'))
  nFea <- ncol(dataset)
  if(missing(max.steps)) max.steps <- switch(selection, backward = nFea - 1, forward = nFea)
  if(missing(nbfuns)) nbfuns <- min(100, nrow(dataset))
  
  inds.fea <- 1:nFea
  fea.subset.cur <- switch(selection, forward = numeric(0), backward = 1:nFea)
  
  if(selection == 'forward') {
    ## there's no '0th' Step in forward selection
    mi.max.prev <- mi.max.prev0 <- 0 
    mword <- 'Included '
    
  } else if(selection == 'backward') {
    ## including all features as a baseline
    mi.max.prev <- mi.max.prev0 <-
      lsmi.extra(dataset %>% select(1:nFea) %>% data.matrix() %>% split(rep(1:nrow(.), nFea)),
                 factor(classes), nbfuns = nbfuns, method.nbfuns = method.nbfuns, method = method)
    mword <- 'Removed '
  }
  
  ## names of the chosen features in chronological order
  mi.selection <- numeric(max.steps)
  for(i in 1:max.steps) {
    mi.step <- numeric(nFea)
    names(mi.step) <- colnames(dataset)
    
    for(i.fea in inds.fea) {
      sbs.step <- switch(selection, forward = c(fea.subset.cur, i.fea), backward = setdiff(fea.subset.cur, i.fea))
      mi.step[i.fea] <-
        lsmi.extra(dataset %>% select(sbs.step) %>% data.matrix() %>% split(rep(1:nrow(.), length(sbs.step))),
                   factor(classes), nbfuns = nbfuns, method.nbfuns = method.nbfuns, method = method)
      
      if(selection == 'forward' & i == 1) mi.importance <- mi.step
    }
    
    
    ind.best <- which.max(mi.step)
    fea.subset.cur <- switch(selection, forward = c(fea.subset.cur, ind.best), backward = setdiff(fea.subset.cur, ind.best))
    inds.fea %<>% setdiff(ind.best)
    
    ## saving MI and feature names on this Step ##
    
    mi.selection[i] <- max(mi.step)
    names(mi.selection)[i] <- colnames(dataset)[ind.best]
    
    if(verbose) {
      message(str_c('[ Step ', i, ' ]'))
      message(str_c('Current MI: ', round(max(mi.step), 3), '; previous MI: ', round(mi.max.prev, 3)))
      message(str_c('Information gain = ', round(max(mi.step) - mi.max.prev, 3)))
      message(str_c(mword, "'", names(mi.selection)[i], "'", collapse = ''))
      message('_______________________________________')
    }
    
    mi.max.prev <- max(mi.step)
  }
  
  if(selection == 'backward') {
    if(verbose) message('Calculating marginal importances...')
    mi.importance <- numeric(length(mi.selection))
    names(mi.importance) <- names(mi.selection)
    
    for(fname in names(mi.importance)) {
      mi.importance[fname] <- 
        lsmi.extra(dataset %>% select(which(colnames(dataset) == fname)) %>% unlist(use.names = FALSE),
                   factor(classes), nbfuns = nbfuns, method.nbfuns = method.nbfuns, method = method)
    }
    
    if(verbose) message('Done!')
  }
  
  if(selection == 'forward') {
    return(data_frame(miStep = mi.selection, 
                      feaStep = names(mi.selection), 
                      step = 1:max.steps, 
                      importance = mi.importance[names(mi.selection)]))
  } else if(selection == 'backward') {
    return(data_frame(miStep = c(mi.max.prev0, mi.selection), 
                      feaStep = c('', names(mi.selection)), 
                      step = 0:max.steps, 
                      importance = c(0.01, mi.importance)))
      }
}

###
# Function for plotting results obtained with the previous function
###


###
## Input args:
###
# miFeaData: data_frame as produced by miFeaSelect()
# allFeatures: character vector with the full set of features from data
# nfeat: how much features to leave in the final subset?
# 
###
## Output values:
###
#
# Returns a ggplot2's plot which can be extended
#

miFlowPlot <- function(miFeaData, allFeatures, nfeat) {
  require(ggplot2)
  require(stringr)
  require(magrittr)
  
  ## determining if it's forward or backward selection
  slctBckwd <- miFeaData$feaStep[1] == ''
  
  if(missing(allFeatures)) allFeatures <- character(0)
  if(!missing(nfeat)) {
    if(!slctBckwd) {
      chosenFeat <- miFeaData$step <= nfeat
    } else {chosenFeat <- max(miFeaData$step) - miFeaData$step <= nfeat}
      miFeaData$`Chosen features` <- factor(ifelse(chosenFeat, 'included', 'excluded'))
  }
  
  minMI <- round(min(miFeaData$miStep), 3)
  maxMI <- round(max(miFeaData$miStep), 3)
  
  feasLeftOut <- setdiff(allFeatures, miFeaData$feaStep)

  stringy <- ifelse(slctBckwd, '-', '+')
  miFeaData$feaStep[miFeaData$feaStep != ''] %<>% str_c(stringy, .)
  ## processing feature names, providing two rows for those with long names
  miFeaData$feaStep %<>% 
    str_replace(' ', '_') %>%
    str_replace(' ', '\n') %>%
    str_replace('_', ' ')
  
  figure <-
    ggplot(miFeaData, aes(x = step, y = miStep, size = importance, label = feaStep)) +
    geom_line(color = 'lightcyan3', alpha = 0.8) +
    geom_point(shape = 18) +
    geom_text(color = 'wheat4', size = 3, vjust = 2) +
    ylab('Mutual Information') +
    scale_x_continuous(breaks = 0:nrow(miFeaData), limits = c(0, max(miFeaData$step)), minor_breaks = NULL) +
    scale_y_continuous(breaks = round(seq(minMI - 0.05, maxMI + 0.05, length.out = 9), 3), 
                       limits = c(minMI - 0.05, maxMI + 0.05)) +
    scale_size_area(limits = c(0.01, maxMI))
  
  ## mentioning features that aren't represented as dots
  if(!missing(allFeatures) & length(feasLeftOut) != 0) {
    # cutting missing features in case there's just too many (like in gene expression data)
    if(length(feasLeftOut) > 4) feasLeftOut %<>% str_c(.[1:4], '... <truncated>')
    figure %<>% add(labs(caption = str_c('Features included on step ', max(miFeaData$step), ': ', str_c(feasLeftOut, collapse = ',  '))))
    }
    
  ## showing which features have been chosen
  if(!missing(nfeat)) {
    return(figure + 
             aes(color = `Chosen features`) + 
             guides(size = guide_legend(title = 'Importance\n(marginal MI)',
                                        ncol = 2, order = 1),
                    color = guide_legend(title = 'Chosen\nfeature set', 
                                         override.aes = list(size = 4),
                                         reverse = TRUE, order = 2)))
    }
  figure
}
