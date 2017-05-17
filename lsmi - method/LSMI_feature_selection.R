###
## implementation of LSMI-based feature selection
###

miFeaSelect <- function(dataset, classes, nbfuns, selection = c('backward', 'forward'), verbose = TRUE,
                        method = 'suzuki', method.nbfuns = 'balanced', max.iter, marginal.importance = TRUE) {
  require(magrittr)
  require(dplyr)
  require(stringr)
  
  selection %<>% match.arg(c('backward', 'forward'))
  nFea <- ncol(dataset)
  if(missing(max.iter)) max.iter <- nFea - 1
  if(missing(nbfuns)) nbfuns <- min(100, nrow(dataset))
  
  inds.fea <- 1:nFea
  fea.subset.cur <- switch(selection, forward = numeric(0), backward = 1:nFea)
  
  if(selection == 'forward') {
    ## there's no '0th' iteration in forward selection
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
  mi.selection <- numeric(max.iter)
  for(i in 1:max.iter) {
    mi.iter <- numeric(nFea)
    names(mi.iter) <- colnames(dataset)
    
    for(i.fea in inds.fea) {
      sbs.iter <- switch(selection, forward = c(fea.subset.cur, i.fea), backward = setdiff(fea.subset.cur, i.fea))
      mi.iter[i.fea] <-
        lsmi.extra(dataset %>% select(sbs.iter) %>% data.matrix() %>% split(rep(1:nrow(.), length(sbs.iter))),
                   factor(classes), nbfuns = nbfuns, method.nbfuns = method.nbfuns, method = method)
      
      if(marginal.importance & selection == 'forward' & i == 1) mi.importance <- mi.iter
    }
    
    
    ind.best <- which.max(mi.iter)
    fea.subset.cur <- switch(selection, forward = c(fea.subset.cur, ind.best), backward = setdiff(fea.subset.cur, ind.best))
    inds.fea %<>% setdiff(ind.best)
    
    ## saving MI and feature names on this iteration ##
    
    mi.selection[i] <- max(mi.iter)
    names(mi.selection)[i] <- colnames(dataset)[ind.best]
    
    if(verbose) {
      message(str_c('[ Iteration ', i, ' ]'))
      message(str_c('Current MI: ', round(max(mi.iter), 3), '; previous MI: ', round(mi.max.prev, 3)))
      message(str_c('Information gain = ', round(max(mi.iter) - mi.max.prev, 3)))
      message(str_c(mword, "'", names(mi.selection)[i], "'", collapse = ''))
      message('_______________________________________')
    }
    
    mi.max.prev <- max(mi.iter)
  }
  
  if(marginal.importance & selection == 'backward') {
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
    return(data_frame(miIter = mi.selection, 
                      feaIter = names(mi.selection), 
                      iteration = 1:max.iter, 
                      importance = mi.importance[names(mi.selection)]))
  } else if(selection == 'backward') {
    return(data_frame(miIter = c(mi.max.prev0, mi.selection), 
                      feaIter = c('', names(mi.selection)), 
                      iteration = 0:max.iter, 
                      importance = c(0.01, mi.importance)))
      }
}

