
##==================================##
##                                  ##
## Functions to process each sample ## 
##                                  ##
##==================================##

### for testing
# setwd("/Volumes/proj-tracerx-lung/tctProjects/frankella/R_packages/eclipse")
# setwd("/camp/project/proj-tracerx-lung/tctProjects/frankella/R_packages/eclipse")
# data <- data.table::fread("data/example_datacsv")
# tree <- data.table::fread('data/tree.csv')
# 
# normSD.for.1.sup.mut.clones <- 0.77
# lod_correction <- TRUE
# test_subsequent_CIN <- TRUE
# 
# #load libraries
# library(data.table)
# library(diptest)
# 
# hard_filters <- c( "primer_abundance_filter", "primer_strand_bias" )
# mrd_filters <- c( "tnc_error_rate" )
# 
# out <- clonal_deconvolution( data = data, 
#                              tree = tree, 
#                              hard_filters = hard_filters, 
#                              mrd_filters = mrd_filters, 
#                              normSD.for.1.sup.mut.clones = normSD.for.1.sup.mut.clones  )
# 

make.random.lognormal <- function(mean, sd, n){
  if(!mean==0){
    location <- log(mean^2 / sqrt(sd^2 + mean^2))
    shape <- sqrt(log(1 + (sd^2 / mean^2)))
    return(rlnorm(n = n, location, shape))
  } else {
    return(rep(0,n))
  }
}


model.subclone.DAOs <- function(clonalpopmean, CCF, clone_name.table){
  VAF.CCF <-  clonalpopmean * CCF
  CCF.expected.subclone <- make.random.lognormal(mean=VAF.CCF,sd=0.75*VAF.CCF,n = nrow(clone_name.table))  * clone_name.table$depth
  max.read.support.to.test <- round(mean(CCF.expected.subclone)*length(CCF.expected.subclone))
  CCF.expected.subclone.rounded <- round(dpois(0:max.read.support.to.test, lambda = mean(CCF.expected.subclone)) * length(CCF.expected.subclone))
  CCF.expected.subclone.unrounded <- dpois(0:max.read.support.to.test, lambda = mean(CCF.expected.subclone)) * length(CCF.expected.subclone)
  if(!sum(CCF.expected.subclone.rounded)==round(sum(CCF.expected.subclone.unrounded))){
    diff <- round(sum(CCF.expected.subclone.unrounded)) - sum(CCF.expected.subclone.rounded)
    round.diffs <- CCF.expected.subclone.rounded - CCF.expected.subclone.unrounded
    names(round.diffs) <- 1:length(round.diffs)
    if(diff < 0){
      neg.indexes <- names(round.diffs[order(round.diffs,decreasing = T)])[1:abs(diff)]
      CCF.expected.subclone.rounded[as.numeric(neg.indexes)] <- CCF.expected.subclone.rounded[as.numeric(neg.indexes)] - 1
    } else {
      pos.indexes <- names(round.diffs[order(round.diffs,decreasing = F)])[1:abs(diff)]
      CCF.expected.subclone.rounded[as.numeric(pos.indexes)] <- CCF.expected.subclone.rounded[as.numeric(pos.indexes)] + 1
    }
  }
  CCF.expected.subclone <- CCF.expected.subclone.rounded
  names(CCF.expected.subclone) <- 0:max.read.support.to.test
  CCF.expected.subclone <- CCF.expected.subclone[!CCF.expected.subclone==0]
  CCF.expected.subclone <- as.numeric(unlist(lapply(1:length(CCF.expected.subclone), function(i) rep(names(CCF.expected.subclone)[i],CCF.expected.subclone[i]))))
  return(CCF.expected.subclone)
}

# Fxn to estimate average ofCN.adj.VAF in a clone accounting for LODs and when only tail of distrebution is observed 
mean.estimate.LOD.correction <- function(CN.adj.VAF, MatchedLOD, normSD.est = NA, 
                                         use.vars = "all", normSD.search.space = c(0.2, 1.2), 
                                         normSD.search.incements = 0.001, Remove.negative = FALSE){
  
  if(all(use.vars=="all")){
    use.vars <- rep(TRUE,length(CN.adj.VAF))
  }
  
  #record original mutorder for later
  originalmutcells <- CN.adj.VAF
  
  #only consider muts previous slected based on LOD & that are good quality
  CN.adj.VAF <- CN.adj.VAF[MatchedLOD & use.vars]
  
  #remove negattive if mutations loss leading in high VAF case leading to 0s
  if(Remove.negative == TRUE){
    CN.adj.VAF <- CN.adj.VAF[!CN.adj.VAF==0]
  }
  
  #calcualate Z-score for each detected mutation (only using those which do not casue LOD bias)
  #add percentages to muts
  CN.adj.VAF <- CN.adj.VAF[order(CN.adj.VAF)]
  Percentages <- seq(0, 1 ,length.out = length(CN.adj.VAF) + 2 )
  Percentages <- Percentages[!Percentages %in% c(0,1)]
  Zscores <- qnorm(Percentages)
  
  Zscores <- Zscores[CN.adj.VAF>0 & !is.na(CN.adj.VAF)]
  CN.adj.VAF <- CN.adj.VAF[CN.adj.VAF>0 & !is.na(CN.adj.VAF)]
  
  # Vary the normSD until to find the best fit (smallest correlation). Only if enouch muts (can't get a 0.5 correlation with just two.. require 4)
  if( length(CN.adj.VAF)>1 ){
    
    #corlelations <- c()
    variation <- c()
    
    for ( normSD in seq( min(normSD.search.space), max(normSD.search.space), normSD.search.incements ) ){
      
      logmeanests <- log10(CN.adj.VAF)/(Zscores * log10(normSD) + 1)
      meanests <- 10 ^ logmeanests
      
      #plot(CN.adj.VAF,meanests)
      
      variation <- c(variation, sd(meanests)/mean(meanests))
      
    } 
    
    onlyonemut <- FALSE
    nosignal <- FALSE
    
  } else { 
    
    if(length(CN.adj.VAF) == 1){
      onlyonemut <- TRUE
      nosignal <- FALSE
      
      #can't calculate CIs with only one data piont and can't estimate normSD as with mulitple muts 
      #hence just calculate LOD (no CIs) based on average SD of other data (specified in  function)
      normSD <- normSD.est
      
      logmeanests <- log10(CN.adj.VAF)/(Zscores * log10(normSD) + 1)
      logpopmean <- mean(logmeanests)
      logsd <- log10(normSD)
      n <- length(CN.adj.VAF)
      logUpCI <- logpopmean - ( qnorm(0.975) * logsd / sqrt(n) )
      logLwCI <- logpopmean + ( qnorm(0.975) * logsd / sqrt(n) )
      
      meanests <- 10 ^ logmeanests
      popmean <- 10 ^ logpopmean
      UpCI <- 10 ^ logUpCI
      LwCI <- 10 ^ logLwCI  
      
      allmutn <- length(originalmutcells)
      
      onlyonemut <- rep(onlyonemut,allmutn)
      normSD.var.result <- rep("too.few.muts",allmutn)
      normSDsolution <- rep(normSD,allmutn)
      popmean <- rep(10 ^ logpopmean,allmutn)
      UpCI <- rep(10 ^ logUpCI,allmutn)
      LwCI <- rep(10 ^ logLwCI,allmutn)
      
      #add back in NA Zscores and meanests for muts we couldn't use (negative mut cell or not LOD match)
      tmp <- rep(NA,allmutn)
      tmp[match(CN.adj.VAF,originalmutcells)] <- Zscores
      Zscores <- tmp
      tmp <- rep(NA,allmutn)
      tmp[match(CN.adj.VAF,originalmutcells)] <- meanests
      meanests <- tmp
      
      output <- as.dataframe(cbind(Zscores, onlyonemut, normSDsolution, normSD.var.result, meanests, popmean, UpCI, LwCI))
      
      names(output) <- c('zscores', 'is_one_mut', 'normalised_sd', 'mean_est_variance', 
                         'mean_estimate', 'ctDNA_fraction', 'ctDNA_fraction_UpCI', 
                         'ctDNA_fraction_LwCI' )
      
      return( output )
      
    } else {
      onlyonemut <- TRUE
      nosignal <- TRUE
    }
    
  }
  
  if(nosignal == FALSE){
    normSDs <- seq( min(normSD.search.space), max(normSD.search.space), normSD.search.incements )
    normSDs <- normSDs[ !variation == Inf & !variation %in% NaN & !is.na(variation)]
    variation <- variation[ !variation == Inf & !variation %in% NaN & !is.na(variation)]
    #plot(variation,normSDs) ## for testing
    solution_index <- which(variation == min(variation) )

    testedaround <- solution_index/length(variation) > 0.1 & solution_index/length(variation) < 0.95 & solution_index > 5 & length(variation) - solution_index > 5
    
    if(testedaround == TRUE){ 
      #positive correlation for above? and negative for below? tests for nice U-quadratic dist 
      lower_test <- solution_index - 50
      if( lower_test  < 0 )  lower_test <- 0  
      upper_test <- solution_index + 50
      if( upper_test  > length(variation) )  upper_test <- length(variation)
                              
      Uquadtest <- cor( lower_test:solution_index, variation[ lower_test:solution_index ] ) < -0.9 & +
        cor(  solution_index:upper_test , variation[ solution_index:upper_test ] ) > 0.9
      
    } else {
      Uquadtest <- FALSE
    }
  } else {
    Uquadtest <- FALSE
  }
  
  #if we have sufficiently searched the normSD space - ie if variations fit u shape
  if(Uquadtest == TRUE){
    
    normSD <- normSDs[ solution_index ]
    
    logmeanests <- log10(CN.adj.VAF)/(Zscores * log10(normSD) + 1)
    logpopmean <- mean(logmeanests)
    logsd <- log10(normSD)
    n <- length(CN.adj.VAF)
    logUpCI <- logpopmean + ( qnorm(0.975) * logsd / sqrt(n) )
    logLwCI <- logpopmean - ( qnorm(0.975) * logsd / sqrt(n) )
    
    meanests <- 10 ^ logmeanests
    popmean <- 10 ^ logpopmean
    UpCI <- 10 ^ logUpCI
    LwCI <- 10 ^ logLwCI  
    
    allmutn <- length(originalmutcells)
    
    onlyonemut <- rep(onlyonemut,allmutn)
    normSD.var.result <- rep(min(variation), allmutn )
    normSDsolution <- rep(normSD,allmutn)
    popmean <- rep(10 ^ logpopmean,allmutn)
    UpCI <- rep(10 ^ logUpCI,allmutn)
    LwCI <- rep(10 ^ logLwCI,allmutn)
    
    #add back in NA Zscores and meanests for muts we couldn't use (negative mut cell or not LOD match)
    tmp <- rep(NA,allmutn)
    tmp[match(CN.adj.VAF,originalmutcells)] <- Zscores
    Zscores <- tmp
    tmp <- rep(NA,allmutn)
    tmp[match(CN.adj.VAF,originalmutcells)] <- meanests
    meanests <- tmp
    
    output <- as.dataframe(cbind(Zscores, onlyonemut, normSDsolution, normSD.var.result, meanests, popmean, UpCI, LwCI))
    
    names(output) <- c('zscores', 'is_one_mut', 'normalised_sd', 'mean_est_variance', 'mean_estimate', 'ctDNA_fraction', 'ctDNA_fraction_UpCI', 'ctDNA_fraction_LwCI' )
    
    return( output )
    
  } else {
    
    allmutn <- length(originalmutcells)
    
    if(nosignal == TRUE){
      
      onlyonemut <- rep(NA,allmutn)
      normSD.var.result <- rep(NA,allmutn)
      meanests <- rep(0,allmutn)
      popmean <- rep(0,allmutn)
      UpCI <- rep(0,allmutn)
      LwCI <- rep(0,allmutn)
      normSDsolution <- rep(NA,allmutn)
      
      if( all(is.na(originalmutcells)) ){
        
        onlyonemut <- rep(NA,allmutn)
        normSD.var.result <- rep(NA,allmutn)
        meanests <- rep(NA,allmutn)
        popmean <- rep(NA,allmutn)
        UpCI <- rep(NA,allmutn)
        LwCI <- rep(NA,allmutn)
        
      }
      
      #add back in NA Zscores and meanests for muts we couldn't use (negative mut cell or not LOD match)
      tmp <- rep(NA,allmutn)
      tmp[match(CN.adj.VAF,originalmutcells)] <- Zscores
      Zscores <- tmp
      tmp <- rep(NA,allmutn)
      tmp[match(CN.adj.VAF,originalmutcells)] <- meanests
      meanests <- tmp
      
    } else {
      
      onlyonemut <- rep(NA,allmutn)
      normSD.var.result <- rep("Poor_SD_Solution",allmutn)
      meanests <- rep("Poor_SD_Solution",allmutn)
      popmean <- rep("Poor_SD_Solution",allmutn)
      UpCI <- rep("Poor_SD_Solution",allmutn)
      LwCI <- rep("Poor_SD_Solution",allmutn)
      normSDsolution <- rep("Poor_SD_Solution",allmutn)
      Zscores <- rep("Poor_SD_Solution",allmutn)
      meanests <- rep("Poor_SD_Solution",allmutn)
      
    }
    
    output <- as.dataframe(cbind(Zscores, onlyonemut, normSDsolution, normSD.var.result, meanests, popmean, UpCI, LwCI))
    
    names(output) <- c('zscores', 'is_one_mut', 'normalised_sd', 'mean_est_variance', 'mean_estimate', 'ctDNA_fraction', 'ctDNA_fraction_UpCI', 'ctDNA_fraction_LwCI' )
    
    return( output )    
  }
}


#' Function to simulate a log normal distribution while defining true mean and sd
#' 
#' 
#'  
#'    
#' @export
random_lognormal<- function(mean, sd, n){
  
  if( !mean == 0 ){
    location <- log( mean ^ 2 / sqrt( sd ^ 2 + mean ^ 2 ) )
    shape <- sqrt( log( 1 + ( sd ^ 2 / mean ^ 2) ) )
    return( rlnorm( n = n, location, shape ) )
  } else {
    return( rep( 0, n ) )
  }
  
}

#' Function to simulate expected VAFs of subclone at a specific CCF
#' 
#' 
#'  
#'    
#' `@export 
model_subclone <- function(clonalpopmean, CCF, clone_table){
  VAF.CCF <-  clonalpopmean * CCF
  CCF.expected.subclone <- make.random.lognormal(mean=VAF.CCF,sd=0.75*VAF.CCF,n = nrow(clone_table))  * clone_table$Deep_Depth
  max.read.support.to.test <- round(mean(CCF.expected.subclone)*length(CCF.expected.subclone))
  CCF.expected.subclone.rounded <- round(dpois(0:max.read.support.to.test, lambda = mean(CCF.expected.subclone)) * length(CCF.expected.subclone))
  CCF.expected.subclone.unrounded <- dpois(0:max.read.support.to.test, lambda = mean(CCF.expected.subclone)) * length(CCF.expected.subclone)
  if(!sum(CCF.expected.subclone.rounded)==round(sum(CCF.expected.subclone.unrounded))){
    diff <- round(sum(CCF.expected.subclone.unrounded)) - sum(CCF.expected.subclone.rounded)
    round.diffs <- CCF.expected.subclone.rounded - CCF.expected.subclone.unrounded
    names(round.diffs) <- 1:length(round.diffs)
    if(diff < 0){
      neg.indexes <- names(round.diffs[order(round.diffs,decreasing = T)])[1:abs(diff)]
      CCF.expected.subclone.rounded[as.numeric(neg.indexes)] <- CCF.expected.subclone.rounded[as.numeric(neg.indexes)] - 1
    } else {
      pos.indexes <- names(round.diffs[order(round.diffs,decreasing = F)])[1:abs(diff)]
      CCF.expected.subclone.rounded[as.numeric(pos.indexes)] <- CCF.expected.subclone.rounded[as.numeric(pos.indexes)] + 1
    }
  }
  CCF.expected.subclone <- CCF.expected.subclone.rounded
  names(CCF.expected.subclone) <- 0:max.read.support.to.test
  CCF.expected.subclone <- CCF.expected.subclone[!CCF.expected.subclone==0]
  CCF.expected.subclone <- as.numeric(unlist(lapply(1:length(CCF.expected.subclone), function(i) rep(names(CCF.expected.subclone)[i],CCF.expected.subclone[i]))))
  return(CCF.expected.subclone)
}


#' Adjust VAF for CN status
#' 
#'  
#'    
#' @export
#Calculate the amount of WT DNA (VAF) is actually coming from the tumour or normal cells
cn_adjust <- function(supporting_reads, depth, total_cn, multiplicity, ccf, is.Clonal = NA, 
                      return.all.intermediate.calculations = FALSE, cnEVO = FALSE ){
  
  # CCFs can be sometimes estimated as > 1, impossible hence limit to 1
  ccf[ ccf > 1 ] <- 1
  
  # estimate the average number of WT copies per cell across the tumour
  Extra.WT.copies <- total_cn - (multiplicity * ccf)
  Extra.WT.copies[ Extra.WT.copies < 0 ] <- 0
  
  
  if(cnEVO == FALSE){
    
    reads.per.totalMutCPN <- supporting_reads / (multiplicity * ccf)
    
  } else {
    
    reads.per.totalMutCPN <- multiplicity * depth
    
  }
  
  WT.tumour.reads <- Extra.WT.copies * reads.per.totalMutCPN
  
  #don't allow more WT reads than depth-(mutCN * varount) - this suggests you 
  #are underestimating the number of MT copies (has there been an amplification?)
  WT.overest <- WT.tumour.reads > (depth - (supporting_reads / multiplicity)) & !is.na(WT.tumour.reads)
  WT.tumour.reads[ WT.overest ] <- (depth - (supporting_reads / multiplicity))[ WT.overest ]
  poss_amplification <- WT.overest
  
  CN.adj.vaf.no.WT.Tumour.DNA.correction <- (supporting_reads / multiplicity) / (depth + supporting_reads)
  
  #calculate "cellularity" by getting VAF per cell (/Mtcopies) and deviding this by VAF signal from normal cells
  Corrected.Depth.2tumour.copies <- (depth - (reads.per.totalMutCPN * (2 - multiplicity)) - WT.tumour.reads)
  Corrected.Depth.1tumour.copy <- (depth - (reads.per.totalMutCPN * (1 - multiplicity)) - WT.tumour.reads)
  Corrected.Depth.for.weight <- (depth - (reads.per.totalMutCPN * (total_cn - multiplicity)) - WT.tumour.reads)
  
  Cellularity.eqivilent <- 2 * (supporting_reads / multiplicity) / Corrected.Depth.2tumour.copies
  CN.adj.VAF <- (supporting_reads / multiplicity) / Corrected.Depth.1tumour.copy
  Tumour.DNA.weight.fraction <- total_cn * (supporting_reads / multiplicity) / Corrected.Depth.for.weight
  
  if(return.all.intermediate.calculations == TRUE){
    
    return(cbind(CN.adj.VAF,
                 Cellularity.eqivilent,
                 Tumour.DNA.weight.fraction,
                 Extra.WT.copies,
                 reads.per.totalMutCPN,
                 WT.tumour.reads,
                 Corrected.Depth.2tumour.copies,
                 Corrected.Depth.1tumour.copy,
                 CN.adj.vaf.no.WT.Tumour.DNA.correction,
                 poss_amplification))
    
  } else {
    
    return(CN.adj.VAF)
    
  }
}

sumlog <- function(p) {
  keep <- (p > 0) & (p <= 1)
  lnp <- log(p[keep])
  chisq <- (-2) * sum(lnp)
  df <- 2 * length(lnp)
  if(sum(1L * keep) < 2)
    stop("Must have at least two valid p values")
  if(length(lnp) != length(p)) {
    warning("Some studies omitted")
  }
  res <- list(chisq = chisq, df = df,
              p = pchisq(chisq, df, lower.tail = FALSE), validp = p[keep])
  class(res) <- c("sumlog", "metap")
  res
}

correct.specific.ps <- function(table,field.to.deduplicate=NA,field.to.correct, indices.to.correct = "all",method ="BH"){
  
  table <- as.dataframe(table)
  
  if(all(indices.to.correct =="all")){
    indices.to.correct <- rep(TRUE,nrow(table))
  }
  
  #check classes are correct
  table[,field.to.deduplicate] <- as.character(table[,field.to.deduplicate])
  table[,field.to.correct] <- as.numeric(table[,field.to.correct])
  
  #need to ensure all dedup fields ordered together (makes code a lot faster) then can restore order at the end
  order.change <- order(paste(table[,field.to.deduplicate]))
  table <- table[order.change,]
  indices.to.correct <- indices.to.correct[order.change]
  
  if(!is.na(field.to.deduplicate)){
    table.filtered <- table[indices.to.correct,]
    ps.to.correct <- table.filtered[!duplicated(table.filtered[,field.to.deduplicate]),field.to.correct]
    unique.dedups <- table.filtered[!duplicated(table.filtered[,field.to.deduplicate]),field.to.deduplicate]
  } else {
    ps.to.correct <- table[indices.to.correct,field.to.correct]
  }
  
  if(grepl("\\.p$|_p$",field.to.correct)) new.field.name <- gsub("\\.p$|_p$",".q",field.to.correct) else new.field.name <- paste0(field.to.correct,".q")
  
  table[,new.field.name] <- NA
  
  if( length(ps.to.correct) > 0){
  ps.corrected <- p.adjust(ps.to.correct, n = length(ps.to.correct), method = method)
  
  if(!is.na(field.to.deduplicate)){
    lengths <- sapply(1:length(ps.corrected), function(i){
      if(!is.na(unique.dedups[i])) return( sum((table[,field.to.deduplicate][indices.to.correct] == unique.dedups[i]) %in% 'TRUE' ) )
      if( is.na(unique.dedups[i])) return( sum(is.na(table[,field.to.deduplicate][indices.to.correct]) ) )
    })
    
    table[,new.field.name][indices.to.correct] <- round(unlist(lapply(1:length(ps.corrected), function(i) rep(ps.corrected[i], lengths[i]))),digits=4)
  } else {
    table[,new.field.name][indices.to.correct] <- ps.corrected
  }
  
  }
  #restore order
  table <- table[match(1:nrow(table),order.change),]
  
  #make the q column next to the p column
  field.to.correct.col <- which(names(table)==field.to.correct)
  table <- table[,c(1:field.to.correct.col,ncol(table),(field.to.correct.col+1):(ncol(table)-1))]
  
  return(table)
  
}

#' Function to correct data when we see evidence of subsequent CIN
#' 
#' 
#' 
#' 
mulitmodal.corrections <- function(data){
  
 
  data$is.multimodal <- data$diptestresult.p < 0.05
  
  suppressPackageStartupMessages(library(mclust))
  
  ###for those positive samples, correct mutation CN and remove mutations which are no longer prersent from further analysis, indcate
  
  mulit.model.clones <- unique(data[data$is.multimodal == TRUE & !is.na(data$is.multimodal), "clone" ])

  data$old_mutcpn <- NA
  data$new_mutcpn <- NA
  data$cn_change <- NA
  data$new_clone_ccf <- NA

  if( length( mulit.model.clones ) == 0 ) return( data )
  
    data <- do.call(rbind, lapply( unique(data$clone), function(clone_name){
      
      clone_table <- data[data$clone %in% clone_name,]
      
      if( ! clone_name %in% mulit.model.clones ) return( clone_table )
      
      VAFs <- clone_table$VAF
      VAFs[VAFs==0] <- clone_table[clone_table$VAF==0,"LOD"]
      
      oldMutCN <- clone_table$multiplicity
      oldMutCN[oldMutCN>1.5] <- 2
      oldMutCN[oldMutCN<1.5] <- 1
      
      BIC <- mclustBIC(VAFs)
      mod1 <- Mclust(VAFs, x = BIC)
      classification <- mod1$classification
      classes <- unique(classification)
      class.VAF <- sapply(classes,function(x) 10^mean(log10(VAFs[classification==x])))
      names(class.VAF) <- classes
      
      #are either classes not signifcantly different from backgroound?
      class.positive <-  sapply(classes,function(class){
        supporting_reads <- clone_table$supporting_reads[classification==class]
        Expected <- clone_table$background_reads[classification==class]
        result <- poisson.test(sum(supporting_reads),sum(Expected),alternative = "greater")$p.value
        out <- result < 0.01
        return(out)
      })
      names(class.positive) <- classes
      
      newMutCN <- oldMutCN
      ##if cluster = background then CN of these mutations must now be 0
      newMutCN[classification==classes[!class.positive]] <- 0
      
      #now ensure that its not an incompleteCIN event where mutations that have been lost are 
      if(sum(!class.positive)==0){
        nmuts <- sapply(classes,function(class) sum( classification==class ) )
        class.1mutcpn.perc <- sapply(classes,function(class) sum(oldMutCN[classification==class]==1,na.rm = T)/sum(classification==class))
        #cluster of 1 copy should be one of the top 2 clusers with the most muts and the one of these with the most 1 cp original muts
        classcp1 <- classes[ order(nmuts,decreasing=T)[1:2] ]
        classcp1 <- classcp1[class.1mutcpn.perc == max(class.1mutcpn.perc[ order(nmuts,decreasing=T)[1:2] ])]
        #if they both have the same number of 1 muts (no WGD eg) should be detelions and cp1 hence 1cp is the higher VAF cluster
        if(length(classcp1)>1){
          classcp1 <- classcp1[class.VAF[classes %in% classcp1] %in% max(class.VAF[unique(mod1$classification) %in% classcp1])]
          newMutCN[classification==classes[classes==classcp1]] <- 1
        }
        
        #now we've got CN 1 then CN 2 must be one above in VAF
        classes.vaf.ordered <- classes [ order(class.VAF, decreasing = F) ]
        classcp2 <- classes.vaf.ordered[ which( classes.vaf.ordered == classcp1 ) + 1 ]
        if( length(classcp2) > 0 ){
          newMutCN[classification == classcp2 ] <- 2
          
          if(!classcp2 == classes.vaf.ordered[ length(classes.vaf.ordered) ]){
            classcpamp <- classes.vaf.ordered[ which( classes.vaf.ordered == classcp2 ):length(classes.vaf.ordered) ]
            newMutCN[classification %in% classcpamp ] <- "gain"
          }
        }
        
        if(!classcp1 == classes.vaf.ordered[1]){
          classcp0 <- classes.vaf.ordered[ 1:which( classes.vaf.ordered == classcp1 ) - 1 ]
          
          newMutCN[classification %in% classcp0 ] <- 0
          
          #if we have mutations at 0 mutCN but with signal above background then we can estimate the CCF of the new clone (!!) & must have been CIN not pre-MRCA mets
          Copy1VAF.new <- as.numeric(class.VAF[names(class.VAF) %in% classcp1] - class.VAF[names(class.VAF) %in% classcp0])
          Copy1VAF.old <- 10^mean(log10(VAFs[classification==classcp0] / oldMutCN [classification==classcp0])) ##might need to remove background niose from 0
          new.clone.CCF <- 1 - (Copy1VAF.old / Copy1VAF.new)
          clone_table$new.clone.CCF <- new.clone.CCF
          
        }
        
      } else {
        
        positive.classes <- classes[class.positive]
        positive.classes <- positive.classes[order(positive.classes)]
        positive.classes.newCN <- 1:length(positive.classes)
        positive.classes.newCN[positive.classes.newCN>2] <- "amp"
        newMutCN[classification %in% positive.classes] <- positive.classes.newCN[match(classification[classification %in% positive.classes],positive.classes)]
        
      }
      
      newMutCN. <- newMutCN
      newMutCN.[newMutCN.=="gain"] <- 3
      CN.change <- as.numeric(newMutCN.) - as.numeric(oldMutCN)
      
      clone_table$old_mutcpn<- oldMutCN
      clone_table$new_mutcpn <- newMutCN
      clone_table$cn_change <- CN.change
      
      return( clone_table )
      
    }) )
    
    sample_data <- as.data.table(data)
    #need approximate CCFs to calculate WT tumour DNA contribution (will iterate to correct CCFs)
    #use meanVAF for this initially
    CCFs <- sample_data[, .( mutation_id, clone, supporting_reads, depth)]
    clonal.clust <- sample_data[ is_clonal == TRUE , unique(clone) ]
    CCFs[, clone_mean_vaf := mean( supporting_reads / depth ), by = clone]
    CCFs[, CCF := clone_mean_vaf / CCFs[ clone == clonal.clust, unique(clone_mean_vaf) ] ]
    
    #make all clonal CCFs == 1 (this shouldn't change through iterations)
    CCFs[ clone==clonal.clust & !is.na(clone), CCF  := 1 ] 
    #don't allow CCF > 1
    CCFs[ CCF > 1 & !is.na( CCF ), CCF := 1 ]
    #for this initial CCF don't allow <5% (causes WT tumour DNA estimate to be overetimated)
    CCFs[ CCF < 0.05 & !is.na( CCF ) , CCF := 0.05] 
    
    
    #for some reason we have slightly different combpyCCFs calcualted within the same clones.. TODO .. 
    iteration <- 1
    mean.Perc.change.next <- c()
    mean.Perc.change.record <-  c()
    
    sample_dataconstant <- sample_data
    i <- 0
    
    root <- find_root( tree )
    tree <- logically.order.tree( as.matrix( tree ) )
    
    while( ( (mean.Perc.change.next > 0.01 || is.null(mean.Perc.change.next)) & any( CCFs[, !is.na(CCF) ] ) ) || i == 0 ) {
      
      i <- i + 1
      
      sample_data <- sample_dataconstant
      
      CN.adj.VAF.outputs <- cn_adjust( supporting_reads = sample_data$reads_no_background, 
                                       depth = sample_data$depth, 
                                       total_cn = sample_data$total_cpn,
                                       multiplicity = sample_data$multiplicity,
                                       ccf = CCFs[ match( CCFs$mutation_id, sample_data$mutation_id ), CCF ],
                                       return.all.intermediate.calculations = TRUE )
      
      cols.order <- names(sample_data)
      
      sample_data[, (colnames(CN.adj.VAF.outputs)) := NULL]
      
      sample_data <- cbind(sample_data,CN.adj.VAF.outputs)
      
      #estimate averageCN.adj.VAF for each clone using tail method
      
      #determine tumour fraction first for each sample 
      clone_table.clonal <- sample_data[ is_clonal == TRUE & !is.na(is_clonal),]
      clonalpopmean <- mean.estimate.LOD.correction(CN.adj.VAF = clone_table.clonal$CN.adj.VAF, 
                                                    MatchedLOD = clone_table.clonal$matched_lod,
                                                    normSD.est = normSD.for.1.sup.mut.clones, 
                                                    use.vars = !clone_table.clonal$hard_filtered)
      
      clonalpopmean <- as.numeric(as.character(unique(clonalpopmean$popmean)))
      
      if(is.na(clonalpopmean) || length(clonalpopmean)==0){
        clonalpopmean <- 0
        sample_data$Tumour.fraction <- NA
      } else {
        sample_data$Tumour.fraction <- clonalpopmean
      }
      
      sample_data <- do.call( rbind, lapply( sample_data[, unique(clone) ], function( clone_name ){ 
        
        # message( clone_name )
        
        clone_table <- sample_data[ clone %in% clone_name ]
        
        clone_name_orig <- clone_table[, unique(clone_orig)]
        
        clone_table <- cbind(clone_table, mean.estimate.LOD.correction(CN.adj.VAF = clone_table$CN.adj.VAF, 
                                                                       MatchedLOD = clone_table$matched_lod,
                                                                       normSD.est = normSD.for.1.sup.mut.clones, 
                                                                       use.vars = !clone_table$hard_filtered) )
        
        #has there been further evolution and CIN? test whether distbution is mulitmodal (use this later to correct p values)
        
        if( length(clone_table[ VAF > 0 & hard_filtered == FALSE, VAF]) > 0 & !is.na( clone_name_orig ) ){
          VAFs <- clone_table$VAF
          #make 0 the LOD
          VAFs[VAFs==0] <- clone_table[ VAF == 0, LOD ]
          # suppress warnings for regularise ties
          dip.result.CN1 <- suppressWarnings( dip.test( log10( VAFs[ clone_table$hard_filtered == FALSE & 
                                                                       clone_table$multiplicity < 1.5 ] ) )$p.value )
          if( dip.result.CN1 == 0 ) dip.result.CN1 <- 0.000000000000001 #stated min P-value by function
          if( sum(clone_table$hard_filtered == FALSE & (clone_table$multiplicity > 1.5) %in% TRUE ) > 1){
            # suppress warnings for regularise ties
            dip.result.CN2 <- suppressWarnings( dip.test(log10(VAFs[ clone_table$hard_filtered == FALSE & 
                                                                       clone_table$multiplicity > 1.5]))$p.value )
            if(dip.result.CN2==0) dip.result.CN2 <- 0.000000000000001 #stated min P-value by function
            combined.p <- sumlog(c(dip.result.CN1,dip.result.CN2))$p
          } else {
            combined.p <- dip.result.CN1
          }
          clone_table$diptestresult.p <- combined.p
          
        } else {
          clone_table$diptestresult.p <- NA
        }
        
        #calculate approximate ctDNA fraction for each clone and SD for ctDNA fraction
        clone_table$median_cn_adj_vaf <- median(clone_table$CN.adj.VAF[!clone_table$CN.adj.VAF %in% Inf],na.rm=T)
        clone_table$frac_muts_unsupported <- sum(clone_table$reads_no_background == 0) / nrow(clone_table)
        clone_table$clone.fraction.perc.at.least.two.support <- sum(clone_table$reads_no_background >= 2) / nrow(clone_table)
        clone_table$muts_followed_in_clone <- nrow(clone_table)
        clone_table$unfiltered.muts.followed.in.clone <- sum(clone_table$Target_Variant_Deep_Error_Filtered == "false")
        clone_table$no.supported.muts <- sum(clone_table$reads_no_background>0)
        clone_table$vaf_no_background <-  clone_table$reads_no_background / clone_table$depth
        clone_table$mean_clone_vaf_no_background <- mean(clone_table$VAF.no.background)
        clone_table$mean_clone_vaf <- mean(clone_table$VAF)
        
        MRD.vars <- clone_table[, hard_filtered == FALSE & mrd_filtered == FALSE | (is.na(hard_filtered) & is.na(mrd_filtered)) ]
        observed <- sum(clone_table$supporting_reads[MRD.vars])
        expected <- round(sum(clone_table$background_reads[MRD.vars]))
        if( !is.na(observed) & !is.na(expected) ){
          clone_table$Clone.detected.p <-  poisson.test(c(observed, expected), alternative =  "greater")$p.value
        } else {
          clone_table$Clone.detected.p <- NA
        }
        clone_table$Clone.meanLOD <- mean(clone_table$LOD,na.rm=T)
        
        
        if(all(clone_table$is_clonal %in% FALSE & !is.na(clone_name_orig)) & (!is.na(clonalpopmean) && length(clonalpopmean)>0)){
          
          clonal.muts <- sample_data$is_clonal == TRUE & !is.na(sample_data$clone)
          
          CCF10.expected.subclone <- model.subclone.DAOs(clonalpopmean,CCF=0.1, clone_table)
          CCF1.expected.subclone <- model.subclone.DAOs(clonalpopmean,0.01, clone_table)
          CCF50.expected.subclone <- model.subclone.DAOs(clonalpopmean,0.5, clone_table)
          
          background <- as.numeric(clone_table$background_reads)
          
          if(all(sample_data[clonal.muts,supporting_reads]==0)) no.signal <- TRUE else no.signal <- FALSE
          
          clone_table$Power.to.detect.10CCF.p <- poisson.test(sum(CCF10.expected.subclone, na.rm = T),sum(background, na.rm = T),alternative = "greater")$p.value
          clone_table$Power.to.detect.1CCF.p <- poisson.test(sum(CCF1.expected.subclone, na.rm = T), sum(background, na.rm = T), alternative = "greater")$p.value
          clone_table$Power.to.detect.50CCF.p <- poisson.test(sum(CCF50.expected.subclone, na.rm = T), sum(background, na.rm = T), alternative = "greater")$p.value
          
          CCF100.expected <- model.subclone.DAOs(clonalpopmean,1,clone_table)
          
          CCF_background_niose <- mean(background) / mean(CCF100.expected)
          if(CCF_background_niose>1 || is.na(CCF_background_niose)) CCF_background_niose <-  1
          clone_table$CCF_background_niose <- CCF_background_niose
          
          #do we have power to detect the CCF in the primary tumour at baseline?
          MseqCCF.expected <- model.subclone.DAOs(clonalpopmean, CCF = mean(clone_table$tumour_ccf), clone_name.table = clone_table)
          clone_table$Power.to.detect.MseqCCF.p <- poisson.test(sum(MseqCCF.expected, na.rm = T), sum(background, na.rm = T), alternative= "greater")$p.value
          
        } else {
          clone_table$Power.to.detect.10CCF.p <- NA
          clone_table$Power.to.detect.1CCF.p <- NA
          clone_table$Power.to.detect.50CCF.p <- NA
          clone_table$CCF_background_niose <- NA
          clone_table$Power.to.detect.MseqCCF.p <- NA
          no.signal <- FALSE
        }
        
        if(no.signal){
          clone_table$Power.to.detect.10CCF.p <- 1
          clone_table$Power.to.detect.1CCF.p <- 1
          clone_table$Power.to.detect.50CCF.p <- 1
          clone_table$Power.to.detect.MseqCCF.p <- 1
          clone_table$CCF_background_niose <- 1
        }
        
        #test whether this clone is signifcantly different to its parent
        
        parent <- as.numeric( tree[ tree[,2] %in% clone_name_orig, 1 ] )
        # if we're looking at the clonal cluster
        if( length( parent ) == 0 | all( is.na( parent )) ) parent <- NA
        
        if( !is.na( parent) ){
          
          parent.CNadjVAF <- as.numeric(sample_data[sample_data$clone == parent, CN.adj.VAF ])
          parent.CNadjVAF <- parent.CNadjVAF[!is.na(parent.CNadjVAF)]
          clone.CNadjVAF <- as.numeric(clone_table$CN.adj.VAF)
          clone.CNadjVAF <- clone.CNadjVAF[!is.na(clone.CNadjVAF)]
          
          if(length(parent.CNadjVAF) > 1 & length(clone.CNadjVAF) > 1){
            # supress cannot compute exact p-value with ties warning
            suppressWarnings( clone_table[, different.to.parent.p := wilcox.test(parent.CNadjVAF,clone.CNadjVAF)$p.value ] )
            clone_table[, different.to.parent.OR := mean(clone.CNadjVAF) / mean(parent.CNadjVAF) ]
          } else {
            
            clone_table[, different.to.parent.p := 1 ]
            clone_table[, different.to.parent.OR := mean(clone.CNadjVAF) / mean(parent.CNadjVAF) ]
          }
        } else {
          clone_table[, different.to.parent.p := NA ]
          clone_table[, different.to.parent.OR := NA ]
        }
        
        return( clone_table )
        
      } ))
      
      MRD.vars <- sample_data[, (hard_filtered == FALSE & mrd_filtered == FALSE | 
                                   (is.na(hard_filtered) & is.na(mrd_filtered)) ) & 
                                !is.na(supporting_reads) & !is.na(background_reads) ]
      
      sample_data$MRD.call.sample.p <- poisson.test(sum(sample_data$supporting_reads[MRD.vars]),sum(sample_data$background_reads[MRD.vars]), alternative = "greater")$p.value
      
      #if no clonal vars set to NA (ie if been unable to run pyclone)
      if(any(sample_data$is_clonal== TRUE)){
        clonalpopmean <- as.numeric(as.character(unique(sample_data[sample_data$is_clonal== TRUE, ctDNA_fraction ])))
      } else {
        clonalpopmean <- NA
      }
      #supress warnings when SD solution is not found (characters)
      
      sample_data$plasma_CCF <- suppressWarnings( as.numeric(as.character(sample_data$ctDNA_fraction)) /  clonalpopmean )
      sample_data$plasma_CCF_UpCI <- suppressWarnings( as.numeric(as.character(sample_data$ctDNA_fraction_UpCI)) /  clonalpopmean )
      sample_data$plasma_CCF_LwCI <- suppressWarnings( as.numeric(as.character(sample_data$ctDNA_fraction_LwCI)) /  clonalpopmean )
      
      ##quantify how much popmean changes from last iteration with popmean table
      if(iteration == 1){
        
        meanpoptable <- sample_data[!duplicated(sample_data$clone), c("clone","is_clonal")]
        # supress warning when popmean is actually NA (not CN etc information)
        meanpoptable <- suppressWarnings( cbind(meanpoptable, signif( as.numeric(as.character( sample_data[!duplicated(clone), ctDNA_fraction] )), 5) ) )
        names(meanpoptable)[ncol(meanpoptable)] <- paste0("popmean_iteration_",iteration)
        
      } else {
        #suppress warnings when SD solution is not found (characters)
        meanpoptable <- suppressWarnings( cbind( meanpoptable, signif(as.numeric(as.character(sample_data[!duplicated(clone), ctDNA_fraction])),5) ) )
        names(meanpoptable)[ncol(meanpoptable)] <- paste0("popmean_iteration_",iteration)
        mean.Perc.change.next <- mean(abs( as.numeric(meanpoptable[[ (ncol(meanpoptable) - 1) ]]) - as.numeric(meanpoptable[[ ncol(meanpoptable) ]]) ) / as.numeric(meanpoptable[[ (ncol(meanpoptable) - 1) ]]), na.rm=T )
        mean.Perc.change.record <- c(mean.Perc.change.record,mean.Perc.change.next)  #to record changes over iteration when testing
      }
      
      iteration <- iteration + 1
      
      CCFs <- sample_data[,c("mutation_id","clone","plasma_CCF")]
      
      if(iteration==20){
        cat(paste(sample,"could not find CCF solution for all subclones "))
        badclones <- meanpoptable$clone[abs( as.numeric(meanpoptable[,(ncol(meanpoptable) - 1)]) - as.numeric(meanpoptable[,ncol(meanpoptable)]) ) / as.numeric(meanpoptable[,ncol(meanpoptable) - 1])>0.01 & !is.na(meanpoptable[,ncol(meanpoptable)])]
        CCFs[sample_data$clone %in% badclones,"plasma_CCF"] <- NA
        break
      }
      
      names(CCFs)[3] <- "CCF"
    }
    
    
    # call mrd for this sample
    #For samples without +MRD set the CCFs and popmeans etc to 0
    out_cols <- c("ctDNA_fraction","ctDNA_fraction_UpCI","ctDNA_fraction_LwCI","plasma_CCF","plasma_CCF_UpCI","plasma_CCF_LwCI")
    suppressWarnings( sample_data[, (out_cols) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = out_cols] )
    sample_data[!sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p),c("ctDNA_fraction","ctDNA_fraction_UpCI","ctDNA_fraction_LwCI","plasma_CCF","plasma_CCF_UpCI","plasma_CCF_LwCI")] <- 0
    
    ##correct all p values
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Clone.detected.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p),method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Power.to.detect.10CCF.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Power.to.detect.1CCF.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Power.to.detect.50CCF.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Power.to.detect.MseqCCF.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="different.to.parent.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    
      
  return(sample_data)
  
}


#' Function to remove the number of estimated supporting reads which derive from
#' background sequencing niose for groups of mutations (in this case clones)
#' @export
remove_background_reads <- function( data., niose_col = 'background_error',
                                     filter_col = 'hard_filtered', support_col = 'supporting_reads',
                                     depth_col = 'depth', group_max = 4, mutid_col = 'mutation_id'){
  
  # Order by this - then don't need randomisation in sampling and this is deterministic & reproducible
  data. <- data.[ order(get(mutid_col)) ] 
  
  # mutations to ignore throughout
  ignore_i <- data.[, get(filter_col) ]
  
  # If all mutations are invaluable then return with NA
  if(all(ignore_i)){
    data.[, `:=`(group = NA,
                 supporting_reads_no_niose = NA) ]
    return(data.)
  } 
  
  # Determine the number of groups so that each group has enough mutations
  # that if the niose was distrebuted evenly each group would have at least 1
  # read of niose - max of 4 groups
  total_niose <- data.[ !ignore_i , floor( sum( get(niose_col) * get(depth_col) ) )]
  group_num = total_niose
  if( group_num > group_max ) group_num = group_max
  if( group_num == 0 ) group_num = 1
  
  # Now cluster mutations based on background noise 
  data.[ !ignore_i, group := ifelse( .N > 1, kmeans(data.[ !ignore_i, get(niose_col)], group_num)$cluster, 1) ]
  groups <- data.[ !is.na(group), unique(group)]
  
  group_read_to_remove <- lapply(groups, function(group_name){
    background.lambda <- data.[ group == group_name, mean(get(niose_col) * get(depth_col)) ]
    max.read.support.to.test <- round( background.lambda * 100 )
    reads.to.remove <- round( dpois(1:max.read.support.to.test, lambda = background.lambda) * data.[ group == group_name, .N] )
    names(reads.to.remove) <- 1:max.read.support.to.test
    out <- reads.to.remove[ !reads.to.remove == 0 ] 
    if(length(out) == 0) out <- NA
    return( out )
  } )
  
  names(group_read_to_remove) <- groups

  subtract_reads <- function( reads, reads.to.remove ){
    #if no reads to remove return read unaltered
    if( all(is.na(reads.to.remove)) ) return(reads)
    
    # record which mutations have already had niose removed 
    # (only do this once to each mutation)
    already_subtracted = FALSE
    reads <- data.table( reads, already_subtracted)
    
    #if no reads to remove return read unaltered
    if( all(is.na(reads.to.remove)) ) return(reads)
    
    # Loop around each read number which needs removing
    for( read_num in rev( names(reads.to.remove) ) ){

      if( all(reads$reads == 0) ) return(reads$reads)
      
      # which mutations have enough reads that we could remove the desired number
      reads[, is_eligible := reads >= read_num ]
      posiblei <- reads[, which(is_eligible & !already_subtracted) ]
      
      # If few mutations available then need to remove and spread the reads to lower read numbers
      # to lower numbers
      if(length(posiblei) < reads.to.remove[read_num] ){
        
        # We still have excess at 0 this means there are more background estimated reads than
        # observed reads - return no support for all mutations
        if(read_num == '1') return( rep(0, nrow(reads)) )
        
        # work out how we distribute the excess - use the next number below first then add any remainder
        num_excess <- reads.to.remove[read_num] - length(posiblei)
        reads_excess <- num_excess * as.numeric(read_num)
        read_num_below <- as.numeric(read_num) - 1 
        read_num_below_add <- floor( reads_excess / read_num_below)
        remainder <- reads_excess %% read_num_below
        if( any(names(reads.to.remove) == as.character(read_num_below)) ){
          reads.to.remove[[ as.character(read_num_below) ]] <- reads.to.remove[[ as.character(read_num_below) ]] + read_num_below_add
        } else {
          reads.to.remove[[ as.character(read_num_below) ]] <- read_num_below_add
          reads.to.remove <- reads.to.remove[ order( as.numeric(names(reads.to.remove)))]
        }
        if( remainder > 0 ){
          if( any(names(reads.to.remove) == as.character(remainder)) ){
            reads.to.remove[[ as.character(remainder) ]] <- reads.to.remove[[ as.character(remainder) ]] + 1
          } else {
            reads.to.remove[[ as.character(remainder) ]] <- 1
            reads.to.remove <- reads.to.remove[ order( as.numeric(names(reads.to.remove)))]
          }
        }
        
        # Now we've redistributed the excess remove the excess
        reads.to.remove[read_num] = length(posiblei)
      }
      #if no possible mutations to add just skip to the next level
      if( length(posiblei) == 0 ) next
      
      # Randomly select the mutations for which read are to be removed 
      # and mark so this only occurs once per mutation
      # don't use sample/runif - just select the first X - as we order by sample then mutation id this
      # should make it deterministic
      selectedi <-posiblei[ 1:reads.to.remove[read_num] ]
      reads[selectedi, already_subtracted := TRUE ]
      reads[selectedi, reads := reads - as.numeric(read_num) ]
    }
    return( reads$reads )
  }
  
  data.[ !ignore_i, 
         supporting_reads_no_niose := subtract_reads( get(support_col), group_read_to_remove[[.GRP]]), 
         by = group ]
  
  return(data.)
  
}

#' Function to annotate in each clone mutations which mutations are unsupported
#' because of a low limit of detection (LOD) ie the depth is too low or the mutation 
#' copy is lower than other mutations etc
#' This is used when modelling the distbution including the unobserved variants
#' @export
annotate_low_LOD <- function( data., niose_col = 'background_error', LOD_col = 'LOD',
                              filter_col = 'hard_filtered', support_col = 'supporting_reads',
                              depth_col = 'depth', mutid_col = 'mutation_id'){
  
  # mutations to ignore throughout
  ignore_i <- data.[, get(filter_col) | is.na(get(LOD_col)) ]
  
  if( data.[ !ignore_i, all(get(support_col) > 0) ] ) return(data.)
  
  # get the median LOD for mutations that have read support
  medianLOD_supported <- data.[ !ignore_i & get(support_col) > 0, median(get(LOD_col)) ]
  
  # get the median of the unsupported variants while removing one at a time
  unsupported_lod <- data.[ !ignore_i & get(support_col) == 0, get(LOD_col) ]
  unsupported_lod <- unsupported_lod[ order(unsupported_lod, decreasing = TRUE) ]
  median_lod_unsupported_range <- sapply( 1:length(unsupported_lod), function(i) median(unsupported_lod[i:length(unsupported_lod)]) )
  
  # how many unsupported variants to discount
  diff <- abs(median_lod_unsupported_range - unsupported_lod )
  num_supported_to_remove <- which( diff == min(diff) )
  mutids_to_remove <- data.[ get(support_col) == 0 ][ order(get(LOD_col), decreasing = TRUE), get(mutid_col) ][ 1:num_supported_to_remove ]
  
  # assign this to LOD unmatched
  data.[, matched_lod := TRUE ]
  data.[ get(mutid_col) %in% mutids_to_remove, matched_lod := FALSE ]
  
  return(data.)
}

#' Function to calculate CCFs from vaf, purity and cn as in ABSOLUTE
#' @export
calculate_ccfs <- function( data., purity, varcount_col = 'supporting_reads', depth_col = 'depth', 
                           tumour_totalCN_col = 'total_cpn', normal_totalCN_col = normal_totalCN ){
  
  calculate_ccf <- function(varcount, depth, purity, tumour_totalCN, normal_totalCN ){
  
    ccf_est <- function (vaf, purity, tumour_totalCN, normal_totalCN){
    return(min(1,c((purity*vaf) / (normal_totalCN*(1-purity) + purity*tumour_totalCN))))
  }
  
  test_vafs <- c( seq(0.0001, 0.001, 0.0001), seq(0.001, 0.01, 0.001), seq(0.01, 1, 0.01) )
  binom_dist    <- dbinom(n.alt,depth, prob=sapply(test_vafs, ccf_est, purity, local.copy.number, normal.copy.number ))
  names(x) <- test_vafs
  
  if(min(binom_dist)==0) x[length(x)] <- 1
  
  
  sub.cint <- function(x, prob = 0.95,n.alt,depth) {
    xnorm   <- x/sum(x)
    xsort   <- sort(xnorm, decreasing = TRUE)
    xcumLik <- cumsum(xsort)
    n = sum(xcumLik < prob) + 1
    LikThresh <- xsort[n]
    cint  <- x[xnorm >= LikThresh]
    all   <- as.numeric(names(x))
    cellu <- as.numeric(names(cint))
    l.t   <- cellu[1]
    r.t   <- cellu[length(cellu)]
    m     <- cellu[which.max(cint)]
    
    prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
    prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val
    
    data.frame(left = l.t, est = m, right = r.t,prob.subclonal=prob.subclonal,prob.clonal=prob.clonal)
  }
  
  return(sub.cint(x,n.alt=n.alt,depth=depth))
  }
}

#' Clonal deconvolution function
#' 
#' 
#' 
#' @export
clonal_deconvolution <- function(data, hard_filters = NA, mrd_filters = NA, tree = NA,
                                 lod_correction = TRUE, test_subsequent_CIN = TRUE, min_muts_clone = 2, 
                                 normSD.for.1.sup.mut.clones){
  
  class_origin <- class( data )
  data <- as.data.table( data )
  
  
  # make a mutation id
  data[, mutation_id := paste( sample_id, chr, pos, alt, sep = ":" ) ]
  
  # allow to run for multiple samples if inputted as 1 data table
  data[, clone_orig := clone ]
  data[, clone := paste( sample_id, clone, sep = "_") ]
  
  # work out which mutations are good for MRD and which should be removed altogether
  data[, hard_filtered := grepl( paste( hard_filters, collapse = "|" ), error_filter ) | is.na(background_error) ]
  data[, mrd_filtered := grepl( paste( mrd_filters, collapse = "|" ), error_filter ) | is.na(background_error) ]
  if( all( is.na(hard_filters) ) )  data[, hard_filtered := FALSE ]
  if( all( is.na(mrd_filters) ) )  data[, mrd_filtered := FALSE ]
  
  # if no background info hard filter
  data[ is.na(background_error), hard_filtered := TRUE ]
  
  # filter background niose reads
  data <- rbindlist( lapply(data[, unique(clone) ], function(clone_name) remove_background_reads(data[ clone == clone_name ]) ) )
  
  
  # calculate VAF
  data[, VAF := supporting_reads / depth]

  # calculate raw mutant CN, equation taken from NEJM TRACERx 100 paper (this is an 
  # average mut cn per cell across the tumour)
  data[, mutCPN := (tumour_vaf * (1 / tumour_cellularity)) * ((tumour_cellularity * total_cpn) + 2 * (1 - tumour_cellularity)) ]
  
  # calculate number of mutant copies per mutated cell (i/e. multiplicity of mutation)
  data[, multiplicity := mutCPN / tumour_ccf ]
  
  # multiplicity < 1 is not possible (causes by noise in data) - correct this
  data[ multiplicity < 1, multiplicity := 1 ]
  
  # assign all mutations as matched LOD for now - will use this later
  data$matched_lod <- TRUE
  
  # calculate the limit of detection for each variant based on cellularity that = 1 supporting read
  # will only compare LOD within a clone so CCF can be any value here
  data[, LOD := cn_adjust( supporting_reads = 1 , depth, total_cpn, multiplicity, ccf = 1 ) ]
  
  # Annotate which variants in each clones that appear to be unsupported becuase of
  # a low LOD
  if(lod_correction == TRUE){
    data <- rbindlist( lapply(data[, unique(clone) ], function(clone_name) annotate_low_LOD(data[ clone == clone_name ]) ) )
  }
  
  ##### Calculate ctDNA cellularity accounting for mutCPN, WT cancer copies #####

  
  data <- do.call( rbind, lapply( data[, unique( sample_id) ], function( sample ){
    
    message( sample )
    
    sample_data <- data[ sample_id == sample ]
    
    #need approximate CCFs to calculate WT tumour DNA contribution (will iterate to correct CCFs)
    #use meanVAF for this initially
    CCFs <- sample_data[, .( mutation_id, clone, supporting_reads, depth)]
    clonal.clust <- sample_data[ is_clonal == TRUE , unique(clone) ]
    CCFs[, clone_mean_vaf := mean( supporting_reads / depth ), by = clone]
    CCFs[, CCF := clone_mean_vaf / CCFs[ clone == clonal.clust, unique(clone_mean_vaf) ] ]
    
    #make all clonal CCFs == 1 (this shouldn't change through iterations)
    CCFs[ clone==clonal.clust & !is.na(clone), CCF  := 1 ] 
    #don't allow CCF > 1
    CCFs[ CCF > 1 & !is.na( CCF ), CCF := 1 ]
    #for this initial CCF don't allow <5% (causes WT tumour DNA estimate to be overetimated)
    CCFs[ CCF < 0.05 & !is.na( CCF ) , CCF := 0.05] 
    
    
    #for some reason we have slightly different combpyCCFs calcualted within the same clones.. TODO .. 
    iteration <- 1
    mean.Perc.change.next <- c()
    mean.Perc.change.record <-  c()
    
    sample_dataconstant <- sample_data
    i <- 0
    
    root <- find_root( tree )
    tree <- logically.order.tree( as.matrix( tree ) )
    
    while( ( (mean.Perc.change.next > 0.01 || is.null(mean.Perc.change.next)) & any( CCFs[, !is.na(CCF) ] ) ) || i == 0 ) {
      
      i <- i + 1
      
      sample_data <- sample_dataconstant
      
      CN.adj.VAF.outputs <- cn_adjust( supporting_reads = sample_data$reads_no_background, 
                                       depth = sample_data$depth, 
                                       total_cn = sample_data$total_cpn,
                                       multiplicity = sample_data$multiplicity,
                                       ccf = CCFs[ match( CCFs$mutation_id, sample_data$mutation_id ), CCF ],
                                       return.all.intermediate.calculations = TRUE )
      
      sample_data <- cbind(sample_data,CN.adj.VAF.outputs)
      
      #estimate averageCN.adj.VAF for each clone using tail method
      
      #determine tumour fraction first for each sample 
      clone_table.clonal <- sample_data[ is_clonal == TRUE & !is.na(is_clonal),]
      clonalpopmean <- mean.estimate.LOD.correction(CN.adj.VAF = clone_table.clonal$CN.adj.VAF, 
                                                    MatchedLOD = clone_table.clonal$matched_lod,
                                                    normSD.est = normSD.for.1.sup.mut.clones, 
                                                    use.vars = !clone_table.clonal$hard_filtered)
      
      clonalpopmean <- as.numeric(as.character(unique(clonalpopmean$ctDNA_fraction)))
      
      if(is.na(clonalpopmean) || length(clonalpopmean)==0){
        clonalpopmean <- 0
        sample_data$Tumour.fraction <- NA
      } else {
        sample_data$Tumour.fraction <- clonalpopmean
      }
      
      sample_data <- do.call( rbind, lapply( sample_data[, unique(clone) ], function( clone_name ){ 
       
        message( clone_name )
        
        clone_table <- sample_data[ clone %in% clone_name ]
        
        clone_name_orig <- clone_table[, unique(clone_orig)]
        
        clone_table <- cbind(clone_table, mean.estimate.LOD.correction(CN.adj.VAF = clone_table$CN.adj.VAF, 
                                                                       MatchedLOD = clone_table$matched_lod,
                                                                       normSD.est = normSD.for.1.sup.mut.clones, 
                                                                       use.vars = !clone_table$hard_filtered) )
        
        #has there been further evolution and CIN? test whether distbution is mulitmodal (use this later to correct p values)
        
        if( length(clone_table[ VAF > 0 & hard_filtered == FALSE, VAF]) > 0 & !is.na( clone_name_orig ) ){
          VAFs <- clone_table$VAF
          #make 0 the LOD
          VAFs[ VAFs==0 & !is.na(VAFs)] <- clone_table[ VAF == 0, LOD ]
          # suppress warnings for regularise ties
          dip.result.CN1 <- suppressWarnings( dip.test( log10( VAFs[ clone_table$hard_filtered == FALSE & 
                                                                       clone_table$multiplicity < 1.5 ] ) )$p.value )
          if( dip.result.CN1 == 0 ) dip.result.CN1 <- 0.000000000000001 #stated min P-value by function
          if( sum(clone_table$hard_filtered == FALSE & (clone_table$multiplicity > 1.5) %in% TRUE ) > 1){
            # suppress warnings for regularise ties
            dip.result.CN2 <- suppressWarnings( dip.test(log10(VAFs[ clone_table$hard_filtered == FALSE & 
                                                                       clone_table$multiplicity > 1.5]))$p.value )
            if(dip.result.CN2==0) dip.result.CN2 <- 0.000000000000001 #stated min P-value by function
            combined.p <- sumlog(c(dip.result.CN1,dip.result.CN2))$p
          } else {
            combined.p <- dip.result.CN1
          }
          clone_table$diptestresult.p <- combined.p
          
        } else {
          clone_table$diptestresult.p <- NA
        }
        
        #calculate approximate ctDNA fraction for each clone and SD for ctDNA fraction
        clone_table$median_cn_adj_vaf <- median(clone_table$CN.adj.VAF[!clone_table$CN.adj.VAF %in% Inf],na.rm=T)
        clone_table$frac_muts_unsupported <- sum(clone_table$reads_no_background == 0) / nrow(clone_table)
        clone_table$clone.fraction.perc.at.least.two.support <- sum(clone_table$reads_no_background >= 2) / nrow(clone_table)
        clone_table$muts_followed_in_clone <- nrow(clone_table)
        clone_table$unfiltered.muts.followed.in.clone <- sum(clone_table$Target_Variant_Deep_Error_Filtered == "false")
        clone_table$no.supported.muts <- sum(clone_table$reads_no_background>0)
        clone_table$vaf_no_background <-  clone_table$reads_no_background / clone_table$depth
        clone_table$mean_clone_vaf_no_background <- mean(clone_table$vaf_no_background)
        clone_table$mean_clone_vaf <- mean(clone_table$VAF)
        
        MRD.vars <- clone_table[, hard_filtered == FALSE & mrd_filtered == FALSE | (is.na(hard_filtered) & is.na(mrd_filtered)) ]
        observed <- sum(clone_table$supporting_reads[MRD.vars])
        expected <- round(sum(clone_table$background_reads[MRD.vars]))
        if( !is.na(observed) & !is.na(expected) ){
          clone_table$Clone.detected.p <-  poisson.test(c(observed, expected), alternative =  "greater")$p.value
        } else {
          clone_table$Clone.detected.p <- NA
        }
        clone_table$Clone.meanLOD <- mean(clone_table$LOD,na.rm=T)
        
        
        if(all(clone_table$is_clonal %in% FALSE & !is.na(clone_name_orig)) & (!is.na(clonalpopmean) && length(clonalpopmean)>0)){
          
          clonal.muts <- sample_data$is_clonal == TRUE & !is.na(sample_data$clone)
          
          CCF10.expected.subclone <- model.subclone.DAOs(clonalpopmean,CCF=0.1, clone_table)
          CCF1.expected.subclone <- model.subclone.DAOs(clonalpopmean,0.01, clone_table)
          CCF50.expected.subclone <- model.subclone.DAOs(clonalpopmean,0.5, clone_table)
          
          background <- as.numeric(clone_table$background_reads)
          
          if(all(sample_data[clonal.muts,supporting_reads]==0)) no.signal <- TRUE else no.signal <- FALSE
          
          clone_table$Power.to.detect.10CCF.p <- poisson.test(sum(CCF10.expected.subclone, na.rm = T),sum(background, na.rm = T),alternative = "greater")$p.value
          clone_table$Power.to.detect.1CCF.p <- poisson.test(sum(CCF1.expected.subclone, na.rm = T), sum(background, na.rm = T), alternative = "greater")$p.value
          clone_table$Power.to.detect.50CCF.p <- poisson.test(sum(CCF50.expected.subclone, na.rm = T), sum(background, na.rm = T), alternative = "greater")$p.value
          
          CCF100.expected <- model.subclone.DAOs(clonalpopmean,1,clone_table)
          
          CCF_background_niose <- mean(background) / mean(CCF100.expected)
          if(CCF_background_niose>1 || is.na(CCF_background_niose)) CCF_background_niose <-  1
          clone_table$CCF_background_niose <- CCF_background_niose
          
          #do we have power to detect the CCF in the primary tumour at baseline?
          MseqCCF.expected <- model.subclone.DAOs(clonalpopmean, CCF = mean(clone_table$tumour_ccf), clone_name.table = clone_table)
          clone_table$Power.to.detect.MseqCCF.p <- poisson.test(sum(MseqCCF.expected, na.rm = T), sum(background, na.rm = T), alternative= "greater")$p.value
          
        } else {
          clone_table$Power.to.detect.10CCF.p <- NA
          clone_table$Power.to.detect.1CCF.p <- NA
          clone_table$Power.to.detect.50CCF.p <- NA
          clone_table$CCF_background_niose <- NA
          clone_table$Power.to.detect.MseqCCF.p <- NA
          no.signal <- FALSE
        }
        
        if(no.signal){
          clone_table$Power.to.detect.10CCF.p <- 1
          clone_table$Power.to.detect.1CCF.p <- 1
          clone_table$Power.to.detect.50CCF.p <- 1
          clone_table$Power.to.detect.MseqCCF.p <- 1
          clone_table$CCF_background_niose <- 1
        }
        
        #test whether this clone is signifcantly different to its parent
        
        parent <- as.numeric( tree[ tree[,2] %in% clone_name_orig, 1 ] )
        # if we're looking at the clonal cluster
        if( length( parent ) == 0 | all( is.na( parent )) ) parent <- NA
        
        if( !is.na( parent) ){
          
          parent.CNadjVAF <- as.numeric(sample_data[sample_data$clone_orig == parent, CN.adj.VAF ])
          parent.CNadjVAF <- parent.CNadjVAF[!is.na(parent.CNadjVAF)]
          clone.CNadjVAF <- as.numeric(clone_table$CN.adj.VAF)
          clone.CNadjVAF <- clone.CNadjVAF[!is.na(clone.CNadjVAF)]
          
          if(length(parent.CNadjVAF) > 1 & length(clone.CNadjVAF) > 1){
            # supress cannot compute exact p-value with ties warning
            suppressWarnings( clone_table[, different.to.parent.p := wilcox.test(parent.CNadjVAF,clone.CNadjVAF)$p.value ] )
            clone_table[, different.to.parent.OR := mean(clone.CNadjVAF) / mean(parent.CNadjVAF) ]
          } else {
            
            clone_table[, different.to.parent.p := 1 ]
            clone_table[, different.to.parent.OR := mean(clone.CNadjVAF) / mean(parent.CNadjVAF) ]
          }
        } else {
          clone_table[, different.to.parent.p := NA ]
          clone_table[, different.to.parent.OR := NA ]
        }
        
        return( clone_table )
        
      } ))
      
      MRD.vars <- sample_data[, (hard_filtered == FALSE & mrd_filtered == FALSE | 
                                   (is.na(hard_filtered) & is.na(mrd_filtered)) ) & 
                                !is.na(supporting_reads) & !is.na(background_reads) ]
      
      sample_data$MRD.call.sample.p <- poisson.test(sum(sample_data$supporting_reads[MRD.vars]),sum(sample_data$background_reads[MRD.vars]), alternative = "greater")$p.value
      
      #if no clonal vars set to NA (ie if been unable to run pyclone)
      if(any(sample_data$is_clonal== TRUE)){
        clonalpopmean <- as.numeric(as.character(unique(sample_data[sample_data$is_clonal== TRUE, ctDNA_fraction ])))
      } else {
        clonalpopmean <- NA
      }
      #supress warnings when SD solution is not found (characters)
      
      sample_data$plasma_CCF <- suppressWarnings( as.numeric(as.character(sample_data$ctDNA_fraction)) /  clonalpopmean )
      sample_data$plasma_CCF_UpCI <- suppressWarnings( as.numeric(as.character(sample_data$ctDNA_fraction_UpCI)) /  clonalpopmean )
      sample_data$plasma_CCF_LwCI <- suppressWarnings( as.numeric(as.character(sample_data$ctDNA_fraction_LwCI)) /  clonalpopmean )
      
      ##quantify how much popmean changes from last iteration with popmean table
      if(iteration == 1){
        
        meanpoptable <- sample_data[!duplicated(sample_data$clone), c("clone","is_clonal")]
        # supress warning when popmean is actually NA (not CN etc information)
        meanpoptable <- suppressWarnings( cbind(meanpoptable, signif( as.numeric(as.character( sample_data[!duplicated(clone), ctDNA_fraction] )), 5) ) )
        names(meanpoptable)[ncol(meanpoptable)] <- paste0("popmean_iteration_",iteration)
        
      } else {
        #suppress warnings when SD solution is not found (characters)
        meanpoptable <- suppressWarnings( cbind( meanpoptable, signif(as.numeric(as.character(sample_data[!duplicated(clone),ctDNA_fraction])),5) ) )
        names(meanpoptable)[ncol(meanpoptable)] <- paste0("popmean_iteration_",iteration)
        mean.Perc.change.next <- mean(abs( as.numeric(meanpoptable[[ (ncol(meanpoptable) - 1) ]]) - as.numeric(meanpoptable[[ ncol(meanpoptable) ]]) ) / as.numeric(meanpoptable[[ (ncol(meanpoptable) - 1) ]]), na.rm=T )
        mean.Perc.change.record <- c(mean.Perc.change.record,mean.Perc.change.next)  #to record changes over iteration when testing
      }
      
      iteration <- iteration + 1
      
      CCFs <- sample_data[,c("mutation_id","clone","plasma_CCF")]
      
      if(iteration==20){
        cat(paste(sample,"could not find CCF solution for all subclones "))
        badclones <- meanpoptable$clone[abs( as.numeric(meanpoptable[,(ncol(meanpoptable) - 1)]) - as.numeric(meanpoptable[,ncol(meanpoptable)]) ) / as.numeric(meanpoptable[,ncol(meanpoptable) - 1])>0.01 & !is.na(meanpoptable[,ncol(meanpoptable)])]
        CCFs[sample_data$clone %in% badclones,"plasma_CCF"] <- NA
        break
      }
      names(CCFs)[3] <- "CCF"
    }
    
    
    # call mrd for this sample
    #For samples without +MRD set the CCFs and popmeans etc to 0
    out_cols <- c("ctDNA_fraction","ctDNA_fraction_UpCI","ctDNA_fraction_LwCI","plasma_CCF","plasma_CCF_UpCI","plasma_CCF_LwCI")
    suppressWarnings( sample_data[, (out_cols) := lapply(.SD,function(x) as.numeric(as.character(x))), .SDcols = out_cols] )
    sample_data[!sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p),
                c("ctDNA_fraction","ctDNA_fraction_UpCI","ctDNA_fraction_LwCI","plasma_CCF","plasma_CCF_UpCI","plasma_CCF_LwCI")] <- 0
    
    ##correct all p values
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Clone.detected.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p),method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Power.to.detect.10CCF.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Power.to.detect.1CCF.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Power.to.detect.50CCF.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="Power.to.detect.MseqCCF.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    sample_data <- correct.specific.ps(sample_data,field.to.deduplicate="clone",field.to.correct="different.to.parent.p", indices.to.correct= sample_data$MRD.call.sample.p < 0.05 & !is.na(sample_data$MRD.call.sample.p) & sample_data$Tumour.fraction > 0.0001 & !is.na(sample_data$Tumour.fraction), method ="BH")
    
    #Look at cases with possible CIN evoluation (signifcant dip.tests)- confirm with genomic localisation of losses (TO DO - tried & this is v difficult). (need to do outside loop to generate q-values)
    #then assign hom and het deletions and possible amplifications, recalculate popmean and CCFs based on muts incldue the correct CNAs
    # if(test_subsequent_CIN == TRUE){
    #   sample_data <- mulitmodal.corrections(sample_data)
    # }
    
    # don't include those with NA ctDNA fractions as these are poor quality clones without a cohesive VAF distribution
    # presence / absence of clones based on very few mutations can be unreliable (eg can't assess distrebution well) - parameter here to
    # only consider mutations with at least a certain number of mutations followed in ctDNA - default is 2
    detected <- sample_data[sample_data$Clone.detected.q<0.1 & !is.na(sample_data$ctDNA_fraction) & muts_followed_in_clone >= min_muts_clone,]
    
    if(any(detected$Clone.detected.q < 0.1 & !is.na(detected$Clone.detected.q))){
      
      clones.detected <- unique(detected$clone_orig)[ !is.na(unique(detected$clone_orig)) ]
      
      #if all CN fail can have no PyClone
      if(!length(clones.detected)==0){
        
        detected.tree <- remove.clones.on.tree( tree = tree, clones.to.keep = clones.detected )
        #does the tree branch?
        branching <- any(duplicated(detected.tree[,1]))
        
        if(all(clones.detected %in% root )){
          out <- "Clonal.only"
        } else {
          if(branching==TRUE){
            out <- "Polyphyletic"
          } else {
            if(any(detected$different.to.parent.q < 0.1 & 
                   !is.na(detected$different.to.parent.q) &  
                   detected$different.to.parent.OR <1)) {
              out <- "Polyclonal"
            } else {
              out <- "Monoclonal"
            }
          }
        }
      } else {
        out <- NA
      }
    } else {
      out <- NA
    }
    
    sample_data$ctDNA_clonal_structure <- out
    
    return( sample_data )
    
  }))
  
  data <- as.data.table(data)
  
  data[, clone := clone_orig ]
  data[, clone_orig := NULL ]
  
  return(data)
  
}

#=====#
# END #
#=====#

