
##==================================##
##                                  ##
## Functions to process each sample ## 
##                                  ##
##==================================##

#' Function to simulate a log normal distribution while defining true mean and sd
#' 
#' 
#'  
#'    
#' 
random_lognormal <- function(mean, sd, n){
  
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
#' 
model_subclone <- function(clonalpopmean, CCF, sampleclone.table){
  VAF.CCF <-  clonalpopmean * CCF
  CCF.expected.subclone <- make.random.lognormal(mean=VAF.CCF,sd=0.75*VAF.CCF,n = nrow(sampleclone.table))  * sampleclone.table$Deep_Depth
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
#' 
#Calculate the amount of WT DNA (VAF) is actually coming from the tumour or normal cells
cn_adjust <- function(Varcount, Depth, TotalCN, cloneMutCPN, CCF, is.Clonal = NA, 
                      return.all.intermediate.calculations = FALSE, cnEVO = FALSE, clone.MutCell = NA){
  
  # plasma CCFs can be sometimes estimated as > 1, impossible hence limit tp 1
  CCF[ CCF > 1 ] <- 1
  
  Extra.WT.copies <- TotalCN - (cloneMutCPN * CCF)
  Extra.WT.copies[Extra.WT.copies<0] <- 0
  
  if(cnEVO == FALSE){
    
    reads.per.totalMutCPN <- Varcount / (cloneMutCPN * CCF)
    
  } else {
    
    reads.per.totalMutCPN <- clone.MutCell * Depth
    
  }
  
  # #reads.per.mutCN should be equal accross clones - average accross (allow better calculations when looks like subsequent amplification has occured)
  # if(!is.na(is.Clonal)){
  # reads.per.totalMutCPN <- sapply(1:length(reads.per.totalMutCPN), function(i) median(reads.per.totalMutCPN[is.Clonal],na.rm = T))
  # }
  
  WT.tumour.reads <- Extra.WT.copies * reads.per.totalMutCPN
  
  #don't allow moree WT reads than depth-(mutCN * varount) - this suggests you are underestimating the number of MT copies (has theere been an amplification?)
  WT.overest <- WT.tumour.reads>(Depth - (Varcount / cloneMutCPN)) & !is.na(WT.tumour.reads)
  WT.tumour.reads[WT.overest] <- (Depth - (Varcount / cloneMutCPN))[WT.overest]
  
  mut.Cell.no.WT.Tumour.DNA.correction <- (Varcount / cloneMutCPN) / (Depth + Varcount)
  
  #calculate "cellularity" by getting VAF per cell (/Mtcopies) and deviding this by VAF signal from normal cells
  Corrected.Depth.2tumour.copies <- (Depth - (reads.per.totalMutCPN*(2-cloneMutCPN)) - WT.tumour.reads)
  Corrected.Depth.1tumour.copy <- (Depth - (reads.per.totalMutCPN*(1-cloneMutCPN)) - WT.tumour.reads)
  Corrected.Depth.for.weight <- (Depth - (reads.per.totalMutCPN*(TotalCN-cloneMutCPN)) - WT.tumour.reads)
  
  Cellularity.eqivilent <- 2 * (Varcount / cloneMutCPN) / Corrected.Depth.2tumour.copies
  CN.adj.VAF <- (Varcount / cloneMutCPN) / Corrected.Depth.1tumour.copy
  Tumour.DNA.weight.fraction <- TotalCN * (Varcount / cloneMutCPN) / Corrected.Depth.for.weight
  
  if(return.all.intermediate.calculations == TRUE){
    
    return(cbind(CN.adj.VAF,
                 Cellularity.eqivilent,
                 Tumour.DNA.weight.fraction,
                 Extra.WT.copies,
                 reads.per.totalMutCPN,
                 WT.tumour.reads,
                 Corrected.Depth.2tumour.copies,
                 Corrected.Depth.1tumour.copy,
                 mut.Cell.no.WT.Tumour.DNA.correction))
    
  } else {
    
    return(CN.adj.VAF)
    
  }
}


#' Function to correct data when we see evidence of subsequent CIN
#' 
#' 
#' 
#' 
mulitmodal.corrections <- function(archer.table){
  
  archer.table$diptestresult.q <- NA
  clones.to.test <- as.character(unique(archer.table[archer.table$muts.followed.in.clone>20 & 
                                                       archer.table$call=="RED" & !is.na(archer.table$call) &
                                                       archer.table$Clone.detected.q<0.1 &
                                                       archer.table$Clone.sum.DAOs.nobackground > 40 &
                                                       !is.na(archer.table$PyCloneCluster) &
                                                       !archer.table$Day_post_op=="0","sample.clone"]))
  clones.to.test <- clones.to.test[!is.na(clones.to.test)]
  
  diptests.to.correct <- sapply(clones.to.test, function(clone) unique(archer.table[archer.table$sample.clone == clone,"diptestresult.p"]))                                                   
  diptests.corrected <- p.adjust(diptests.to.correct, n = length(diptests.to.correct), method= "fdr")
  names(diptests.corrected) <- clones.to.test
  
  for(clone in clones.to.test){
    archer.table[archer.table$sample.clone == clone, "diptestresult.q"] <- diptests.corrected[names(diptests.corrected)==clone]
  }
  
  archer.table$is.multimodal <- archer.table$diptestresult.q < 0.1
  
  ##also look for whether we can see signal when we combine p-values accross muliple samples from the same patient?
  patients.to.test <- sapply(strsplit(clones.to.test,split=" "),"[[",1)
  barcodes.to.test <- patients.to.test
  #remove baselines
  baseline.index <- sapply(patients.to.test, function(sample) !all(as.numeric(archer.table[archer.table$Barcode==sample,"Day_post_op"]) <= 0)) 
  patients.to.test <- patients.to.test[baseline.index]
  patients.to.test <- sapply(patients.to.test, function(sample) unique(archer.table[archer.table$Barcode==sample,"SampleID"]))
  clones.to.test <- sapply(strsplit(clones.to.test[baseline.index],split=" "),"[[",2,USE.NAMES = F)
  no.dups.index <- !duplicated(paste(patients.to.test,clones.to.test))
  patients.to.test <- patients.to.test[no.dups.index]
  clones.to.test <- clones.to.test[no.dups.index]
  pat.clones <- paste(patients.to.test,clones.to.test)
  
  patient.clone.ps.combined <- sapply(pat.clones, function(pat.clone){
    pat <- strsplit(pat.clone,split=" ")[[1]][1]
    clone <- strsplit(pat.clone,split=" ")[[1]][2]
    p.values <- archer.table[archer.table$SampleID==pat & 
                               archer.table$PyCloneCluster==clone & 
                               archer.table$Barcode %in% barcodes.to.test,] 
    p.values <- p.values[!duplicated(p.values$Barcode),"diptestresult.p"] 
    p.values <- p.values[!is.na(p.values)]
    if(length(p.values)>1){
      return(sumlog(p.values)$p)
    } else {
      return(NULL)
    }
  })
  names(patient.clone.ps.combined) <- pat.clones
  patient.clone.ps.combined <- patient.clone.ps.combined[!(is.null(patient.clone.ps.combined) | sapply(patient.clone.ps.combined,length)==0)]
  
  archer.table$combined.diptest.p <- NA
  for(pat.clone in names(patient.clone.ps.combined)){
    pat <- strsplit(pat.clone,split=" ")[[1]][1]
    clone <- strsplit(pat.clone,split=" ")[[1]][2]
    archer.table[archer.table$SampleID==pat & 
                   archer.table$PyCloneCluster==clone &
                   archer.table$diptestresult.p < 0.05 &
                   !is.na(archer.table$diptestresult.p), "combined.diptest.p"] <- patient.clone.ps.combined[names(patient.clone.ps.combined)==pat.clone]
  }
  
  #correct
  patient.clone.qs.combined <- p.adjust(patient.clone.ps.combined,n=length(patient.clone.ps.combined),method = "fdr")
  
  archer.table$combined.diptest.q <- NA
  for(pat.clone in names(patient.clone.qs.combined)){
    pat <- strsplit(pat.clone,split=" ")[[1]][1]
    clone <- strsplit(pat.clone,split=" ")[[1]][2]
    archer.table[archer.table$SampleID==pat & 
                   archer.table$PyCloneCluster==clone &
                   archer.table$diptestresult.p < 0.05 &
                   !is.na(archer.table$diptestresult.p), "combined.diptest.q"] <- patient.clone.qs.combined[names(patient.clone.qs.combined)==pat.clone]
  }
  
  archer.table[archer.table$combined.diptest.q<0.1 & !is.na(archer.table$combined.diptest.q),"is.multimodal"] <- "TRUE"
  
  #save(archer.table,file=paste0(Output,"temp_archer.table.file.rda"))  
  #load(file=paste0(Output,"temp_archer.table.file.rda"))  #archer.table
  
  suppressPackageStartupMessages(library(mclust))
  
  ###for those positive samples, correct mutation CN and remove mutations which are no longer prersent from further analysis, indcate
  
  mulit.model.clones <- unique(archer.table[archer.table$is.multimodal == TRUE & !is.na(archer.table$is.multimodal),"sample.clone"])
  mulit.model.samples <- unique(archer.table[archer.table$is.multimodal == TRUE & !is.na(archer.table$is.multimodal),"Barcode"])
  
  archer.table$oldMutCN <- NA
  archer.table$newMutCN <- NA
  archer.table$CN.change <- NA
  archer.table$new.clone.CCF <- NA
  
  for(sample in mulit.model.samples){
    print(sample)
    
    sample.table <- archer.table[archer.table$Barcode==sample,]
    sample.clones <- mulit.model.clones[grepl(sample,mulit.model.clones)]
    
    for(sampleclone in sample.clones){
      
      sampleclone.table <- archer.table[archer.table$sample.clone==sampleclone,]
      
      VAFs <- sampleclone.table$VAF
      VAFs[VAFs==0] <- sampleclone.table[sampleclone.table$VAF==0,"LOD"]
      
      oldMutCN <- sampleclone.table$cellMutCPN
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
        DAOs <- sampleclone.table$DAOs[classification==class]
        Expected <- sampleclone.table$Expected_DAOS[classification==class]
        result <- poisson.test(sum(DAOs),sum(Expected),alternative = "greater")$p.value
        out <- result < 0.01
        return(out)
      })
      names(class.positive) <- classes
      
      
      newMutCN <- oldMutCN
      ##if cluster = background then CN of these mutations must now be 0
      newMutCN[classification==classes[!class.positive]] <- 0
      
      #now ensure that its not an incompleteCIN event where mutations that have been lost are 
      if(sum(!class.positive)==0){
        nmuts <- table(classification)
        class.1mutcpn.perc <- sapply(classes,function(class) sum(oldMutCN[classification==class]==1,na.rm = T)/sum(classification==class))
        #cluster of 1 copy should be one of the top 2 clusers with the most muts and the one of these with the most 1 cp original muts
        classcp1 <- classes[order(nmuts,decreasing=T) <= 2]
        classcp1 <- classcp1[class.1mutcpn.perc == max(class.1mutcpn.perc[order(nmuts,decreasing=T) <= 2])]
        #if they both have the same number of 1 muts (no WGD eg) should be detelions and cp1 hence 1cp is the higher VAF cluster
        if(length(classcp1)>1){
          classcp1 <- classcp1[class.VAF[classes %in% classcp1] %in% max(class.VAF[unique(mod1$classification) %in% classcp1])]
          newMutCN[classification==classes[classes==classcp1]] <- 1
          
        }
        
        classcp2 <- classcp1 + 1
        if(any(classes==classcp2)){
          newMutCN[classification==classes[classes==classcp2]] <- 2
          
          if(!classcp2 == max(classes)){
            classcpamp <- classcp2 + 1
            newMutCN[classification==classes[classes==classcpamp]] <- "amp"
          }
        }
        
        if(!classcp1 == min(classes)){
          classcp0 <- classcp1 - 1
          newMutCN[classification==classes[classes==classcp0]] <- 0
          
          #if we have mutations at 0 mutCN but with signal above background then we can estimate the CCF of the new clone (!!) & must have been CIN not pre-MRCA mets
          Copy1VAF.new <- as.numeric(class.VAF[names(class.VAF)==classcp1] - class.VAF[names(class.VAF)==classcp0])
          Copy1VAF.old <- 10^mean(log10(VAFs[classification==classcp0] / oldMutCN [classification==classcp0])) ##might need to remove background niose from 0
          new.clone.CCF <- 1 - (Copy1VAF.old / Copy1VAF.new)
          sampleclone.table$new.clone.CCF <- new.clone.CCF
          
        }
        
      } else {
        
        positive.classes <- classes[class.positive]
        positive.classes <- positive.classes[order(positive.classes)]
        positive.classes.newCN <- 1:length(positive.classes)
        positive.classes.newCN[positive.classes.newCN>2] <- "amp"
        newMutCN[classification %in% positive.classes] <- positive.classes.newCN[match(classification[classification %in% positive.classes],positive.classes)]
        
      }
      
      newMutCN. <- newMutCN
      newMutCN.[newMutCN.=="amp"] <- 3
      CN.change <- as.numeric(newMutCN.) - as.numeric(oldMutCN)
      
      sampleclone.table$oldMutCN <- oldMutCN
      sampleclone.table$newMutCN <- newMutCN
      sampleclone.table$CN.change <- CN.change
      
      sample.table[sample.table$sample.clone == sampleclone,] <- sampleclone.table
      
    }
    
    CCFs <-  sample.table[,c("mutation_id","combPyCCF")]
    #make all clonal CCFs == 1 (this shouldn't change through iterations)
    CCFs$combPyCCF[sample.table$PyCloneClonal %in% "C"] <- 1
    
    #for some reason we have slightly different combpyCCFs calcualted within the same clones.. TODO .. 
    
    iteration <- 1
    mean.Perc.change.next <- c()
    mean.Perc.change.record <-  c()
    
    sample.table.constant <- sample.table
    i <- 0
    
    while( (mean.Perc.change.next > 0.01 || is.null(mean.Perc.change.next)) & (any(!is.na(CCFs[,2]) == TRUE) || i == 0) ){
      i <- i + 1
      #cat(paste0(i," "))
      sample.table <- sample.table.constant
      
      CN.adjust.VAF <- CN.adjust( Varcount = sampleclone.table$DAOs.no.background, 
                                  Depth = sampleclone.table$Deep_Depth, 
                                  TotalCN = sampleclone.table$panTotalCN,
                                  cloneMutCPN = sampleclone.table$cellMutCPN,
                                  CCF = CCFs[match(sampleclone.table$mutation_id,CCFs$mutation_id),2],
                                  return.all.intermediate.calculations = TRUE )
      
      sampleclone.table[,names(sampleclone.table) %in% colnames(CN.adjust.VAF)] <- CN.adjust.VAF
      
      #determine tumour fraction first
      sampleclone.table.clonal <- sample.table[sample.table$PyCloneClonal == "C" & !is.na(sample.table$PyCloneClonal),]
      clonalpopmean <- mean.estimate.LOD.correction(CN.adj.VAF = sampleclone.table.clonal$CN.adj.VAF, MatchedLOD = sampleclone.table.clonal$Matched_LOD,normSD.est = normSD.for.1.sup.mut.clones, use.vars = sampleclone.table.clonal$Use.in.Phylogenetics)
      clonalpopmean <- as.numeric(as.character(unique(clonalpopmean$popmean)))
      if(is.na(clonalpopmean) || length(clonalpopmean)==0){
        clonalpopmean <- 0
      }
      
      for(sampleclone in sample.clones){
        
        modelling.outputs <- mean.estimate.LOD.correction(CN.adj.VAF = sampleclone.table$CN.adj.VAF, MatchedLOD = rep(TRUE,nrow(sampleclone.table)),
                                                          normSD.est = normSD.for.1.sup.mut.clones, 
                                                          use.vars = sampleclone.table$Use.in.Phylogenetics & !sampleclone.table$newMutCN == "amp", 
                                                          Remove.negative = TRUE)
        
        if(any(sampleclone.table$PyCloneClonal) %in% "C"){
          clonalpopmean <- as.numeric(as.character(unique(sampleclone.table$popmean)))
          if(is.na(clonalpopmean) || length(clonalpopmean)==0){
            clonalpopmean <- 0
            sample.table$Tumour.fraction <- NA
          } else {
            sample.table$Tumour.fraction <- clonalpopmean
          }
        }
        
        sampleclone.table[,names(sampleclone.table) %in% names(modelling.outputs)] <- modelling.outputs
        
        
        #calulate approximate ctDNA fraction for each clone and SD for ctDNA fraction
        sampleclone.table$mutcell.normSD <- sd(sampleclone.table$CN.adj.VAF)/10^(mean(log10(sampleclone.table$CN.adj.VAF),na.rm=T))
        sampleclone.table$mutcell.log.normSD <- sd(log10(sampleclone.table$CN.adj.VAF))/mean(log10(sampleclone.table$CN.adj.VAF),na.rm=T)
        sampleclone.table$mutcell.mean <- 10^(mean(log10(sampleclone.table$CN.adj.VAF),na.rm=T))
        sampleclone.table$Clone.mean.cellularity.equiv <-mean(sampleclone.table$Cellularity.eqivilent,na.rm=T)
        sampleclone.table$Clone.mean.Tumour.DNA.weight.perc <-mean(sampleclone.table$Tumour.DNA.weight.fraction,na.rm=T)
        
        #calculate power to detect clone at 10% CCF 1% CCF 50% CCF and then the CCF equivilent to doubel the background noise
        if(all(sampleclone.table$PyCloneClonal=="S" & !is.na(sampleclone.table$PyCloneClonal)) & (!is.na(clonalpopmean) && length(clonalpopmean)>0)){
          
          clonal.muts <- sample.table$PyCloneClonal=="C" & !is.na(sample.table$PyCloneClonal)
          
          CCF10.expected.subclone <- model.subclone.DAOs(clonalpopmean,CCF=0.1,sampleclone.table)
          CCF1.expected.subclone <- model.subclone.DAOs(clonalpopmean,0.01,sampleclone.table)
          CCF50.expected.subclone <- model.subclone.DAOs(clonalpopmean,0.5,sampleclone.table)
          
          background <- as.numeric(sampleclone.table$Expected_DAOS) 
          
          if(all(sample.table[clonal.muts,"DAOs"]==0)) no.signal <- TRUE else no.signal <- FALSE
          
          sampleclone.table$Power.to.detect.10CCF.p <- poisson.test(sum(CCF10.expected.subclone),sum(background),alternative = "greater")$p.value
          sampleclone.table$Power.to.detect.1CCF.p <- poisson.test(sum(CCF1.expected.subclone),sum(background), alternative = "greater")$p.value
          sampleclone.table$Power.to.detect.50CCF.p <- poisson.test(sum(CCF50.expected.subclone),sum(background), alternative = "greater")$p.value
          
          CCF100.expected <- model.subclone.DAOs(clonalpopmean,1,sampleclone.table)
          CCF_background_niose <- mean(background) / mean(CCF100.expected)
          if(CCF_background_niose>1 || is.na(CCF_background_niose)) CCF_background_niose <-  1
          sampleclone.table$CCF_background_niose <- CCF_background_niose
          
          #do we have power to detect the CCF in the primary tumour at baseline?
          MseqCCF.expected <- model.subclone.DAOs(clonalpopmean,unique(sampleclone.table$combPyCCF),sampleclone.table)
          sampleclone.table$Power.to.detect.MseqCCF.p <- poisson.test(sum(MseqCCF.expected), sum(background), alternative= "greater")$p.value
          
          
        } else {
          sampleclone.table$Power.to.detect.10CCF.p <- NA
          sampleclone.table$Power.to.detect.1CCF.p <- NA
          sampleclone.table$Power.to.detect.50CCF.p <- NA
          sampleclone.table$CCF_background_niose <- NA
          sampleclone.table$Power.to.detect.MseqCCF.p <- NA
          no.signal <- FALSE
        }
        
        if(no.signal){
          sampleclone.table$Power.to.detect.10CCF.p <- 1
          sampleclone.table$Power.to.detect.1CCF.p <- 1
          sampleclone.table$Power.to.detect.50CCF.p <- 1
          sampleclone.table$Power.to.detect.MseqCCF.p <- 1
          sampleclone.table$CCF_background_niose <- 1
        }
        
        sample.table[sample.table$sample.clone == sampleclone,] <- sampleclone.table
      }
      
      
      #if no clonal vars set to NA (ie if been unable to run pyclone)
      if(any(sample.table$PyCloneClonal %in% "C")){
        clonalpopmean <- as.numeric(as.character(unique(sample.table[sample.table$PyCloneClonal %in% "C","popmean"])))
      } else {
        clonalpopmean <- NA
      }
      
      sample.table$plasma_CCF <- as.numeric(as.character(sample.table$popmean)) /  clonalpopmean
      sample.table$plasma_CCF_UpCI <- as.numeric(as.character(sample.table$UpCI)) /  clonalpopmean
      sample.table$plasma_CCF_LwCI <- as.numeric(as.character(sample.table$LwCI)) /  clonalpopmean
      
      ##quantify how much popmean changes from last iteration with popmean table
      if(iteration == 1){
        meanpoptable <- sample.table[!duplicated(sample.table$PyCloneCluster), c("PyCloneCluster","PyCloneClonal")]
        meanpoptable <- cbind(meanpoptable, signif(as.numeric(as.character(sample.table[!duplicated(sample.table$PyCloneCluster),"popmean"])),5))
        names(meanpoptable)[ncol(meanpoptable)] <- paste0("popmean_iteration_",iteration)
        
      } else {
        meanpoptable <- cbind(meanpoptable, signif(as.numeric(as.character(sample.table[!duplicated(sample.table$PyCloneCluster),"popmean"])),5))
        names(meanpoptable)[ncol(meanpoptable)] <- paste0("popmean_iteration_",iteration)
        mean.Perc.change.next <- mean(abs( as.numeric(meanpoptable[,(ncol(meanpoptable) - 1)]) - as.numeric(meanpoptable[,ncol(meanpoptable)]) ) / as.numeric(meanpoptable[,ncol(meanpoptable) - 1]),na.rm=T)
        mean.Perc.change.record <- c(mean.Perc.change.record,mean.Perc.change.next)  #to record changes over iteration when testing
      }
      
      iteration <- iteration + 1
      
      CCFs <- sample.table[,c("mutation_id","plasma_CCF")]
      
      if(iteration==100){
        cat(paste(sample,"could not find CCF solution for all subclones "))
        badclones <- meanpoptable$PyCloneCluster[abs( as.numeric(meanpoptable[,(ncol(meanpoptable) - 1)]) - as.numeric(meanpoptable[,ncol(meanpoptable)]) ) / as.numeric(meanpoptable[,ncol(meanpoptable) - 1])>0.01 & !is.na(meanpoptable[,ncol(meanpoptable)])]
        CCFs[sample.table$PyCloneCluster %in% badclones,"plasma_CCF"] <- NA
      }
      
    }
    
    archer.table[archer.table$Barcode==sample,] <- sample.table
    
  }
  
  return(archer.table)
  
}

#' Clonal desconvoltuion function
#' 
#' 
#' 
#' 
Calculate.CCFs <- function(Archermuts = Archermuts, normSD.for.1.sup.mut.clones = normSD.for.1.sup.mut.clones, LODcorrection = TRUE, test.subsequent.CIN.evo = TRUE){
  
  # make a field for wherther a mutation can be used in the phylogenetic analysis
  Archermuts$Use.in.Phylogenetics <- TRUE
  
  ### mutations cannot be used if they have been hardfiltered by Archer (eg if they have heavy strand bais) - NAs for when no support
  Archermuts[!Archermuts$Target_Variant_Hard_Filtered %in% c("False",NA), "Use.in.Phylogenetics"] <- FALSE
  
  ### mutations cannot be use if they don't have CN and mutation info overlaid
  Archermuts[ is.na(Archermuts$cellMutCPN) | is.na(Archermuts$combPyCCF) | is.na(Archermuts$panTotalCN) | Archermuts$panTotalCN %in% NaN, "Use.in.Phylogenetics"] <- FALSE
  
  #calculate the List of detection for each variant based onCN.adj.VAF that = 1 supporting read
  Archermuts$LOD <- CN.adjust( Varcount=1 , Depth=Archermuts$Deep_Depth, TotalCN=Archermuts$panTotalCN, cloneMutCPN=Archermuts$cellMutCPN, CCF=Archermuts$combPyCCF)
  Archermuts$Support <- Archermuts$DAOs > 0
  Archermuts$Matched_LOD <- TRUE
  Archermuts$DAOs.no.background <- Archermuts$DAOs
  
  #Calculate normal VAF to compare
  Archermuts$VAF.no.background <- Archermuts$DAOs.no.background / Archermuts$Deep_Depth
  Archermuts$VAF <- Archermuts$DAOs / Archermuts$Deep_Depth
  
  Archermutslist <- list()
  
  for(sampleclone in unique(paste(Archermuts$Barcode, Archermuts$PyCloneCluster))){
    
    cat(paste0(which(unique(paste(Archermuts$Barcode, Archermuts$PyCloneCluster))==sampleclone))," ")
    
    sampleclone.table <- Archermuts[paste(Archermuts$Barcode,Archermuts$PyCloneCluster) %in% sampleclone,]
    
    if(grepl("NA",sampleclone)){
      Archermutslist[[which(unique(paste(Archermuts$Barcode, Archermuts$PyCloneCluster)) %in% sampleclone)]] <- sampleclone.table
      next
    }
    
    ##identify mutations to remove to account for background noise when calculating CCFs
    ##old method - random hence not consistant or reproducable
    #reads.to.remove <- rpois(nrow(sampleclone.table), lambda = mean(sampleclone.table$Expected_DAOS))
    #reads.to.remove <- table(reads.to.remove[!reads.to.remove==0])
    ##new method
    ###do this for each D group
    for(Dgroup in unique(sampleclone.table[sampleclone.table$Use.in.Phylogenetics==TRUE,"Dgroup"])){
      cat(Dgroup)
      background.lambda <- mean(sampleclone.table[sampleclone.table$Dgroup %in% Dgroup & sampleclone.table$Use.in.Phylogenetics==TRUE,"Expected_DAOS"])
      max.read.support.to.test <- round(background.lambda*100)
      reads.to.remove <- round(dpois(1:max.read.support.to.test, lambda = background.lambda)*sum(sampleclone.table$Dgroup %in% Dgroup))
      names(reads.to.remove) <- 1:max.read.support.to.test
      reads.to.remove <- factor(reads.to.remove[!reads.to.remove==0], levels = 1:sum(sampleclone.table$Dgroup %in% Dgroup))
      
      if(length(reads.to.remove)>0){
        #only needed for old method
        #needs to speciify 0s (eg 0 cases of 1 supporting read to remove) so it can be added to if required ## needed for old method
        #if(length(reads.to.remove)>0){
        #  if(any(!seq(1:max(as.numeric(names(reads.to.remove)))) == as.numeric(names(reads.to.remove)))){
        #    toremovefull <- table(seq(1:max(as.numeric(names(reads.to.remove)))))
        #    toremovefull[seq(1:max(as.numeric(names(reads.to.remove))))] <- 0
        #    toremovefull[as.numeric(names(reads.to.remove))] <- reads.to.remove
        #    reads.to.remove <- toremovefull
        #  }
        #}
        
        mutidsnotcorrected <- sampleclone.table[!sampleclone.table$DAOs == 0 & sampleclone.table$Use.in.Phylogenetics==TRUE,"mutation_id"]
        
        for(background in rev(names(reads.to.remove))){
          eliablemuts <- sampleclone.table[sampleclone.table$Dgroup %in% Dgroup &
                                             sampleclone.table$DAOs >= as.numeric(background) &
                                             sampleclone.table$mutation_id %in% mutidsnotcorrected &
                                             sampleclone.table$Use.in.Phylogenetics==TRUE,"mutation_id"]
          
          if(length(eliablemuts) == 0){
            if(background>1){
              reads.to.add <- as.numeric(background) * as.numeric(reads.to.remove[names(reads.to.remove)==background])
              leveldown <- names(reads.to.remove)[which(names(reads.to.remove)== background)-1]
              reads.to.remove[leveldown] <- as.numeric(reads.to.remove[leveldown]) + floor(reads.to.add/as.numeric(names(reads.to.remove[leveldown])))
              reads.to.carry.over <- reads.to.add %% as.numeric(names(reads.to.remove[leveldown]))
              for(i in rev(names(reads.to.remove)[1:(length(reads.to.remove))])){
                if(reads.to.carry.over == 0){
                  break
                }
                reads.to.remove[i] <- reads.to.remove[i] + floor(reads.to.carry.over/as.numeric(names(reads.to.remove[i])))
                reads.to.carry.over <- reads.to.carry.over %% as.numeric(names(reads.to.remove[i]))
              }
            }
            next
          }
          if(length(eliablemuts) >= as.numeric(reads.to.remove[names(reads.to.remove) %in% background])){
            chosen <- eliablemuts[runif(reads.to.remove[names(reads.to.remove) %in% background], 1, length(eliablemuts))]
            mutidsnotcorrected <- mutidsnotcorrected[!mutidsnotcorrected %in% chosen]
            sampleclone.table[sampleclone.table$mutation_id %in% chosen,"DAOs.no.background"] <- sampleclone.table[sampleclone.table$mutation_id %in% chosen,"DAOs"] - as.numeric(background)
            
          } else {
            
            mutidsnotcorrected <- mutidsnotcorrected[!mutidsnotcorrected %in% eliablemuts]
            sampleclone.table[sampleclone.table$mutation_id %in% eliablemuts,"DAOs.no.background"] <- sampleclone.table[sampleclone.table$mutation_id %in% eliablemuts,"DAOs"] - as.numeric(background)
            
            if(background>1){
              reads.to.add <- as.numeric(background) * (as.numeric(reads.to.remove[names(reads.to.remove)==background]) - length(eliablemuts))
              leveldown <- names(reads.to.remove)[which(names(reads.to.remove)== background)-1]
              reads.to.remove[leveldown] <- reads.to.remove[leveldown] + floor(reads.to.add/as.numeric(names(reads.to.remove[leveldown])))
              reads.to.carry.over <- reads.to.add %% as.numeric(names(reads.to.remove[leveldown]))
              for(i in rev(names(reads.to.remove)[1:(length(reads.to.remove))])){
                if(reads.to.carry.over == 0){
                  break
                }
                reads.to.remove[i] <- reads.to.remove[i] + floor(reads.to.carry.over/as.numeric(names(reads.to.remove[i])))
                reads.to.carry.over <- reads.to.carry.over %% as.numeric(names(reads.to.remove[i]))
              }
            }
          }
        }
      }
    }
    
    if(LODcorrection == TRUE){
      #Highlight unsupported muts which can be removed in further analysis so they are not biased by LOD
      if(sum(!sampleclone.table$Support & !is.na(sampleclone.table$LOD) & sampleclone.table$Use.in.Phylogenetics==TRUE) > 0 & sum(sampleclone.table$Support & !is.na(sampleclone.table$LOD) & sampleclone.table$Use.in.Phylogenetics==TRUE) > 0){
        medianSuportedLOD <- median(sampleclone.table[sampleclone.table$Support & sampleclone.table$Use.in.Phylogenetics==TRUE,"LOD"],na.rm = T)
        iteration <- 1
        repeat{
          if(!iteration == 1){
            oldLODsdiff <- LODsdiff
          }
          
          sampleclone.table.unsupported <- sampleclone.table[sampleclone.table$Matched_LOD %in% 'TRUE' & !is.na(sampleclone.table$LOD) & !sampleclone.table$Support & sampleclone.table$Use.in.Phylogenetics==TRUE,]
          maxLODunsupportedmut <- sampleclone.table.unsupported[sampleclone.table.unsupported$LOD == max(sampleclone.table.unsupported$LOD[sampleclone.table.unsupported$Use.in.Phylogenetics==TRUE],na.rm = T),"mutation_id"]
          
          sampleclone.table[sampleclone.table$mutation_id %in% maxLODunsupportedmut,"Matched_LOD"] <- FALSE
          
          medianUnsuportedLOD <- median(sampleclone.table[!sampleclone.table$Support & sampleclone.table$Matched_LOD %in% 'TRUE' & !is.na(sampleclone.table$LOD) & sampleclone.table$Use.in.Phylogenetics==TRUE,"LOD"],na.rm = T)
          
          LODsdiff <- medianUnsuportedLOD - medianSuportedLOD
          
          iteration <- iteration + 1
          
          #when you've removed enough unsupported muts so LODunsupported is less than LODsupported then check whether they would be more similar
          #had you kept the previous mut and if so restore it
          #also if no more unsupported muts left
          if(!any(!sampleclone.table$Support & sampleclone.table$Matched_LOD & !is.na(sampleclone.table$LOD) & sampleclone.table$Use.in.Phylogenetics==TRUE) | LODsdiff<0){
            break
            if(abs(LODsdiff)<abs(oldLODsdiff) || LODsdiff %in% NA){
              break
            } else {
              #sampleclone.table[sampleclone.table$mutation_id %in% maxLODunsupportedmut,"Matched_LOD"] <- TRUE
              break
            }
          }
        }
      }
    }
    Archermutslist[[which(unique(paste(Archermuts$Barcode, Archermuts$PyCloneCluster)) %in% sampleclone)]] <- sampleclone.table
  }
  Archermuts <- do.call(rbind, Archermutslist)
  
  #save(Archermuts,file = paste0(Input,date,"tmp_Archermuts_additional_anotations.RData"))
  #load(file = paste0(Input,"/20200326tmp_Archermuts_additional_anotations.RData"))#Archermuts
  
  ##### Calculate ctDNA cellularity accounting for mutCPN, WT cancer copies #####
  
  #loop over each sample improving CCF accuracy with each iteration
  Archermutslist <- list()
  
  for(sample in unique(Archermuts$Barcode)){
    cat(paste0(sample," "))
    
    sample.table <- Archermuts[Archermuts$Barcode == sample,]
    #sample.table <- sample.table[,1:96]
    #need approximate CCFs to calculate WT tumour DNA contrebution (will iterate to correct CCFs)
    #use meanVAF for this initially
    CCFs <- sample.table[,c("mutation_id","Barcode","PyCloneCluster")]
    clonal.clust <- unique(sample.table[sample.table$PyCloneClonal %in% "C","PyCloneCluster"])
    CCFs$CCF <- sapply(paste(CCFs$Barcode,CCFs$PyCloneCluster), function(sample.clone){
      sample <- strsplit(sample.clone,split=" ")[[1]][1]
      CCF <- mean(sample.table[paste(sample.table$Barcode,sample.table$PyCloneCluster)==sample.clone,"DAOs"],na.rm = T) /
        mean(sample.table[sample.table$Barcode==sample & sample.table$PyCloneCluster==clonal.clust,"DAOs"],na.rm = T)
      return(CCF)
    })
    #make all clonal CCFs == 1 (this shouldn't change through iterations)
    CCFs[CCFs$PyCloneCluster==clonal.clust & !is.na(CCFs$PyCloneCluster),"CCF"] <- 1
    #don't allow CCF > 1
    CCFs[CCFs$CCF>1 &!is.na(CCFs$CCF),"CCF"] <- 1
    #for this initial CCF don't allow <1% (causes WT tumour DNA estimate to be overetimated)
    CCFs[CCFs$CCF<0.01 &!is.na(CCFs$CCF),"CCF"] <- 0.01
    
    
    #for some reason we have slightly different combpyCCFs calcualted within the same clones.. TODO .. 
    iteration <- 1
    mean.Perc.change.next <- c()
    mean.Perc.change.record <-  c()
    
    sample.table.constant <- sample.table
    i <- 0
    
    while( (mean.Perc.change.next > 0.01 || is.null(mean.Perc.change.next)) & (any(!is.na(CCFs[,"CCF"]) == TRUE) || i == 0) ){
      i <- i + 1
      cat(paste0("iteration: ",i," "))
      sample.table <- sample.table.constant
      
      CN.adj.VAF.outputs <- CN.adjust( Varcount = sample.table$DAOs.no.background, 
                                       Depth = sample.table$Deep_Depth, 
                                       TotalCN = sample.table$panTotalCN,
                                       cloneMutCPN = sample.table$cellMutCPN,
                                       #is.Clonal = sample.table$PyCloneClonal == "C",
                                       CCF = CCFs[match(paste(CCFs$mutation_id,CCFs$Barcode),paste(sample.table$mutation_id,sample.table$Barcode)),"CCF"],
                                       return.all.intermediate.calculations = TRUE )
      
      sample.table <- cbind(sample.table,CN.adj.VAF.outputs)
      
      sample.table.list <- list()
      
      #estimate averageCN.adj.VAF for each clone using tail method
      
      #determine tumour fraction first
      sampleclone.table.clonal <- sample.table[sample.table$PyCloneClonal == "C" & !is.na(sample.table$PyCloneClonal),]
      clonalpopmean <- mean.estimate.LOD.correction(CN.adj.VAF = sampleclone.table.clonal$CN.adj.VAF, MatchedLOD = sampleclone.table.clonal$Matched_LOD,normSD.est = normSD.for.1.sup.mut.clones, use.vars = sampleclone.table.clonal$Use.in.Phylogenetics)
      clonalpopmean <- as.numeric(as.character(unique(clonalpopmean$popmean)))
      if(is.na(clonalpopmean) || length(clonalpopmean)==0){
        clonalpopmean <- 0
        sample.table$Tumour.fraction <- NA
      } else {
        sample.table$Tumour.fraction <- clonalpopmean
      }
      
      #extract tree so you can test whether daughter >/< parent
      sample.tree <- extract.tree.mutTable(mutTable=sample.table, clone.field="PyCloneCluster", parent.field="Tree.parent")
      
      for(clone in unique(sample.table$PyCloneCluster)){
        cat(paste0(clone," "))
        
        sampleclone.table <- sample.table[sample.table$PyCloneCluster %in% clone,]
        
        sampleclone.table <- cbind(sampleclone.table, mean.estimate.LOD.correction(CN.adj.VAF = sampleclone.table$CN.adj.VAF, MatchedLOD = sampleclone.table$Matched_LOD,normSD.est = normSD.for.1.sup.mut.clones, use.vars = sampleclone.table$Use.in.Phylogenetics))
        
        #has there been further evolution and CIN? test whether distbution is mulitmodal (use this later to correct p values)
        library(diptest)
        if(length(sampleclone.table$VAF[sampleclone.table$VAF>0 & sampleclone.table$Use.in.Phylogenetics==TRUE])>0){
          VAFs <- sampleclone.table$VAF
          #make 0 the LOD
          VAFs[VAFs==0] <- sampleclone.table[sampleclone.table$VAF==0,"LOD"]
          dip.result.CN1 <- dip.test(log10(VAFs[ sampleclone.table$Use.in.Phylogenetics==TRUE & sampleclone.table$cellMutCPN<1.5]))$p.value
          if(dip.result.CN1==0) dip.result.CN1 <- 0.000000000000001 #stated min P-value by function
          if(sum(sampleclone.table$Use.in.Phylogenetics==TRUE & sampleclone.table$cellMutCPN>1.5)>1){
            dip.result.CN2 <- dip.test(log10(VAFs[ sampleclone.table$Use.in.Phylogenetics==TRUE & sampleclone.table$cellMutCPN>1.5]))$p.value
            if(dip.result.CN2==0) dip.result.CN2 <- 0.000000000000001 #stated min P-value by function
            combined.p <- sumlog(c(dip.result.CN1,dip.result.CN2))$p
          } else {
            combined.p <- dip.result.CN1
          }
          sampleclone.table$diptestresult.p <- combined.p
          
        } else {
          sampleclone.table$diptestresult.p <- NA
        }
        
        #calulate approximate ctDNA fraction for each clone and SD for ctDNA fraction
        sampleclone.table$clone.fraction.median <- median(sampleclone.table$CN.adj.VAF[!sampleclone.table$CN.adj.VAF %in% Inf],na.rm=T)
        sampleclone.table$clone.fraction.perc.no.support <- sum(sampleclone.table$DAOs.no.background == 0) / nrow(sampleclone.table)
        sampleclone.table$clone.fraction.perc.at.least.two.support <- sum(sampleclone.table$DAOs.no.background >= 2) / nrow(sampleclone.table)
        sampleclone.table$muts.followed.in.clone <- nrow(sampleclone.table)
        sampleclone.table$unfiltered.muts.followed.in.clone <- sum(sampleclone.table$Target_Variant_Deep_Error_Filtered == "false")
        sampleclone.table$no.supported.muts <- sum(sampleclone.table$DAOs.no.background>0)
        sampleclone.table$mutcell.normSD <- sd(sampleclone.table$CN.adj.VAF)/10^(mean(log10(sampleclone.table$CN.adj.VAF),na.rm=T))
        sampleclone.table$mutcell.log.normSD <- sd(log10(sampleclone.table$CN.adj.VAF))/mean(log10(sampleclone.table$CN.adj.VAF),na.rm=T)
        sampleclone.table$mutcell.mean <- 10^(mean(log10(sampleclone.table$CN.adj.VAF),na.rm=T))
        sampleclone.table$Clone.VAF.nobackground <- mean(sampleclone.table$VAF.no.background)
        sampleclone.table$Clone.VAF <- mean(sampleclone.table$VAF)
        sampleclone.table$Clone.sum.DAOs.nobackground <- sum(sampleclone.table$DAOs.no.background)
        sampleclone.table$Clone.sum.DAOs <-sum(sampleclone.table$DAOs)
        sampleclone.table$Clone.mean.cellularity.equiv <-mean(sampleclone.table$Cellularity.eqivilent,na.rm=T)
        sampleclone.table$Clone.mean.Tumour.DNA.weight.perc <-mean(sampleclone.table$Tumour.DNA.weight.fraction,na.rm=T)
        
        MRD.vars <- sampleclone.table$Target_Variant_Hard_Filtered == "False" & sampleclone.table$Target_Variant_Deep_Error_Filtered == "False" | (is.na(sampleclone.table$Target_Variant_Hard_Filtered) & is.na(sampleclone.table$Target_Variant_Deep_Error_Filtered))
        sampleclone.table$Clone.detected.p <-  poisson.test(c(sum(sampleclone.table$DAOs[MRD.vars]),round(sum(sampleclone.table$Expected_DAOS[MRD.vars]))), alternative =  "greater")$p.value
        sampleclone.table$Clone.meanLOD <- mean(sampleclone.table$LOD,na.rm=T)
        
        #test whether this clone is signifcantly different to its parent
        sampleclone.table$different.to.parent.p <- NA
        sampleclone.table$different.to.parent.OR <- NA
        if(!all(sampleclone.table$Tree.parent %in% c(NA,"root"))){
          parent.CNadjVAF <- as.numeric(sample.table[sample.table$PyCloneCluster == unique(sampleclone.table$Tree.parent),"CN.adj.VAF"])
          parent.CNadjVAF <- parent.CNadjVAF[!is.na(parent.CNadjVAF)]
          clone.CNadjVAF <- as.numeric(sampleclone.table$CN.adj.VAF)
          clone.CNadjVAF <- clone.CNadjVAF[!is.na(clone.CNadjVAF)]
          
          if(length(parent.CNadjVAF) > 1 & length(clone.CNadjVAF) > 1){
            sampleclone.table$different.to.parent.p <- wilcox.test(parent.CNadjVAF,clone.CNadjVAF)$p.value
            sampleclone.table$different.to.parent.OR <- mean(clone.CNadjVAF) / mean(parent.CNadjVAF) 
          } 
          
        }
        
        if(all(sampleclone.table$PyCloneClonal=="S" & !is.na(sampleclone.table$PyCloneClonal)) & (!is.na(clonalpopmean) && length(clonalpopmean)>0)){
          
          clonal.muts <- sample.table$PyCloneClonal=="C" & !is.na(sample.table$PyCloneClonal)
          
          CCF10.expected.subclone <- model.subclone.DAOs(clonalpopmean,CCF=0.1,sampleclone.table)
          CCF1.expected.subclone <- model.subclone.DAOs(clonalpopmean,0.01,sampleclone.table)
          CCF50.expected.subclone <- model.subclone.DAOs(clonalpopmean,0.5,sampleclone.table)
          
          background <- as.numeric(sampleclone.table$Expected_DAOS)
          
          if(all(sample.table[clonal.muts,"DAOs"]==0)) no.signal <- TRUE else no.signal <- FALSE
          
          sampleclone.table$Power.to.detect.10CCF.p <- poisson.test(sum(CCF10.expected.subclone),sum(background),alternative = "greater")$p.value
          sampleclone.table$Power.to.detect.1CCF.p <- poisson.test(sum(CCF1.expected.subclone), sum(background), alternative = "greater")$p.value
          sampleclone.table$Power.to.detect.50CCF.p <- poisson.test(sum(CCF50.expected.subclone), sum(background), alternative = "greater")$p.value
          
          CCF100.expected <- model.subclone.DAOs(clonalpopmean,1,sampleclone.table)
          CCF_background_niose <- mean(background) / mean(CCF100.expected)
          if(CCF_background_niose>1 || is.na(CCF_background_niose)) CCF_background_niose <-  1
          sampleclone.table$CCF_background_niose <- CCF_background_niose
          
          #do we have power to detect the CCF in the primary tumour at baseline?
          MseqCCF.expected <- model.subclone.DAOs(clonalpopmean,unique(sampleclone.table$combPyCCF),sampleclone.table)
          sampleclone.table$Power.to.detect.MseqCCF.p <- poisson.test(sum(MseqCCF.expected), sum(background), alternative= "greater")$p.value
          
        } else {
          sampleclone.table$Power.to.detect.10CCF.p <- NA
          sampleclone.table$Power.to.detect.1CCF.p <- NA
          sampleclone.table$Power.to.detect.50CCF.p <- NA
          sampleclone.table$CCF_background_niose <- NA
          sampleclone.table$Power.to.detect.MseqCCF.p <- NA
          no.signal <- FALSE
        }
        
        if(no.signal){
          sampleclone.table$Power.to.detect.10CCF.p <- 1
          sampleclone.table$Power.to.detect.1CCF.p <- 1
          sampleclone.table$Power.to.detect.50CCF.p <- 1
          sampleclone.table$Power.to.detect.MseqCCF.p <- 1
          sampleclone.table$CCF_background_niose <- 1
        }
        
        sample.table.list[[which(unique(sample.table$PyCloneCluster) %in% clone)]] <- sampleclone.table
      }
      
      sample.table <- do.call(rbind, sample.table.list)
      
      MRD.vars <- (sample.table$Target_Variant_Hard_Filtered == "False" & sample.table$Target_Variant_Deep_Error_Filtered == "False") | (is.na(sample.table$Target_Variant_Hard_Filtered) & is.na(sample.table$Target_Variant_Deep_Error_Filtered))
      
      sample.table$MRD.call.sample.p <- poisson.test(sum(sample.table$DAOs[MRD.vars]),sum(sample.table$Expected_DAOS[MRD.vars]), alternative = "greater")$p.value
      
      #if no clonal vars set to NA (ie if been unable to run pyclone)
      if(any(sample.table$PyCloneClonal %in% "C")){
        clonalpopmean <- as.numeric(as.character(unique(sample.table[sample.table$PyCloneClonal %in% "C","popmean"])))
      } else {
        clonalpopmean <- NA
      }
      
      sample.table$plasma_CCF <- as.numeric(as.character(sample.table$popmean)) /  clonalpopmean
      sample.table$plasma_CCF_UpCI <- as.numeric(as.character(sample.table$UpCI)) /  clonalpopmean
      sample.table$plasma_CCF_LwCI <- as.numeric(as.character(sample.table$LwCI)) /  clonalpopmean
      
      ##quantify how much popmean changes from last iteration with popmean table
      if(iteration == 1){
        meanpoptable <- sample.table[!duplicated(sample.table$PyCloneCluster), c("PyCloneCluster","PyCloneClonal")]
        meanpoptable <- cbind(meanpoptable, signif(as.numeric(as.character(sample.table[!duplicated(sample.table$PyCloneCluster),"popmean"])),5))
        names(meanpoptable)[ncol(meanpoptable)] <- paste0("popmean_iteration_",iteration)
        
      } else {
        meanpoptable <- cbind(meanpoptable, signif(as.numeric(as.character(sample.table[!duplicated(sample.table$PyCloneCluster),"popmean"])),5))
        names(meanpoptable)[ncol(meanpoptable)] <- paste0("popmean_iteration_",iteration)
        mean.Perc.change.next <- mean(abs( as.numeric(meanpoptable[,(ncol(meanpoptable) - 1)]) - as.numeric(meanpoptable[,ncol(meanpoptable)]) ) / as.numeric(meanpoptable[,ncol(meanpoptable) - 1]),na.rm=T)
        mean.Perc.change.record <- c(mean.Perc.change.record,mean.Perc.change.next)  #to record changes over iteration when testing
      }
      
      iteration <- iteration + 1
      
      CCFs <- sample.table[,c("mutation_id","Barcode","PyCloneCluster","plasma_CCF")]
      
      if(iteration==100){
        cat(paste(sample,"could not find CCF solution for all subclones "))
        badclones <- meanpoptable$PyCloneCluster[abs( as.numeric(meanpoptable[,(ncol(meanpoptable) - 1)]) - as.numeric(meanpoptable[,ncol(meanpoptable)]) ) / as.numeric(meanpoptable[,ncol(meanpoptable) - 1])>0.01 & !is.na(meanpoptable[,ncol(meanpoptable)])]
        CCFs[sample.table$PyCloneCluster %in% badclones,"plasma_CCF"] <- NA
      }
      
      names(CCFs)[4] <- "CCF"
    }
    
    Archermutslist[[which(unique(Archermuts$Barcode) %in% sample)]] <- sample.table
    
  }
  Archermuts <- do.call(rbind, Archermutslist)
  
  save(Archermuts,file = paste0(Input,date,"tmp_Archermuts_additional.anotations.RData"))
  #load(file = paste0(Input,"/20200327tmp_Archermuts_additional.anotations.RData"))#Archermuts
  
  Archermuts$sample.clone <- paste(Archermuts$Barcode,Archermuts$PyCloneCluster)
  
  #For samples without +MRD set the CCFs and popmeans etc to 0
  Archermuts[!Archermuts$call == "RED" & !is.na(Archermuts$call),c("popmean","UpCI","LwCI","plasma_CCF","plasma_CCF_UpCI","plasma_CCF_LwCI")] <- 0
  
  ##correct all p values
  Archermuts <- correct.specific.ps(Archermuts,field.to.deduplicate="sample.clone",field.to.correct="Clone.detected.p", indices.to.correct= Archermuts$call == "RED" & !is.na(Archermuts$call),method ="BH")
  Archermuts <- correct.specific.ps(Archermuts,field.to.deduplicate="sample.clone",field.to.correct="Power.to.detect.10CCF.p", indices.to.correct= Archermuts$call == "RED" & !is.na(Archermuts$call) & Archermuts$Tumour.fraction > 0.0001 & !is.na(Archermuts$Tumour.fraction),method ="BH")
  Archermuts <- correct.specific.ps(Archermuts,field.to.deduplicate="sample.clone",field.to.correct="Power.to.detect.1CCF.p", indices.to.correct= Archermuts$call == "RED" & !is.na(Archermuts$call) & Archermuts$Tumour.fraction > 0.0001 & !is.na(Archermuts$Tumour.fraction),method ="BH")
  Archermuts <- correct.specific.ps(Archermuts,field.to.deduplicate="sample.clone",field.to.correct="Power.to.detect.50CCF.p", indices.to.correct= Archermuts$call == "RED" & !is.na(Archermuts$call) & Archermuts$Tumour.fraction > 0.0001 & !is.na(Archermuts$Tumour.fraction),method ="BH")
  Archermuts <- correct.specific.ps(Archermuts,field.to.deduplicate="sample.clone",field.to.correct="Power.to.detect.MseqCCF.p", indices.to.correct= Archermuts$call == "RED" & !is.na(Archermuts$call) & Archermuts$Tumour.fraction > 0.0001 & !is.na(Archermuts$Tumour.fraction),method ="BH")
  Archermuts <- correct.specific.ps(Archermuts,field.to.deduplicate="sample.clone",field.to.correct="different.to.parent.p", indices.to.correct= Archermuts$call == "RED" & !is.na(Archermuts$call) & Archermuts$Clone.detected.q < 0.1 & !is.na(Archermuts$Clone.detected.q),method ="BH")
  Archermuts <- correct.specific.ps(Archermuts,field.to.deduplicate="SampleID",field.to.correct="MRD.call.sample.p",method ="BH")
  
  #Look at cases with possible CIN evoluation (signifcant dip.tests)- confirm with genomic localisation of losses (TO DO - tried & this is v difficult). (need to do outside loop to generate q-values)
  #then assign hom and het deletions and possible amplifications, recalculate popmean and CCFs based on muts incldue the correct CNAs
  if(test.subsequent.CIN.evo == TRUE){
    Archermuts <- mulitmodal.corrections(Archermuts)
  }
  
  Clone.detected.ps <- Archermuts[!duplicated(Archermuts$sample.clone) & Archermuts$call == "RED" & !is.na(Archermuts$call),"Clone.detected.p"]
  sampleclones <- Archermuts[!duplicated(Archermuts$sample.clone) & Archermuts$call == "RED" & !is.na(Archermuts$call),"sample.clone"]
  Clone.detected.qs <- p.adjust(Clone.detected.ps, n = length(Clone.detected.ps), method = "fdr")
  Archermuts$Clone.detected.q <- 1
  Archermuts$Clone.detected.q[Archermuts$call == "RED" & !is.na(Archermuts$call)] <- round(unlist(lapply(1:length(Clone.detected.qs), function(i) rep(Clone.detected.qs[i],sum(Archermuts$sample.clone %in% sampleclones[i])))),digits=4)
  
  
  #determine on sample level the clonality type (monoclonal,polyclonal, polyphyletic etc) -> use this later to detemine metastasis type on patient level
  Archermuts$Sample.clonality.type <- unlist(lapply(unique(Archermuts$Barcode), function(sample){
    print(sample)
    sample.index <- Archermuts$Barcode == sample
    
    clones.detected.q <- Archermuts[sample.index, "Clone.detected.q"]
    if(any(clones.detected.q<0.1 & !is.na(clones.detected.q))){
      sample.table <- Archermuts[sample.index,]
      
      detected <- sample.table[sample.table$Clone.detected.q<0.1,]
      clones.detected <- unique(detected$PyCloneCluster)[ !is.na(unique(detected$PyCloneCluster)) ]
      
      #if all CN fail can have no PyClone
      if(!length(clones.detected)==0){
        
        sample.tree <- extract.tree.mutTable(mutTable=sample.table, clone.field = "PyCloneCluster",parent.field = "Tree.parent")
        detected.tree <- remove.clones.on.tree(tree = sample.tree, clones.to.keep = clones.detected)
        #does the tree branch?
        branching <- any(duplicated(detected.tree[,1]))
        
        parents.of.detected <- sample.table[sample.table$PyCloneCluster %in% clones.detected,"Tree.parent"]
        parents.of.detected <- parents.of.detected[!is.na(parents.of.detected)]
        
        if(all(parents.of.detected=="root")){
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
    
    out <- rep(out,sum(sample.index))
    return(out)
  }))
  
  #make heirchary of options
  hierarchy <- data.frame(hierarchy=1:4, option=c("Clonal.only","Monoclonal","Polyclonal","Polyphyletic"),stringsAsFactors = F)
  
  ##now on a patient level for metastases
  Archermuts$patient.metastasis.type <- unlist(lapply(unique(Archermuts$SampleID), function(pat){
    
    pat.index <-Archermuts$SampleID==pat
    pat.clonal.index <- Archermuts$PyCloneClonal == "C" & Archermuts$SampleID==pat
    
    clones.detected.q <- Archermuts[pat.clonal.index, "Clone.detected.q"]
    
    if(any(clones.detected.q<0.1 & !is.na(clones.detected.q))){
      pat.table <- Archermuts[pat.index,]
      
      detected <- pat.table[pat.table$Clone.detected.q<0.1,]
      clones.detected <- unique(detected$PyCloneCluster)[ !is.na(unique(detected$PyCloneCluster)) ]
      
      #if all CN fail can have no PyClone
      if(!length(clones.detected)==0){
        
        pat.tree <- extract.tree.mutTable(pat.table,clone.field = "PyCloneCluster",parent.field = "Tree.parent")
        detected.tree <- remove.clones.on.tree(pat.tree, clones.to.keep = clones.detected)
        #does the tree branch?
        branching <- any(duplicated(detected.tree[,1]))
        
        parents.of.detected <- pat.table[pat.table$PyCloneCluster %in% clones.detected,"Tree.parent"]
        parents.of.detected <- parents.of.detected[!is.na(parents.of.detected)]
        
        if(all(parents.of.detected=="root")){
          out <- "Clonal.only"
        } else {
          if(branching==TRUE){
            out <- "Polyphyletic"
          } else {
            if(any(detected$different.to.parent.q < 0.1 & !is.na(detected$different.to.parent.q) & detected$different.to.parent.OR <1)){
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
    out <- rep(out,sum(pat.index))
    return(out)
  }))
  
  ##look for when tree is being violated
  Archermuts$Tree.violation <- Archermuts$different.to.parent.q < 0.1 & !is.na(Archermuts$different.to.parent.q) & Archermuts$different.to.parent.OR > 1.2
  
  #How many regions do we have for each case? (how well sampled are they?)
  library(stringr)
  Archermuts$RegionNumber <- str_count(Archermuts$RegionSum,":")
  
  return(Archermuts)
  
}

