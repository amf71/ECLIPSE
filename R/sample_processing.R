
##==================================##
##                                  ##
## Functions to process each sample ## 
##                                  ##
##==================================##

#' Function to remove the number of estimated supporting reads which derive from
#' background sequencing niose for groups of mutations (in this case clones)
#' @export
remove_background_reads <- function( data., niose_col = 'background_error',
                                     filter_col = 'hard_filtered', varcount_col = 'supporting_reads',
                                     depth_col = 'depth', group_max = 4, mutid_col = 'mutation_id'){
  
  # Order by mut id - then don't need randomisation in sampling and the background 
  # noise removal is deterministic & reproducible
  data. <- data.[ order(get(mutid_col)) ] 
  
  # mutations to ignore throughout
  ignore_i <- data.[, get(filter_col) ]
  
  # If all mutations are not evaluable then return with NA
  if(all(ignore_i)){
    data.[, `:=`(group = NA,
                 supporting_reads_no_niose = NA) ]
    return(data.)
  } 
  
  # Determine the number of groups so that each group has enough mutations
  # that if the noise was distributed evenly each group would have at least 1
  # read of noise - max of 4 groups - also make sure there are fewer groups than
  # mutations (and obviously at least 1 group)
  group_num <- data.[ !ignore_i , floor( sum( get(niose_col) * get(depth_col) ) )]
  if( group_num >= data.[ !ignore_i, length(unique(get(niose_col)))] ) group_num = data.[ !ignore_i, length(unique(get(niose_col))) - 1]
  if( group_num > group_max ) group_num = group_max
  if( group_num == 0 ) group_num = 1
  
  # Now cluster mutations based on background noise 
  data.[ !ignore_i, group := ifelse( .N > 1, kmeans(data.[ !ignore_i, get(niose_col)], group_num)$cluster, 1) ]
  groups <- data.[ !is.na(group), unique(group)]
  
  # Work out hwo mnay reads need removing from each group
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
  
  # Function to remove these supporting reads which we attribute to sequencing noise 
  # Sequences noise comes as a possion distribution - ie we would expect this many mutations
  # with 3 supporting reads, this number with 2 supporting reads and this number with 1
  # supporting reads due to noise. This function will try to find mutations have the correct
  # number of supporting reads to remove and if none exist then still remove enough reads
  # from other mutations (eg there are no mutations with 3 supporting reads we will instead 
  # remove the same number of supporting reads from mutations with 1 / 2 supporting reads). 
  # It also makes sure a given mutation can only have reads removed from it once. 
  subtract_reads <- function( reads, reads.to.remove ){
    
    #if no reads to remove return read unaltered
    if( all(is.na(reads.to.remove)) ) return(reads)
    
    # record which mutations have already had noise removed 
    # (only remove reads once to any mutation)
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
      
      # If too few mutations available where enough reads can be removed then 
      # need to take these reads and spread them to lower read numbers
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
      # don't use sample/runif - just select the first X - as we order by sample then mutation id 
      # should make it deterministic and repsroducible
      selectedi <-posiblei[ 1:reads.to.remove[read_num] ]
      reads[selectedi, already_subtracted := TRUE ]
      reads[selectedi, reads := reads - as.numeric(read_num) ]
    }
    return( reads$reads )
  }
  
  # apply this function for each group of mutations and the corresponding number
  # of supporting reads we estimate is due to background niose and hence should be
  # removed
  data.[ !ignore_i, 
         supporting_reads_no_niose := subtract_reads( get(varcount_col), group_read_to_remove[[.GRP]]), 
         by = group ]
  
  return(data.)
  
}

# Function to combine p values with fisher method
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

#' Function to detect if a large number of mutations in a clone have undergone
#' copy number changes so that our ccf estimates would no longer be valid
#' This function will try to estimate the new CN states for each mutation
#' Sometimes, if mutation losses have occured but not in all cells this 
#' function can also estimate the ccf of the new clone which has emerged 
#' with this CIN.
#' @export
correct_new_CIN <- function( data., varcount_col = 'supporting_reads', depth_col = 'depth', vaf_col = 'vaf_nobackground',
                             tumour_totalCN_col = 'total_cpn', multimodal_p_thes = 0.05,
                             background_reads_col = 'background_error', filter_col = 'hard_filtered', multiplicity_col = 'int_multiplicity',
                             normal_totcn_col = 'normal_cpn'){
  
  # mutations to ignore throughout
  ignore_i <- data.[, get(filter_col) | is.na(get(tumour_totalCN_col)) ]
  
  # if no good quality mutations then simple add NA columns
  if( all(ignore_i) ){
    data.[, `:=`(multi_modal_p = NA,
                 subsequecent_cin = FALSE, 
                 int_multiplicity_preCIN = NA,
                 subsequent_amplfiication = NA, 
                 cn_change = NA,
                 new_clone_ccf = NA)]
    return( data. )
  } 
  
  ## For mutations at CN == 1 or CN == 2 test if the VAF distrebution is mulitmodal
  ## If so that suggests that the CN states have changed

  # when vaf is 0 add lod (ie 1/depth) so the we don't get multimodal vaf due to simply 0s and 1s in varcount
  data.[, vaf_lod := get(vaf_col) ]
  data.[ get(vaf_col) == 0, vaf_lod := 1/get(depth_col) ]
  
  # Only do this if you have a reasonable number rof mutation (20) at either CN == 1 or 2 
  # If you have 20 mutations in both categories then combine the pvalues. 
  # p value can be outputted as 0 by the dip test but actually the limit is actually e-14 
  # in the documentation - can't be 0 if we are to accurately combine p values
  multi_modal_p_cn1 <- ifelse(data.[!ignore_i & get(multiplicity_col) < 1.5, .N] < 20, 
                              NA, data.[ !ignore_i & get(multiplicity_col) < 1.5, max(1e-14, dip.test(vaf_lod)$p.value) ] )
  multi_modal_p_cn2 <- ifelse(data.[!ignore_i & get(multiplicity_col) > 1.5, .N] < 20, 
                              NA, data.[ !ignore_i & get(multiplicity_col) > 1.5, max(1e-14, dip.test(vaf_lod)$p.value) ] )
  p_values <- c(multi_modal_p_cn1, multi_modal_p_cn2)
  p_values <- p_values[ !is.na(p_values) ]
  if(length(p_values) == 0) p_values <- NA
  multi_modal_p <- ifelse( length(p_values) > 1, sumlog(c(p_values[1], p_values[2]))$p, p_values )
  data.[, multi_modal_p := multi_modal_p ]
  
  is_multimodel <- multi_modal_p < multimodal_p_thes & !is.na(multi_modal_p)

  # remove vaf_lod col
  data.[, vaf_lod := NULL ]
  
  # If not multumodal then assign subsequecent_cin to F and add NA columns for other fields
  if( !is_multimodel ){
    data.[, `:=`(subsequecent_cin = FALSE, 
                 int_multiplicity_preCIN = NA,
                 subsequent_amplfiication = NA, 
                 cn_change = NA,
                 new_clone_ccf = NA)]
    return( data. )
  } 
  
  # If the data is multimodel this suggest the CN of many mutations has changed
  # This must be because a new clone has come to dominate the sample and this
  # clone must have had many CN changes
  
  # Replace unobserved variants with the LOD (ie 1/depth) for log transformation
  # and clustering
  vafs_lod <- data.[ !(ignore_i), ifelse(get(vaf_col) > 0, get(vaf_col), 1/get(depth_col) ) ]
  int_multip_old <- data.[ !(ignore_i), get(multiplicity_col) ]
  
  # Clsuter the mutations into thier modes using mcClust
  # capture.output to remove printed messages
  capture.output( BIC <- mclustBIC(vafs_lod) )
  mode1 <- Mclust(vafs_lod, x = BIC)
  classification <- mode1$classification
  classes <- unique(classification)
  # for each mode (class) get the mean vaf so we can order them by CN
  class_vaf <- sapply(classes,function(x) 10^mean(log10(vafs_lod[classification==x])))
  names(class_vaf) <- classes
  
  # Now work out what the new CN is for each mode (class)
  int_multip_new <- int_multip_old
  
  # are any classes not significantly different from backgroound?
  # If this is the case these mutations must be at CN = 0
  class.positive <-  sapply(classes, function(class){
    supporting_reads <- data.[ (!ignore_i), get(varcount_col) ][ classification == class ]
    Expected <- data.[ (!ignore_i), get(background_reads_col) ][ classification == class ]
    result <- poisson.test(sum(supporting_reads),sum(Expected),alternative = "greater")$p.value
    out <- result < 0.01
    return(out)
  })
  names(class.positive) <- classes
  
  # If the loweest vaf class is absent (at backgroun niose) and there at CN == 0
  # then the mode above it must CN = 1 and above that CN = 2 etc - easy
  # If none are absent (ie they are all positive) then this is more tricky
  if(sum(!class.positive)==0){
    
    # How many mutations are in ecah mode/class
    nmuts <- sapply(classes,function(class) sum( classification==class ) )
    # What percentage of those mutations were previously at 1 copy
    class.1mutcpn.perc <- sapply(classes,function(class) sum(int_multip_old[classification==class]==1,na.rm = T)/sum(classification==class))
    
    # Cluster of 1 copy should be one of the bottom 2 clusers with the most muts and the one of these with the most 1 cp original muts
    classcp1 <- classes[ order(class_vaf,decreasing=F)[1:2] ]
    classcp1 <- classcp1[class.1mutcpn.perc[ classes %in% classcp1 ] == max(class.1mutcpn.perc[ classes %in% classcp1 ] )]
    
    # If they both have the same number of 1 muts (no WGD) then two classes should be deletions and cp1 hence 1cp is the higher VAF cluster
    # The class/mode can still be CN = 0 even if its positive (ie > background) because it may not be deleted in all cells
    if(length(classcp1)>1){
      classcp1 <- classcp1[class_vaf[classes %in% classcp1] %in% max(class_vaf[unique(mode1$classification) %in% classcp1])][1]
      int_multip_new[classification==classes[classes==classcp1]] <- 1
    }
    
    # Now we've got CN 1 then CN 2 must be one above that one in VAF
    classes.vaf.ordered <- classes [ order(class_vaf, decreasing = F) ]
    classcp2 <- classes.vaf.ordered[ which( classes.vaf.ordered == classcp1 ) + 1 ]
    if( !is.na(classcp2) ){
      int_multip_new[classification == classcp2 ] <- 2
      
      # Any clusters above 2 we can't be certain what the CN is as there is probably too few of mutations
      # for the clustering to be accurate. Just give them NA CN (later we add a column indicating
      # that they have been amplified)
      if(!classcp2 == classes.vaf.ordered[ length(classes.vaf.ordered) ]){
        classcpamp <- classes.vaf.ordered[ which( classes.vaf.ordered == classcp2 ):length(classes.vaf.ordered) ]
        int_multip_new[classification %in% classcpamp ] <- NA
      }
      
    }
    
    # Is there a class / mode below the CN == 1 mode? if so this must be CN = 0 but some cells still retain
    # the mutation (as the class is > background) and we can estimate the % of cell which have lost these
    # mutations from this (this will be the CCF of the new clone which has grown out with CIN)
    if(!classcp1 == classes.vaf.ordered[1]){
      classcp0 <- classes.vaf.ordered[ 1:(which( classes.vaf.ordered == classcp1 ) - 1) ]
      
      int_multip_new[classification %in% classcp0 ] <- 0
      
      cn1vaf_new <- as.numeric(class_vaf[names(class_vaf) %in% classcp1] - class_vaf[names(class_vaf) %in% classcp0])
      cn1vaf_old <- 10^mean(log10(vafs_lod[classification==classcp0] / int_multip_old [classification==classcp0])) 
      new.clone.ccf <- 1 - (cn1vaf_old / cn1vaf_new)
    } else new.clone.ccf = NA
    
  } else {
    
    # if we have a class which is negative (ie not > background niose) simply assign this to 0
    # and then the one above in baf to CnN = 1 and the one above that to 2 etc 
    positive.classes <- classes[class.positive]
    positive.classes <- positive.classes[order(positive.classes)]
    positive.classes.newCN <- 1:length(positive.classes)
    positive.classes.newCN[positive.classes.newCN>2] <- NA
    int_multip_new[classification %in% positive.classes] <- positive.classes.newCN[match(classification[classification %in% positive.classes],positive.classes)]
    int_multip_new[!classification %in% positive.classes] <- 0
    
    # If no remaining evidence of lost mutations then new clone must be
    # completely dominant
    new.clone.ccf = 1
    
  }
  
  # ork out what the CN change has been for each mutations in the new clone
  int_multip_new_copy <- int_multip_new
  int_multip_new_copy[ is.na(int_multip_new_copy) ] <- 3 # could be more this but will calculate the minimum CN gain for amplified mutations
  CN.change <- as.numeric(int_multip_new_copy) - as.numeric(int_multip_old)
  
  # add these info to the data
  data.[ !ignore_i, (multiplicity_col) := int_multip_new ]
  data.[ !ignore_i, `:=`(subsequecent_cin = TRUE, 
                         int_multiplicity_preCIN = int_multip_old,
                         subsequent_amplfiication = is.na(int_multip_new), 
                         cn_change = CN.change,
                         new_clone_ccf = new.clone.ccf )]
  
  # now recalculate ccf based on new CN states
  data.[, ccf := sapply(1:nrow(data.), function(i) calculate_ccf( data.[i, get(vaf_col)], data.[i, purity], 
                                                                         data.[i, get(tumour_totalCN_col)], 
                                                                         data.[i, get(normal_totcn_col)], 
                                                                         data.[i, get(multiplicity_col)] ) )]
  return(data.)
  
}

#' Function to extract expected noramlised SDs using clonal mutations and higher 
#' ctDNA fraction samples
#' @export
extract_normalised_sd <- function(data., hard_filters, tumour_vaf_col = 'mean_tumour_vaf', tumour_purity_col = 'mean_tumour_cellularity',
                                  tumour_totcn_col = 'total_cpn', varcount_col = 'dao', depth_col = 'ddp', normal_totcn_col = 'normal_cpn',
                                  background_col = 'tnc_error_rate', sample_id_col = 'tracerx_id', tumour_ccf_col = 'mean_tumour_ccf',
                                  is_clonal_col = 'is_clonal', filter_col = 'failed filters', quality_signal_niose = 10){
  
  data.[, hard_filtered := grepl( paste( hard_filters, collapse = "|" ), get(filter_col) ) | is.na(get(background_col)) ]
  
  # calculate raw mutant CN, equation taken from NEJM TRACERx 100 paper (this is an 
  # average mut cn per cell across the tumour)
  data.[, mutCPN := (get(tumour_vaf_col) * (1 / get(tumour_purity_col))) * ((get(tumour_purity_col) * get(tumour_totcn_col)) + 2 * (1 - get(tumour_purity_col))) ]
  # calculate number of mutant copies per mutated cell (i/e. multiplicity of mutation)
  data.[, multiplicity := mutCPN / get(tumour_ccf_col) ]
  # multiplicity < 1 is not possible (causes by noise in data) - correct this
  data.[ multiplicity < 1, multiplicity := 1 ]
  data.[, int_multiplicity := round(multiplicity) ]
  
  data.[, vaf := get(varcount_col) / get(depth_col) ]
  
  data.[, high_quality := hard_filtered == FALSE & !is.na(int_multiplicity) & get(depth_col) > 0 ] 
  data.[, purity_est := calculate_purity( vaf_nobackground, 
                                         get(tumour_totcn_col), get(normal_totcn_col),
                                         int_multiplicity  ),  by = seq_len(nrow(data.)) ]
  data.[, purity := mean( purity_est [ (is_clonal & high_quality) ] ), 
       by = get(sample_id_col) ]
  
  # calculate ccf using NEJM formula
  data.[, ccf := calculate_ccf( vaf,  purity,
                               get(tumour_totcn_col), get(normal_totcn_col),
                               get(multiplicity_col)  ),  by = seq_len(nrow(data)) ]
  # get the clonal purity (/ CN adjusted vaf for each mutation)
  data.[, clonal_purity_mut := ccf * purity ]
  
  data.[ , mutation_present := ppois(get(varcount_col) - 1, (get(background_col)*get(depth_col)), lower.tail = FALSE) < 0.01 ]
  
  # Use the samples where the signal is >10x the background noise (default - can change this)
  clones <- data.[ (high_quality) & get(is_clonal_col) == TRUE,
                   .(ctDNA_frac = mean(clonal_purity_mut, na.rm = TRUE),
                     background = mean( get(background_col) ),
                     sd = sd(clonal_purity_mut, na.rm = TRUE),
                     perc_called_muts = sum(mutation_present)/.N ),
                   by = get(sample_id_col) ][, `:=`(signal = ctDNA_frac / background,
                                                    normsd = sd / ctDNA_frac ) ]
  
  
  
  # restrict to clones where all mutations are observed if there enough
  if( clones[ perc_called_muts == 1, .N > 4 ]) clones <- clones[ perc_called_muts == 1 ]
  
  # Make sure no correlation between signal and mornsd at this threshold - if so suggests LOD is 
  # effecting the normalised SD hence repeatedly increase the theshold by 20% to remove any such cases
  # and there is no correlation
  is_corrleated <- clones[ signal > quality_signal_niose, cor.test(log10(signal), normsd)$p.value < 0.05 ]
  while( is_corrleated ){
    quality_signal_niose <- quality_signal_niose * 1.2
    is_corrleated <- clones[ signal > quality_signal_niose, cor.test(log10(signal), normsd)$p.value < 0.05 ]
  }
  
  # calculate confidence intervals for expected normalised SD for clones in this data
  # hci95 should be used as normalisedSD_max in the outlier_test function
  hci95 <- clones[ signal > quality_signal_niose, mean(normsd) + (1.96 * sd(normsd)) ]
  lci95 <- clones[ signal > quality_signal_niose, mean(normsd) + (-1.96 * sd(normsd)) ]
  normsd <- clones[ signal > quality_signal_niose, mean(normsd) ]
  
  out <- c(normsd, lci95, hci95, clones[ signal > quality_signal_niose, .N], quality_signal_niose)
  names(out) <- c('mean_normsd', 'lci95', 'hci95', 'num_qual_clones', 'signal_to_niose_qual_theshold')
  
  return(out)
  
}

#' Function to identify mutations which have no read support because they have a lower
#' limit of detection (eg not deeply sequenced enough or the multiplicity is too low)
#' @export
identify_low_LOD <- function( data., data_col = 'ccf', LOD_col = 'ccf_LOD',
                          filter_col = 'hard_filtered', mutid_col = 'mutation_id'){
  
  
  # mutations to ignore throughout
  ignore_i <- data.[, get(filter_col) | is.na(get(LOD_col)) ]
  
  ### First work out which variants are unsupported due to low LOD ###
  # for now assign all mutations as good quality LOD and alter this 
  # as we go along
  data.[, matched_lod := TRUE ]
  
  # if no quality mutations or if all mutations are either supported or unsupported
  # all ccfs will be NA is purity = 0
  #  just return data with matched lod = TRUE
  if( data.[ (!ignore_i & !is.na(get(data_col))), all(get(data_col) > 0) | all(get(data_col) == 0)  ] ){
    return( data. )
  } 
  
  # get the median LOD for mutations that have read support
  medianLOD_supported <- data.[ !ignore_i & get(data_col) > 0, median(get(LOD_col)) ]
  
  # get the median of the unsupported variants while removing one at a time
  unsupported_lod <- data.[ !ignore_i & get(data_col) == 0, get(LOD_col) ]
  unsupported_lod <- unsupported_lod[ order(unsupported_lod, decreasing = TRUE) ]
  median_lod_unsupported_range <- sapply( 1:length(unsupported_lod), function(i) median(unsupported_lod[i:length(unsupported_lod)]) )
  
  # Calculate how many unsupported variants to discount due to low LOD
  # This is done so that if the LOD medians are as similar as possible between
  # observed and unobserved variants
  diff <- median_lod_unsupported_range - medianLOD_supported 
  diff_abs <- abs(diff)
  is_best_pos <- diff[ which( diff_abs == min(diff_abs) )[1]] > 0
  
  # if min diff is negative then should keep this otherwise remove 
  num_supported_to_remove <- which( diff_abs == min(diff_abs) )[1] - as.numeric(!is_best_pos)
  mutids_to_remove <- data.[ get(data_col) == 0 ][ order(get(LOD_col), decreasing = TRUE), get(mutid_col) ][ 1:num_supported_to_remove ]
  
  # assign these mutations to LOD unmatched
  data.[ get(mutid_col) %in% mutids_to_remove, matched_lod := FALSE ]
  
  return(data.)
  
}


#' Function to identify outliers until SD is acceptable and identify clones which 
#' have an incoherent vaf distribution (ie probably not a true clone) because the normalised SD is still
#' too high after removing 4 outliers (or if >50% of mutations are removed as outliers)
#' The normalisedSD_max is claculated by applying the extract_normalised_sd to the whole cohort
#' @export
outlier_test <- function( data., data_col = 'cn_adj_vaf', LOD_col = 'cn_adj_vaf_LOD', 
                          normalisedSD_max = 1.33, filter_col = 'hard_filtered', ## TODO eventually should remove default in normalisedSD_max
                          outlier_perc_limit = 0.5, outlier_num_limit = 4 ){
  
  
  # mutations to ignore throughout
  ignore_i <- data.[, get(filter_col) | is.na(get(LOD_col)) ]
  
  # if no good mutations just return NA columns
  if( all(ignore_i) ){
    data.[, `:=`(zscore = NA,
                 is_outlier = NA,
                 clone_normsd = NA,
                 poor_quality_clone = NA) ]
    return( data. )
  } 
  
  ### First work out which variants are unsupported due to low LOD ###
  
  ## Assign Z scores to quality mutations to identify outliers
  # if mutations is unobserved then use the LOD instead to avoid 0s being called as outliers when 
  # they should not be
  data.[ (!ignore_i & matched_lod), data_lod := ifelse(get(data_col) > 0, get(data_col), get(LOD_col) ) ]
  
  data.[ (!ignore_i & matched_lod), 
         zscore := (data_lod - mean(data_lod)) / sd(data_lod) ]
  
  # fucntion to test if the normalises SD is acceptable and if not remove outliers until it is
  outlier_remove <- function( zscores, values, sd_limit = normalisedSD_max ){
    outlier = rep(FALSE, length(zscores))
    sdnorm <- sd(values[!outlier]) / mean(values[!outlier])
    # NA when only one mutation - let this pass through
    while( sdnorm > sd_limit & !is.na(sdnorm) ){
      outlier[!outlier][ which(abs(zscores[!outlier]) == max(abs(zscores[!outlier])))[1] ] <- TRUE
      sdnorm <- sd(values[!outlier]) / mean(values[!outlier])
    }
    return( list(outlier, sdnorm) )
  }
  outlier_test <-  outlier_remove(  data.[ (!ignore_i & matched_lod), zscore], 
                                  data.[ (!ignore_i & matched_lod), data_lod] )
  data.[ (!ignore_i & matched_lod), is_outlier := outlier_test[[1]] ]
  data.[, clone_normsd := outlier_test[[2]] ]
  
  # estimate how many mutations were identified as outliers - if too many then assign this clone as
  # bad quality
  data.[, perc_outlier := sum(is_outlier[(!ignore_i & matched_lod)])/sum((!ignore_i & matched_lod)) ]
  data.[, num_outlier := sum(is_outlier[(!ignore_i & matched_lod)])]
  
  data.[, poor_quality_clone := perc_outlier > outlier_perc_limit | num_outlier > outlier_num_limit ]
  data.[, clone_normsd := ifelse(poor_quality_clone, NA, clone_normsd) ]
  data.[, is_outlier := ifelse(poor_quality_clone, NA, is_outlier) ]
  
  # Remove some data columns we don't need
  data.[, `:=`(perc_outlier = NULL,
               num_outlier = NULL,
               data_lod = NULL )]
  
  return(data.)
  
}

#' Function to estimate our power to significantly detect either sample purity (type = 'power_purity')
#' clone ccf (type = 'power_clone_ccf') or the ccf of the average subclone in a sample (type = 'power_sample_ccf')
#' based on the background noise, copy number status, depth of sequencing etc. Default is for detection at
#' p = 0.01
#' @export
power_calc <- function( data., type, niose_col = 'background_error',
                            filter_col = 'mrd_filtered', normal_cn_col = 'normal_cn', purity_col = 'purity',
                            multiplicity_col = 'int_multiplicity', tumour_totcn_col = 'total_cpn', depth_col = 'depth',
                            clonal_col = 'is_clonal', p_theshold = 0.05, num_subclonal_mutations_sim = 5 ){
  
  # work out which mutations we want depending on the type of power analysis we're doing
  if( type == 'power_purity'){ subset <- data.[, get(clonal_col) ] ; num_muts_correction = 1}
  if( type == 'power_clone_ccf'){ subset <- rep(TRUE, data.[, .N]) ; purity = data.[, unique(get(purity_col)) ] ; num_muts_correction = 1}
  if( type == 'power_sample_ccf'){ 
    subset <- rep(TRUE, data.[, .N])
    purity = data.[, unique(get(purity_col)) ]
    num_muts_correction = data.[ (get(filter_col) & get(clonal_col)), num_subclonal_mutations_sim/.N ]
  }
  
  # if not good mutations jsut reutn the output col with NAs
  if( data.[ (get(filter_col) & subset), .N == 0] ){
    data.[, (type) := NA ]
    return(data.)
  } 
  
  # calculate CCF for this vaf
  # Allow na.rm here for CIN clones which have NA multiplicity for amps
  # and also when at sample level to remove NA clones without CN info
  if( !type == 'power_purity' ){
    
    # calculate vaf using number of observations required for p = 0.01
    lambda <- data.[ (get(filter_col) & subset), sum(get(niose_col)*get(depth_col))*num_muts_correction ]
    signifcant_vaf_95 <- data.[ (get(filter_col) & subset), qpois((1 - p_theshold), lambda = lambda) / (sum(get(depth_col))*num_muts_correction) ]
    
    power_95 <- calculate_ccf(signifcant_vaf_95, purity, data.[(get(filter_col) & subset), mean(get(tumour_totcn_col), na.rm = T)],
                              data.[(get(filter_col) & subset), mean(get(normal_cn_col))], 
                              data.[(get(filter_col) & subset), mean(get(multiplicity_col), na.rm = T)])
    
  } 
  
  # calculate the purity for this vaf as in calc_purity function
  if( type == 'power_purity' ){
    
    # calculate vaf using number of observations required for p = 0.01
    lambda <- data.[ (get(filter_col) & get(clonal_col)), sum(get(niose_col)*get(depth_col))*num_muts_correction ]
    signifcant_vaf_95 <- data.[ (get(filter_col) & get(clonal_col)), qpois((1 - p_theshold), lambda = lambda) / (sum(get(depth_col))*num_muts_correction) ]
    
    power_95 <- data[ (get(filter_col) & get(clonal_col)) , calculate_purity( signifcant_vaf_95, 
                                           mean(get(tumour_totcn_col)) , mean(get(normal_totcn_col)),
                                           mean(get(multiplicity_col)) ) ]

  }
  
  # ccf or purity can't be > 1 - limit as in exome pipeline
  data.[, (type) := ifelse(power_95 > 1, 1, power_95) ]
  
  return(data.)
  
}

#' Function to calculate CCFs from vaf, purity and cn
#' Equation from NEJM 2017 paper
#' @export
calculate_ccf <- function(vaf, purity, tumour_totalCN, normal_totalCN, multiplicity ){
  
  # Sometimes after subsequent CN change a mutation can have multiplicity of 0 - can't calculate ccf for this
  multiplicity <- ifelse( multiplicity == 0, NA, multiplicity)
  
  # if we're missing any of this info (or if the purity is 0) then ccf can't be calculated
  if( any(is.na(c(vaf, purity, tumour_totalCN, normal_totalCN, multiplicity))) | purity == 0 ){
    return( as.numeric(NA) ) 
  }
  
  return( (vaf * 1/purity)*((purity*tumour_totalCN)+normal_totalCN*(1-purity)) / multiplicity )
  
}

#' Function to assess whether subclones are in 100% of cells or not in a given sample
#' ie are th ccfs significantly different from the ccfs of the clonal mtuations
#' @export
call_clonal <- function(data., data_col = 'ccf', 
                        filter_col = 'hard_filtered', 
                        clone_col = 'clone', clonal_col = 'is_clonal',
                        p_thsehold = 0.05 ){
  
  # if not good mutations jsut reutn the output col with NAs
  if( data.[ (!get(filter_col) & !is.na(get(data_col)) & get(clonal_col)), .N == 0 ] ){
    data.[, `:=`(is_subclonal_sample_p = NA,
                 is_subclonal_sample = NA ) ]
    return(data.)
  } 
  
  clonal_ccfs <- data.[ (!get(filter_col) & get(clonal_col)), get(data_col) ]
  #for CIn samples where ccf still unknown for some
  clonal_ccfs <- clonal_ccfs[ !is.na(clonal_ccfs) ]
  
  # suppress warning on ties
  data.[, is_subclonal_sample_p := ifelse( !all(is.na(  get(data_col)[ (!get(filter_col)) ] )), 
                                            suppressWarnings( wilcox.test( get(data_col)[ (!get(filter_col)) ],
                                                                           clonal_ccfs, alternative = 'less')$p.value ), 
                                            as.numeric(NA) ),
        by = eval(clone_col) ]
  
  data.[, is_subclonal_sample := is_subclonal_sample_p < p_thsehold ]
  
  return(data.)
}


#' function to estimate the purity using ccf equation for clonal variants
#' This uses a rearrangement of the equation to calculate ccfs which is validate
#' for mutations we know have a ccf of 1 (those we think are clonal)
#' @export
calculate_purity <- function(vaf, tumour_totalCN, normal_totalCN, multiplicity ){
  
  # Sometimes after subsequent CN change a mutation can have multiplicity of 0 - can't calculate ccf for this
  multiplicity <- ifelse( multiplicity == 0, NA, multiplicity)
  
  # if we're missing any of this info (or if the purity is 0) then ccf can't be calculated
  if( any(is.na(c(vaf, tumour_totalCN, normal_totalCN, multiplicity))) ){
    return( as.numeric(NA) ) 
  }
  
  return( normal_totalCN / ( (multiplicity/vaf) - tumour_totalCN + normal_totalCN ) )
  
}


#' Clonal deconvolution function
#' 
#' 
#' 
#' @export
clonal_deconvolution <- function(data, normalisedSD_max = 0.56, sample_id_col = 'sample_id', niose_col = 'background_error', 
                                 chromosome_col = 'chr', position_col = 'pos', alt_base_col = 'alt', 
                                 varcount_col = 'supporting_reads', depth_col = 'depth', clone_col = 'clone', 
                                 is_clonal_col = 'is_clonal', tumour_vaf_col = 'tumour_vaf', 
                                 tumour_totcn_col = 'total_cpn', normal_totcn_col = 'normal_cpn',
                                 tumour_purity = 'tumour_cellularity', tumour_ccf_col = 'tumour_ccf', 
                                 hard_filtered_col = NA, mrd_filtered_col = NA, multiplicity_col = NA,
                                 testing = FALSE ){
  
  class_origin <- class( data )
  data <- as.data.table( data )
  
  # make a mutation id
  data[, mutation_id := paste( get(sample_id_col), get(chromosome_col), get(position_col), get(alt_base_col), sep = ":" ) ]
  
  # allow to run for multiple samples if inputted as 1 data table
  data[, sample_clone := paste( get(sample_id_col), get(clone_col), sep = "_") ]
  
  # if filter cols are left NA then just use all mutations
  if( all( is.na(hard_filtered_col) ) ){  data[, hard_filtered := FALSE ] ; hard_filtered_col = 'hard_filtered' }
  if( all( is.na(mrd_filtered_col) ) ){  data[, mrd_filtered := FALSE ] ; mrd_filtered_col = 'mrd_filtered' }
  
  message( 'removing background niose..')
  
  # if no background niose info then hard filter
  data[ is.na(get(niose_col)), (hard_filtered_col) := TRUE ]
  
  # filter background noise reads
  data <- rbindlist( lapply(data[, unique(sample_clone) ], function(clone_name){ 
    if( testing ) print(clone_name) 
    remove_background_reads(data[ sample_clone == clone_name ],
                            niose_col = niose_col,
                            filter_col = hard_filtered_col, 
                            varcount_col = varcount_col,
                            depth_col = depth_col, group_max = 4, 
                            mutid_col = 'mutation_id') 
    } ) )
  
  # calculate raw VAF and VAF after removal of background noise
  data[, vaf := get(varcount_col) / get(depth_col) ]
  data[, vaf_nobackground := supporting_reads_no_niose / get(depth_col) ]
  
  # calculate multiplicity using the mutant cn (mut cn averaged accross all tumour cells)
  # and the ccf info from the primary tumour. Calculate the mutant cn using the equation in
  # the NEJM paper. 
  # Option to supply the mulitplciity directly to tool but if multiplicity_col = NA 
  # tool will calclate it as above
  if( is.na(multiplicity_col) ){
  data[, mutCPN := (get(tumour_vaf_col) * (1 / get(tumour_purity))) * ((get(tumour_purity) * get(tumour_totcn_col)) + get(normal_totcn_col) * (1 - get(tumour_purity))) ]
  # calculate number of mutant copies per mutated cell (i/e. multiplicity of mutation)
  data[, multiplicity := mutCPN / get(tumour_ccf_col) ]
  multiplicity_col = 'multiplicity'
  }
  # multiplicity < 1 is not possible (caused by noise in data) - correct this
  # asit will cause issues downstream
  # Average to integer multiplicity (as must be in reality) as done in Tx exome pipeline
  data[  get(multiplicity_col) < 1, (multiplicity_col) := 1 ]
  data[, int_multiplicity := round(get(multiplicity_col)) ]
  
  message( 'calculating ccfs..')
  
  # make a first attempt at purity estimation before we identify outliers & 
  # low LOD variants after which we will recalculate it
  # int_multiplicity will be NA if we don't have CN data for this mutation - exclude
  data[, high_quality := get(hard_filtered_col) == FALSE & !is.na(int_multiplicity) ] 
  data[, purity_est := calculate_purity( vaf_nobackground, 
                               get(tumour_totcn_col), get(normal_totcn_col),
                               get(multiplicity_col)  ),  by = seq_len(nrow(data)) ]
  data[, purity := mean( purity_est [ (is_clonal & high_quality) ] ), 
       by = get(sample_id_col) ]
  
  # calculate ccf using NEJM formula
  data[, ccf := calculate_ccf( vaf_nobackground,  purity,
                        get(tumour_totcn_col), get(normal_totcn_col),
                        get(multiplicity_col)  ),  by = seq_len(nrow(data)) ]
  # get the clonal purity (/ CN adjusted vaf for each mutation)
  data[, clonal_purity_mut := ccf * purity ]
  
  message( 'correcting for subsequent CIN..')
  
  # correct for any subsequent CIN where it is very obvious by detecting
  # multimodal distributions in mutations that should be the same CN and 
  # assigning the mutations each mode thier new CN state
  data <- rbindlist( lapply(data[, unique(sample_clone) ], function(clone_name){ 
    if( testing ) print(clone_name)
    correct_new_CIN(data[ sample_clone == clone_name ],
                    varcount_col = varcount_col, 
                    depth_col = depth_col, 
                    vaf_col = 'vaf_nobackground',
                    tumour_totalCN_col = tumour_totcn_col, 
                    multimodal_p_thes = 0.05,
                    background_reads_col = niose_col, 
                    filter_col = hard_filtered_col, 
                    multiplicity_col = multiplicity_col,
                    normal_totcn_col = normal_totcn_col) 
    } ) )

  message( 'estimating limit of detection..')
  
  # Now calculate LOD for ccf in each mutation - ie what is the ccf if only 1 read was observed?
  data[, ccf_LOD := calculate_ccf( 1/get(depth_col),  purity,
                                   get(tumour_totcn_col), get(normal_totcn_col),
                                   get(multiplicity_col)  ),  by = seq_len(nrow(data)) ]
  data[, clonal_purity_mut_LOD := ccf_LOD * purity ]
  
  message( 'determining poor quality mutations and clones..')
  
  # Annotate which variants in each clones that appear to be unsupported because of
  # a low LOD
  data <- rbindlist( lapply(data[, unique(sample_clone) ], function(clone_name){ 
    if( testing ) print(clone_name)
    identify_low_LOD(data[ sample_clone == clone_name ],
                     data_col = 'ccf', 
                     LOD_col = 'ccf_LOD',
                     mutid_col = 'mutation_id',
                     filter_col = hard_filtered_col) 
    } ) )
  
  # Annotate outliers for each clone (they may have been incorrectly assigned to this cluster or the CN may have changed) 
  # and identify clones with bad ccf distributions (probably not real clones eg could be a clone from neutral tail)
  # normalisedSD_max is calculated using the upper 95% CI for normalised SD in clonal clones of high vaf 
  # (see the extract_normalised_sd function which should be applied to data from a whole cohort)
  data <- rbindlist( lapply(data[, unique(sample_clone) ], function(clone_name){ 
    if( testing ) print(clone_name)
    outlier_test( data[ sample_clone == clone_name ], 
                 data_col = 'clonal_purity_mut', 
                 LOD_col = 'clonal_purity_mut_LOD', 
                 normalisedSD_max = normalisedSD_max, 
                 filter_col = hard_filtered_col ) 
    } ) )
  
  message( 'determining clone ccfs..')
  
  # now account for all the poor wuality mutations we have identified in the past
  # section
  # int_multiplicity can be 0 or NA (when amplified to an unknown amount) if we've detected 
  # subsequent cin - can't use these mutations
  data[, high_quality := matched_lod == TRUE & (is_outlier == FALSE | is.na(is_outlier)) & 
         get(hard_filtered_col) == FALSE & !int_multiplicity == 0 & !is.na(int_multiplicity) ] 
  
  # Now work out the purity again using onlny good quality mutations
  data[, purity_est := calculate_purity( vaf_nobackground, 
                                         get(tumour_totcn_col), get(normal_totcn_col),
                                         get(multiplicity_col)  ),  by = seq_len(nrow(data)) ]
  data[, purity := mean( purity_est [ (is_clonal & high_quality) ] ), 
       by = get(sample_id_col) ]
  
  # and now recalculate the ccfs etc based on this more accurate purity
  # calculate ccf using NEJM formula
  data[, ccf := calculate_ccf( vaf_nobackground,  purity,
                               get(tumour_totcn_col), get(normal_totcn_col),
                               get(multiplicity_col)  ),  by = seq_len(nrow(data)) ]
  # get the clonal purity (/ CN adjusted vaf) for each mutation
  data[, clonal_purity_mut := ccf * purity ]
  
  # Now work out the clone CCFs , 1.96 = Z score for 95% CIs
  data[, clone_ccf := mean(ccf[ (high_quality) ] ), by = sample_clone ]
  data[, clone_ccf_UpCI := mean(ccf[ (high_quality) ] ) + (1.96 * sd(ccf[ (high_quality) ])), by = sample_clone ]
  data[, clone_ccf_LwCI := mean(ccf[ (high_quality) ] )  + (-1.96 * sd(ccf[ (high_quality) ])), by = sample_clone ]

  # limit CCFs > 1 for clones as in the exome pipeline
  data[, clone_ccf := ifelse(clone_ccf > 1, 1, clone_ccf) ]
  data[, clone_ccf_UpCI := ifelse(clone_ccf_UpCI > 1, 1, clone_ccf_UpCI) ]
  data[, clone_ccf_LwCI := ifelse(clone_ccf_LwCI > 1, 1, clone_ccf_LwCI) ]
  
  # clone purity should be proportional to the number of cancer cells present in
  # the body for each clone - use for fishplots - open to suggestions for name! (CN adjusted VAF?)
  data[, clone_purity := clone_ccf * purity ]
  data[, clone_purity_UpCI := clone_ccf_UpCI * purity ]
  data[, clone_purity_LwCI := clone_ccf_LwCI * purity ]
  

  ## Now do a MRD style test for each clone to determine whether it is present or absent from a sample
  ## -1 as ppios returns the pvalue for greater than not greater than or equal to
  data[, clone_present_p := ppois(sum(get(varcount_col)[!get(mrd_filtered_col)]) - 1, 
                                  sum((get(niose_col)*get(depth_col))[!get(mrd_filtered_col)]), lower.tail = FALSE), 
       by = sample_clone]
  
  ## Also do this for all mutations in the sample - mrd test
  data[, mrd_p :=  ppois(sum(get(varcount_col)[!get(mrd_filtered_col)]) - 1, 
                                  sum((get(niose_col)*get(depth_col))[!get(mrd_filtered_col)]), lower.tail = FALSE), 
       by = get(sample_id_col)]
  
  ## Also do this at the mutation level
  data[ (!get(hard_filtered_col)), mutation_present_p := ppois(get(varcount_col) - 1, 
                                     (get(niose_col)*get(depth_col)), lower.tail = FALSE) ]
  
  message( 'performing power calculations..')
  
  data[, high_quality := get(hard_filtered_col) == FALSE & get(mrd_filtered_col) == FALSE & 
                         !int_multiplicity == 0 & !is.na(int_multiplicity) & !is.na(get(niose_col))] 
  
  ## Now do power calculations for each sample clone for what CCF would be detectable given the LOD
  ## Then calculate the equivalent for purity 
  data <- rbindlist( lapply(data[, unique(sample_clone) ], function(clone_name){ 
    if( testing ) print(clone_name)
    power_calc( data[ sample_clone == clone_name ], 
                type = 'power_clone_ccf',
                niose_col = niose_col, 
                filter_col = 'high_quality', 
                normal_cn_col = normal_totcn_col, 
                purity_col = 'purity',
                multiplicity_col = 'int_multiplicity', 
                tumour_totcn_col = tumour_totcn_col, 
                depth_col = depth_col,
                clonal_col = is_clonal_col, 
                p_theshold = 0.01 ) 
    } ) )
  
  data <- rbindlist( lapply(data[, unique(get(sample_id_col)) ], function(sample){ 
    if( testing ) print(sample)
    power_calc( data[ get(sample_id_col) == sample ], 
                type = 'power_sample_ccf',
                niose_col = niose_col, 
                filter_col = 'high_quality', 
                normal_cn_col = normal_totcn_col, 
                purity_col = 'purity',
                multiplicity_col = 'int_multiplicity', 
                tumour_totcn_col = tumour_totcn_col, 
                depth_col = depth_col,
                clonal_col = is_clonal_col,  
                p_theshold = 0.01, 
                num_subclonal_mutations_sim = 5 ) 
    } ) )
  
  data <- rbindlist( lapply(data[, unique(get(sample_id_col)) ], function(sample){ 
    if( testing ) print(sample)
    power_calc( data[ get(sample_id_col) == sample ], 
                type = 'power_purity',
                niose_col = niose_col, 
                filter_col = 'high_quality', 
                normal_cn_col = normal_totcn_col, 
                purity_col = 'purity',
                multiplicity_col = 'int_multiplicity', 
                tumour_totcn_col = tumour_totcn_col, 
                depth_col = depth_col,
                clonal_col = is_clonal_col, 
                p_theshold = 0.01 ) 
    } ) )

  message( 'estimating sample-level clonalty..')

  ## Finally check for each subclone whether the CCFs are significantly different to the clonal cluster in a given sample
  ## If not this indicates this subclone may have become clonal as far as we can detect (though there may still be 
  ## small # of cells not from this subclone that are undetectable but grow back at later time points)
  data <- rbindlist( lapply(data[, unique(get(sample_id_col)) ], function(sample){ 
    if( testing ) print(sample)
    call_clonal(data[ get(sample_id_col) == sample ],
                data_col = 'ccf', 
                filter_col = hard_filtered_col, 
                clone_col = clone_col, 
                clonal_col = is_clonal_col,
                p_thsehold = 0.05 ) 
    } ) )
  
  # remove sample_clone col we made at the start & quality col
  data[, `:=`(sample_clone = NULL,
              high_quality = NULL,
              purity_est = NULL) ]
  
  # if data.frame was inputted then convert back
  if(!any(class_origin == 'data.table')) data <- as.data.frame(data)
  
  return(data)
  
}
 
#=====#
# END #
#=====#

# normalisedSD_max = 0.56; sample_id_col = 'tracerx_id'; niose_col = 'tnc_error_rate';
# chromosome_col = 'chromosome'; position_col = 'position'; alt_base_col = 'alternate';
# varcount_col = 'dao'; depth_col = 'ddp'; clone_col = paste("PyCloneCluster", pyclone_type, sep = '_');
# is_clonal_col = 'is_clonal'; tumour_vaf_col = 'mean_tumour_vaf';
# tumour_totcn_col = 'tumour_total_cpn'; normal_totcn_col = 'normal_total_cpn';
# tumour_purity = 'mean_tumour_cellularity'; tumour_ccf_col = 'mean_tumour_ccf';
# hard_filtered_col = 'hard_filtered'; mrd_filtered_col = 'mrd_filtered'
# multiplicity_col = NA; testing = T