
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
         supporting_reads_no_niose := subtract_reads( get(varcount_col), group_read_to_remove[[.GRP]]), 
         by = group ]
  
  return(data.)
  
}

# Function to combine p value with fisher method
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

#' Function to calculate CCFs from vaf, purity and cn as in ABSOLUTE
#' Allow either a purity or a purity column as input and same for varcount
#' @export
correct_new_CIN <- function( data., varcount_col = 'supporting_reads', depth_col = 'depth', vaf_col = 'vaf_nobackground',
                             tumour_totalCN_col = 'total_cpn', multimodal_p_thes = 0.05,
                             background_reads_col = 'background_error', filter_col = 'hard_filtered', multiplicity_col = 'int_multiplicity',
                             normal_totcn_col = 'normal_cpn'){
  
  # mutations to ignore throughout
  ignore_i <- data.[, get(filter_col) | is.na(get(tumour_totalCN_col)) ]
  
  # p value can be outputted as 0 but actually the limit is actually e-14
  mulit_modal_p_cn1 <- ifelse(data.[!ignore_i & multiplicity < 1.5, .N] < 20, 
                              NA, data.[ !ignore_i & multiplicity < 1.5, max(1e-14, dip.test(get(vaf_col))$p.value) ] )
  mulit_modal_p_cn2 <- ifelse(data.[!ignore_i & multiplicity > 1.5, .N] < 20, 
                              NA, data.[ !ignore_i & multiplicity > 1.5, max(1e-14, dip.test(get(vaf_col))$p.value) ] )
  p_values <- c(mulit_modal_p_cn1, mulit_modal_p_cn2)
  p_values <- p_values[ !is.na(p_values) ]
  if(length(p_values) == 0) p_values <- NA
  mulit_modal_p <- ifelse( length(p_values) > 1, sumlog(c(p_values[1], p_values[2]))$p, p_values )
  data.[, mulit_modal_p := mulit_modal_p ]
  
  is_multimodel <- mulit_modal_p < multimodal_p_thes & !is.na(mulit_modal_p)

  if( !is_multimodel ){
    data.[ !ignore_i, `:=`(subsequecent_cin = FALSE, 
                           int_multiplicity_preCIN = NA,
                           subsequent_amplfiication = NA, 
                           cn_change = NA,
                           new_clone_ccf = NA)]
    return( data. )
  } 
  
  # If the data is multimodel this suggest the CN of many mutations has changed
  # This must be because a new clone has come to dominate the sample and this
  # clone must have had many CN changes
  
  vafs_lod <- data.[ !(ignore_i), ifelse(get(vaf_col) > 0, get(vaf_col), 1/get(depth_col) ) ]
  int_multip_old <- data.[ !(ignore_i), get(multiplicity_col) ]
  
  capture.output( BIC <- mclustBIC(vafs_lod) )
  mode1 <- Mclust(vafs_lod, x = BIC)
  classification <- mode1$classification
  classes <- unique(classification)
  class_vaf <- sapply(classes,function(x) 10^mean(log10(vafs_lod[classification==x])))
  names(class_vaf) <- classes
  
  # are either classes not signifcantly different from backgroound?
  class.positive <-  sapply(classes, function(class){
    supporting_reads <- data.[ (!ignore_i), get(varcount_col) ][ classification == class ]
    Expected <- data.[ (!ignore_i), get(background_reads_col) ][ classification == class ]
    result <- poisson.test(sum(supporting_reads),sum(Expected),alternative = "greater")$p.value
    out <- result < 0.01
    return(out)
  })
  names(class.positive) <- classes
  
  # Now work out what the new CN is
  int_multip_new <- int_multip_old
  
  if(sum(!class.positive)==0){
    
    nmuts <- sapply(classes,function(class) sum( classification==class ) )
    class.1mutcpn.perc <- sapply(classes,function(class) sum(int_multip_old[classification==class]==1,na.rm = T)/sum(classification==class))
    
    #cluster of 1 copy should be one of the top 2 clusers with the most muts and the one of these with the most 1 cp original muts
    classcp1 <- classes[ order(nmuts,decreasing=T)[1:2] ]
    classcp1 <- classcp1[class.1mutcpn.perc[ classes %in% classcp1 ] == max(class.1mutcpn.perc[ classes %in% classcp1 ] )]
    
    #if they both have the same number of 1 muts (no WGD eg) should be deletions and cp1 hence 1cp is the higher VAF cluster
    if(length(classcp1)>1){
      classcp1 <- classcp1[class_vaf[classes %in% classcp1] %in% max(class_vaf[unique(mode1$classification) %in% classcp1])]
      int_multip_new[classification==classes[classes==classcp1]] <- 1
    }
    
    #now we've got CN 1 then CN 2 must be one above in VAF
    classes.vaf.ordered <- classes [ order(class_vaf, decreasing = F) ]
    classcp2 <- classes.vaf.ordered[ which( classes.vaf.ordered == classcp1 ) + 1 ]
    if( !is.na(classcp2) ){
      int_multip_new[classification == classcp2 ] <- 2
      
      if(!classcp2 == classes.vaf.ordered[ length(classes.vaf.ordered) ]){
        classcpamp <- classes.vaf.ordered[ which( classes.vaf.ordered == classcp2 ):length(classes.vaf.ordered) ]
        int_multip_new[classification %in% classcpamp ] <- NA
      }
      
    }
    
    # Now get the mutations that have been lost
    if(!classcp1 == classes.vaf.ordered[1]){
      classcp0 <- classes.vaf.ordered[ 1:which( classes.vaf.ordered == classcp1 ) - 1 ]
      
      int_multip_new[classification %in% classcp0 ] <- 0
      
      cn1vaf_new <- as.numeric(class_vaf[names(class_vaf) %in% classcp1] - class_vaf[names(class_vaf) %in% classcp0])
      cn1vaf_old <- 10^mean(log10(vafs_lod[classification==classcp0] / int_multip_old [classification==classcp0])) ##might need to remove background niose from 0
      new.clone.ccf <- 1 - (cn1vaf_old / cn1vaf_new)
    } else new.clone.ccf = NA
    
  } else {
    
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
  
  int_multip_new_copy <- int_multip_new
  int_multip_new_copy[ is.na(int_multip_new_copy) ] <- 3 # could be more this but will calculate the minimum CN gain for amplified mutations
  CN.change <- as.numeric(int_multip_new_copy) - as.numeric(int_multip_old)
  
  data.[ !ignore_i, `:=`(int_multiplicity = int_multip_new,
                         subsequecent_cin = TRUE, 
                         int_multiplicity_preCIN = int_multip_old,
                         subsequent_amplfiication = is.na(int_multip_new), 
                         cn_change = CN.change,
                         new_clone_ccf = new.clone.ccf )]
  
  # now recalculate ccf based on new CN
  data.[, cn_adj_vaf := sapply(1:nrow(data.), function(i) calculate_ccf( data.[i, get(vaf_col)], 1, 
                                                                         data.[i, get(tumour_totcn_col)], get(normal_totcn_col), 
                                                                         data.[i, int_multiplicity] ) )]
  return(data.)
  
}

#' Function to extract expected noramlised SDs using clonal mutations and higher 
#' ctDNA fraction samples
#' @export
extract_normalised_sd <- function(data., tumour_vaf_col = 'mean_tumour_vaf', tumour_purity_col = 'mean_tumour_cellularity',
                                  tumour_totcn_col = 'total_cpn', varcount_col = 'dao', depth_col = 'ddp', normal_totcn_col = 'normal_cpn',
                                  background_col = 'tnc_error_rate', sample_id_col = 'tracerx_id', tumour_ccf_col = 'mean_tumour_ccf',
                                  is_clonal_col = 'is_clonal', filter_col = 'failed filters', quality_signal_niose = 10){
  
  data[, hard_filtered := grepl( paste( hard_filters, collapse = "|" ), get(filter_col) ) | is.na(get(background_col)) ]
  
  # calculate raw mutant CN, equation taken from NEJM TRACERx 100 paper (this is an 
  # average mut cn per cell across the tumour)
  data.[, mutCPN := (get(tumour_vaf_col) * (1 / get(tumour_purity_col))) * ((get(tumour_purity_col) * get(tumour_totcn_col)) + 2 * (1 - get(tumour_purity_col))) ]
  # calculate number of mutant copies per mutated cell (i/e. multiplicity of mutation)
  data.[, multiplicity := mutCPN / get(tumour_ccf_col) ]
  # multiplicity < 1 is not possible (causes by noise in data) - correct this
  data.[ multiplicity < 1, multiplicity := 1 ]
  data.[, int_multiplicity := round(multiplicity) ]
  
  data.[, purity := 1] ## REMOVE
  
  data.[, vaf := get(varcount_col) / get(depth_col) ]
  
  data.[, cn_adj_vaf := sapply(1:nrow(data.), function(i) calculate_ccf( data.[i, vaf], data.[i, purity], 
                                                                         data.[i, get(tumour_totcn_col)], data.[i, get(normal_totcn_col)], 
                                                                         data.[i, int_multiplicity] ) )]
  
  # Use the samples where the signal is >5x the background noise
  clones <- data.[ get(filter_col) == '' & get(is_clonal_col) == TRUE,
                   .(ctDNA_frac = mean(cn_adj_vaf, na.rm = TRUE),
                     background = mean( get(background_col) ),
                     sd = sd(cn_adj_vaf, na.rm = TRUE)),
                   by = get(sample_id_col) ][, `:=`(signal = ctDNA_frac / background,
                                                    normsd = sd / ctDNA_frac ) ]
  
  # make sure no ocrrelation between signal and mornsd at this theshod - if so increase
  is_corrleated <- clones[ signal > quality_signal_niose, cor.test(log10(signal), normsd)$p.value < 0.05 ]
  while( is_corrleated ){
    quality_signal_niose <- quality_signal_niose * 1.2
    is_corrleated <- clones[ signal > quality_signal_niose, cor.test(log10(signal), normsd)$p.value < 0.05 ]
  }
  
  hci95 <- clones[ signal > quality_signal_niose, mean(normsd) + (1.96 * sd(normsd)) ]
  lci95 <- clones[ signal > quality_signal_niose, mean(normsd) + (-1.96 * sd(normsd)) ]
  normsd <- clones[ signal > quality_signal_niose, mean(normsd) ]
  
  out <- c(normsd, lci95, hci95, clones[ signal > quality_signal_niose, .N], quality_signal_niose)
  names(out) <- c('mean_normsd', 'lci95', 'hci95', 'num_qual_clones', 'signal_to_niose_qual_theshold')
  
  return(out)
  
}

#' Function to identy mutations which have no read support beause they have a lower
#' limit of detection (eg not deeply sequenced enough or mutCN too low)
#' @export
identify_low_LOD <- function( data., data_col = 'cn_adj_vaf', LOD_col = 'cn_adj_vaf_LOD',
                          filter_col = 'hard_filtered'){
  
  
  # mutations to ignore throughout
  ignore_i <- data.[, get(filter_col) | is.na(get(LOD_col)) ]
  
  ### First work out which variants are unsupported due to low LOD ###
  data.[, matched_lod := TRUE ]
  
  if( data.[ !ignore_i, all(get(data_col) > 0) | all(get(data_col) == 0) ] ){
    return( data. )
  } 
  
  # get the median LOD for mutations that have read support
  medianLOD_supported <- data.[ !ignore_i & get(data_col) > 0, median(get(LOD_col)) ]
  
  # get the median of the unsupported variants while removing one at a time
  unsupported_lod <- data.[ !ignore_i & get(data_col) == 0, get(LOD_col) ]
  unsupported_lod <- unsupported_lod[ order(unsupported_lod, decreasing = TRUE) ]
  median_lod_unsupported_range <- sapply( 1:length(unsupported_lod), function(i) median(unsupported_lod[i:length(unsupported_lod)]) )
  
  # how many unsupported variants to discount
  diff <- abs(median_lod_unsupported_range - medianLOD_supported )
  num_supported_to_remove <- which( diff == min(diff) )[1]
  mutids_to_remove <- data.[ get(data_col) == 0 ][ order(get(LOD_col), decreasing = TRUE), get(mutid_col) ][ 1:num_supported_to_remove ]
  
  # assign this to LOD unmatched
  data.[ get(mutid_col) %in% mutids_to_remove, matched_lod := FALSE ]
  
  return(data.)
  
}


#' Function to identy outliers until SD is acceptable and identify clones which 
#' have an incoherent vaf distrebution (ie probably not a true clone) becasue you ahve to 
#' do this too much
#' @export
outlier_test <- function( data., data_col = 'cn_adj_vaf', LOD_col = 'cn_adj_vaf_LOD', 
                          normalisedSD_max = 1.33, filter_col = 'hard_filtered', ## eventually shoudl remove default in ormalisedSD_max
                          outlier_perc_limit = 0.5, outlier_num_limit = 4 ){
  
  
  # mutations to ignore throughout
  ignore_i <- data.[, get(filter_col) | is.na(get(LOD_col)) ]
  
  ### First work out which variants are unsupported due to low LOD ###
  
  if( all(ignore_i) ){
    data.[, `:=`(zscore = NA,
                 is_outlier = NA,
                 clone_normsd = NA,
                 poor_quality_clone = NA) ]
    return( data. )
  } 
  ## do I need to check if normally dist and if not transform?
  
  ## Now assign Z scores to quality mutations to identify outliers
  data.[ (!ignore_i & matched_lod), data_lod := ifelse(get(data_col) > 0, get(data_col), get(LOD_col) ) ]
  
  data.[ (!ignore_i & matched_lod), 
         zscore := (data_lod - mean(data_lod)) / sd(data_lod) ]
  
  outlier_remove <- function( zscores, values, sd_limit = normalisedSD_max ){
    outlier = rep(FALSE, length(zscores))
    sdnorm <- sd(values[!outlier]) / mean(values[!outlier])
    #Na when only one mutation - let this pass through
    while( sdnorm > sd_limit & !is.na(sdnorm) ){
      outlier[ which(abs(zscores) == max(abs(zscores)))[1] ] <- TRUE
    }
    return( list(outlier, sdnorm) )
  }
  
  data.[ (!ignore_i & matched_lod), is_outlier := outlier_remove( zscore, data_lod )[[1]] ]
  data.[, clone_normsd := outlier_remove( zscore[!ignore_i & matched_lod], 
                                          data_lod[!ignore_i & matched_lod] )[[2]] ]
  
  data.[, perc_outlier := sum(is_outlier[(!ignore_i & matched_lod)])/sum((!ignore_i & matched_lod)) ]
  data.[, num_outlier := sum(is_outlier[(!ignore_i & matched_lod)])]
  
  data.[, poor_quality_clone := perc_outlier > outlier_perc_limit | num_outlier > outlier_num_limit ]
  data.[, clone_normsd := ifelse(poor_quality_clone, NA, clone_normsd) ]
  data.[, is_outlier := ifelse(poor_quality_clone, NA, clone_normsd) ]
  
  # Remove some data we don't need
  data.[, `:=`(perc_outlier = NULL,
               num_outlier = NULL,
               data_lod = NULL )]
  
  return(data.)
  
}

#' Function to estimate our power to signifcantly () detect either sample purity (type = power_purity)
#' clone ccf (type = power_clone_ccf) or the ccf of the average subclone in a sample (type = power_sample_ccf)
#' based on the background niose, copy number status, depth of sequencing etc. Default is for detection at
#' p = 0.05
#' @export
power_calc <- function( data., type, niose_col = 'background_error', 
                            filter_col = 'mrd_filtered', normal_cn_col = 'normal_cn', purity_col = 'purity',
                            multiplicity_col = 'int_multiplicity', tumour_totcn_col = 'total_cpn', depth_col = 'depth',
                            clonal_col = 'is_clonal', p_theshold = 0.05, num_subclonal_mutations_sim = 5 ){
  
  if( data.[ (!get(filter_col)), .N == 0] ){
    data.[, (type) := NA ]
    return(data.)
  } 
  
  if( type == 'power_purity'){ subset <- data.[, get(clonal_col) ] ; purity = 1 ; num_muts_correction = 1}
  if( type == 'power_clone_ccf'){ subset <- rep(TRUE, data.[, .N]) ; purity = data.[, unique(get(purity_col)) ] ; num_muts_correction = 1}
  if( type == 'power_sample_ccf'){ 
    subset <- data.[, !get(clonal_col) ] 
    purity = data.[, unique(get(purity_col)) ] 
    num_muts_correction = data.[ (!get(filter_col) & subset), num_subclonal_mutations_sim/.N ]
  }
  
  # calculate vaf using number of observations required for p = 0.05
  lambda <- data.[ (!get(filter_col) & subset), sum(get(niose_col)*get(depth_col))*num_muts_correction ]
  signifcant_vaf_95 <- data.[ (!get(filter_col) & subset), qpois((1 - p_theshold), lambda = lambda) / (sum(get(depth_col))*num_muts_correction) ]
  
  # calculate CCF for this vaf
  # Allow na.rm here for CIN clones which have NA multiplicity for amps
  # and also when at sample level to remove NA clones without CN info
  power_95 <- calculate_ccf(signifcant_vaf_95, purity, data.[(!get(filter_col) & subset), mean(get(tumour_totcn_col), na.rm = T)],
                                data.[(!get(filter_col) & subset), mean(get(normal_cn_col))], 
                                data.[(!get(filter_col) & subset), mean(get(multiplicity_col), na.rm = T)])
  data.[, (type) := ifelse(power_95 > 1, 1, power_95) ]
  
  return(data.)
  
}

#' Function to assess whether subclones are in 100% of cells or not in a sample
#' @export
call_clonal <- function(data., data_col = 'cn_adj_vaf', 
                        filter_col = 'hard_filtered', 
                        clone_col = 'clone', clonal_col = 'is_clonal',
                        p_thsehold = 0.05 ){
  
  clonal_cn_adj_vafs <- data.[ (!get(filter_col) & get(clonal_col)), get(data_col) ]
  #for CIn samples where cn_ajd vaf still unknown for some
  clonal_cn_adj_vafs <- clonal_cn_adj_vafs[ !is.na(clonal_cn_adj_vafs) ]
  
  data.[, is_subclonal_sample_p := ifelse( !all(is.na(  get(data_col)[ (!get(filter_col)) ] )), 
                                            wilcox.test( get(data_col)[ (!get(filter_col)) ],
                                                         clonal_cn_adj_vafs, alternative = 'less')$p.value, 
                                            NA),
        by = get(clone_col) ]
  
  data.[, is_subclonal_sample := is_not_clonal_sample_p < p_thsehold ]
  
  return(data.)
}


#' Clonal deconvolution function
#' 
#' 
#' 
#' @export
clonal_deconvolution <- function(data, sample_id_col = 'sample_id', niose_col = 'background_error', 
                                 chromosome_col = 'chr', position_col = 'pos', alt_base_col = 'alt', 
                                 varcount_col = 'supporting_reads', depth_col = 'depth', clone_col = 'clone', 
                                 is_clonal_col = 'is_clonal', tumour_vaf_col = 'tumour_vaf', 
                                 tumour_totcn_col = 'total_cpn', normal_totcn_col = 'normal_cpn',
                                 tumour_purity = 'tumour_cellularity', tumour_ccf_col = 'tumour_ccf', 
                                 hard_filtered_col = NA, mrd_filtered_col = NA ){
  
  class_origin <- class( data )
  data <- as.data.table( data )
  
  # make a mutation id
  data[, mutation_id := paste( get(sample_id_col), get(chromosome_col), get(position_col), get(alt_base_col), sep = ":" ) ]
  
  # allow to run for multiple samples if inputted as 1 data table
  data[, clone_orig := get(clone_col) ]
  data[, clone := paste( get(sample_id_col), get(clone_col), sep = "_") ]
  
  # if filter cols are left NA then just use all mutations
  if( all( is.na(hard_filtered_col) ) ){  data[, hard_filtered := FALSE ] ; hard_filtered_col = 'hard_filtered' }
  if( all( is.na(mrd_filtered_col) ) ){  data[, mrd_filtered := FALSE ] ; mrd_filtered_col = 'mrd_filtered' }
  
  # if no background info hard filter
  data[ is.na(get(niose_col)), (hard_filtered_col) := TRUE ]
  
  # filter background noise reads
  data <- rbindlist( lapply(data[, unique(clone) ], function(clone_name) remove_background_reads(data[ clone == clone_name ],
                                                                                                 niose_col = niose_col,
                                                                                                 filter_col = hard_filtered_col, 
                                                                                                 varcount_col = varcount_col,
                                                                                                 depth_col = depth_col, group_max = 4, 
                                                                                                 mutid_col = 'mutation_id') ) )
  
  # calculate VAF
  data[, vaf := get(varcount_col) / get(depth_col) ]
  data[, vaf_nobackground := supporting_reads_no_niose / get(depth_col) ]
  
  # calculate raw mutant CN, equation taken from NEJM TRACERx 100 paper (this is an 
  # average mut cn per cell across the tumour)
  data[, mutCPN := (get(tumour_vaf_col) * (1 / get(tumour_purity))) * ((get(tumour_purity) * get(tumour_totcn_col)) + get(normal_totcn_col) * (1 - get(tumour_purity))) ]
  # calculate number of mutant copies per mutated cell (i/e. multiplicity of mutation)
  data[, multiplicity := mutCPN / get(tumour_ccf_col) ]
  # multiplicity < 1 is not possible (causes by noise in data) - correct this
  data[  multiplicity < 1, multiplicity := 1 ]
  data[, int_multiplicity := round(multiplicity) ]
  
  # calculate cn_adj_vaf (ccf but if purity was 1)
  data[, cn_adj_vaf := sapply(1:nrow(data), function(i) calculate_ccf( data[i, vaf_nobackground], 1, 
                                                                       data[i, get(tumour_totcn_col)], data[i, get(normal_totcn_col)], 
                                                                       data[i, int_multiplicity] ) )]
  
  # correct for any subsequent CIN where it is very obvious
  data <- rbindlist( lapply(data[, unique(clone) ], function(clone_name) correct_new_CIN(data[ clone == clone_name ],
                                                                                         varcount_col = varcount_col, 
                                                                                         depth_col = depth_col, 
                                                                                         vaf_col = 'vaf_nobackground',
                                                                                         tumour_totalCN_col = tumour_totcn_col, 
                                                                                         multimodal_p_thes = 0.05,
                                                                                         background_reads_col = niose_col, 
                                                                                         filter_col = hard_filtered_col, 
                                                                                         multiplicity_col = 'int_multiplicity',
                                                                                         normal_totcn_col = normal_totcn_col) ) )

  # Now calculate LOD for this value in each mutation - ie what is ony 1 read was observed?
  data[, cn_adj_vaf_LOD := sapply(1:nrow(data), function(i) calculate_ccf( data[i, 1/get(depth_col)], purity = 1, 
                                                                           data[i, get(tumour_totcn_col)], 
                                                                           data[i, get(normal_totcn_col)], 
                                                                           data[i, int_multiplicity] ) )]
  
  # Annotate which variants in each clones that appear to be unsupported because of
  # a low LOD, also note outliers for each clone (they may have been incorrectly assigned to this cluster or have the wrong CN) 
  # and clones with bad vaf distributions (probably not real clones eg could be from neutral tail)
  data <- rbindlist( lapply(data[, unique(clone) ], function(clone_name) identify_low_LOD(data[ clone == clone_name ],
                                                                                          data_col = 'cn_adj_vaf', 
                                                                                          LOD_col = 'cn_adj_vaf_LOD',
                                                                                          filter_col = hard_filtered_col) ) )
  
  data <- rbindlist( lapply(data[, unique(clone) ], function(clone_name) outlier_test(data[ clone == clone_name ], 
                                                                                          data_col = 'cn_adj_vaf', 
                                                                                          LOD_col = 'cn_adj_vaf_LOD', 
                                                                                          normalisedSD_max = 1.33, 
                                                                                          filter_col = hard_filtered_col) ) )
  
  # Now work out the cellularity
  # int_multiplicity can be 0 or NA (when amplified to an unknown amount) if we've detected subsequent cin
  data[, purity := mean(cn_adj_vaf[ get(is_clonal_col) == TRUE & matched_lod == TRUE & 
                                      outlier == FALSE & get(hard_filtered_col) == FALSE & 
                                      !int_multiplicity == 0 & !is.na(int_multiplicity) ]), by = get(sample_id_col) ]

  # Now work out the clone CCFs
  data[, clone_cn_adj_vaf := mean(cn_adj_vaf[ matched_lod == TRUE & 
  outlier == FALSE & get(hard_filtered_col) == FALSE & 
    !int_multiplicity == 0 & !is.na(int_multiplicity)  ]), by = get(clone_col)]
  data[, ccf := cn_adj_vaf / ifelse(purity == 0, NA, purity) ]
  data[, clone_ccf := mean(ccf[ matched_lod == TRUE & 
                                            outlier == FALSE & get(hard_filtered_col) == FALSE & 
                                            !int_multiplicity == 0 & !is.na(int_multiplicity)  ]), by = get(clone_col)]
  data[, clone_ccf := ifelse(all(is.na(ccf)), NA, clone_ccf), by = get(clone_col)]
  data[, clone_ccf := ifelse(clone_ccf > 1, 1, clone_ccf) ]
  
  ## Now do a MRD style test for each clone to determine whether it is present or absent from a sample
  ## -1 as ppiose returns the pvalue for greater than not greater than or equal to
  data[, clone_present_p := ppois(sum(get(varcount_col)[!mrd_filtered_col]) - 1, 
                                  sum((get(background_col)*get(depth_col))[!get(mrd_filtered_col)]), lower.tail = FALSE), 
       by = get(clone_col)]
  
  ## Also do this for all mutations in the sample - mrd test
  data[, mrd_p :=  ppois(sum(get(varcount_col)[!mrd_filtered]) - 1, 
                                  sum((get(background_col)*get(depth_col))[!get(mrd_filtered_col)]), lower.tail = FALSE), 
       by = get(sample_id_col)]
  
  ## Also do this at the mutation level
  data[ (!get(hard_filtered_col)), mutation_present_p := ppois(get(varcount_col) - 1, 
                                     (get(background_col)*get(depth_col)), lower.tail = FALSE) ]
  
  ## Now do power calculations for each sam clone for what CCF would be detectable given the LOD
  ## Then calculate the equivalent for purity 
  data <- rbindlist( lapply(data[, unique(clone) ], function(clone_name) power_calc( data[ clone == clone_name ], 
                                                                                     type = 'power_clone_ccf',
                                                                                     niose_col = niose_col, 
                                                                                     filter_col = mrd_filtered_col, 
                                                                                     normal_cn_col = normal_cn_col, 
                                                                                     purity_col = 'purity',
                                                                                     multiplicity_col = 'int_multiplicity', 
                                                                                     tumour_totcn_col = tumour_totcn_col, 
                                                                                     depth_col = depth_col,
                                                                                     clonal_col = is_clonal_col, 
                                                                                     p_theshold = 0.05, 
                                                                                     num_subclonal_mutations_sim = 5 ) ) )
  
  data <- rbindlist( lapply(data[, unique(sample_id) ], function(sample) power_calc( data[ sample_id == sample ], 
                                                                                     type = 'power_sample_ccf',
                                                                                     niose_col = niose_col, 
                                                                                     filter_col = mrd_filtered_col, 
                                                                                     normal_cn_col = normal_cn_col, 
                                                                                     purity_col = 'purity',
                                                                                     multiplicity_col = 'int_multiplicity', 
                                                                                     tumour_totcn_col = tumour_totcn_col, 
                                                                                     depth_col = depth_col,
                                                                                     clonal_col = is_clonal_col,  
                                                                                     p_theshold = 0.05, 
                                                                                     num_subclonal_mutations_sim = 5 ) ) )
  
  data <- rbindlist( lapply(data[, unique(sample_id) ], function(sample) power_calc( data[ sample_id == sample ], 
                                                                                     type = 'power_purity',
                                                                                     niose_col = niose_col, 
                                                                                     filter_col = mrd_filtered_col, 
                                                                                     normal_cn_col = normal_cn_col, 
                                                                                     purity_col = 'purity',
                                                                                     multiplicity_col = 'int_multiplicity', 
                                                                                     tumour_totcn_col = tumour_totcn_col, 
                                                                                     depth_col = depth_col,
                                                                                     clonal_col = is_clonal_col, 
                                                                                     p_theshold = 0.05, 
                                                                                     num_subclonal_mutations_sim = 5 ) ) )
  
  ## Finally check for each subclone whether the CCFs are significantly different to the clonal cluster
  ## If not this indicates this subclone has become clonal as far as we can detect (though there may still be 
  ## small # of cells from this subclone that are undetectable but grow back at later time points)
  data <- rbindlist( lapply(data[, unique(sample_id) ], function(sample) call_clonal(data[ sample_id == sample ],
                                                                                     data_col = 'cn_adj_vaf', 
                                                                                     filter_col = hard_filtered_col, 
                                                                                     clone_col = clone_col, 
                                                                                     clonal_col = is_clonal_col,
                                                                                     p_thsehold = 0.05) ) )
  
  # correct clone col
  data[, clone := clone_orig ]
  data[, clone_orig := NULL ]
  
  # if data.frame was inputted they convert back
  if(!any(class_origin == 'data.table')) data <- as.data.frame(data)
  
  return(data)
  
}
 
#=====#
# END #
#=====#

