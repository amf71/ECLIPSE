##======================================================##
##                                                      ##
## Functions to process each all samples from a patient ## 
##                                                      ##
##======================================================##

# Author : Alexander Frankell
# Date: 14-09-2020


#' Function to plot Fish plots from the processed data across a patient 
#' 
#' 
#' @export
draw_fish_plot <- function( data, clinical_data = NA, tree, min_clone_muts = 2 ){
  
  data <- as.data.frame( data )
  tree <- as.matrix( tree )
  if(!all(is.na(clinical_data))) clinical_data <- as.data.frame(clinical_data)
  # limit to 1 row per clone (only need clone level info)
  data <- data[ !duplicated(paste(data$clone, data$sample_id)), ]
  
  data <- data[ data$muts_followed_in_clone >= min_clone_muts, ]
  accepted.clones <- unique( data$clone[data$clone %in% as.character(unique(c(tree[,1],tree[,2])))] )
  
  data <- data[ data$clone %in% accepted.clones,]
  data$ctDNA_fraction[is.na(data$ctDNA_fraction)] <- 0
  
  ##adjust by whether clones/MRD is detected
  data$ctDNA_fraction[ data$Clone.detected.p > 0.05 | is.na(data$Clone.detected.p) ] <- 0
  
  #subset the clones that we can track - remove warnings where no clones are being removed
  suppressWarnings( tree <- remove.clones.on.tree(tree, clones.to.keep = as.character(accepted.clones) ) )
  
  #can't plot if only clonal looks real
  root <- find_root( tree )
  if(all(tree == root)){
    stop("only clonal clone accepted")
  }
  
  clones <- as.character( unique(c(tree[,1],tree[,2])) )
  clones <- clones[clones %in% as.character( accepted.clones) ]
  
  ##make objects for fish plotting
  days <- unique(data$days_post_surgery)[order(unique(data$days_post_surgery),decreasing = F)]
  days <- days[!is.na(days)]
  timepoints <- 1:length(days)
  
  frac.table <- do.call(cbind, lapply(days, function(time){
    out <- data[data$days_post_surgery==time,]
    out <- out[match(clones,out$clone),"ctDNA_fraction"]
    return(out)
  }))
  
  if(all(frac.table==0 | is.na(frac.table))){
    print("all samples negative")
    return(NULL)
  }
  
  frac.table[is.na(frac.table)] <- 0
  
  options(scipen = 999)
  #also add clones from primary tumour data at the begining
  frac.table <- cbind(sapply(clones, function(clone) median(data[data$clone==clone,"tumour_ccf"]))*max(frac.table),frac.table)
  
  parent.order <- unique(tree[,1]) 
  is.daughter <- parent.order %in% unique(tree[,2])
  is.parent <- unique(tree[,2]) %in%  parent.order
  parent.order <- rev(c(parent.order[!is.daughter],
                        unique(tree[,2])[is.parent]))
  
  tree <- logically.order.tree(tree)

  parent.order <- rev(unique(tree[,1]))
  
  for(col in 1:ncol(frac.table)){
    #don't let subclones push up the total ctDNA fraction - should be determined by the clonal muts
    clonal.frac <- frac.table[clones == root,col]
    frac.table[frac.table[,col]>clonal.frac,col] <- clonal.frac
    
    #don't let daughters have higher CCF than parents (if so make parents the size of daughters)
    #for each parent each the are >/equal to all daughters
    for(parent in parent.order){
      parentrow <- which(clones==parent)
      daughters <- tree[tree[,1]==parent,2]
      daughterrows <- which(clones%in%daughters)
      parent.ctDNA_fraction <- frac.table[parentrow,col]
      daughterctDNA_fraction <-  sum(frac.table[daughterrows,col])
      if(parent.ctDNA_fraction < daughterctDNA_fraction*1.002){
        frac.table[parentrow,col] <- daughterctDNA_fraction*1.002
      }
    }
  }
  
  #adjust the M-seq back to max plasma fraction
  frac.table[,1] <- frac.table[,1] * (max(frac.table[,2:ncol(frac.table)])/max(frac.table[,1]))
  
  #add a correction to see the scale properly
  frac.table <- frac.table/(max(frac.table)/0.95)*100
  
  #remove clones not present (even in primary...)
  clones <- clones[!rowSums(frac.table) == 0]
  frac.table <- frac.table[!rowSums(frac.table) == 0,]
  
  parents <- as.numeric(sapply(clones, function(clone){
    if(all(!tree[,2]==clone)){
      return(0)
    } else {
      return(which(clones==tree[tree[,2]==clone,1]))
    }
  }))
  
  #doesn't allow present -> absent -> present - make absent vv small 0.00001 so it can still come back when goes below LOD
  
  n <- nrow(frac.table)
  suppressMessages(cols <- brewer.pal(n, "Paired"))
  if(length(cols)<n){
    diff <- n - length(cols)
    cols <- c(cols,cols[1:diff])
  }
  if(length(cols)>n){
    cols <- cols[1:n]
  }
  
  days.orig <- days
  ##need to sapce out days to be visualisable
  days.gap <- sapply(2:length(days), function(i) days[i]-days[i-1])
  if(any(days.gap < max(days)/50)){
    for(day in days[c(FALSE,days.gap < max(days)/50)]){
      diff <- max(days)/50 - days.gap[days[2:length(days)] == day] 
      days[which(days ==day):length(days)] <- days[which(days ==day):length(days)] + diff
    }
  }
  
  days.scale.factor <- max(days)/length(days)
  
  fish = createFishObject(frac.table,parents,timepoints=c(0,days+days.scale.factor),fix.missing.clones=TRUE,col = cols)
  
  #calculate the layout of the drawing
  fish = layoutClones(fish)
  
  #draw the plot, using the splining method (recommended)
  #and providing both timepoints to label and a plot title
  fishPlot(fish = fish,shape="spline",title=unique(data$patient_id),
           cex.title=2, vlines=c(0,days+days.scale.factor), col.vline="grey",
           vlab=c("Tumour M-seq",days.orig),bg.type = "solid",bg.col = "white",
           Clinical.data = clinical_data, days.scale.factor = days.scale.factor)
  
  # fish
  
  
}


### modified version of fishPlot function from fishPlot R package (source = CRAN)
fishPlot <- function(fish,shape="polygon", vlines=NULL, col.vline="#FFFFFF99", vlab=NULL,
border=0.5, col.border="#777777", pad.left=0.2, ramp.angle=0.5,
title=NULL, title.btm=NULL, cex.title=NULL, cex.vlab=0.7,
bg.type="gradient", bg.col=c("bisque","darkgoldenrod1","darkorange3"),
Clinical.data = NA,days.scale.factor=NA, angle.vlab = 0, adjust.vlab = c(0,1)){
  
  #make sure we have the right number of colors
  checkCol(fish)
  
  pad = (max(fish@timepoints)-min(fish@timepoints))*pad.left;
  
  ##record the plot as an object so I can add tto it
  
  #set up the plot
  plot(-100,-100,col="white",
       ylim=c(0,100),
       xlim=c(min(fish@timepoints)-pad, max(fish@timepoints)),
       yaxt="n", xaxt="n",
       bty="n", xlab="", ylab="")
  
  lim=par()
  bckImage = png::readPNG(createBackgroundImage(bg.col))
  ##create raster background image for smooth gradient
  if(bg.type=="gradient"){
    rasterImage(bckImage, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
  }
  ##add background color to plot
  if(bg.type=="solid"){
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bg.col)
  }
  #(if neither is set, bg will just be white)
  
  ##draw the clusters one at a time, being sure that parents go before children
  parentsList = 0
  while(length(parentsList) > 0){
    parent = parentsList[[1]]
    children =  which(fish@parents==parent)
    parentsList = parentsList[-1]
    parentsList = c(parentsList, children)
    numChildren = length(children)
    for(i in children){
      pad.left=pad
      if(parent>0){
        pad.left=pad*0.4
      }
      
      if(shape=="bezier"){
        drawClustBezier(fish@xpos[[i]], fish@ytop[[i]], fish@ybtm[[i]],
                        fish@col[i], fish@nest.level[i],
                        pad.left=pad.left, border=border, col.border=col.border,
                        annot = fish@clone.annots[i], annot.angle=fish@clone.annots.angle,
                        annot.col=fish@clone.annots.col)
      } else {
        if(shape=="spline"){
          drawClustSpline(fish@xpos[[i]], fish@ytop[[i]], fish@ybtm[[i]],
                          fish@col[i], fish@nest.level[i],
                          pad.left=pad.left, border=border, col.border=col.border,
                          annot = fish@clone.annots[i], annot.angle=fish@clone.annots.angle,
                          annot.col=fish@clone.annots.col)
        } else {
          if(!shape=="polygon"){
            print(paste("unknown shape \"",shape,"\". Using polygon representation"))
          }
          drawClustPolygon(fish@xpos[[i]], fish@ytop[[i]], fish@ybtm[[i]],
                           fish@col[i], fish@nest.level[i], ramp.angle=ramp.angle,
                           pad.left=pad.left, border=border, col.border=col.border,
                           annot = fish@clone.annots[i], annot.angle=fish@clone.annots.angle,
                           annot.col=fish@clone.annots.col)
        }
      }
    }
  }
  #draw timepoint labels/lines
  if(!is.null(vlines)){
    abline(v=vlines,col=col.vline,xpd=F)
    
    if(!is.null(vlab)){
      text(vlines,103,vlab,pos=3,cex=cex.vlab,col="grey20",
           xpd=NA, adj = adjust.vlab, srt = angle.vlab)
    }
  }
  
  if(!is.null(title)){
    #get the center
    xmax = tail(fish@timepoints,n=1)
    cent = (xmax/2)-(pad/2)
    text(cent,112,title,pos=3,cex=cex.title,xpd=T)
  }
  
  
  if(!is.null(title.btm)){
    text(min(fish@timepoints)-(pad*1.2),2,title.btm,pos=4,cex=cex.title)
  }
  
  if(!all(is.na(Clinical.data)) & nrow(Clinical.data)>0){
    
    Clinical.data[,c(3:6,8)] <- apply(Clinical.data[,c(3:6,8)],2,as.numeric)
    
    Clinical.data[,c(3:5,8)] <- Clinical.data[,c(3:5,8)] + days.scale.factor
    
    relapse.days <- unique(Clinical.data$DFS)
    
    if(!is.na(relapse.days)){
      points(x = relapse.days, y = 50,pch = 23)
    }
    
    treatment <- Clinical.data[!is.na(Clinical.data$Start_date) & !Clinical.data$Description == 'sample',]
    
    #treatment.cols <- 
    treatment.type.old <- ""
    offset <- 0
    
    if(sum(!is.na(treatment$SampleID))>0){
      for(i in 1:nrow(treatment)){
        
        treatment.type <- treatment$Description[i]
        rect(xleft = treatment$Start_date[i], ybottom = -10, xright = treatment$End_date[i], ytop = 110,
             col = "#FFFF0080", border = NA)
        if(!treatment.type==treatment.type.old){
          text(treatment$Start_date[i],-12 + offset,treatment.type,pos=3,cex=1.5,col="grey20",xpd=NA)
          offset <- offset - 4
        }
        
        treatment.type.old <- treatment.type
      }
    }
    
  }
  
  
}

plot.full.data <- function(data, pateint){
  
  #CNadjVAF NA when no signal or clone == NA
  data <- data[!is.na(data$CN.adj.VAF),]
  
  #remove samples wheere no MRD
  data <- data[!data$MRD.call.sample.p < 0.05,]
  
  if(nrow(data)==0){
    print("no MRD positive samples")
    return(NULL)
  }
  
  data$mut.class <- "PASS"
  data[data$matched_lod==FALSE,"mut.class"] <- "removed.to.match.LOD"
  data[!data$hard_filtered %in% FALSE,"mut.class"] <- "Library_Qual_filtered"
  
  data[data$CN.adj.VAF==0,"CN.adj.VAF"] <- 0.000001
  data[data$ctDNA_fraction==0 & !is.na(data$ctDNA_fraction),"ctDNA_fraction"] <- 0.000001
  
  #missing.viol <- setdiff(unique(data$sample.clone),unique(modeleddists$sample.clone)) 
  
  #for now remove negative CNadjVAF
  data <- data[!(is.na(data$CN.adj.VAF) | data$CN.adj.VAF<0),]
  
  data[, sample_clone := paste0(sample_id, clone) ]
  
  modeleddists <- do.call(rbind,lapply(unique(data$sample.clone), function(sample.clone){
    sample <- strsplit(as.character(sample.clone),split=" ")[[1]][1]
    days <- unique(data[data$sample_id==sample ,days_post_surgery])
    clone <- strsplit(as.character(sample.clone),split=" ")[[1]][2]
    
    m <- unique(as.numeric(data[data$sample.clone==sample.clone,"ctDNA_fraction"]))
    
    if(!m==0 & !is.na(m)){

      s <- as.numeric(unique(data[data$sample.clone==sample.clone,"ctDNA_fraction"]))*as.numeric(unique(data[data$sample.clone==sample.clone,"normSDsolution"]))
      location <- log(m)
      shape <- sqrt(log(1 + (s^2 / m^2)))
      modeldist <- rlnorm(n = 10000, location, shape)
      data.frame(sample.clone = sample.clone, sample=sample, clone=clone, days_post_surgery=days, modeldist = modeldist, stringsAsFactors = F)
    } else {
      return(NULL)
    }
  }))
  
  data$CN.adj.VAF <- as.numeric(data$CN.adj.VAF)
  data$ctDNA_fraction <- as.numeric(data$ctDNA_fraction)
  modeleddists$modeldist <- as.numeric(modeleddists$modeldist)
  data$clone <- as.factor(as.character(data$clone))
  data$days_post_surgery <- as.factor(as.character(data$days_post_surgery))
  modeleddists$clone <- as.factor(as.character(modeleddists$clone))
  modeleddists$days_post_surgery <- as.factor(as.character(modeleddists$days_post_surgery))
  
  #add in fake data above the yaxis limits (not shown) to ensure the factors all line up
  modeleddists.dup <- modeleddists
  data.dup <- data
  data.dup2 <- data[!duplicated(data$sample.clone),]
  data.dup2 <- rbind(data.dup2,data.dup2,data.dup2)
  data.dup2$mut.class <- c(rep("PASS",nrow(data.dup2)/3),
                                  rep("removed.to.match.LOD",nrow(data.dup2)/3),
                                  rep("Library_Qual_filtered",nrow(data.dup2)/3))
  data.dup2$CN.adj.VAF <- 10
  all.sample.clones <- data$sample.clone[!duplicated(data$sample.clone)]
  missing.sample.clones <- setdiff(all.sample.clones,unique(modeleddists.dup$sample.clone))
  days <- data$days_post_surgery[!duplicated(data$sample.clone)]
  missing.days <- days[all.sample.clones %in% missing.sample.clones]
  
  modeleddists.dup2 <- do.call(rbind,lapply(1:length(missing.sample.clones), function(i){
    data.frame(sample.clone = missing.sample.clones[i],
               sample_id = sapply(strsplit(missing.sample.clones[i],split=" "),"[[", 1),
               clone = sapply(strsplit(missing.sample.clones[i],split=" "),"[[", 2),
               days_post_surgery = missing.days[i], 
               modeldist = rep(0.000001,10),
               #modeldist =  make.random.lognormal(mean =0.000001,sd = 0.000001,n=10),
               stringsAsFactors = F)}))
  modeleddists.dup <- rbind(modeleddists.dup,modeleddists.dup2)
  data.dup <- rbind(data.dup,data.dup2)
  
  modeleddists.dup$days_post_surgery <- factor(modeleddists.dup$days_post_surgery,levels=unique(modeleddists.dup$days_post_surgery)[order(as.numeric(as.character(unique(modeleddists.dup$days_post_surgery))),decreasing = F)])
  data.dup$days_post_surgery <- factor(data.dup$days_post_surgery,levels=unique(data.dup$days_post_surgery)[order(as.numeric(as.character(unique(data.dup$days_post_surgery))),decreasing = F)])
  
  modeleddists.dup$clone <- factor(modeleddists.dup$clone,levels=unique(modeleddists.dup$clone)[order(as.numeric(as.character(unique(modeleddists.dup$clone))),decreasing = F)])
  data.dup$clone <- factor(data.dup$clone,levels=unique(data.dup$clone)[order(as.numeric(as.character(unique(data.dup$clone))),decreasing = F)])
  
  ggplot(data.dup,aes(x = days_post_surgery, y = CN.adj.VAF, fill = clone)) + 
    geom_violin(data=modeleddists.dup,aes(x = days_post_surgery, fill = clone, y = modeldist),alpha = 1,scale = "width") +
    geom_jitter(data=data.dup[data.dup$mut.class=="PASS",], aes(x = days_post_surgery, fill = clone),alpha = 0.5,position=position_jitterdodge(dodge.width=0.9,jitter.width=0.4),colour="blue",shape = 1) +
    geom_jitter(data=data.dup[data.dup$mut.class=="removed.to.match.LOD",], aes(x = days_post_surgery, fill = clone),position=position_jitterdodge(dodge.width=0.9,jitter.width=0.4),colour="red",shape = 1) +
    geom_jitter(data=data.dup[data.dup$mut.class=="Library_Qual_filtered",], aes(x = days_post_surgery, fill = clone),position=position_jitterdodge(dodge.width=0.9,jitter.width=0.4),colour="green",shape = 1) +
    geom_point(aes(x = days_post_surgery, y=ctDNA_fraction, fill = clone),colour="red",shape = 1,position=position_dodge(width = 0.9),size=3) +
    #scale_fill_manual(values = rep("grey",length(unique(modeleddists$clone))))+
    scale_y_continuous(trans = "log10",breaks =10^(-6:0),labels = c(0,10^(-5:0)),limits = c(10^-6,10^0)) +
    ggtitle(pat)+
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5))
  
  
}

#=======#
#  END  #
#=======#