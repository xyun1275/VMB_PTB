setwd('/Desktop/1_KH/VMB_PTB/')

##### functions #####
library(readxl)
library(ggplot2)
library(stringr)
library(parallel)
library(dplyr)
library(devtools)

library(vegan) # for diversity
library(stringr) 
library(ggplot2)
library(ALDEx2) # for diffential abundance
library(GUniFrac) # for 'Rarefy'
library(tidyr) # for 'gather'
library(reshape2) # for 'melt' in bar plot
library(parallel)
library(compositions) # for 'clr' function

get_abundance_table = function(reads_table) {
  reads_table_abundance = sweep(reads_table,2,colSums(reads_table),"/")
  return(reads_table_abundance)
}

prepare_reads_table = function(reads_table, metadata, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 1) {
  # total reads threshold
  keep <- colSums(reads_table) >= total_reads_threshold 
  reads_table <- reads_table[,keep]
  metadata <- metadata[keep,]
  
  # get abundance table
  reads_table_abundance = sweep(reads_table,2,colSums(reads_table),"/")
  
  # species threshold
  keep <- rep(T, nrow(reads_table_abundance))
  
  trials = c(1: nrow(reads_table_abundance))
  func_1 = function(trial) {
    c = sum(reads_table_abundance[trial,] >= species_threshold) / ncol(reads_table_abundance) >= 0.05        # input
    d = sum(reads_table_abundance[trial,] >= species_threshold/10) / ncol(reads_table_abundance) >= 0.15        # input
    keep = c|d
    return(keep)
  }
  keep = mclapply(trials, func_1, mc.cores = mc.cores)
  keep = unlist(keep)
  
  reads_table <- reads_table[keep,]
  
  # total reads threshold
  keep <- colSums(reads_table) >= total_reads_threshold 
  reads_table <- reads_table[,keep]
  metadata <- metadata[keep,]
  
  return(c(list(reads_table = reads_table),list(metadata = metadata)))
}

prepare_reads_table_2 = function(reads_table, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 1) {
  # total reads threshold
  keep <- colSums(reads_table) >= total_reads_threshold 
  reads_table <- reads_table[,keep]
  
  # get abundance table
  reads_table_abundance = sweep(reads_table,2,colSums(reads_table),"/")
  
  # species threshold
  keep <- rep(T, nrow(reads_table_abundance))
  
  trials = c(1: nrow(reads_table_abundance))
  func_1 = function(trial) {
    c = sum(reads_table_abundance[trial,] >= species_threshold) / ncol(reads_table_abundance) >= 0.05        # input
    d = sum(reads_table_abundance[trial,] >= species_threshold/10) / ncol(reads_table_abundance) >= 0.15        # input
    keep = c|d
    return(keep)
  }
  keep = mclapply(trials, func_1, mc.cores = mc.cores)
  keep = unlist(keep)
  
  reads_table <- reads_table[keep,]
  
  # total reads threshold
  keep <- colSums(reads_table) >= total_reads_threshold 
  reads_table <- reads_table[,keep]
  
  return(reads_table)
}

alpha_diversity = function(reads_table, metadata = NA, factor_name = NA, paired = F, order = NA, rarefy_to = NA) {
  
  if (is.na(metadata)[1]) {
    print('no metadata')
    return(NA)
  }
  
  if (is.na(factor_name)) {
    factor_name = ''
  }
  
  # rarefy to normalize data
  reads_table = t(reads_table)
  
  if (is.na(rarefy_to)) {
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  } else {
    reads_table = Rarefy(reads_table, depth = rarefy_to)
  }
  
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  # calculate diversity
  alpha.shannon_diversity <- data.frame(diversity(reads_table))
  alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table))
  alpha.ovserved_OTU <- data.frame(colSums(t(reads_table) != 0))
  
  alpha = as.data.frame(matrix(data = NA,ncol=3,nrow = nrow(reads_table)))
  colnames(alpha) = c('alpha.shannon','alpha.evenness','alpha.ovserved_OTU')
  
  alpha$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table.
  alpha$alpha.evenness <- alpha.evenness$diversity.reads_table.
  alpha$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table.....0.
  
  metadata = cbind(metadata, alpha)
  
  colnames(metadata)[1] = 'factor'
  
  metadata = metadata[order(metadata$factor),]
  
  if (is.na(order)[1] ) {
    metadata$factor <- factor(metadata$factor , levels = unique(metadata$factor))
  } else {
    metadata$factor <- factor(metadata$factor , levels = order)
  }
  
  alpha.shannon = ggplot(metadata, aes(x=factor, y=alpha.shannon)) + geom_violin(trim=T, aes(fill=factor), lwd = 0.5)+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1, alpha= 0.1)+ theme_bw()+
    labs(x = NULL, y = "Shannon index", fill=factor_name)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  
  
  alpha.evenness = ggplot(metadata, aes(x=factor, y=alpha.evenness)) + geom_violin(trim=T, aes(fill=factor), lwd = 0.5)+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1, alpha= 0.1)+theme_bw()+
    labs(x = NULL, y = "Evenness", fill=factor_name)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  
  alpha.ovserved_OTU = ggplot(metadata, aes(x=factor, y=alpha.ovserved_OTU)) +geom_violin(trim=T, aes(fill=factor), lwd = 0.5)+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1, alpha= 0.1)+theme_bw()+
    labs(x = NULL, y = "Number of observed taxa", fill=factor_name)+
    geom_jitter(size = 0.1, alpha= 0.1)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  
  # calculate significance
  factor_levels = unique(metadata$factor)
  n = length(factor_levels)
  
  Shannon_sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
  colnames(Shannon_sig) = factor_levels
  row.names(Shannon_sig) = factor_levels
  
  Evenness_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(Evenness_sig) = factor_levels
  row.names(Evenness_sig) = factor_levels
  
  OTU_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(OTU_sig) = factor_levels
  row.names(OTU_sig) = factor_levels
  
  for (a in 1:(n-1)) {
    for (b in (a+1) : n) {
      factor_level1 <- subset(metadata,  factor == factor_levels[a],
                              drop = TRUE)
      factor_level2 <- subset(metadata,  factor == factor_levels[b],
                              drop = TRUE)
      
      Shannon_sig[a,b] <- wilcox.test(factor_level1$alpha.shannon, 
                                      factor_level2$alpha.shannon, paired = paired)$p.value
      Evenness_sig[a,b] <- wilcox.test(factor_level1$alpha.evenness, 
                                       factor_level2$alpha.evenness, paired = paired)$p.value
      OTU_sig[a,b] <- wilcox.test(factor_level1$alpha.ovserved_OTU, 
                                  factor_level2$alpha.ovserved_OTU, paired = paired)$p.value
      
    }
  }
  output = c(list(alpha = alpha), list(shannon = alpha.shannon), 
             list(evenness =alpha.evenness) , list(ovserved_OTU =alpha.ovserved_OTU),
             list(sig_Shannon = Shannon_sig),list(sig_Evenness = Evenness_sig),
             list(sig_OTU = OTU_sig))
  
  return(output)
}

beta_diversity = function(reads_table, metadata = NA, factor_name = NA, 
                          order = NA, NMDS_skip = T, ref_group = NA, 
                          rarefy_to = NA, pheatmap_fontsize = 50,
                          treeheight = 10, pheatmap_y = T) {
  
  if (is.na(metadata)[1]) {
    print('no metadata')
    return(NA)
  }
  
  if (is.na(factor_name)[1]) {
    factor_name = ''
  }
  
  # rarefy to normalize data
  reads_table = as.data.frame(t(reads_table))
  
  if (is.na(rarefy_to)) {
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  } else {
    reads_table = Rarefy(reads_table, depth = rarefy_to)
  }
  
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  metadata=as.matrix(metadata)
  
  # Bray_Curtis
  Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
  Bray_Curtis <- as.data.frame(Bray_Curtis)
  
  Bray_Curtis_2 = Bray_Curtis
  #    Bray_Curtis_2[row(Bray_Curtis_2) <= col(Bray_Curtis_2)] =NA
  
  # within sample distance
  group_dis = gather(Bray_Curtis_2)
  group_dis$key2 = rep(row.names(Bray_Curtis_2),ncol(Bray_Curtis_2))
  
  Source = matrix(data = NA, ncol = length(metadata), nrow = length(metadata))
  
  for (a in 1:length(metadata)) {
    Source[a,] = metadata
  }
  Source = gather(as.data.frame(Source))
  group_dis$Source = Source$value
  group_dis$Target = rep(metadata,length(metadata))
  
  group_dis = group_dis[!is.na(group_dis$value),]
  group_dis$Source <- as.factor(group_dis$Source)
  #  group_dis$value = as.numeric(as.character(group_dis$value))
  
  keep = group_dis$Source == group_dis$Target
  within_dis = group_dis[keep,]
  keep = within_dis$key != within_dis$key2
  within_dis = within_dis[keep,]
  #  within_dis$value = as.numeric(as.character(within_dis$value))
  
  if (!is.na(order)) {
    within_dis$Source = as.factor(within_dis$Source)
    within_dis$Source = factor(within_dis$Source, levels= order)
  }
  
  within_dis_p = ggplot(within_dis, aes(x=Source, y=value)) + geom_violin(trim=T, aes(fill=Source))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1, alpha= 0.1)+ theme_bw()+
    labs(x = NULL, y = "Within sample distance", fill=factor_name)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  within_dis_p
  ggsave("within_group_distance.pdf", width=3, height=2)
  
  
  # significance among within sample distance, Wilcoxon test
  group_level = unique(within_dis$Source)
  n = length(group_level)
  
  within_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
  colnames(within_dis_sig) = group_level
  row.names(within_dis_sig) = group_level
  for (a in 1:(n-1)) {
    for (b in (a+1): n) {
      keep = metadata == group_level[a] | metadata == group_level[b]
      metadata_2 = metadata[keep]
      reads_table_2 = reads_table[keep,]
      data = mrpp(reads_table_2, as.matrix(metadata_2), permutations = 999, distance = "bray",
                  weight.type = 1, strata = NULL, parallel = getOption("mc.cores"))
      within_dis_sig[a,b] <- data$Pvalue
      
    }
  }
  
  
  write.csv(within_dis_sig,'within_group_distance.csv', row.names = T, quote = F)
  
  # distance among groups
  {
    group_level = unique(group_dis$Source)
    n = length(group_level)
    distance_median = matrix(data=NA, nrow = n, ncol =n)
    colnames(distance_median) = group_level
    row.names(distance_median) = group_level
    
    group_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
    colnames(group_dis_sig) = group_level
    row.names(group_dis_sig) = group_level
    
    group_betadisper_sig <- matrix(data = NA, ncol=n, nrow = n)
    colnames(group_betadisper_sig) = group_level
    row.names(group_betadisper_sig) = group_level
    
    for (a in 1:n) {
      for (b in a:n) {
        distance_data = group_dis$value[group_dis$Source == row.names(distance_median)[a] & group_dis$Target == colnames(distance_median)[b]]
        distance_median[a,b] <- median(distance_data)
        distance_median[b,a] <- distance_median[a,b]
        
        if (a != b) {
          keep = metadata == group_level[a] | metadata == group_level[b]
          metadata_2 = as.character(metadata[keep])
          reads_table_2 = reads_table[keep,]
          
          metadata_2 = as.data.frame(metadata_2)
          pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")[1,5] 
          group_dis_sig[b,a] = pvalue
          group_dis_sig[a,b] = pvalue
          
          d = vegdist(reads_table_2)
          perm.eg.betadisper <- betadisper(d, group = as.factor(metadata_2$metadata_2), type = "centroid")
          pvalue <- adonis2(dist(perm.eg.betadisper$distances) ~ as.factor(metadata_2$metadata_2))[1,5] 
          
          group_betadisper_sig[b,a] = pvalue
          group_betadisper_sig[a,b] = pvalue
        }
        
      }
    }
    write.csv(group_dis_sig,'between_group_distance.csv', row.names = T, quote = F)
    
    if (pheatmap_y == F & n > 2) {
      group_dis_2_p = corrplot(distance_median,p.mat = group_dis_sig, method = 'shade', diag = F, type="upper",
                               sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
                               is.corr=FALSE, col=colorRampPalette(c("white","red"))(100)[c(10:100)],
                               order="original",tl.col = "black")
      
      pdf("between_group_distance_2.pdf")
      corrplot(distance_median,p.mat = group_dis_sig, method = 'shade', diag = F, type="upper",
               sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
               is.corr=FALSE, col=colorRampPalette(c("white","red"))(100)[c(10:100)],
               order="original",tl.col = "black")
      dev.off()
      
    } else {
      group_dis_sig_2 = group_dis_sig
      group_dis_sig_2[is.na(group_dis_sig)] = ''
      group_dis_sig_2[group_dis_sig > 0.05] = ''
      group_dis_sig_2[group_dis_sig <= 0.05 & group_dis_sig > 0.01] = '*'
      group_dis_sig_2[group_dis_sig <= 0.01 & group_dis_sig > 0.001] = '**'
      group_dis_sig_2[group_dis_sig <= 0.001] = '***'
      
      save_heatmap_pdf <- function(x, filename, width=8, height=8) {
        stopifnot(!missing(x))
        stopifnot(!missing(filename))
        pdf(filename, width=width, height=height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
      }
      
      group_dis_2_p = pheatmap(distance_median, cluster_rows=TRUE, show_rownames=TRUE, 
                               cluster_cols=T, show_colnames=T, 
                               color=colorRampPalette(c("white","red"))(100),
                               fontsize = pheatmap_fontsize, display_numbers = group_dis_sig_2, 
                               treeheight_row = treeheight, treeheight_col = treeheight)
      save_heatmap_pdf(group_dis_2_p, "between_group_distance_2.pdf", width=2, height=2)
    }
    
  }
  
  
  
  # Running Nonmetric Multidimensional Scaling (NMDS) Ordination
  if (NMDS_skip == T) {
    # output
    output = c(Bray_Curtis = list(Bray_Curtis), within_dis_p = list(within_dis_p),
               within_dis_sig = list(within_dis_sig), group_dis_2_p = list(group_dis_2_p),
               group_dis_sig = list(group_dis_sig))
    return(output)
    
  } else {
    
    colnames(metadata)[1] = 'factor'
    
    NMDS <-
      metaMDS(Bray_Curtis,
              distance = "bray",
              k = 2,
              maxit = 999, 
              trymax = 20,
              wascores = TRUE)
    
    mds_data <- as.data.frame(NMDS$points)
    mds_data$factor <- metadata
    
    if (!is.na(order)[1]) {
      mds_data$factor = as.factor(mds_data$factor)
      mds_data$factor = factor(mds_data$factor, levels= order)
    }
    
    NMDS = ggplot(mds_data, aes(x = MDS1, y = MDS2, color = factor)) +
      geom_point(size = 0.3)+
      scale_colour_discrete(factor_name)+
      labs(x = 'NMDS1', y = "NMDS2")+ 
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    NMDS+theme_bw()
    ggsave("NMDS.pdf", width=4, height=3)
    
    NMDS_2 = ggplot(mds_data, aes(x = MDS1, y = MDS2, color = factor)) +
      geom_point(size = 0.3)+
      scale_colour_discrete(factor_name)+
      stat_ellipse(type = "t")+
      labs(x = 'NMDS1', y = "NMDS2")+ 
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    NMDS_2+theme_bw()
    ggsave("NMDS_2.pdf", width=4, height=3)
    
    # output
    output = c(Bray_Curtis = list(Bray_Curtis), within_dis_p = list(within_dis_p), 
               within_dis_sig = list(within_dis_sig),group_dis_2_p = list(group_dis_2_p),
               group_dis_sig = list(group_dis_sig), group_betadisper_sig = list(group_betadisper_sig),
               NMDS =list(NMDS), NMDS_2 =list(NMDS_2))
    return(output)
  }
  
}

dif_abundance <- function(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, paired_test = F, order_reverse = F, style = 1, order = NA) {
  #   style =1
  #   reads_table= reads_table
  #   metadata= metadata$Flag
  #   paired_test = F
  #   order_reverse = F
  #   style =2 
  #   order = c('TB','PTB')
  #   fold_change_th = 1
  #   pvalue_th = 0.05
  
  if (!is.na(order)[1]) {
    metadata = factor(metadata, levels = order)
  }
  
  conds <- metadata
  
  x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
  
  # paired Wilcoxon Rank Sum test and Welch's t-test
  x.tt <- aldex.ttest(x, paired.test= paired_test)
  
  x.effect <- aldex.effect(x)
  
  x.all <- data.frame(cbind(x.tt,x.effect))
  if (sum(is.nan(x.all$wi.eBH)) > 1) {
    x.all$wi.eBH = x.all$we.eBH
    x.all$we.eBH = 'we.eBH used in this test'
  }
  
  x.all$`-Log10(adj-pvalue)` <- -log10(x.all$wi.eBH)
  x.all$Taxa = row.names(x.all)
  
  if (max(abs(x.all$diff.btw)) < fold_change_th | min(x.tt$wi.eBH) > 0.05) {
    print('No taxon has significant abundance change')
    return(list(data = x.all))
    
  }
  
  if (order_reverse == T) {
    x.all$diff.btw = -x.all$diff.btw
  }
  
  # draw figure
  das <- x.all[(x.all$`-Log10(adj-pvalue)` >= -log10(pvalue_th) & (x.all$diff.btw >=fold_change_th | x.all$diff.btw <=-fold_change_th)),]
  
  if (nrow(das)==0) {
    print('No taxon has significant abundance change')
    return(list(data = x.all))
  }
  
  das$Species <- row.names(das)
  das <- das[order(das$diff.btw),] 
  
  metadata = as.factor(metadata)
  lev = levels(metadata)
  
  if (order_reverse == T) {
    das$Color <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
  } else {
    das$Color <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
    
  }
  
  das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
  
  das$diff.btw[das$diff.btw == Inf] = 10
  das$diff.btw[das$diff.btw == -Inf] = -10
  das$diff.btw[das$diff.btw <= -10] = -10
  das$diff.btw[das$diff.btw >= 10] = 10
  
  
  if (style == 1) {
    theme_set(theme_bw())  
    
    p <- ggplot(das, aes(Species, diff.btw)) + 
      geom_point(aes(col=Color, size=`-Log10(adj-pvalue)`)) + 
      coord_flip() +          # convert x y axis
      labs(x = 'Taxa', y = "Median difference in clr values")+ 
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))  
    
    
  } else {
    taxa_list = row.names(das)
    
    keep = which(row.names(reads_table) %in% taxa_list)
    
    reads_table_abundance = sweep(reads_table,2,colSums(reads_table),"/")
    reads_table_2 <- reads_table_abundance[keep,]
    reads_table_2 <- as.data.frame(t(reads_table_2))
    
    #   colnames(reads_table_2)=taxa_list
    reads_table_3 = gather(reads_table_2)
    reads_table_3$Type = rep(conds, length(taxa_list))
    colnames(reads_table_3) = c('Taxa', 'Abundance','Type')
    
    reads_table_3$Taxa = str_replace_all(reads_table_3$Taxa,'.*__','')
    
    if (!is.na(order)[1]) {
      reads_table_3$Type = factor(reads_table_3$Type, levels = order)
    }
    
    p <- ggplot(reads_table_3, aes(x=Taxa, y=Abundance,fill=Type)) + geom_violin(trim=T, aes(fill=Type))+
      geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
      geom_jitter(size = 0.1, alpha= 0.1)+theme_bw()+
      ylab("Abundance (%)")+
      theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))
  }
  
  return(c(list(p = p), list(data = x.all)))
  
}




##### get reads table and metadata #####
setwd('/Desktop/1_KH/VMB_PTB/data')
{
  reads_table_all = read.csv('reads_table.csv', row.names = 1)

  # metadata cleansing
  {
    metadata_all = read.csv('/Desktop/1_KH/VMB_PTB/data/metadata_all.csv', row.names = 1)

    metadata_all[metadata_all ==''] = NA; metadata_all[metadata_all =="Not Reported"] = NA
    metadata_all[metadata_all =="not applicable"] = NA; metadata_all[metadata_all =="Decline"] = NA; metadata_all[metadata_all =="Other"] = NA

    metadata_all <- metadata_all %>% mutate(env_biome = recode(env_biome, 
                                                                  `Cervicovaginal swab` = 'vagina',
                                                                  `human vaginal metagenome` = 'vagina',
                                                                  `laboratory` = 'environment',
                                                                  `vaginal metagenome` = 'vagina',
                                                                  `Vaginal mucose` = 'vagina',
                                                                  `vaginal swab eluate` = 'vagina'))
    metadata_all$env_biome[str_detect(metadata_all$env_biome,'kit')] = 'environment'
    
    metadata_all$NugentScore = as.numeric(as.character(metadata_all$NugentScore))
    metadata_all$pH = as.numeric(as.character(metadata_all$pH))
    
    metadata_all <- metadata_all %>% mutate(Outcome = recode(Outcome, 
                                                               `MiPTB` = 'PTB',
                                                               `preterm` = 'PTB',
                                                               `Preterm` = 'PTB',
                                                               `PRETERM` = 'PTB',
                                                               `SpPTB` = 'PTB', 
                                                             `SPTB` = 'PTB', 
                                                             `term` = 'Term',
                                                             `TERM` = 'Term',
                                                             `FetalLoss` = 'Fetal Loss'))
    
    metadata_all <- metadata_all %>% mutate(Race = recode(Race, 
                                                             `AmericanIndian` = 'American Indian',
                                                             `Black` = 'B/AA',
                                                             `white` = 'White'))
    metadata_all$seq_methods[str_detect(metadata_all$seq_methods,'Roche 454')] = 'Roche 454'
    metadata_all <- metadata_all %>% mutate(seq_methods = recode(seq_methods, `Illumina HiSeq 2500` = 'Illumina HiSeq'))
    
    table(metadata_all$Outcome)
  }

  metadata_all$Run = row.names(metadata_all)
  metadata_all = metadata_all[metadata_all$Run %in% colnames(reads_table_all),]
  metadata_all[metadata_all == 'NA'] = NA
  reads_table_all = reads_table_all[,colnames(reads_table_all) %in% metadata_all$Run]
  reads_table_all = reads_table_all[metadata_all$Run]
  metadata_all$gest_day_collection = as.numeric(as.character(metadata_all$gest_day_collection))
  metadata_all$gest_day_delivery = as.numeric(as.character(metadata_all$gest_day_delivery))
  metadata_all$Age = as.numeric(as.character(metadata_all$Age))
  
  keep = metadata_all$env_biome == 'vagina' & metadata_all$Outcome %in% c('PTB', 'Term')
  metadata_all=metadata_all[keep,]
  reads_table_all = reads_table_all[,keep]
  
  x = data.frame(total_reads = colSums(reads_table_all))
  x$total_reads = log10(x$total_reads)
  ggplot(x, aes(x=total_reads)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
    xlab('Sample total reads (Log10)') + theme_bw()
  ggsave('/Desktop/1_KH/VMB_PTB/results/overall_profiles/Sample total reads_ori.pdf', width = 3, height = 2)
  
  # colors
  {
    library(RColorBrewer)
    x <- brewer.pal(n = 12, name = "Set3")
    x = x[c(1:11)]
    y <- brewer.pal(n = 9, name = "Set1")
    y = y[c(1,2,3,4,7)]
    brewer_colors = c(x,y)
    
    mycolors = list()
    for (a in 1:ncol(metadata_all)) {
      data <- brewer_colors[c(1:length(unique(metadata_all[,a])))]
      names(data) <- sort(unique(metadata_all[,a]))
      mycolors[[a]] = data
    }
    names(mycolors) = colnames(metadata_all)
  }
  
  # metadata profiles
  {
    x = colnames(metadata_all)
    for (a in 1:length(x)) {
      if (length(unique(metadata_all[,a])) >8) {
        if (sum(is.numeric(metadata_all[,a]) & !is.na(metadata_all[,a])) < 20) {
          next
        } else {
          y = as.data.frame(metadata_all[,a])
          colnames(y) = 'aaa'
          ggplot(y, aes(x= aaa)) +
            geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
            xlab(x[a]) + theme_bw() + labs(title = paste0('n = ', sum(!is.na(y$aaa))))
          ggsave(paste0('/Desktop/1_KH/VMB_PTB/results/Metadata_profile/Metadata_profile_', x[a],'.pdf'), width = 3, height = 2)
        }
      } else {
        data <- as.data.frame(table(metadata_all[,a]))
        data = data[order(data$Var1),]
        ggplot(data, aes(x="", y=Freq, fill=Var1)) +
          geom_bar(stat="identity", width=1, color="white") + 
          coord_polar("y", start=0) +
          theme_void() + 
          scale_fill_manual(values= mycolors[[a]], name = paste0('Total number = ', sum(data$Freq), '\n',x[a])) 
        ggsave(paste0('/Desktop/1_KH/VMB_PTB/results/Metadata_profile/Metadata_profile_', x[a],'.pdf'), width = 6, height = 3.5)
      }
    }
  }
}
rm(x,a, taxonomy, data,data_2, file_list,keep,n,taxa_list,reads_table,y)




##### reads and metadata profiles #################
setwd('/Desktop/1_KH/VMB_PTB/results/overall_profiles')
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata_all = metadata_all[keep,]
  table(metadata_all$BioProject)
  reads_table_all = reads_table_all[,keep]
  
  
  # alpha rarefaction
  {
    for (b in 1:30) {
      reads_table_2 = reads_table_all
      reads_table_2 = Rarefy(reads_table_2, depth = b*200)
      reads_table_2 <- reads_table_2$otu.tab.rff
      reads_table_2 <- as.data.frame(reads_table_2)
      
      alpha.ovserved_OTU <- data.frame(colSums(t(reads_table_2) != 0))
      alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table_2.....0.
      
      if (b == 1) {
        output_2 = data.frame(Depth = b*200, taxa = alpha.ovserved_OTU)
        output = output_2
      } else {
        output_2 = data.frame(Depth = b*200, taxa = alpha.ovserved_OTU)
        output = rbind(output, output_2)
      }
    }
    
    wilcox.test(as.numeric(as.character(output$taxa[output$Depth == 1600])), 
                as.numeric(as.character(output$taxa[output$Depth == 1800])))
    
    wilcox.test(as.numeric(as.character(output$taxa[output$Depth == 1800])), 
                as.numeric(as.character(output$taxa[output$Depth == 2000])))
    
    wilcox.test(as.numeric(as.character(output$taxa[output$Depth == 2000])), 
                as.numeric(as.character(output$taxa[output$Depth == 2200])))
    
    output$Depth = as.factor(output$Depth)
    ggplot(output, aes(x=Depth, y=taxa)) + 
      geom_boxplot(outlier.shape=NA, width=0.5, fill = '#69b3a2') + theme_bw()+
      labs( y = "Number of observed taxa")+
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave('alpha_rarefaction.pdf', width = 5, height = 2)
  }
  keep = colSums(reads_table_all) >= 5000; sum(keep)
  reads_table = reads_table_all[,keep]
  metadata = metadata_all[keep,]
  output = prepare_reads_table(reads_table, metadata, total_reads_threshold = 5000, species_threshold = 0.0001, mc.cores = 8)
  reads_table = output$reads_table; metadata = output$metadata
  reads_table_all = reads_table;  metadata_all = metadata
  
  x = data.frame(total_reads = colSums(reads_table_all))
  x$total_reads = log10(x$total_reads)
  ggplot(x, aes(x=total_reads)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
    xlab('Sample total reads (Log10)') + theme_bw()
  ggsave('/Desktop/1_KH/VMB_PTB/results/overall_profiles/Sample total reads.pdf', width = 3, height = 2)
  # colors
  {
    library(RColorBrewer)
    x <- brewer.pal(n = 12, name = "Set3")
    x = x[c(1:11)]
    y <- brewer.pal(n = 9, name = "Set1")
    y = y[c(1,2,3,4,7)]
    brewer_colors = c(x,y)
    
    mycolors = list()
    for (a in 1:ncol(metadata_all)) {
      data <- brewer_colors[c(1:length(unique(metadata_all[,a])))]
      names(data) <- sort(unique(metadata_all[,a]))
      mycolors[[a]] = data
    }
    names(mycolors) = colnames(metadata_all)
  }
  # metadata profiles
  {
    x = colnames(metadata_all)
    for (a in 1:length(x)) {
      if (length(unique(metadata_all[,a])) <=1) {next}
      
      if (length(unique(metadata_all[,a])) >8) {
        if (sum(is.numeric(metadata_all[,a]) & !is.na(metadata_all[,a])) < 8) {
          next
        } else {
          y = as.data.frame(metadata_all[,a])
          colnames(y) = 'aaa'
          ggplot(y, aes(x= aaa)) +
            geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
            xlab(x[a]) + theme_bw() + labs(title = paste0('n = ', sum(!is.na(y$aaa))))
          ggsave(paste0('/Desktop/1_KH/VMB_PTB/results/Metadata_profile/Metadata_profile_', x[a],'.pdf'), width = 3, height = 2)
        }
      } else {
        data <- as.data.frame(table(metadata_all[,a]))
        data = data[order(data$Var1),]
        ggplot(data, aes(x="", y=Freq, fill=Var1)) +
          geom_bar(stat="identity", width=1, color="white") + 
          coord_polar("y", start=0) +
          theme_void() + 
          scale_fill_manual(values= mycolors[[a]], name = paste0('Total number = ', sum(data$Freq), '\n',x[a])) 
        ggsave(paste0('/Desktop/1_KH/VMB_PTB/results/Metadata_profile/Metadata_profile_', x[a],'.pdf'), width = 6, height = 3.5)
      }
    }
  }
  rm(data, metadata, output, reads_table,y, a, keep, x, output_2, reads_table_2, alpha.ovserved_OTU, b)
}






##### batch effect removal #####
{
  library(sva) 
  # 1 normalization by clr
  reads_table = as.matrix(decostand(as.matrix(t(reads_table_all)),method = 'rclr'))
  reads_table = as.data.frame(t(reads_table))
  
  # 2 batch removal by ComBat
  otu_corrected <- ComBat(dat = reads_table, 
                          batch = as.factor(metadata_all$BioProject),
                          mod = model.matrix(~ as.factor(metadata_all$Outcome), data=reads_table))
  otu_corrected = as.data.frame(otu_corrected)
  
  # 3 convert clr values back to reads table
  x = exp(t(otu_corrected))
  x = x/rowSums(x)
  x = x*10000
  x = round(x) 
  x = as.data.frame(t(x))
  reads_table_all_batch_rm = x

  # 4 compare by t-SNE
  #tsne <- Rtsne(as.matrix(t(reads_table_all)), dims = 2, perplexity=20, verbose=TRUE, max_iter = 2000)
  x = prcomp(as.data.frame((reads_table_all)), center = TRUE, scale. = TRUE)
  pic <- x$rotation[,c(1,2)]
  #pic <- tsne$Y
  pic <- data.frame(pic,metadata_all)
  colnames(pic)[1:2] <- c('X1','X2')
  pic$Outcome
  ggplot(pic, aes(X1, X2,color = BioProject))  +
    geom_point(size=0.5) + 
    scale_color_manual(values= mycolors[[1]]) +
    xlab(paste0("PCoA1")) +
    ylab(paste0("PCoA2"))+
    stat_ellipse(type = "t") + 
    coord_fixed()+ 
    theme(
      axis.title.x = element_text( size=7),
      axis.title.y = element_text( size=7),
      legend.text = element_text(size=7),
      legend.title = element_text(size=7),
      plot.title = element_text(hjust = 0.5, size = 7)
    ) + theme_bw()
  ggsave('tsne_ori_BioProject.pdf', width=5, height=5)
  
  #tsne <- Rtsne(as.matrix(t(reads_table_all_batch_rm)), dims = 2, perplexity=20, verbose=TRUE, max_iter = 2000)
  #pic <- tsne$Y
  
  x = prcomp(as.data.frame((reads_table_all_batch_rm)), center = TRUE, scale. = TRUE)
  pic <- x$rotation[,c(1,2)]
  pic <- data.frame(pic,metadata_all)
  colnames(pic)[1:2] <- c('X1','X2')
  pic$Outcome
  ggplot(pic, aes(X1, X2,color = BioProject))  +
    geom_point(size=0.5) + 
    scale_color_manual(values= mycolors[[1]]) +
    xlab(paste0("t-SNE1")) +
    ylab(paste0("t-SNE2"))+
    stat_ellipse(type = "t") + 
    coord_fixed()+ 
    theme(
      axis.title.x = element_text( size=7),
      axis.title.y = element_text( size=7),
      legend.text = element_text(size=7),
      legend.title = element_text(size=7),
      plot.title = element_text(hjust = 0.5, size = 7)
    ) + theme_bw()
  ggsave('tsne_BioProject.pdf', width=5, height=5)
}
rm(adjusted_reads, metadata, otu_corrected, pca_results,pic, reads_table, x, tsne)
##### overall profiles #####
setwd('/Desktop/1_KH/VMB_PTB/results/overall_profiles/')

# Differential abundance
{
  metadata = metadata_all$Outcome
  reads_table = as.data.frame(t(reads_table_all_batch_rm))
  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(t(reads_table))
  
  pvalue_th = 0.05; fold_change_th = 0.5; paired_test = F; order_reverse = F; style = 1; order = NA
  {
    if (!is.na(order)[1]) {
      metadata = factor(metadata, levels = order)
    }
    
    conds <- metadata
    
    x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
    
    # paired Wilcoxon Rank Sum test and Welch's t-test
    x.tt <- aldex.ttest(x, paired.test= paired_test)
    
    x.effect <- aldex.effect(x)
    
    x.all <- data.frame(cbind(x.tt,x.effect))
    if (sum(is.nan(x.all$wi.eBH)) > 1) {
      x.all$wi.eBH = x.all$we.eBH
      x.all$we.eBH = 'we.eBH used in this test'
    }
    
    x.all$`-Log10(FDR)` <- -log10(x.all$wi.eBH)
    x.all$Taxa = row.names(x.all)
    
    write.csv(x.all,'Differential_abundance.csv')
    
    # draw figure
    das <- x.all[x.all$`-Log10(FDR)` >= -log10(pvalue_th),]
    
    das$Species <- row.names(das)
    das <- das[order(das$diff.btw),] 
    
    metadata = as.factor(metadata)
    lev = levels(metadata)
    
    if (order_reverse == T) {
      das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
    } else {
      das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
      
    }
    das$Species = str_replace_all(das$Species, '_',' ')
    das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
    
    das$diff.btw[das$diff.btw == Inf] = 10
    das$diff.btw[das$diff.btw == -Inf] = -10
    das$diff.btw[das$diff.btw <= -10] = -10
    das$diff.btw[das$diff.btw >= 10] = 10
    
    
    if (style == 1) {
      theme_set(theme_bw())  
      
      
      ggplot(das, aes(Species, diff.btw)) +
        geom_line() +  
        geom_hline(yintercept=1, linetype='dotted', col = 'grey') +
        geom_hline(yintercept=-1, linetype='dotted', col = 'grey') +
        geom_point(aes(col=Outcome, size=`-Log10(FDR)`)) +
        scale_size_continuous(range  = c(0.1, 5), limits = c(0, 25), 
                              breaks = c(1, 6, 11, 16, 21))+ 
        scale_color_brewer(palette="Set1")+
        coord_flip() +          # convert x y axis
        labs(x = 'Taxa', y = "Median difference in clr values")+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))  
      ggsave('Differential_abundance.pdf', width=4, height=2.7)
    } 
    
  }
}

# PCoA
{
  Padonis = adonis2(as.matrix(t(reads_table_all_batch_rm)) ~ Outcome, data = metadata_all,method = "bray", parallel = 8)
  Padonis$`Pr(>F)`
  
  #tsne <- Rtsne(as.matrix(t(reads_table_all_batch_rm)), dims = 2, perplexity=20, verbose=TRUE, max_iter = 5000)
  #x <- vegdist(as.matrix(t(reads_table_all_batch_rm)), method = "bray")
  #x <- cmdscale(bray_curtis_dist, eig = TRUE, k = 3)  # k = number of dimensions
  #x <- as.data.frame(x$points)

  #reads_table = as.data.frame(t(reads_table_all_batch_rm))
  #reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  #reads_table <- reads_table$otu.tab.rff
  #reads_table <- as.data.frame(reads_table)
  #dbRDA = capscale(reads_table ~ Outcome, data = metadata_all, distance = "bray",parallel = 8)

  reads_table = as.matrix(decostand(as.data.frame(t(reads_table_all_batch_rm)),method = 'rclr'))
  x = prcomp(as.data.frame(t(reads_table)), center = TRUE, scale. = TRUE)
  
  pic <- x$rotation
  pic <- data.frame(pic,metadata_all)
  
  colnames(pic)[1:2] <- c('X1','X2')
  pic$Outcome
  ggplot(pic, aes(X1, X2,color = Outcome))  +
    geom_point(size=0.5) +
    xlab(paste0("PCoA1")) +
    ylab(paste0("PCoA2"))+
    stat_ellipse(type = "t") + 
    coord_fixed()+ 
    scale_color_brewer(palette="Set1")+
    theme(
      axis.title.x = element_text( size=7),
      axis.title.y = element_text( size=7),
      legend.text = element_text(size=7),
      legend.title = element_text(size=7),
      plot.title = element_text(hjust = 0.5, size = 7)
    ) 
  ggsave('PCoA_Outcome.pdf', width=4, height=3)
}

# alpha
{
  reads_table = reads_table_all_batch_rm
  metadata = metadata_all$Outcome
  
  # alpha diversity 
  reads_table = t(reads_table)
  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  alpha <- data.frame(diversity(reads_table))
  alpha$Evenness <- alpha$diversity.reads_table./log(specnumber(reads_table))
  alpha$Ovserved_phage <- colSums(t(reads_table) != 0)
  colnames(alpha) = c('Shannon','Evenness','Ovserved_taxa')
  alpha$Outcome = metadata
  x = kruskal.test(Shannon ~ Outcome,data = alpha)
  alpha$Outcome = metadata
  P = formatC(x$p.value, format = "e", digits = 3)
  
  ggplot(data=alpha, aes(x=Outcome, y=Shannon)) +geom_violin(trim=T, aes(fill=Outcome))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1, alpha= 0.05)+theme_bw()+
    theme_classic()+
    scale_fill_brewer(palette="Set1")+
    labs(y = 'Shannon index')+
    ggtitle(paste0('P = ',P)) +
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          plot.title = element_text(size = 7))
  ggsave('Shannon.pdf',width=2.5, height=2)
}

# Adonis marginal
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata = metadata_all[keep, colnames(metadata_all) %in% c('Outcome', 'Age', 'Race', 'gest_day_collection')]
  reads_table = reads_table_all_batch_rm[,keep]
  
  reads_table2 = as.data.frame(t(reads_table))
  reads_table2 = Rarefy(reads_table2, depth = min(rowSums(reads_table2)))
  reads_table2 <- reads_table2$otu.tab.rff
  reads_table2 <- as.data.frame(reads_table2)
  
  x = metadata$Age
  y = quantile(x, probs = seq(0, 1, by = 1/3))
  x <- cut(x, breaks = quantile(x, probs = seq(0, 1, by = 1/3), na.rm = TRUE), 
                     labels = c("Low", "Medium", "High"), 
                     include.lowest = TRUE)
  ggplot(metadata, aes(x= Age)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
    xlab('Age') + 
    geom_vline(xintercept = y, 
               linetype = "dashed", # Dashed lines
               size = 0.3) +          # Line thickness
    theme_bw() 
  ggsave('Age_classify.pdf',height = 2, width = 3)
  
  metadata$Age = as.character(x)
  
  x = metadata$gest_day_collection
  y = quantile(x, probs = seq(0, 1, by = 1/3))
  x <- cut(x, breaks = quantile(x, probs = seq(0, 1, by = 1/3), na.rm = TRUE), 
           labels = c("Low", "Medium", "High"), 
           include.lowest = TRUE)
  ggplot(metadata, aes(x= gest_day_collection)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
    xlab('Gestational age at sample collection') + 
    geom_vline(xintercept = y, 
               linetype = "dashed", # Dashed lines
               size = 0.3) +          # Line thickness
    theme_bw() 
  ggsave('gest_day_collection_classify.pdf',height = 2, width = 3)
  metadata$gest_day_collection = as.character(x)

  pvalue <- adonis2(reads_table2 ~ ., data = metadata,method = "bray", parallel = 8, by = 'margin')
  x = as.data.frame(pvalue)
  write.csv(x,'Adonis_marginal.csv')
  
  metadata <- as.data.frame(apply(metadata, c(1, 2), function(x) paste0(": ",x)))
  colnames(metadata)[2] = 'Gestational_age'
  
  dbRDA = capscale(reads_table2 ~ ., data = metadata, distance = "bray",parallel = 8)
  {
    ## extract % explained by the first 2 axes
    perc <- round(100*(summary(dbRDA)$cont$importance[2, 1:2]), 2)
    
    ## extract scores - these are coordinates in the RDA space
    sc_si <- scores(dbRDA, display="sites", choices=c(1,2), scaling=1)
    sc_sp <- scores(dbRDA, display="species", choices=c(1,2), scaling=1)
    sc_sp = sc_sp/2
    sc_bp <- scores(dbRDA, display="bp", choices=c(1, 2), scaling=1)
    sc_bp = sc_bp*55
    
    #      n = which(row.names(sc_sp) %in% c('Streptococcus_mitis','Streptococcus_infantis',
    #                                        'Escherichia_coli'))
    #      sc_sp[n,] = sc_sp[n,]/2
    
    rownames(sc_sp) = str_replace_all(rownames(sc_sp),'_'," ")
    
    arrow_color <- rgb(1, 0, 0, alpha = 0.5)
    {
      # Set up a blank plot with scaling, axes, and labels
      pdf("dbRDA.pdf",width = 12, height = 12)
      plot(dbRDA,
           scaling = 1, # set scaling type 
           type = "none", # this excludes the plotting of any points from the results
           frame = FALSE,
           xlim = c(-13,17), 
           ylim = c(-10,8),
           # label the plot (title, and axes)
           main = "Triplot RDA - scaling 1",
           xlab = paste0("RDA1 (", perc[1], "%)"), 
           ylab = paste0("RDA2 (", perc[2], "%)"))
      # add points for site scores
      points(sc_si, 
             pch = 21, # set shape (here, circle with a fill colour)
             col = "black", # outline colour
             bg = "steelblue", # fill colour
             cex = 1) # size
      # add points for species scores
      points(sc_sp, 
             pch = 22, # set shape (here, square with a fill colour)
             col = "black",
             bg = "#f2bd33", 
             cex = 1)
      # add text labels for species abbreviations
      text(sc_sp + c(0.2, -0.2), # adjust text coordinates to avoid overlap with points 
           labels = rownames(sc_sp), 
           col = "grey40", 
           font = 2, # bold
           cex = 0.8)
      # add arrows for effects of the expanatory variables
      arrows(0,0, # start them from (0,0)
             sc_bp[,1], sc_bp[,2], # end them at the score value
             col = arrow_color, 
             lwd = 3)
      # add text labels for arrows
      text(x = sc_bp[,1] +1.5, # adjust text coordinate to avoid overlap with arrow tip
           y = sc_bp[,2] -0.2, 
           labels = rownames(sc_bp), 
           col = "#2596be", 
           cex = 1, 
           font = 2)
      dev.off()
    }
    
  }
}

# Taxa exclusive index
{
  Predominant_taxa = vagitype(reads_table_all_batch_rm, th = 0.5)
  
  reads_table = t(reads_table_all_batch_rm)
  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  # calculate diversity
  alpha.shannon_diversity <- data.frame(diversity(reads_table))
  alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table))
  alpha.ovserved_OTU <- data.frame(colSums(t(reads_table) != 0))
  
  alpha = as.data.frame(matrix(data = NA,ncol=3,nrow = nrow(reads_table)))
  colnames(alpha) = c('alpha.shannon','alpha.evenness','alpha.ovserved_OTU')
  
  alpha$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table.
  alpha$alpha.evenness <- alpha.evenness$diversity.reads_table.
  alpha$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table.....0.
  
  data = data.frame(Predominant_taxa = as.character(Predominant_taxa), Shannno_index = alpha$alpha.shannon, Observed_taxa = alpha$alpha.ovserved_OTU, Evenness = alpha$alpha.evenness)
  data_2 = as.data.frame(table(data$Predominant_taxa))
  data = data[data$Predominant_taxa %in% data_2$Var1[data_2$Freq>10],]
  data = data[str_detect(data$Predominant_taxa,'Lactobacillus'),]
  data = data[!str_detect(data$Predominant_taxa,'iners'),]
  
  ggplot(data, aes(x=Predominant_taxa, y=Shannno_index)) + geom_violin(trim=T, aes(fill=Predominant_taxa), lwd = 0.5)+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1, alpha= 0.3)+ theme_bw()+
    labs(x = NULL, y = "Shannon index", fill=Predominant_taxa)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))+scale_x_discrete(guide = guide_axis(n.dodge = 2))
  ggsave('Predominant_taxa.pdf',width = 4, height = 2.5)

  wilcox.test(data$Shannno_index[data$Predominant_taxa =='Lactobacillus_crispatus'], 
              data$Shannno_index[data$Predominant_taxa =='Lactobacillus_gasseri'], paired = F)$p.value
  
  wilcox.test(data$Shannno_index[data$Predominant_taxa =='Lactobacillus_crispatus'], 
              data$Shannno_index[data$Predominant_taxa =='Lactobacillus_jensenii'], paired = F)$p.value
  
  wilcox.test(data$Shannno_index[data$Predominant_taxa =='Lactobacillus_gasseri'], 
              data$Shannno_index[data$Predominant_taxa =='Lactobacillus_jensenii'], paired = F)$p.value
      

}

# case match
{
  metadata_x = metadata_all; reads_table_x = reads_table_all_batch_rm
  
  # case match
  {
    {
      metadata_PTB = metadata_x[metadata_x$Outcome == "PTB",]
      metadata_Term = metadata_PTB; metadata_Term[1:nrow(metadata_Term),1:ncol(metadata_Term)] = NA
      
      if (nrow(metadata_PTB) < 6) {next}
      
      metadata_2 = metadata_x[metadata_x$Outcome == "Term",]
      
      for (a in 1: nrow(metadata_PTB)) {
        
        n =1
        keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
          metadata_2$Age <= metadata_PTB$Age[a] +n & 
          metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
          metadata_2$BioProject == metadata_PTB$BioProject[a] & 
          metadata_2$Race == metadata_PTB$Race[a] & 
          metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
          metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
        while (sum(keep[!is.na(keep)]) == 0 & n < 3) {
          n = n+1 
          keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
            metadata_2$Age <= metadata_PTB$Age[a] +n & 
            metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
            metadata_2$BioProject == metadata_PTB$BioProject[a] & 
            metadata_2$Race == metadata_PTB$Race[a] & 
            metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
            metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
        }
        
        if (n < 3) {
          metadata_Term[a,] = metadata_2[which(keep)[1],]
          metadata_2 = metadata_2[-which(keep)[1],]
        }
        
      }
      
      sum(is.na(metadata_Term$ParticipantID))
      
      keep = !is.na(metadata_Term$Run)
      metadata_Term = metadata_Term[keep,]; metadata_PTB = metadata_PTB[keep,]
      
      data = data.frame(Age = c(metadata_PTB$Age,metadata_Term$Age), 
                        Outcome = c(rep('PTB',length(metadata_PTB$Age)), rep('Term',length(metadata_Term$Age))))
      
      if (nrow(data) < 6) {next}
      
      P = kruskal.test(Age ~ Outcome,data = data)
      if (P$p.value < 0.001) {
        P = formatC(P$p.value, format = "e", digits = 3)
      } else {
        P = formatC(P$p.value, format = "f", digits = 3)
      }
      
      ggplot(data, aes(x=Outcome, y=Age)) + geom_violin(trim=T, aes(fill = Outcome))+
        geom_boxplot(fill='white', color="black", width=0.1) + 
        geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
        scale_fill_brewer(palette="Set1")+
        ggtitle(paste0(bb,'\n','P = ',P)) +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              plot.title = element_text(size = 7))
      
      ggsave(paste0("Race_",str_replace(bb,'\\/',"_"),"_",'Age.pdf'), height = 2,width = 2)
      
      data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                        Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
      
      P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
      if (P$p.value < 0.001) {
        P = formatC(P$p.value, format = "e", digits = 3)
      } else {
        P = formatC(P$p.value, format = "f", digits = 3)
      }
      ggplot(data, aes(x=Outcome, y=Gestational_age_in_sample_collection)) + geom_violin(trim=T, aes(fill = Outcome))+
        geom_boxplot(fill='white', color="black", width=0.1) + 
        geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
        labs(y = 'Gestational age at sample collection')+
        scale_fill_brewer(palette="Set1")+
        ggtitle(paste0(bb,'\n','P = ',P)) +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              plot.title = element_text(size = 7))
      ggsave(paste0("Race_",str_replace(bb,'\\/',"_"),"_",'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
      
      metadata = data.frame(ParticipantID = c(metadata_PTB$Patient_ID,metadata_Term$Patient_ID), 
                            BioProject = c(metadata_PTB$BioProject,metadata_Term$BioProject), 
                            Run = c(metadata_PTB$Run,metadata_Term$Run), 
                            gest_day_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection),
                            Age = c(metadata_PTB$Age,metadata_Term$Age), 
                            Race = c(metadata_PTB$Race,metadata_Term$Race), 
                            Outcome = c(metadata_PTB$Outcome,metadata_Term$Outcome))
      
      reads_table = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata$Run]
      reads_table = reads_table[metadata$Run]
      
      metadata_x = metadata; reads_table_x = reads_table
      rm(data, metadata,metadata_2, metadata_PTB, metadata_Term, reads_table)
    }
    table(metadata_x$Race)
    
    # Differential abundance
    {
      metadata = metadata_x$Outcome
      reads_table = as.data.frame(t(reads_table_x))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table <- reads_table$otu.tab.rff
      reads_table <- as.data.frame(t(reads_table))
      
      pvalue_th = 0.05; fold_change_th = 0.5; paired_test = F; order_reverse = F; style = 1; order = NA
      {
        if (!is.na(order)[1]) {
          metadata = factor(metadata, levels = order)
        }
        
        conds <- metadata
        
        x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
        
        # paired Wilcoxon Rank Sum test and Welch's t-test
        x.tt <- aldex.ttest(x, paired.test= paired_test)
        
        x.effect <- aldex.effect(x)
        
        x.all <- data.frame(cbind(x.tt,x.effect))
        if (sum(is.nan(x.all$wi.eBH)) > 1) {
          x.all$wi.eBH = x.all$we.eBH
          x.all$we.eBH = 'we.eBH used in this test'
        }
        
        x.all$`-Log10(FDR)` <- -log10(x.all$wi.eBH)
        x.all$Taxa = row.names(x.all)
        
        write.csv(x.all,'Differential_abundance.csv')
        
        # draw figure
        das <- x.all[x.all$`-Log10(FDR)` >= -log10(pvalue_th),]
        
        das$Species <- row.names(das)
        das <- das[order(das$diff.btw),] 
        
        metadata = as.factor(metadata)
        lev = levels(metadata)
        
        if (order_reverse == T) {
          das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
        } else {
          das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
          
        }
        das$Species = str_replace_all(das$Species, '_',' ')
        das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
        
        das$diff.btw[das$diff.btw == Inf] = 10
        das$diff.btw[das$diff.btw == -Inf] = -10
        das$diff.btw[das$diff.btw <= -10] = -10
        das$diff.btw[das$diff.btw >= 10] = 10
        
          theme_set(theme_bw())  
          
          
          ggplot(das, aes(Species, diff.btw)) +
            geom_line() +  
            geom_hline(yintercept=1, linetype='dotted', col = 'grey') +
            geom_hline(yintercept=-1, linetype='dotted', col = 'grey') +
            geom_point(aes(col=Outcome, size=`-Log10(FDR)`)) +
            scale_size_continuous(range  = c(0.1, 5), limits = c(0, 25), 
                                  breaks = c(1, 6, 11, 16, 21))+ 
            scale_color_brewer(palette="Set1")+
            coord_flip() +          # convert x y axis
            labs(x = 'Taxa', y = "Median difference in clr values")+ 
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7))  
          ggsave('Differential_abundance_case_match.pdf', width=4, height=2.7)

        
      }
    }
    
    # PCoA
    {
      Padonis = adonis2(as.matrix(t(reads_table_x)) ~ Outcome, data = metadata_x,method = "bray", parallel = 8)
      Padonis$`Pr(>F)`
      
      reads_table = as.matrix(decostand(as.data.frame(t(reads_table_all_batch_rm)),method = 'rclr'))
      x = prcomp(as.data.frame(t(reads_table)), center = TRUE, scale. = TRUE)
      
      pic <- x$rotation
      pic <- data.frame(pic,metadata_all)
      
      colnames(pic)[1:2] <- c('X1','X2')
      pic$Outcome
      ggplot(pic, aes(X1, X2,color = Outcome))  +
        geom_point(size=0.5) +
        xlab(paste0("PCoA1")) +
        ylab(paste0("PCoA2"))+
        stat_ellipse(type = "t") + 
        coord_fixed()+ 
        scale_color_brewer(palette="Set1")+
        theme(
          axis.title.x = element_text( size=7),
          axis.title.y = element_text( size=7),
          legend.text = element_text(size=7),
          legend.title = element_text(size=7),
          plot.title = element_text(hjust = 0.5, size = 7)
        ) 
      ggsave('PCoA_Outcome_case_match.pdf', width=4, height=3)
    }
    
    # alpha
    {
      reads_table = reads_table_x
      metadata = metadata_x$Outcome
      
      # alpha diversity 
      reads_table = t(reads_table)
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table <- reads_table$otu.tab.rff
      reads_table <- as.data.frame(reads_table)
      
      alpha <- data.frame(diversity(reads_table))
      alpha$Evenness <- alpha$diversity.reads_table./log(specnumber(reads_table))
      alpha$Ovserved_phage <- colSums(t(reads_table) != 0)
      colnames(alpha) = c('Shannon','Evenness','Ovserved_taxa')
      alpha$Outcome = metadata
      x = kruskal.test(Shannon ~ Outcome,data = alpha)
      alpha$Outcome = metadata
      P = formatC(x$p.value, format = "e", digits = 3)
      
      ggplot(data=alpha, aes(x=Outcome, y=Shannon)) +geom_violin(trim=T, aes(fill=Outcome))+
        geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
        geom_jitter(size = 0.1, alpha= 0.05)+theme_bw()+
        theme_classic()+
        scale_fill_brewer(palette="Set1")+
        labs(y = 'Shannon index')+
        ggtitle(paste0('P = ',P)) +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              plot.title = element_text(size = 7))
      ggsave('Shannon_case_match.pdf',width=2.5, height=2)
    }
  }
}

##### case match #####
setwd('/Desktop/1_KH/VMB_PTB/results/case_match/')
load("~/Desktop/1_KH/VMB_PTB/others/Untitled.RData")
# analysis_all_dif_race
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata_all_3 = metadata_all[keep,]
  table(metadata_all_3$BioProject)
  reads_table_all_3 = reads_table_all_batch_rm[,keep]
  
  cohort = unique(metadata_all_3$Race)
  for (bb in cohort) {
    keep = metadata_all_3$Race == bb; sum(keep)
    
    metadata_all_2 = metadata_all_3[keep,]
    reads_table_all_2 = reads_table_all_3[,keep]

    metadata_x = metadata_all_2; reads_table_x = reads_table_all_2
    
    # case match
    {
      # case matched by the same bioproject, race, similar age, ga_at_collection, and different participants
      {
        metadata_PTB = metadata_x[metadata_x$Outcome == "PTB",]
        metadata_Term = metadata_PTB; metadata_Term[1:nrow(metadata_Term),1:ncol(metadata_Term)] = NA
        
        if (nrow(metadata_PTB) < 6) {next}
        
        metadata_2 = metadata_x[metadata_x$Outcome == "Term",]
        
        for (a in 1: nrow(metadata_PTB)) {
          
          n =1
          keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
            metadata_2$Age <= metadata_PTB$Age[a] +n & 
            metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
            metadata_2$BioProject == metadata_PTB$BioProject[a] & 
            metadata_2$Race == metadata_PTB$Race[a] & 
            metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
            metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
          while (sum(keep[!is.na(keep)]) == 0 & n < 3) {
            n = n+1 
            keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
              metadata_2$Age <= metadata_PTB$Age[a] +n & 
              metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
              metadata_2$BioProject == metadata_PTB$BioProject[a] & 
              metadata_2$Race == metadata_PTB$Race[a] & 
              metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
              metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
          }
          
          if (n < 3) {
            metadata_Term[a,] = metadata_2[which(keep)[1],]
            metadata_2 = metadata_2[-which(keep)[1],]
          }
          
        }
        
        sum(is.na(metadata_Term$ParticipantID))
        
        keep = !is.na(metadata_Term$Run)
        metadata_Term = metadata_Term[keep,]; metadata_PTB = metadata_PTB[keep,]
        
        data = data.frame(Age = c(metadata_PTB$Age,metadata_Term$Age), 
                          Outcome = c(rep('PTB',length(metadata_PTB$Age)), rep('Term',length(metadata_Term$Age))))
        
        if (nrow(data) < 6) {next}
        
        P = kruskal.test(Age ~ Outcome,data = data)
        if (P$p.value < 0.001) {
          P = formatC(P$p.value, format = "e", digits = 3)
        } else {
          P = formatC(P$p.value, format = "f", digits = 3)
        }
        
        ggplot(data, aes(x=Outcome, y=Age)) + geom_violin(trim=T, aes(fill = Outcome))+
          geom_boxplot(fill='white', color="black", width=0.1) + 
          geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
          scale_fill_brewer(palette="Set1")+
          ggtitle(paste0(bb,'\n','P = ',P)) +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                plot.title = element_text(size = 7))
        
        ggsave(paste0("Race_",str_replace(bb,'\\/',"_"),"_",'Age.pdf'), height = 2,width = 2)
        
        data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                          Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
        
        P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
        if (P$p.value < 0.001) {
          P = formatC(P$p.value, format = "e", digits = 3)
        } else {
          P = formatC(P$p.value, format = "f", digits = 3)
        }
        ggplot(data, aes(x=Outcome, y=Gestational_age_in_sample_collection)) + geom_violin(trim=T, aes(fill = Outcome))+
          geom_boxplot(fill='white', color="black", width=0.1) + 
          geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
          labs(y = 'Gestational age at sample collection')+
          scale_fill_brewer(palette="Set1")+
          ggtitle(paste0(bb,'\n','P = ',P)) +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                plot.title = element_text(size = 7))
        ggsave(paste0("Race_",str_replace(bb,'\\/',"_"),"_",'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
        
        metadata = data.frame(ParticipantID = c(metadata_PTB$Patient_ID,metadata_Term$Patient_ID), 
                              BioProject = c(metadata_PTB$BioProject,metadata_Term$BioProject), 
                              Run = c(metadata_PTB$Run,metadata_Term$Run), 
                              gest_day_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection),
                              Age = c(metadata_PTB$Age,metadata_Term$Age), 
                              Race = c(metadata_PTB$Race,metadata_Term$Race), 
                              Outcome = c(metadata_PTB$Outcome,metadata_Term$Outcome))
        
        reads_table = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata$Run]
        reads_table = reads_table[metadata$Run]
        
        metadata_x = metadata; reads_table_x = reads_table
        rm(data, metadata,metadata_2, metadata_PTB, metadata_Term, reads_table)
      }
    }
    table(metadata_x$Race,metadata_x$Outcome)
    table(metadata_x$BioProject,metadata_x$Outcome)
    
    {
      # Differential abundance
      {
        metadata = metadata_x$Outcome
        reads_table = as.data.frame(t(reads_table_x))
        reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
        reads_table <- reads_table$otu.tab.rff
        reads_table <- as.data.frame(t(reads_table))
        
        pvalue_th = 0.05; fold_change_th = 0.5; paired_test = F; order_reverse = F; style = 1; order = NA
        
        if (!is.na(order)[1]) {
          metadata = factor(metadata, levels = order)
        }
        
        conds <- metadata
        
        x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
        
        # paired Wilcoxon Rank Sum test and Welch's t-test
        x.tt <- aldex.ttest(x, paired.test= paired_test)
        
        x.effect <- aldex.effect(x)
        
        x.all <- data.frame(cbind(x.tt,x.effect))
        if (sum(is.nan(x.all$wi.eBH)) > 1) {
          x.all$wi.eBH = x.all$we.eBH
          x.all$we.eBH = 'we.eBH used in this test'
        }
        
        x.all$`-Log10(FDR)` <- -log10(x.all$wi.eBH)
        x.all$Taxa = row.names(x.all)
        x.all = x.all[order(x.all$wi.ep),]
        write.csv(x.all,paste0("Race_",str_replace(bb,'\\/',"_"),"_",'Differential_abundance.csv'))
        
        # draw figure
        das <- x.all[x.all$`-Log10(FDR)` >= -log10(pvalue_th),]
        
        das$Species <- row.names(das)
        das <- das[order(das$diff.btw),] 
        
        metadata = as.factor(metadata)
        lev = levels(metadata)
        
        if (order_reverse == T) {
          das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
        } else {
          das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
          
        }
        
        das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
        
        das$diff.btw[das$diff.btw == Inf] = 10
        das$diff.btw[das$diff.btw == -Inf] = -10
        das$diff.btw[das$diff.btw <= -10] = -10
        das$diff.btw[das$diff.btw >= 10] = 10
        
        
        theme_set(theme_bw())  
        
        ggplot(das, aes(Species, diff.btw)) +
          geom_line() +  
          geom_hline(yintercept=1, linetype='dotted', col = 'grey') +
          geom_hline(yintercept=-1, linetype='dotted', col = 'grey') +
          geom_point(aes(col=Outcome, size=`-Log10(FDR)`)) +
          scale_size_continuous(range  = c(0.1, 5), limits = c(0, 25), 
                                breaks = c(2, 4, 6, 8, 10))+ 
          coord_flip() +          # convert x y axis
          labs(x = 'Taxa', y = "Median difference in clr values")+ 
          ggtitle(paste0(bb)) +
          scale_color_brewer(palette="Set1")+
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                plot.title = element_text(hjust = 0.5, size = 7))
        ggsave(paste0("Race_",str_replace(bb,'\\/',"_"),"_",'Differential_abundance.pdf'), width=4, height=1+nrow(das)*0.1)
        
        if (bb =='White') {
          x <- brewer.pal(n = 1, name = "Set1")
          ggplot(das, aes(Species, diff.btw)) +
            geom_line() +  
            geom_hline(yintercept=1, linetype='dotted', col = 'grey') +
            geom_hline(yintercept=-1, linetype='dotted', col = 'grey') +
            geom_point(aes(col=Outcome, size=`-Log10(FDR)`),color=x[2]) +
            scale_size_continuous(range  = c(0.1, 5), limits = c(0, 25), 
                                  breaks = c(2, 4, 6, 8, 10))+ 
            coord_flip() +          # convert x y axis
            labs(x = 'Taxa', y = "Median difference in clr values")+ 
            ggtitle(paste0(bb)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(hjust = 0.5, size = 7))
          ggsave(paste0("Race_",str_replace(bb,'\\/',"_"),"_",'Differential_abundance.pdf'), width=4, height=1+nrow(das)*0.1)
          
        }
        
      }
      
      # t-SNE
      {
        Padonis = adonis2(as.matrix(t(reads_table_x)) ~ Outcome, data = metadata_x,method = "bray", parallel = 8)
        P = Padonis$`Pr(>F)`; P = P[1]
        if (P < 0.001) {
          P = formatC(P, format = "e", digits = 3)
        } else {
          P = formatC(P, format = "f", digits = 3)
        }
        
        reads_table = as.matrix(decostand(as.data.frame(t(reads_table_x)),method = 'rclr'))
        x = prcomp(as.data.frame(t(reads_table)), center = TRUE, scale. = TRUE)
        
        pic <- x$rotation
        pic <- data.frame(pic,metadata_x)
        
        colnames(pic)[1:2] <- c('X1','X2')
        pic$Outcome
        ggplot(pic, aes(X1, X2,color = Outcome))  +
          geom_point(size=0.5) +
          xlab(paste0("PCoA1")) +
          ylab(paste0("PCoA2"))+
          stat_ellipse(type = "t") + 
          coord_fixed()+ 
          scale_color_brewer(palette="Set1")+
          ggtitle(paste0(bb,'\n','P = ',P)) +
          theme(
            axis.title.x = element_text( size=12),
            axis.title.y = element_text( size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=12),
            plot.title = element_text(hjust = 0.5, size = 12)
          ) 
        
        ggsave(paste0("Race_",str_replace(bb,'\\/',"_"),"_", 'tsne_Outcome.pdf'), width=4, height=6)
      }
      
      # alpha
      {
        reads_table = reads_table_x
        metadata = metadata_x$Outcome
        
        # alpha diversity 
        reads_table = t(reads_table)
        reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
        reads_table <- reads_table$otu.tab.rff
        reads_table <- as.data.frame(reads_table)
        
        alpha <- data.frame(diversity(reads_table))
        alpha$Evenness <- alpha$diversity.reads_table./log(specnumber(reads_table))
        alpha$Ovserved_phage <- colSums(t(reads_table) != 0)
        colnames(alpha) = c('Shannon','Evenness','Ovserved_taxa')
        alpha$Outcome = metadata
        x = kruskal.test(Shannon ~ Outcome,data = alpha)
        P=x$p.value
        if (P < 0.001) {
          P = formatC(P, format = "e", digits = 3)
        } else {
          P = formatC(P, format = "f", digits = 3)
        }
        
        ggplot(data=alpha, aes(x=Outcome, y=Shannon)) +geom_violin(trim=T, aes(fill=Outcome))+
          geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
          geom_jitter(size = 0.1, alpha= 0.25)+theme_bw()+
          theme_classic()+
          scale_fill_brewer(palette="Set1")+
          labs(y = 'Shannon index')+
          ggtitle(paste0(bb,'\n', "n = ", nrow(alpha)/2, " × 2", '\n', 'P = ',P)) +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                plot.title = element_text(size = 7))
        ggsave(paste0("Race_",str_replace(bb,'\\/',"_"),"_",'Shannon.pdf'),width=2.5, height=2)
      }
      
    }
  }
  
}

# analysis_all_dif_race_dif_cohort
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata_all_3 = metadata_all[keep,]
  table(metadata_all_3$BioProject)
  reads_table_all_3 = reads_table_all[,keep]
  
  race = unique(metadata_all_3$Race)
  cohort = unique(metadata_all_3$BioProject)
  for (bb in cohort) {
    keep = metadata_all_3$BioProject == bb; sum(keep)
    
    metadata_all_2 = metadata_all_3[keep,]
    reads_table_all_2 = reads_table_all_3[,keep]
    
    for (aa in race) {
      keep = metadata_all_2$Race == aa; sum(keep)
      
      metadata_x = metadata_all_2[keep,]
      reads_table_x = reads_table_all_2[,keep]
      
      # case match
      {
        # case matched by the same bioproject, race, similar age, ga_at_collection, and different participants
        {
          metadata_PTB = metadata_x[metadata_x$Outcome == "PTB",]
          metadata_Term = metadata_PTB; metadata_Term[1:nrow(metadata_Term),1:ncol(metadata_Term)] = NA
          
          if (nrow(metadata_PTB) < 10) {next}
          
          metadata_2 = metadata_x[metadata_x$Outcome == "Term",]
          
          for (a in 1: nrow(metadata_PTB)) {
            
            n =1
            keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
              metadata_2$Age <= metadata_PTB$Age[a] +n & 
              metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
              metadata_2$BioProject == metadata_PTB$BioProject[a] & 
              metadata_2$Race == metadata_PTB$Race[a] & 
              metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
              metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
            while (sum(keep[!is.na(keep)]) == 0 & n < 3) {
              n = n+1 
              keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
                metadata_2$Age <= metadata_PTB$Age[a] +n & 
                metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
                metadata_2$BioProject == metadata_PTB$BioProject[a] & 
                metadata_2$Race == metadata_PTB$Race[a] & 
                metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
                metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
            }
            
            if (n < 3) {
              metadata_Term[a,] = metadata_2[which(keep)[1],]
              metadata_2 = metadata_2[-which(keep)[1],]
            }
            
          }
          
          sum(is.na(metadata_Term$ParticipantID))
          
          keep = !is.na(metadata_Term$Run)
          metadata_Term = metadata_Term[keep,]; metadata_PTB = metadata_PTB[keep,]
          
          if (nrow(metadata_Term) < 6) {next}
          
          data = data.frame(Age = c(metadata_PTB$Age,metadata_Term$Age), 
                            Outcome = c(rep('PTB',length(metadata_PTB$Age)), rep('Term',length(metadata_Term$Age))))
          
          P = kruskal.test(Age ~ Outcome,data = data)
          if (P$p.value < 0.001) {
            P = formatC(P$p.value, format = "e", digits = 3)
          } else {
            P = formatC(P$p.value, format = "f", digits = 3)
          }
          
          ggplot(data, aes(x=Outcome, y=Age)) + geom_violin(trim=T, aes(fill = Outcome))+
            geom_boxplot(fill='white', color="black", width=0.1) + 
            geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
            scale_fill_brewer(palette="Set1")+
            ggtitle(paste0(bb,'\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          
          ggsave(paste0("Race_",bb,"_",str_replace(aa,'\\/','_'),'Age.pdf'), height = 2,width = 2)
          
          data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                            Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
          
          P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
          if (P$p.value < 0.001) {
            P = formatC(P$p.value, format = "e", digits = 3)
          } else {
            P = formatC(P$p.value, format = "f", digits = 3)
          }
          ggplot(data, aes(x=Outcome, y=Gestational_age_in_sample_collection)) + geom_violin(trim=T, aes(fill = Outcome))+
            geom_boxplot(fill='white', color="black", width=0.1) + 
            geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
            labs(y = 'Gestational age at sample collection')+
            scale_fill_brewer(palette="Set1")+
            ggtitle(paste0(bb,'\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          ggsave(paste0("Race_",bb,"_",str_replace(aa,'\\/','_'),'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
          
          metadata = data.frame(ParticipantID = c(metadata_PTB$Patient_ID,metadata_Term$Patient_ID), 
                                BioProject = c(metadata_PTB$BioProject,metadata_Term$BioProject), 
                                Run = c(metadata_PTB$Run,metadata_Term$Run), 
                                gest_day_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection),
                                Age = c(metadata_PTB$Age,metadata_Term$Age), 
                                Race = c(metadata_PTB$Race,metadata_Term$Race), 
                                Outcome = c(metadata_PTB$Outcome,metadata_Term$Outcome))
          
          reads_table = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata$Run]
          reads_table = reads_table[metadata$Run]
          
          metadata_x = metadata; reads_table_x = reads_table
          rm(data, metadata,metadata_2, metadata_PTB, metadata_Term, reads_table)
        }
      }
      table(metadata_x$Race,metadata_x$Outcome)
      table(metadata_x$BioProject,metadata_x$Outcome)
      
      {
        # race
        data <- as.data.frame(table(metadata_x$Race))
        data = data[order(data$Var1),]
        ggplot(data, aes(x="", y=Freq, fill=Var1)) +
          geom_bar(stat="identity", width=1, color="white") + 
          coord_polar("y", start=0) +
          theme_void() + 
          scale_fill_manual(values= mycolors[[11]], name = paste0('Total number = ', sum(data$Freq), '\n','Race')) 
        ggsave(paste0("Race_",bb,"_", str_replace(aa,'\\/','_'),'_race.pdf'), width=4, height=6)
        
        # Differential abundance
        {
          metadata = metadata_x$Outcome
          reads_table = as.data.frame(t(reads_table_x))
          reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
          reads_table <- reads_table$otu.tab.rff
          reads_table <- as.data.frame(t(reads_table))
          
          pvalue_th = 0.05; fold_change_th = 0.5; paired_test = F; order_reverse = F; style = 1; order = NA
          
          if (!is.na(order)[1]) {
            metadata = factor(metadata, levels = order)
          }
          
          conds <- metadata
          
          x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
          
          # paired Wilcoxon Rank Sum test and Welch's t-test
          x.tt <- aldex.ttest(x, paired.test= paired_test)
          
          x.effect <- aldex.effect(x)
          
          x.all <- data.frame(cbind(x.tt,x.effect))
          if (sum(is.nan(x.all$wi.eBH)) > 1) {
            x.all$wi.eBH = x.all$we.eBH
            x.all$we.eBH = 'we.eBH used in this test'
          }
          
          x.all$`-Log10(FDR)` <- -log10(x.all$wi.eBH)
          x.all$Taxa = row.names(x.all)
          x.all = x.all[order(x.all$wi.ep),]
          write.csv(x.all,paste0("Race_",bb,"_",str_replace(aa,'\\/','_'),'Differential_abundance.csv'))
          
          # draw figure
          das <- x.all[x.all$`-Log10(FDR)` >= -log10(pvalue_th),]
          
          das$Species <- row.names(das)
          das <- das[order(das$diff.btw),] 
          
          metadata = as.factor(metadata)
          lev = levels(metadata)
          
          if (order_reverse == T) {
            das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
          } else {
            das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
            
          }
          
          das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
          
          das$diff.btw[das$diff.btw == Inf] = 10
          das$diff.btw[das$diff.btw == -Inf] = -10
          das$diff.btw[das$diff.btw <= -10] = -10
          das$diff.btw[das$diff.btw >= 10] = 10
          
          
          theme_set(theme_bw())  
          
          ggplot(das, aes(Species, diff.btw)) +
            geom_line() +  
            geom_hline(yintercept=1, linetype='dotted', col = 'grey') +
            geom_hline(yintercept=-1, linetype='dotted', col = 'grey') +
            geom_point(aes(col=Outcome, size=`-Log10(FDR)`)) +
            scale_size_continuous(range  = c(0.1, 5), limits = c(0, 25), 
                                  breaks = c(2, 4, 6, 8, 10))+ 
            coord_flip() +          # convert x y axis
            labs(x = 'Taxa', y = "Median difference in clr values")+ 
            ggtitle(paste0(bb)) +
            scale_color_brewer(palette="Set1")+
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(hjust = 0.5, size = 7))
          ggsave(paste0("Race_",bb,"_",str_replace(aa,'\\/','_'),'Differential_abundance.pdf'), width=4, height=1+nrow(das)*0.1)
          
          
        }
        
        # t-SNE
        {
          Padonis = adonis2(as.matrix(t(reads_table_x)) ~ Outcome, data = metadata_x,method = "bray", parallel = 8)
          P = Padonis$`Pr(>F)`; P = P[1]
          if (P < 0.001) {
            P = formatC(P, format = "e", digits = 3)
          } else {
            P = formatC(P, format = "f", digits = 3)
          }
          
          reads_table = as.matrix(decostand(as.data.frame(t(reads_table_x)),method = 'rclr'))
          x = prcomp(as.data.frame(t(reads_table)), center = TRUE, scale. = TRUE)
          
          pic <- x$rotation
          pic <- data.frame(pic,metadata_x)
          
          colnames(pic)[1:2] <- c('X1','X2')
          pic$Outcome
          ggplot(pic, aes(X1, X2,color = Outcome))  +
            geom_point(size=0.5) +
            xlab(paste0("PCoA1")) +
            ylab(paste0("PCoA2"))+
            stat_ellipse(type = "t") + 
            coord_fixed()+ 
            scale_color_brewer(palette="Set1")+
            ggtitle(paste0(bb,'\n','P = ',P)) +
            theme(
              axis.title.x = element_text( size=12),
              axis.title.y = element_text( size=12),
              legend.text = element_text(size=12),
              legend.title = element_text(size=12),
              plot.title = element_text(hjust = 0.5, size = 12)
            ) 
          
          ggsave(paste0("Race_",bb,"_", str_replace(aa,'\\/','_'),'tsne_Outcome.pdf'), width=4, height=6)
        }
        
        # alpha
        {
          reads_table = reads_table_x
          metadata = metadata_x$Outcome
          
          # alpha diversity 
          reads_table = t(reads_table)
          reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
          reads_table <- reads_table$otu.tab.rff
          reads_table <- as.data.frame(reads_table)
          
          alpha <- data.frame(diversity(reads_table))
          alpha$Evenness <- alpha$diversity.reads_table./log(specnumber(reads_table))
          alpha$Ovserved_phage <- colSums(t(reads_table) != 0)
          colnames(alpha) = c('Shannon','Evenness','Ovserved_taxa')
          alpha$Outcome = metadata
          x = kruskal.test(Shannon ~ Outcome,data = alpha)
          P=x$p.value
          if (P < 0.001) {
            P = formatC(P, format = "e", digits = 3)
          } else {
            P = formatC(P, format = "f", digits = 3)
          }
          
          ggplot(data=alpha, aes(x=Outcome, y=Shannon)) +geom_violin(trim=T, aes(fill=Outcome))+
            geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
            geom_jitter(size = 0.1, alpha= 0.25)+theme_bw()+
            theme_classic()+
            scale_fill_brewer(palette="Set1")+
            labs(y = 'Shannon index')+
            ggtitle(paste0(bb,'\n',"n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          ggsave(paste0("Race_",bb,"_",str_replace(aa,'\\/','_'),'Shannon.pdf'),width=2.5, height=2)
        }
        
      }
    }
  }
  
}

# analysis_all_dif_race_dif_ages
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata_all_3 = metadata_all[keep,]
  table(metadata_all_3$BioProject)
  reads_table_all_3 = reads_table_all_batch_rm[,keep]
  
  race = unique(metadata_all_3$Race)
  
  x = metadata_all[metadata_all$Outcome == "PTB" & !is.na(metadata_all$Age),]
  Age_quantile = quantile(x$Age[!is.na(x$Age)], probs = seq(0, 1, length.out = 4)); Age_quantile = as.numeric(Age_quantile)
  Age_quantile[1] = Age_quantile[1]-1
  Age_quantile
  Age_quantile_list = c("(17 - 25]", "(25- 32]", "(32- 42]")

  for (bb in 1:3) {
    keep = metadata_all_3$Age > Age_quantile[bb] & metadata_all_3$Age <= Age_quantile[bb+1]; sum(keep)
    metadata_all_2 = metadata_all_3[keep,]; reads_table_all_2 = reads_table_all_3[,keep]
    
    for (aa in race) {
      keep = metadata_all_2$Race == aa; sum(keep)
      
      metadata_x = metadata_all_2[keep,]
      reads_table_x = reads_table_all_2[,keep]
      
      # case match
      {
        # case matched by the same bioproject, race, similar age, ga_at_collection, and different participants
        {
          metadata_PTB = metadata_x[metadata_x$Outcome == "PTB",]
          metadata_Term = metadata_PTB; metadata_Term[1:nrow(metadata_Term),1:ncol(metadata_Term)] = NA
          
          if (nrow(metadata_PTB) < 10) {next}
          
          metadata_2 = metadata_x[metadata_x$Outcome == "Term",]
          
          for (a in 1: nrow(metadata_PTB)) {
            
            n =1
            keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
              metadata_2$Age <= metadata_PTB$Age[a] +n & 
              metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
              metadata_2$BioProject == metadata_PTB$BioProject[a] & 
              metadata_2$Race == metadata_PTB$Race[a] & 
              metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
              metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
            while (sum(keep[!is.na(keep)]) == 0 & n < 3) {
              n = n+1 
              keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
                metadata_2$Age <= metadata_PTB$Age[a] +n & 
                metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
                metadata_2$BioProject == metadata_PTB$BioProject[a] & 
                metadata_2$Race == metadata_PTB$Race[a] & 
                metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
                metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
            }
            
            if (n < 3) {
              metadata_Term[a,] = metadata_2[which(keep)[1],]
              metadata_2 = metadata_2[-which(keep)[1],]
            }
            
          }
          
          sum(is.na(metadata_Term$ParticipantID))
          
          keep = !is.na(metadata_Term$Run)
          metadata_Term = metadata_Term[keep,]; metadata_PTB = metadata_PTB[keep,]
          
          if (nrow(metadata_Term) < 6) {next}
          
          data = data.frame(Age = c(metadata_PTB$Age,metadata_Term$Age), 
                            Outcome = c(rep('PTB',length(metadata_PTB$Age)), rep('Term',length(metadata_Term$Age))))
          
          P = kruskal.test(Age ~ Outcome,data = data)
          if (P$p.value < 0.001) {
            P = formatC(P$p.value, format = "e", digits = 3)
          } else {
            P = formatC(P$p.value, format = "f", digits = 3)
          }
          
          ggplot(data, aes(x=Outcome, y=Age)) + geom_violin(trim=T, aes(fill = Outcome))+
            geom_boxplot(fill='white', color="black", width=0.1) + 
            geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
            scale_fill_brewer(palette="Set1")+
            ggtitle(paste0("Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          
          ggsave(paste0("Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Age.pdf'), height = 2,width = 2)
          
          data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                            Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
          
          P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
          if (P$p.value < 0.001) {
            P = formatC(P$p.value, format = "e", digits = 3)
          } else {
            P = formatC(P$p.value, format = "f", digits = 3)
          }
          ggplot(data, aes(x=Outcome, y=Gestational_age_in_sample_collection)) + geom_violin(trim=T, aes(fill = Outcome))+
            geom_boxplot(fill='white', color="black", width=0.1) + 
            geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
            labs(y = 'Gestational age at sample collection')+
            scale_fill_brewer(palette="Set1")+
            ggtitle(paste0("Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          ggsave(paste0("Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
          
          metadata = data.frame(ParticipantID = c(metadata_PTB$Patient_ID,metadata_Term$Patient_ID), 
                                BioProject = c(metadata_PTB$BioProject,metadata_Term$BioProject), 
                                Run = c(metadata_PTB$Run,metadata_Term$Run), 
                                gest_day_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection),
                                Age = c(metadata_PTB$Age,metadata_Term$Age), 
                                Race = c(metadata_PTB$Race,metadata_Term$Race), 
                                Outcome = c(metadata_PTB$Outcome,metadata_Term$Outcome))
          
          reads_table = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata$Run]
          reads_table = reads_table[metadata$Run]
          
          metadata_x = metadata; reads_table_x = reads_table
          rm(data, metadata,metadata_2, metadata_PTB, metadata_Term, reads_table)
        }
      }
      table(metadata_x$Race,metadata_x$Outcome)
      table(metadata_x$BioProject,metadata_x$Outcome)
      
      {
        # Differential abundance
        {
          metadata = metadata_x$Outcome
          reads_table = as.data.frame(t(reads_table_x))
          reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
          reads_table <- reads_table$otu.tab.rff
          reads_table <- as.data.frame(t(reads_table))
          
          pvalue_th = 0.05; fold_change_th = 0.5; paired_test = F; order_reverse = F; style = 1; order = NA
          
          if (!is.na(order)[1]) {
            metadata = factor(metadata, levels = order)
          }
          
          conds <- metadata
          
          x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
          
          # paired Wilcoxon Rank Sum test and Welch's t-test
          x.tt <- aldex.ttest(x, paired.test= paired_test)
          
          x.effect <- aldex.effect(x)
          
          x.all <- data.frame(cbind(x.tt,x.effect))
          if (sum(is.nan(x.all$wi.eBH)) > 1) {
            x.all$wi.eBH = x.all$we.eBH
            x.all$we.eBH = 'we.eBH used in this test'
          }
          
          x.all$`-Log10(FDR)` <- -log10(x.all$wi.eBH)
          x.all$Taxa = row.names(x.all)
          x.all = x.all[order(x.all$wi.ep),]
          write.csv(x.all,paste0("Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Differential_abundance.csv'))
          
          # draw figure
          das <- x.all[x.all$`-Log10(FDR)` >= -log10(pvalue_th),]
          
          das$Species <- row.names(das)
          das <- das[order(das$diff.btw),] 
          
          metadata = as.factor(metadata)
          lev = levels(metadata)
          
          if (order_reverse == T) {
            das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
          } else {
            das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
            
          }
          
          das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
          
          das$diff.btw[das$diff.btw == Inf] = 10
          das$diff.btw[das$diff.btw == -Inf] = -10
          das$diff.btw[das$diff.btw <= -10] = -10
          das$diff.btw[das$diff.btw >= 10] = 10
          
          
          theme_set(theme_bw())  
          
          ggplot(das, aes(Species, diff.btw)) +
            geom_line() +  
            geom_hline(yintercept=1, linetype='dotted', col = 'grey') +
            geom_hline(yintercept=-1, linetype='dotted', col = 'grey') +
            geom_point(aes(col=Outcome, size=`-Log10(FDR)`)) +
            scale_size_continuous(range  = c(0.1, 5), limits = c(0, 25), 
                                  breaks = c(2, 4, 6, 8, 10))+ 
            coord_flip() +          # convert x y axis
            labs(x = 'Taxa', y = "Median difference in clr values")+ 
            ggtitle(paste0("Race: ", aa, '\n', "Age: ", Age_quantile_list[bb])) +
            scale_color_brewer(palette="Set1")+
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(hjust = 0.5, size = 7))
          ggsave(paste0("Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Differential_abundance.pdf'), width=4, height=1+nrow(das)*0.1)
          
          
        }
        
        # t-SNE
        {
          Padonis = adonis2(as.matrix(t(reads_table_x)) ~ Outcome, data = metadata_x,method = "bray", parallel = 8)
          P = Padonis$`Pr(>F)`; P = P[1]
          if (P < 0.001) {
            P = formatC(P, format = "e", digits = 3)
          } else {
            P = formatC(P, format = "f", digits = 3)
          }
          
          reads_table = as.matrix(decostand(as.data.frame(t(reads_table_x)),method = 'rclr'))
          x = prcomp(as.data.frame(t(reads_table)), center = TRUE, scale. = TRUE)
          
          pic <- x$rotation
          pic <- data.frame(pic,metadata_x)
          
          colnames(pic)[1:2] <- c('X1','X2')
          pic$Outcome
          ggplot(pic, aes(X1, X2,color = Outcome))  +
            geom_point(size=0.5) +
            xlab(paste0("PCoA1")) +
            ylab(paste0("PCoA2"))+
            stat_ellipse(type = "t") + 
            coord_fixed()+ 
            scale_color_brewer(palette="Set1")+
            ggtitle(paste0("Race: ", aa, '\n', "Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
            theme(
              axis.title.x = element_text( size=12),
              axis.title.y = element_text( size=12),
              legend.text = element_text(size=12),
              legend.title = element_text(size=12),
              plot.title = element_text(hjust = 0.5, size = 12)
            ) 
          
          ggsave(paste0("Race_",Age_quantile_list[bb],"_", str_replace(aa,'\\/','_'),'tsne_Outcome.pdf'), width=4, height=6)
        }
        
        # alpha
        {
          reads_table = reads_table_x
          metadata = metadata_x$Outcome
          
          # alpha diversity 
          reads_table = t(reads_table)
          reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
          reads_table <- reads_table$otu.tab.rff
          reads_table <- as.data.frame(reads_table)
          
          alpha <- data.frame(diversity(reads_table))
          alpha$Evenness <- alpha$diversity.reads_table./log(specnumber(reads_table))
          alpha$Ovserved_phage <- colSums(t(reads_table) != 0)
          colnames(alpha) = c('Shannon','Evenness','Ovserved_taxa')
          alpha$Outcome = metadata
          x = kruskal.test(Shannon ~ Outcome,data = alpha)
          P=x$p.value
          if (P < 0.001) {
            P = formatC(P, format = "e", digits = 3)
          } else {
            P = formatC(P, format = "f", digits = 3)
          }
          
          ggplot(data=alpha, aes(x=Outcome, y=Shannon)) +geom_violin(trim=T, aes(fill=Outcome))+
            geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
            geom_jitter(size = 0.1, alpha= 0.25)+theme_bw()+
            theme_classic()+
            scale_fill_brewer(palette="Set1")+
            labs(y = 'Shannon index')+
            ggtitle(paste0("Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n',"n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          ggsave(paste0("Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Shannon.pdf'),width=2.5, height=2)
        }
        
      }
    }
  }
  
}

# analysis_all_dif_race_dif_age_dif_cohort
{
  x = metadata_all[metadata_all$Outcome == "PTB" & !is.na(metadata_all$Age),]
  Age_quantile = quantile(x$Age[!is.na(x$Age)], probs = seq(0, 1, length.out = 4)); Age_quantile = as.numeric(Age_quantile)
  Age_quantile[1] = Age_quantile[1]-1
  Age_quantile
  Age_quantile_list = c("(17 - 25]", "(25- 32]", "(32- 42]")
  
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) &
    !is.na(metadata_all$gest_day_collection) & metadata_all$Race == 'B/AA'; sum(keep)
  
  metadata_all_3 = metadata_all[keep,]
  table(metadata_all_3$BioProject)
  reads_table_all_3 = reads_table_all[,keep]
  
  cohort = unique(metadata_all_3$BioProject)
  for (bb in cohort) {
    keep = metadata_all_3$BioProject == bb; sum(keep)
    
    metadata_all_2 = metadata_all_3[keep,]
    reads_table_all_2 = reads_table_all_3[,keep]
    
    for (aa in 1: 3) {
      keep = metadata_all_2$Age > Age_quantile[aa] & metadata_all_2$Age <= Age_quantile[aa+1]; sum(keep)
      metadata_x = metadata_all_2[keep,]; reads_table_x = reads_table_all_2[,keep]
      
      # case match
      {
        # case matched by the same bioproject, race, similar age, ga_at_collection, and different participants
        {
          metadata_PTB = metadata_x[metadata_x$Outcome == "PTB",]
          metadata_Term = metadata_PTB; metadata_Term[1:nrow(metadata_Term),1:ncol(metadata_Term)] = NA
          
          if (nrow(metadata_PTB) < 10) {next}
          
          metadata_2 = metadata_x[metadata_x$Outcome == "Term",]
          
          for (a in 1: nrow(metadata_PTB)) {
            
            n =1
            keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
              metadata_2$Age <= metadata_PTB$Age[a] +n & 
              metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
              metadata_2$BioProject == metadata_PTB$BioProject[a] & 
              metadata_2$Race == metadata_PTB$Race[a] & 
              metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
              metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
            while (sum(keep[!is.na(keep)]) == 0 & n < 3) {
              n = n+1 
              keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
                metadata_2$Age <= metadata_PTB$Age[a] +n & 
                metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
                metadata_2$BioProject == metadata_PTB$BioProject[a] & 
                metadata_2$Race == metadata_PTB$Race[a] & 
                metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
                metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
            }
            
            if (n < 3) {
              metadata_Term[a,] = metadata_2[which(keep)[1],]
              metadata_2 = metadata_2[-which(keep)[1],]
            }
            
          }
          
          sum(is.na(metadata_Term$ParticipantID))
          
          keep = !is.na(metadata_Term$Run)
          metadata_Term = metadata_Term[keep,]; metadata_PTB = metadata_PTB[keep,]
          
          data = data.frame(Age = c(metadata_PTB$Age,metadata_Term$Age), 
                            Outcome = c(rep('PTB',length(metadata_PTB$Age)), rep('Term',length(metadata_Term$Age))))
          
          P = kruskal.test(Age ~ Outcome,data = data)
          if (P$p.value < 0.001) {
            P = formatC(P$p.value, format = "e", digits = 3)
          } else {
            P = formatC(P$p.value, format = "f", digits = 3)
          }
          
          ggplot(data, aes(x=Outcome, y=Age)) + geom_violin(trim=T, aes(fill = Outcome))+
            geom_boxplot(fill='white', color="black", width=0.1) + 
            geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
            scale_fill_brewer(palette="Set1")+
            ggtitle(paste0("Race: B/AA", '\n',"Age: ", Age_quantile_list[aa],'\n', bb,"\n","n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          
          ggsave(paste0("Age_",bb,"_",aa,'Age.pdf'), height = 2,width = 2)
          
          data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                            Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
          
          P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
          if (P$p.value < 0.001) {
            P = formatC(P$p.value, format = "e", digits = 3)
          } else {
            P = formatC(P$p.value, format = "f", digits = 3)
          }
          ggplot(data, aes(x=Outcome, y=Gestational_age_in_sample_collection)) + geom_violin(trim=T, aes(fill = Outcome))+
            geom_boxplot(fill='white', color="black", width=0.1) + 
            geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
            labs(y = 'Gestational age at sample collection')+
            scale_fill_brewer(palette="Set1")+
            ggtitle(paste0("Race: B/AA", '\n',"Age: ", Age_quantile_list[aa],'\n', bb,"\n","n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          ggsave(paste0("Age_",bb,"_",aa,'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
          
          metadata = data.frame(ParticipantID = c(metadata_PTB$Patient_ID,metadata_Term$Patient_ID), 
                                BioProject = c(metadata_PTB$BioProject,metadata_Term$BioProject), 
                                Run = c(metadata_PTB$Run,metadata_Term$Run), 
                                gest_day_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection),
                                Age = c(metadata_PTB$Age,metadata_Term$Age), 
                                Race = c(metadata_PTB$Race,metadata_Term$Race), 
                                Outcome = c(metadata_PTB$Outcome,metadata_Term$Outcome))
          
          reads_table = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata$Run]
          reads_table = reads_table[metadata$Run]
          
          metadata_x = metadata; reads_table_x = reads_table
          rm(data, metadata,metadata_2, metadata_PTB, metadata_Term, reads_table)
        }
      }
      table(metadata_x$Race,metadata_x$Outcome)
      table(metadata_x$BioProject,metadata_x$Outcome)
      
      {
        # Differential abundance
        {
          metadata = metadata_x$Outcome
          reads_table = as.data.frame(t(reads_table_x))
          reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
          reads_table <- reads_table$otu.tab.rff
          reads_table <- as.data.frame(t(reads_table))
          
          pvalue_th = 0.05; fold_change_th = 0.5; paired_test = F; order_reverse = F; style = 1; order = NA

            if (!is.na(order)[1]) {
              metadata = factor(metadata, levels = order)
            }
            
            conds <- metadata
            
            x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
            
            # paired Wilcoxon Rank Sum test and Welch's t-test
            x.tt <- aldex.ttest(x, paired.test= paired_test)
            
            x.effect <- aldex.effect(x)
            
            x.all <- data.frame(cbind(x.tt,x.effect))
            if (sum(is.nan(x.all$wi.eBH)) > 1) {
              x.all$wi.eBH = x.all$we.eBH
              x.all$we.eBH = 'we.eBH used in this test'
            }
            
            x.all$`-Log10(FDR)` <- -log10(x.all$wi.eBH)
            x.all$Taxa = row.names(x.all)
            x.all = x.all[order(x.all$wi.ep),]
            write.csv(x.all,paste0("Age_",bb,"_",aa,'Differential_abundance.csv'))
            
            # draw figure
            das <- x.all[x.all$`-Log10(FDR)` >= -log10(pvalue_th),]
            
            das$Species <- row.names(das)
            das <- das[order(das$diff.btw),] 
            
            metadata = as.factor(metadata)
            lev = levels(metadata)
            
            if (order_reverse == T) {
              das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
            } else {
              das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
              
            }
            
            das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
            
            das$diff.btw[das$diff.btw == Inf] = 10
            das$diff.btw[das$diff.btw == -Inf] = -10
            das$diff.btw[das$diff.btw <= -10] = -10
            das$diff.btw[das$diff.btw >= 10] = 10
            

              theme_set(theme_bw())  
              
              ggplot(das, aes(Species, diff.btw)) +
                geom_line() +  
                geom_hline(yintercept=1, linetype='dotted', col = 'grey') +
                geom_hline(yintercept=-1, linetype='dotted', col = 'grey') +
                geom_point(aes(col=Outcome, size=`-Log10(FDR)`)) +
                scale_size_continuous(range  = c(0.1, 5), limits = c(0, 25), 
                                      breaks = c(2, 4, 6, 8, 10))+ 
                coord_flip() +          # convert x y axis
                labs(x = 'Taxa', y = "Median difference in clr values")+ 
                ggtitle(paste0(bb)) +
                scale_color_brewer(palette="Set1")+
                theme(axis.title = element_text(size = 7), 
                      axis.text = element_text(size = 7), 
                      legend.text = element_text(size = 7), 
                      legend.title = element_text(size = 7),
                      plot.title = element_text(hjust = 0.5, size = 7))
              ggsave(paste0("Age_",bb,"_",aa,'Differential_abundance.pdf'), width=4, height=1+nrow(das)*0.1)


        }
        
        # t-SNE
        {
          Padonis = adonis2(as.matrix(t(reads_table_x)) ~ Outcome, data = metadata_x,method = "bray", parallel = 8)
          P = Padonis$`Pr(>F)`; P = P[1]
          if (P < 0.001) {
            P = formatC(P, format = "e", digits = 3)
          } else {
            P = formatC(P, format = "f", digits = 3)
          }
          
          reads_table = as.matrix(decostand(as.data.frame(t(reads_table_x)),method = 'rclr'))
          x = prcomp(as.data.frame(t(reads_table)), center = TRUE, scale. = TRUE)
          
          pic <- x$rotation
          pic <- data.frame(pic,metadata_x)
          
          colnames(pic)[1:2] <- c('X1','X2')
          pic$Outcome
          ggplot(pic, aes(X1, X2,color = Outcome))  +
            geom_point(size=0.5) +
            xlab(paste0("PCoA1")) +
            ylab(paste0("PCoA2"))+
            stat_ellipse(type = "t") + 
            coord_fixed()+ 
            scale_color_brewer(palette="Set1")+
            ggtitle(paste0("Race: B/AA", '\n',"Age: ", Age_quantile_list[aa],'\n', bb,"\n","n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
            theme(
              axis.title.x = element_text( size=12),
              axis.title.y = element_text( size=12),
              legend.text = element_text(size=12),
              legend.title = element_text(size=12),
              plot.title = element_text(hjust = 0.5, size = 12)
            ) 
          
          ggsave(paste0("Age_",bb,"_", aa,'tsne_Outcome.pdf'), width=4, height=6)
        }
        
        # alpha
        {
          reads_table = reads_table_x
          metadata = metadata_x$Outcome
          
          # alpha diversity 
          reads_table = t(reads_table)
          reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
          reads_table <- reads_table$otu.tab.rff
          reads_table <- as.data.frame(reads_table)
          
          alpha <- data.frame(diversity(reads_table))
          alpha$Evenness <- alpha$diversity.reads_table./log(specnumber(reads_table))
          alpha$Ovserved_phage <- colSums(t(reads_table) != 0)
          colnames(alpha) = c('Shannon','Evenness','Ovserved_taxa')
          alpha$Outcome = metadata
          x = kruskal.test(Shannon ~ Outcome,data = alpha)
          P=x$p.value
          if (P < 0.001) {
            P = formatC(P, format = "e", digits = 3)
          } else {
            P = formatC(P, format = "f", digits = 3)
          }
          
          ggplot(data=alpha, aes(x=Outcome, y=Shannon)) +geom_violin(trim=T, aes(fill=Outcome))+
            geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
            geom_jitter(size = 0.1, alpha= 0.25)+theme_bw()+
            theme_classic()+
            scale_fill_brewer(palette="Set1")+
            labs(y = 'Shannon index')+
            ggtitle(paste0("Race: B/AA", '\n',"Age: ", Age_quantile_list[aa],'\n', bb,"\n","n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          ggsave(paste0("Age_",bb,"_",aa,'Shannon.pdf'),width=2.5, height=2)
        }
        
      }
    }
  }
  
}

# analysis_all_dif_race_dif_ga_collection
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata_all_3 = metadata_all[keep,]
  table(metadata_all_3$BioProject)
  reads_table_all_3 = reads_table_all_batch_rm[,keep]
  
  race = unique(metadata_all_3$Race)
  
  x = metadata_all[metadata_all$Outcome == "PTB" & !is.na(metadata_all$Age) & !is.na(metadata_all$gest_day_collection),]
  gest_day_collection_quantile = quantile(x$gest_day_collection, probs = seq(0, 1, length.out = 4)); gest_day_collection_quantile = as.numeric(gest_day_collection_quantile)
  gest_day_collection_quantile[1] = gest_day_collection_quantile[1]-1
  gest_day_collection_quantile
  gest_day_collection_quantile_list = c("(44 - 135.0]", "(135- 179]", "(179- 268]")
  
  for (cc in 1:3) {
    keep = metadata_all_3$gest_day_collection > gest_day_collection_quantile[cc] & metadata_all_3$gest_day_collection <= gest_day_collection_quantile[cc+1]; sum(keep)
    metadata_all_4 = metadata_all_3[keep,]; reads_table_all_4 = reads_table_all_3[,keep]
    
      for (aa in race) {
        keep = metadata_all_2$Race == aa
        if (sum(keep) < 10) {next}
        
        metadata_x = metadata_all_4[keep,]
        reads_table_x = reads_table_all_4[,keep]
        
        # case match
        {
          # case matched by the same bioproject, race, similar age, ga_at_collection, and different participants
          {
            metadata_PTB = metadata_x[metadata_x$Outcome == "PTB",]
            metadata_Term = metadata_PTB; metadata_Term[1:nrow(metadata_Term),1:ncol(metadata_Term)] = NA
            
            if (nrow(metadata_PTB) < 10) {next}
            
            metadata_2 = metadata_x[metadata_x$Outcome == "Term",]
            
            for (a in 1: nrow(metadata_PTB)) {
              
              n =1
              keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
                metadata_2$Age <= metadata_PTB$Age[a] +n & 
                metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
                metadata_2$BioProject == metadata_PTB$BioProject[a] & 
                metadata_2$Race == metadata_PTB$Race[a] & 
                metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
                metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
              while (sum(keep[!is.na(keep)]) == 0 & n < 3) {
                n = n+1 
                keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
                  metadata_2$Age <= metadata_PTB$Age[a] +n & 
                  metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
                  metadata_2$BioProject == metadata_PTB$BioProject[a] & 
                  metadata_2$Race == metadata_PTB$Race[a] & 
                  metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
                  metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
              }
              
              if (n < 3) {
                metadata_Term[a,] = metadata_2[which(keep)[1],]
                metadata_2 = metadata_2[-which(keep)[1],]
              }
              
            }
            
            sum(is.na(metadata_Term$ParticipantID))
            
            keep = !is.na(metadata_Term$Run)
            metadata_Term = metadata_Term[keep,]; metadata_PTB = metadata_PTB[keep,]
            
            if (nrow(metadata_Term) < 6) {next}
            
            data = data.frame(Age = c(metadata_PTB$Age,metadata_Term$Age), 
                              Outcome = c(rep('PTB',length(metadata_PTB$Age)), rep('Term',length(metadata_Term$Age))))
            
            P = kruskal.test(Age ~ Outcome,data = data)
            if (P$p.value < 0.001) {
              P = formatC(P$p.value, format = "e", digits = 3)
            } else {
              P = formatC(P$p.value, format = "f", digits = 3)
            }
            
            ggplot(data, aes(x=Outcome, y=Age)) + geom_violin(trim=T, aes(fill = Outcome))+
              geom_boxplot(fill='white', color="black", width=0.1) + 
              geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
              scale_fill_brewer(palette="Set1")+
              ggtitle(paste0("Race: ", aa, '\n',
                             "Gestational age: ", gest_day_collection_quantile_list[cc],'\n',
                             'P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_",str_replace(aa,'\\/','_'),'Age.pdf'), height = 2,width = 2)
            
            data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                              Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
            
            P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
            if (P$p.value < 0.001) {
              P = formatC(P$p.value, format = "e", digits = 3)
            } else {
              P = formatC(P$p.value, format = "f", digits = 3)
            }
            ggplot(data, aes(x=Outcome, y=Gestational_age_in_sample_collection)) + geom_violin(trim=T, aes(fill = Outcome))+
              geom_boxplot(fill='white', color="black", width=0.1) + 
              geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
              labs(y = 'Gestational age at sample collection')+
              scale_fill_brewer(palette="Set1")+
              ggtitle(paste0("Race: ", aa, '\n',
                             "Gestational age: ", gest_day_collection_quantile_list[cc],'\n',
                             'P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_",str_replace(aa,'\\/','_'),'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
            
            metadata = data.frame(ParticipantID = c(metadata_PTB$Patient_ID,metadata_Term$Patient_ID), 
                                  BioProject = c(metadata_PTB$BioProject,metadata_Term$BioProject), 
                                  Run = c(metadata_PTB$Run,metadata_Term$Run), 
                                  gest_day_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection),
                                  Age = c(metadata_PTB$Age,metadata_Term$Age), 
                                  Race = c(metadata_PTB$Race,metadata_Term$Race), 
                                  Outcome = c(metadata_PTB$Outcome,metadata_Term$Outcome))
            
            reads_table = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata$Run]
            reads_table = reads_table[metadata$Run]
            
            metadata_x = metadata; reads_table_x = reads_table
            rm(data, metadata,metadata_2, metadata_PTB, metadata_Term, reads_table)
          }
        }
        table(metadata_x$Race,metadata_x$Outcome)
        table(metadata_x$BioProject,metadata_x$Outcome)
        
        {
          # Differential abundance
          {
            metadata = metadata_x$Outcome
            reads_table = as.data.frame(t(reads_table_x))
            reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
            reads_table <- reads_table$otu.tab.rff
            reads_table <- as.data.frame(t(reads_table))
            
            pvalue_th = 0.05; fold_change_th = 0.5; paired_test = F; order_reverse = F; style = 1; order = NA
            
            if (!is.na(order)[1]) {
              metadata = factor(metadata, levels = order)
            }
            
            conds <- metadata
            
            x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
            
            # paired Wilcoxon Rank Sum test and Welch's t-test
            x.tt <- aldex.ttest(x, paired.test= paired_test)
            
            x.effect <- aldex.effect(x)
            
            x.all <- data.frame(cbind(x.tt,x.effect))
            if (sum(is.nan(x.all$wi.eBH)) > 1) {
              x.all$wi.eBH = x.all$we.eBH
              x.all$we.eBH = 'we.eBH used in this test'
            }
            
            x.all$`-Log10(FDR)` <- -log10(x.all$wi.eBH)
            x.all$Taxa = row.names(x.all)
            x.all = x.all[order(x.all$wi.ep),]
            write.csv(x.all,paste0("Ga_",gest_day_collection_quantile_list[cc],"_",str_replace(aa,'\\/','_'),'Differential_abundance.csv'))
            
            # draw figure
            das <- x.all[x.all$`-Log10(FDR)` >= -log10(pvalue_th),]
            
            das$Species <- row.names(das)
            das <- das[order(das$diff.btw),] 
            
            metadata = as.factor(metadata)
            lev = levels(metadata)
            
            if (order_reverse == T) {
              das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
            } else {
              das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
              
            }
            
            das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
            
            das$diff.btw[das$diff.btw == Inf] = 10
            das$diff.btw[das$diff.btw == -Inf] = -10
            das$diff.btw[das$diff.btw <= -10] = -10
            das$diff.btw[das$diff.btw >= 10] = 10
            
            
            theme_set(theme_bw())  
            
            ggplot(das, aes(Species, diff.btw)) +
              geom_line() +  
              geom_hline(yintercept=1, linetype='dotted', col = 'grey') +
              geom_hline(yintercept=-1, linetype='dotted', col = 'grey') +
              geom_point(aes(col=Outcome, size=`-Log10(FDR)`)) +
              scale_size_continuous(range  = c(0.1, 5), limits = c(0, 25), 
                                    breaks = c(2, 4, 6, 8, 10))+ 
              coord_flip() +          # convert x y axis
              labs(x = 'Taxa', y = "Median difference in clr values")+ 
              ggtitle(paste0("Race: ", aa, '\n', 
                             "Gestational age: ", gest_day_collection_quantile_list[cc])) +
              scale_color_brewer(palette="Set1")+
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(hjust = 0.5, size = 7))
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_",str_replace(aa,'\\/','_'),'Differential_abundance.pdf'), width=4, height=1+nrow(das)*0.1)
            
            
          }
          
          # t-SNE
          {
            Padonis = adonis2(as.matrix(t(reads_table_x)) ~ Outcome, data = metadata_x,method = "bray", parallel = 8)
            P = Padonis$`Pr(>F)`; P = P[1]
            if (P < 0.001) {
              P = formatC(P, format = "e", digits = 3)
            } else {
              P = formatC(P, format = "f", digits = 3)
            }
            
            reads_table = as.matrix(decostand(as.data.frame(t(reads_table_x)),method = 'rclr'))
            x = prcomp(as.data.frame(t(reads_table)), center = TRUE, scale. = TRUE)
            
            pic <- x$rotation
            pic <- data.frame(pic,metadata_x)
            
            colnames(pic)[1:2] <- c('X1','X2')
            pic$Outcome
            ggplot(pic, aes(X1, X2,color = Outcome))  +
              geom_point(size=0.5) +
              xlab(paste0("PCoA1")) +
              ylab(paste0("PCoA2"))+
              stat_ellipse(type = "t") + 
              coord_fixed()+ 
              scale_color_brewer(palette="Set1")+
              ggtitle(paste0("Race: ", aa, '\n', 
                             "Gestational age: ", gest_day_collection_quantile_list[cc],'\n',
                             'P = ',P)) +
              theme(
                axis.title.x = element_text( size=12),
                axis.title.y = element_text( size=12),
                legend.text = element_text(size=12),
                legend.title = element_text(size=12),
                plot.title = element_text(hjust = 0.5, size = 12)
              ) 
            
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_", str_replace(aa,'\\/','_'),'tsne_Outcome.pdf'), width=4, height=6)
          }
          
          # alpha
          {
            reads_table = reads_table_x
            metadata = metadata_x$Outcome
            
            # alpha diversity 
            reads_table = t(reads_table)
            reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
            reads_table <- reads_table$otu.tab.rff
            reads_table <- as.data.frame(reads_table)
            
            alpha <- data.frame(diversity(reads_table))
            alpha$Evenness <- alpha$diversity.reads_table./log(specnumber(reads_table))
            alpha$Ovserved_phage <- colSums(t(reads_table) != 0)
            colnames(alpha) = c('Shannon','Evenness','Ovserved_taxa')
            alpha$Outcome = metadata
            x = kruskal.test(Shannon ~ Outcome,data = alpha)
            P=x$p.value
            if (P < 0.001) {
              P = formatC(P, format = "e", digits = 3)
            } else {
              P = formatC(P, format = "f", digits = 3)
            }
            
            ggplot(data=alpha, aes(x=Outcome, y=Shannon)) +geom_violin(trim=T, aes(fill=Outcome))+
              geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
              geom_jitter(size = 0.1, alpha= 0.25)+theme_bw()+
              theme_classic()+
              scale_fill_brewer(palette="Set1")+
              labs(y = 'Shannon index')+
              ggtitle(paste0("Race: ", aa, '\n',
                             "Gestational age: ", gest_day_collection_quantile_list[cc],'\n',
                             "n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_",str_replace(aa,'\\/','_'),'Shannon.pdf'),width=2.5, height=2)
          }
          
        }
      }
    }

  
  
}

# analysis_all_dif_race_dif_ages_dif_ga_collection
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata_all_3 = metadata_all[keep,]
  table(metadata_all_3$BioProject)
  reads_table_all_3 = reads_table_all_batch_rm[,keep]
  
  race = unique(metadata_all_3$Race)
  
  x = metadata_all[metadata_all$Outcome == "PTB" & !is.na(metadata_all$Age) & !is.na(metadata_all$gest_day_collection),]
  Age_quantile = quantile(x$Age[!is.na(x$Age)], probs = seq(0, 1, length.out = 4)); Age_quantile = as.numeric(Age_quantile)
  Age_quantile[1] = Age_quantile[1]-1
  Age_quantile
  Age_quantile_list = c("(17 - 25.0]", "(25.0- 31.0]", "(31.0- 42.0]")
  
  gest_day_collection_quantile = quantile(x$gest_day_collection, probs = seq(0, 1, length.out = 4)); gest_day_collection_quantile = as.numeric(gest_day_collection_quantile)
  gest_day_collection_quantile[1] = gest_day_collection_quantile[1]-1
  gest_day_collection_quantile
  gest_day_collection_quantile_list = c("(44 - 135.0]", "(135- 179]", "(179- 268]")
  
  for (cc in 1:3) {
    keep = metadata_all_3$gest_day_collection > gest_day_collection_quantile[cc] & metadata_all_3$gest_day_collection <= gest_day_collection_quantile[cc+1]; sum(keep)
    metadata_all_4 = metadata_all_3[keep,]; reads_table_all_4 = reads_table_all_3[,keep]
    
    for (bb in 1:3) {
      keep = metadata_all_4$Age > Age_quantile[bb] & metadata_all_4$Age <= Age_quantile[bb+1]; sum(keep)
      metadata_all_2 = metadata_all_4[keep,]; reads_table_all_2 = reads_table_all_4[,keep]
      
      for (aa in race) {
        keep = metadata_all_2$Race == aa
        if (sum(keep) < 10) {next}

        metadata_x = metadata_all_2[keep,]
        reads_table_x = reads_table_all_2[,keep]
        
        # case match
        {
          # case matched by the same bioproject, race, similar age, ga_at_collection, and different participants
          {
            metadata_PTB = metadata_x[metadata_x$Outcome == "PTB",]
            metadata_Term = metadata_PTB; metadata_Term[1:nrow(metadata_Term),1:ncol(metadata_Term)] = NA
            
            if (nrow(metadata_PTB) < 10) {next}
            
            metadata_2 = metadata_x[metadata_x$Outcome == "Term",]
            
            for (a in 1: nrow(metadata_PTB)) {
              
              n =1
              keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
                metadata_2$Age <= metadata_PTB$Age[a] +n & 
                metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
                metadata_2$BioProject == metadata_PTB$BioProject[a] & 
                metadata_2$Race == metadata_PTB$Race[a] & 
                metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
                metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
              while (sum(keep[!is.na(keep)]) == 0 & n < 3) {
                n = n+1 
                keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
                  metadata_2$Age <= metadata_PTB$Age[a] +n & 
                  metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
                  metadata_2$BioProject == metadata_PTB$BioProject[a] & 
                  metadata_2$Race == metadata_PTB$Race[a] & 
                  metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
                  metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
              }
              
              if (n < 3) {
                metadata_Term[a,] = metadata_2[which(keep)[1],]
                metadata_2 = metadata_2[-which(keep)[1],]
              }
              
            }
            
            sum(is.na(metadata_Term$ParticipantID))
            
            keep = !is.na(metadata_Term$Run)
            metadata_Term = metadata_Term[keep,]; metadata_PTB = metadata_PTB[keep,]
            
            if (nrow(metadata_Term) < 6) {next}
            
            data = data.frame(Age = c(metadata_PTB$Age,metadata_Term$Age), 
                              Outcome = c(rep('PTB',length(metadata_PTB$Age)), rep('Term',length(metadata_Term$Age))))
            
            P = kruskal.test(Age ~ Outcome,data = data)
            if (P$p.value < 0.001) {
              P = formatC(P$p.value, format = "e", digits = 3)
            } else {
              P = formatC(P$p.value, format = "f", digits = 3)
            }
            
            ggplot(data, aes(x=Outcome, y=Age)) + geom_violin(trim=T, aes(fill = Outcome))+
              geom_boxplot(fill='white', color="black", width=0.1) + 
              geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
              scale_fill_brewer(palette="Set1")+
              ggtitle(paste0("Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n',
                             "Gestational age: ", gest_day_collection_quantile_list[cc],'\n',
                             'P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Age.pdf'), height = 2,width = 2)
            
            data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                              Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
            
            P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
            if (P$p.value < 0.001) {
              P = formatC(P$p.value, format = "e", digits = 3)
            } else {
              P = formatC(P$p.value, format = "f", digits = 3)
            }
            ggplot(data, aes(x=Outcome, y=Gestational_age_in_sample_collection)) + geom_violin(trim=T, aes(fill = Outcome))+
              geom_boxplot(fill='white', color="black", width=0.1) + 
              geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
              labs(y = 'Gestational age at sample collection')+
              scale_fill_brewer(palette="Set1")+
              ggtitle(paste0("Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n',
                             "Gestational age: ", gest_day_collection_quantile_list[cc],'\n',
                             'P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
            
            metadata = data.frame(ParticipantID = c(metadata_PTB$Patient_ID,metadata_Term$Patient_ID), 
                                  BioProject = c(metadata_PTB$BioProject,metadata_Term$BioProject), 
                                  Run = c(metadata_PTB$Run,metadata_Term$Run), 
                                  gest_day_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection),
                                  Age = c(metadata_PTB$Age,metadata_Term$Age), 
                                  Race = c(metadata_PTB$Race,metadata_Term$Race), 
                                  Outcome = c(metadata_PTB$Outcome,metadata_Term$Outcome))
            
            reads_table = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata$Run]
            reads_table = reads_table[metadata$Run]
            
            metadata_x = metadata; reads_table_x = reads_table
            rm(data, metadata,metadata_2, metadata_PTB, metadata_Term, reads_table)
          }
        }
        table(metadata_x$Race,metadata_x$Outcome)
        table(metadata_x$BioProject,metadata_x$Outcome)
        
        {
          # Differential abundance
          {
            metadata = metadata_x$Outcome
            reads_table = as.data.frame(t(reads_table_x))
            reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
            reads_table <- reads_table$otu.tab.rff
            reads_table <- as.data.frame(t(reads_table))
            
            pvalue_th = 0.05; fold_change_th = 0.5; paired_test = F; order_reverse = F; style = 1; order = NA
            
            if (!is.na(order)[1]) {
              metadata = factor(metadata, levels = order)
            }
            
            conds <- metadata
            
            x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
            
            # paired Wilcoxon Rank Sum test and Welch's t-test
            x.tt <- aldex.ttest(x, paired.test= paired_test)
            
            x.effect <- aldex.effect(x)
            
            x.all <- data.frame(cbind(x.tt,x.effect))
            if (sum(is.nan(x.all$wi.eBH)) > 1) {
              x.all$wi.eBH = x.all$we.eBH
              x.all$we.eBH = 'we.eBH used in this test'
            }
            
            x.all$`-Log10(FDR)` <- -log10(x.all$wi.eBH)
            x.all$Taxa = row.names(x.all)
            x.all = x.all[order(x.all$wi.ep),]
            write.csv(x.all,paste0("Ga_",gest_day_collection_quantile_list[cc],"_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Differential_abundance.csv'))
            
            # draw figure
            das <- x.all[x.all$`-Log10(FDR)` >= -log10(pvalue_th),]
            
            das$Species <- row.names(das)
            das <- das[order(das$diff.btw),] 
            
            metadata = as.factor(metadata)
            lev = levels(metadata)
            
            if (order_reverse == T) {
              das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
            } else {
              das$Outcome <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
              
            }
            
            das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
            
            das$diff.btw[das$diff.btw == Inf] = 10
            das$diff.btw[das$diff.btw == -Inf] = -10
            das$diff.btw[das$diff.btw <= -10] = -10
            das$diff.btw[das$diff.btw >= 10] = 10
            
            
            theme_set(theme_bw())  
            
            ggplot(das, aes(Species, diff.btw)) +
              geom_line() +  
              geom_hline(yintercept=1, linetype='dotted', col = 'grey') +
              geom_hline(yintercept=-1, linetype='dotted', col = 'grey') +
              geom_point(aes(col=Outcome, size=`-Log10(FDR)`)) +
              scale_size_continuous(range  = c(0.1, 5), limits = c(0, 25), 
                                    breaks = c(2, 4, 6, 8, 10))+ 
              coord_flip() +          # convert x y axis
              labs(x = 'Taxa', y = "Median difference in clr values")+ 
              ggtitle(paste0("Race: ", aa, '\n', "Age: ", Age_quantile_list[bb],'\n',
                             "Gestational age: ", gest_day_collection_quantile_list[cc])) +
              scale_color_brewer(palette="Set1")+
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(hjust = 0.5, size = 7))
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Differential_abundance.pdf'), width=4, height=1+nrow(das)*0.1)
            
            
          }
          
          # t-SNE
          {
            Padonis = adonis2(as.matrix(t(reads_table_x)) ~ Outcome, data = metadata_x,method = "bray", parallel = 8)
            P = Padonis$`Pr(>F)`; P = P[1]
            if (P < 0.001) {
              P = formatC(P, format = "e", digits = 3)
            } else {
              P = formatC(P, format = "f", digits = 3)
            }
            
            reads_table = as.matrix(decostand(as.data.frame(t(reads_table_x)),method = 'rclr'))
            x = prcomp(as.data.frame(t(reads_table)), center = TRUE, scale. = TRUE)
            
            pic <- x$rotation
            pic <- data.frame(pic,metadata_x)
            
            colnames(pic)[1:2] <- c('X1','X2')
            pic$Outcome
            ggplot(pic, aes(X1, X2,color = Outcome))  +
              geom_point(size=0.5) +
              xlab(paste0("PCoA1")) +
              ylab(paste0("PCoA2"))+
              stat_ellipse(type = "t") + 
              coord_fixed()+ 
              scale_color_brewer(palette="Set1")+
              ggtitle(paste0("Race: ", aa, '\n', "Age: ", Age_quantile_list[bb],'\n',
                             "Gestational age: ", gest_day_collection_quantile_list[cc],'\n',
                             'P = ',P)) +
              theme(
                axis.title.x = element_text( size=12),
                axis.title.y = element_text( size=12),
                legend.text = element_text(size=12),
                legend.title = element_text(size=12),
                plot.title = element_text(hjust = 0.5, size = 12)
              ) 
            
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_",Age_quantile_list[bb],"_", str_replace(aa,'\\/','_'),'tsne_Outcome.pdf'), width=4, height=6)
          }
          
          # alpha
          {
            reads_table = reads_table_x
            metadata = metadata_x$Outcome
            
            # alpha diversity 
            reads_table = t(reads_table)
            reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
            reads_table <- reads_table$otu.tab.rff
            reads_table <- as.data.frame(reads_table)
            
            alpha <- data.frame(diversity(reads_table))
            alpha$Evenness <- alpha$diversity.reads_table./log(specnumber(reads_table))
            alpha$Ovserved_phage <- colSums(t(reads_table) != 0)
            colnames(alpha) = c('Shannon','Evenness','Ovserved_taxa')
            alpha$Outcome = metadata
            x = kruskal.test(Shannon ~ Outcome,data = alpha)
            P=x$p.value
            if (P < 0.001) {
              P = formatC(P, format = "e", digits = 3)
            } else {
              P = formatC(P, format = "f", digits = 3)
            }
            
            ggplot(data=alpha, aes(x=Outcome, y=Shannon)) +geom_violin(trim=T, aes(fill=Outcome))+
              geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
              geom_jitter(size = 0.1, alpha= 0.25)+theme_bw()+
              theme_classic()+
              scale_fill_brewer(palette="Set1")+
              labs(y = 'Shannon index')+
              ggtitle(paste0("Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n',
                             "Gestational age: ", gest_day_collection_quantile_list[cc],'\n',
                             "n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            ggsave(paste0("Ga_",gest_day_collection_quantile_list[cc],"_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Shannon.pdf'),width=2.5, height=2)
          }
          
        }
      }
    }
  }
  
  
}




##### age #####
setwd('/Desktop/1_KH/VMB_PTB/results/age/')
load("~/Desktop/1_KH/VMB_PTB/others/Untitled.RData")

# bar plot
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) & metadata_all$Race %in% c('B/AA', 'White', 'Asian') &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata_all_3 = metadata_all[keep,]
  table(metadata_all_3$BioProject)
  reads_table_all_3 = reads_table_all_batch_rm[,keep]
  
  race = unique(metadata_all_3$Race)
  
  x = metadata_all[metadata_all$Outcome == "PTB" & !is.na(metadata_all$Age),]
  Age_quantile = quantile(x$Age[!is.na(x$Age)], probs = seq(0, 1, length.out = 4)); Age_quantile = as.numeric(Age_quantile)
  Age_quantile[1] = Age_quantile[1]-1
  Age_quantile
  Age_quantile_list = c("(17 - 25]", "(25- 32]", "(32- 42]")
  
  for (aa in race) {
      keep = metadata_all_3$Race == aa; sum(keep)
      
      metadata_all_2 = metadata_all_3[keep,]
      reads_table_all_2 = reads_table_all_3[,keep]
      
    for (bb in 1:3) {
      keep = metadata_all_2$Age > Age_quantile[bb] & metadata_all_2$Age <= Age_quantile[bb+1]; sum(keep)
      metadata_x = metadata_all_2[keep,]; reads_table_x = reads_table_all_2[,keep]
        
      # case match
      {
        # case matched by the same bioproject, race, similar age, ga_at_collection, and different participants
        {
          metadata_PTB = metadata_x[metadata_x$Outcome == "PTB",]
          metadata_Term = metadata_PTB; metadata_Term[1:nrow(metadata_Term),1:ncol(metadata_Term)] = NA
          
          if (nrow(metadata_PTB) < 10) {next}
          
          metadata_2 = metadata_x[metadata_x$Outcome == "Term",]
          
          for (a in 1: nrow(metadata_PTB)) {
            
            n =1
            keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
              metadata_2$Age <= metadata_PTB$Age[a] +n & 
              metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
              metadata_2$BioProject == metadata_PTB$BioProject[a] & 
              metadata_2$Race == metadata_PTB$Race[a] & 
              metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
              metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
            while (sum(keep[!is.na(keep)]) == 0 & n < 3) {
              n = n+1 
              keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
                metadata_2$Age <= metadata_PTB$Age[a] +n & 
                metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
                metadata_2$BioProject == metadata_PTB$BioProject[a] & 
                metadata_2$Race == metadata_PTB$Race[a] & 
                metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
                metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
            }
            
            if (n < 3) {
              metadata_Term[a,] = metadata_2[which(keep)[1],]
              metadata_2 = metadata_2[-which(keep)[1],]
            }
            
          }
          
          sum(is.na(metadata_Term$ParticipantID))
          
          keep = !is.na(metadata_Term$Run)
          metadata_Term = metadata_Term[keep,]; metadata_PTB = metadata_PTB[keep,]
          
          if (nrow(metadata_Term) < 6) {next}
          
          data = data.frame(Age = c(metadata_PTB$Age,metadata_Term$Age), 
                            Outcome = c(rep('PTB',length(metadata_PTB$Age)), rep('Term',length(metadata_Term$Age))))
          
          P = kruskal.test(Age ~ Outcome,data = data)
          if (P$p.value < 0.001) {
            P = formatC(P$p.value, format = "e", digits = 3)
          } else {
            P = formatC(P$p.value, format = "f", digits = 3)
          }
          
          ggplot(data, aes(x=Outcome, y=Age)) + geom_violin(trim=T, aes(fill = Outcome))+
            geom_boxplot(fill='white', color="black", width=0.1) + 
            geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
            scale_fill_brewer(palette="Set1")+
            ggtitle(paste0("Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          
          ggsave(paste0("Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Age.pdf'), height = 2,width = 2)
          
          data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                            Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
          
          P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
          if (P$p.value < 0.001) {
            P = formatC(P$p.value, format = "e", digits = 3)
          } else {
            P = formatC(P$p.value, format = "f", digits = 3)
          }
          ggplot(data, aes(x=Outcome, y=Gestational_age_in_sample_collection)) + geom_violin(trim=T, aes(fill = Outcome))+
            geom_boxplot(fill='white', color="black", width=0.1) + 
            geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
            labs(y = 'Gestational age at sample collection')+
            scale_fill_brewer(palette="Set1")+
            ggtitle(paste0("Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  plot.title = element_text(size = 7))
          ggsave(paste0("Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
          
          metadata = data.frame(ParticipantID = c(metadata_PTB$Patient_ID,metadata_Term$Patient_ID), 
                                BioProject = c(metadata_PTB$BioProject,metadata_Term$BioProject), 
                                Run = c(metadata_PTB$Run,metadata_Term$Run), 
                                gest_day_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection),
                                Age = c(metadata_PTB$Age,metadata_Term$Age), 
                                Race = c(metadata_PTB$Race,metadata_Term$Race), 
                                Outcome = c(metadata_PTB$Outcome,metadata_Term$Outcome))
          
          reads_table = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata$Run]
          reads_table = reads_table[metadata$Run]
          
          metadata_x = metadata; reads_table_x = reads_table
          rm(data, metadata,metadata_2, metadata_PTB, metadata_Term, reads_table)
        }
      }
      table(metadata_x$Race,metadata_x$Outcome)
      table(metadata_x$BioProject,metadata_x$Outcome)
      
      x = data.frame(Outcome = metadata_x$Outcome)
      data = barplot(reads_table_x, x, type_th = 0.1, taxa_num = 9, VMB = T)
      data$p1
      ggsave(paste0("Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'_bar_1.pdf'), height = 3,width = 8)
      
      data$p2
      ggsave(paste0("Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'_bar_2.pdf'), height = 6,width = 6)
      
    }
  }
  
}

# association linear model for taxa and age
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) & metadata_all$Race %in% c('B/AA', 'White', 'Asian') &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata_all_3 = metadata_all[keep,]
  table(metadata_all_3$BioProject)
  reads_table_all_3 = reads_table_all_batch_rm[,keep]
  reads_table = as.data.frame(t(reads_table_all_3))
  reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
  reads_table_all_3 = as.data.frame(t(reads_table)); dim(reads_table_all_3)
    
  metadata_all_3 = data.frame(Age = metadata_all_3$Age) 
  metadata_all_3 =  as.data.frame(t(metadata_all_3)); dim(metadata_all_3)
  colnames(metadata_all_3) = colnames(reads_table_all_3)
  
  reads_table = rbind(metadata_all_3, reads_table_all_3)
  data = newwork_rcorr(reads_table, normalization_method = NA, type = 'spearman', pvalue = 0.05, cor_parameter= 0, 
                            style = 1, bar_max = 2, bar_min = -2, pheatmap_fontsize = 5, treeheight = 50, alpha = 0.05,
                            FDR = T, ABS = F)
  data$p
  
  reads_table_all_3 = reads_table_all_3[row.names(reads_table_all_3) %in% c('Lactobacillus_crispatus','Lactobacillus_gasseri',
                                                                               'Lactobacillus_jensenii', 'Lactobacillus_iners',
                                                                            'Lachnospiraceae_BVAB1','Gardnerella_vaginalis'),]
  for (a in seq(nrow(reads_table_all_3))) {

      data = data.frame(Taxa = as.numeric(reads_table_all_3[a,]), Age = as.numeric(as.character(metadata_all_3)))
      
      x = lm(data); x = summary(x); x = x$coefficients[2,4]
      x = formatC(x, format = "e", digits = 2)
      y = str_replace_all(row.names(reads_table_all_3)[a],'Lactobacillus_', 'L. ')
      y = str_replace_all(row.names(reads_table_all_3)[a],'Lachnospiraceae_BVAB1', 'BVAB1')
      y = str_replace_all(row.names(reads_table_all_3)[a],'Gardnerella_vaginalis', 'G. vaginalis')
      
      ggplot(data, aes(x = Age, y = Taxa)) +
        geom_point(size = 0.3)+
        labs(x = 'Age', y = paste0('CLR value of ',y))+ 
        ggtitle(paste0('P = ',x)) +
        geom_smooth(method='lm', formula= y~x)+ 
        theme(axis.title = element_text(size = 7, face = "bold"), 
              axis.text = element_text(size = 7, face = "bold"), 
              plot.title = element_text(size = 7, face = "bold"), 
              legend.text = element_text(size = 7, face = "bold"))+theme_bw()
      ggsave(paste0('Age_',row.names(reads_table_all_3)[a], '.pdf'), width = 2.5, height = 2.5)

    
    
  }
  
}

# mixed-effects model
{
  keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
    !is.na(metadata_all$Race) & metadata_all$Race %in% c('B/AA', 'White', 'Asian') &
    !is.na(metadata_all$gest_day_collection); sum(keep)
  
  metadata_all_2 = metadata_all[keep,]
  reads_table_all_2 = reads_table_all_batch_rm[,keep]
  
  metadata_x = metadata_all_2; reads_table_x = reads_table_all_2
  
  # case match
  {
    # case matched by the same bioproject, race, similar age, ga_at_collection, and different participants
    {
      metadata_PTB = metadata_x[metadata_x$Outcome == "PTB",]
      metadata_Term = metadata_PTB; metadata_Term[1:nrow(metadata_Term),1:ncol(metadata_Term)] = NA
      
      if (nrow(metadata_PTB) < 6) {next}
      
      metadata_2 = metadata_x[metadata_x$Outcome == "Term",]
      
      for (a in 1: nrow(metadata_PTB)) {
        
        n =1
        keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
          metadata_2$Age <= metadata_PTB$Age[a] +n & 
          metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
          metadata_2$BioProject == metadata_PTB$BioProject[a] & 
          metadata_2$Race == metadata_PTB$Race[a] & 
          metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
          metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
        while (sum(keep[!is.na(keep)]) == 0 & n < 3) {
          n = n+1 
          keep = metadata_2$Age >= metadata_PTB$Age[a] -n & 
            metadata_2$Age <= metadata_PTB$Age[a] +n & 
            metadata_2$Patient_ID != metadata_PTB$Patient_ID[a] & 
            metadata_2$BioProject == metadata_PTB$BioProject[a] & 
            metadata_2$Race == metadata_PTB$Race[a] & 
            metadata_2$gest_day_collection >= metadata_PTB$gest_day_collection[a] -n*1 & 
            metadata_2$gest_day_collection <= metadata_PTB$gest_day_collection[a] +n*1; sum(keep[!is.na(keep)])
        }
        
        if (n < 3) {
          metadata_Term[a,] = metadata_2[which(keep)[1],]
          metadata_2 = metadata_2[-which(keep)[1],]
        }
        
      }
      
      sum(is.na(metadata_Term$ParticipantID))
      
      keep = !is.na(metadata_Term$Run)
      metadata_Term = metadata_Term[keep,]; metadata_PTB = metadata_PTB[keep,]
      
      data = data.frame(Age = c(metadata_PTB$Age,metadata_Term$Age), 
                        Outcome = c(rep('PTB',length(metadata_PTB$Age)), rep('Term',length(metadata_Term$Age))))
      
      if (nrow(data) < 6) {next}
      
      P = kruskal.test(Age ~ Outcome,data = data)
      if (P$p.value < 0.001) {
        P = formatC(P$p.value, format = "e", digits = 3)
      } else {
        P = formatC(P$p.value, format = "f", digits = 3)
      }
      
      ggplot(data, aes(x=Outcome, y=Age)) + geom_violin(trim=T, aes(fill = Outcome))+
        geom_boxplot(fill='white', color="black", width=0.1) + 
        geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
        scale_fill_brewer(palette="Set1")+
        ggtitle(paste0('P = ',P)) +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              plot.title = element_text(size = 7))
      
      ggsave(paste0("age_",'Age.pdf'), height = 2,width = 2)
      
      data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                        Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
      
      P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
      if (P$p.value < 0.001) {
        P = formatC(P$p.value, format = "e", digits = 3)
      } else {
        P = formatC(P$p.value, format = "f", digits = 3)
      }
      ggplot(data, aes(x=Outcome, y=Gestational_age_in_sample_collection)) + geom_violin(trim=T, aes(fill = Outcome))+
        geom_boxplot(fill='white', color="black", width=0.1) + 
        geom_jitter(size = 0.2, alpha = 0.1)+theme_bw()+
        labs(y = 'Gestational age at sample collection')+
        scale_fill_brewer(palette="Set1")+
        ggtitle(paste0('P = ',P)) +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              plot.title = element_text(size = 7))
      ggsave(paste0("age_",'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
      
      metadata = data.frame(ParticipantID = c(metadata_PTB$Patient_ID,metadata_Term$Patient_ID), 
                            BioProject = c(metadata_PTB$BioProject,metadata_Term$BioProject), 
                            Run = c(metadata_PTB$Run,metadata_Term$Run), 
                            gest_day_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection),
                            Age = c(metadata_PTB$Age,metadata_Term$Age), 
                            Race = c(metadata_PTB$Race,metadata_Term$Race), 
                            Outcome = c(metadata_PTB$Outcome,metadata_Term$Outcome))
      
      reads_table = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata$Run]
      reads_table = reads_table[metadata$Run]
      
      metadata_x = metadata; reads_table_x = reads_table
      rm(data, metadata,metadata_2, metadata_PTB, metadata_Term, reads_table)
    }
  }
  table(metadata_x$Race,metadata_x$Outcome)
  table(metadata_x$BioProject,metadata_x$Outcome)

  reads_table = as.data.frame(t(reads_table_x))
  reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
  reads_table_x = as.data.frame(t(reads_table)); dim(reads_table_x)
  
  metadata_x = data.frame(Age = metadata_x$Age, Outcome = metadata_x$Outcome) 
  metadata_x =  as.data.frame(t(metadata_x)); dim(metadata_x)
  colnames(metadata_x) = colnames(reads_table_x)
  metadata_x = as.data.frame(t(metadata_x)); reads_table_x = as.data.frame(t(reads_table_x)); 
  
  library(lme4)     # For mixed-effects modeling
  library(car)      # For ANOVA
  
  microbiota_data = cbind(metadata_x, reads_table_x)
  # keep = microbiota_data$Outcome =='PTB'; microbiota_data$Outcome[keep] = 1;microbiota_data$Outcome[!keep] = 0;
  microbiota_data$Age = as.numeric(as.character(microbiota_data$Age))
  #microbiota_data$Outcome = as.numeric(as.character(microbiota_data$Outcome))
  
  # Initialize a results data frame
  results <- data.frame(
    Bacterium = character(),
    Coef_Age = numeric(),
    Coef_AgePreterm = numeric(),
    Pval_Age = numeric(),
    Pval_Interaction = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Fit models for each bacterium
  for (bacterium in names(microbiota_data)[-c(1,2)]) {
    tryCatch({
      model <- lm(as.formula(paste(bacterium, "~ Age * Outcome")), data = microbiota_data)
    }, error = function(e) {
      cat("Error occurred at iteration - skipping. Error message:", e$message, "\n")
    })

    summary_model <- summary(model)
    
    # Extract coefficients and p-values
    results <- rbind(
      results,
      data.frame(
        Bacterium = bacterium,
        Coef_Age = coef(summary_model)["Age", "Estimate"],
        Coef_AgePreterm = coef(summary_model)["Age:OutcomeTerm", "Estimate"],
        Pval_Age = coef(summary_model)["Age", "Pr(>|t|)"],
        Pval_Interaction = coef(summary_model)["Age:OutcomeTerm", "Pr(>|t|)"]
      )
    )
    
  }
  
  # Apply multiple testing correction
  results$Pval_Age_FDR <- p.adjust(results$Pval_Age, method = "fdr")  # Benjamini-Hochberg
  results$Pval_Interaction_FDR <- p.adjust(results$Pval_Interaction, method = "fdr")
  results = results[!duplicated(results$Bacterium),]
  
  # Filter significant results
  significant_results <- results %>%
    filter(Pval_Age_FDR < 0.05 & Pval_Interaction_FDR < 0.05)

  for (bacterium in names(microbiota_data)[-c(1,2)]) {
    if (bacterium %in% c('Lactobacillus_crispatus','Lactobacillus_gasseri',
                         'Lactobacillus_jensenii', 'Lactobacillus_iners',
                         'Lachnospiraceae_BVAB1','Gardnerella_vaginalis')) {
      
      x = as.numeric(as.character(microbiota_data[,which(colnames(microbiota_data) == bacterium)]))
      data <- data.frame(Taxa_Abundance = x,
                         Age = microbiota_data$Age,
                         Outcome = microbiota_data$Outcome)
      FDR = results$Pval_Interaction_FDR[which(results$Bacterium == bacterium)]; FDR = formatC(FDR, format = "e", digits = 3)
      y = results$Coef_AgePreterm[which(results$Bacterium == bacterium)]; y = formatC(y, format = "f", digits = 3)
      y1 = results$Coef_Age[which(results$Bacterium == bacterium)]; y1 = formatC(y1, format = "f", digits = 3)
      y2 = results$Pval_Age[which(results$Bacterium == bacterium)]; y2 = formatC(y2, format = "e", digits = 3)
      
      ggplot(data, aes(x = Age, y = Taxa_Abundance, color = Outcome)) +
        geom_point(alpha = 0.5, size = 0.3) +  # Scatter plot of raw data
        geom_smooth(method = "lm", se = TRUE) +  # Add regression lines with confidence intervals
        scale_color_brewer(palette="Set1")+
        labs(
          title = paste0(str_replace_all(bacterium,'_',' '), '\n',
          #               'Coefficient Age = ', y1,'\n',
          #               'FDR for Age = ', y2, '\n',
                         'Coefficient Age:Term interaction = ', y,'\n',
                         'FDR for interaction = ', FDR),
          x = "Age",
          y = "CLR value of a bacterial taxon"
        ) +
        theme_minimal()
      ggsave(paste0('Interaction_',bacterium, '.pdf'), width = 4, height = 3)
    }
  }
  write.csv(results,'Age_preterm_interaction.csv')
}



