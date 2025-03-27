library(readxl)
library(ggplot2)
library(stringr)
library(parallel)
library(dplyr)
library(devtools)
library(microbiomeMarker)
library(phyloseq)

library(ipw)
library(survey)
library(tableone)
setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/')

# Megasphaera_sp._oral_taxon_841
##### get reads table and metadata #####
setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/data')
{
  reads_table_all = read.csv('reads_table.csv', row.names = 1)

  # metadata cleansing
  {
    metadata_all = read.csv('/Users/binzhu/Desktop/1_KH/VMB_PTB/data/metadata_all.csv', row.names = 1)

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
  ggsave('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/overall_profiles/Sample total reads_ori.pdf', width = 3, height = 2)
  
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
          ggsave(paste0('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/Metadata_profile/Metadata_profile_', x[a],'.pdf'), width = 3, height = 2)
        }
      } else {
        data <- as.data.frame(table(metadata_all[,a]))
        data = data[order(data$Var1),]
        ggplot(data, aes(x="", y=Freq, fill=Var1)) +
          geom_bar(stat="identity", width=1, color="white") + 
          coord_polar("y", start=0) +
          theme_void() + 
          scale_fill_manual(values= mycolors[[a]], name = paste0('Total number = ', sum(data$Freq), '\n',x[a])) 
        ggsave(paste0('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/Metadata_profile/Metadata_profile_', x[a],'.pdf'), width = 6, height = 3.5)
      }
    }
  }
}
rm(x,a, taxonomy, data,data_2, file_list,keep,n,taxa_list,reads_table,y)




##### reads and metadata profiles #################
setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/overall_profiles')
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
  ggsave('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/overall_profiles/Sample total reads.pdf', width = 3, height = 2)
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
          ggsave(paste0('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/Metadata_profile/Metadata_profile_', x[a],'.pdf'), width = 3, height = 2)
        }
      } else {
        data <- as.data.frame(table(metadata_all[,a]))
        data = data[order(data$Var1),]
        ggplot(data, aes(x="", y=Freq, fill=Var1)) +
          geom_bar(stat="identity", width=1, color="white") + 
          coord_polar("y", start=0) +
          theme_void() + 
          scale_fill_manual(values= mycolors[[a]], name = paste0('Total number = ', sum(data$Freq), '\n',x[a])) 
        ggsave(paste0('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/Metadata_profile/Metadata_profile_', x[a],'.pdf'), width = 6, height = 3.5)
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
setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/overall_profiles/')
load("~/Desktop/1_KH/VMB_PTB/others/Untitled.RData")
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
    sc_si = sc_si *10
    sc_sp <- scores(dbRDA, display="species", choices=c(1,2), scaling=1)
    sc_sp = sc_sp/2
    sc_bp <- scores(dbRDA, display="bp", choices=c(1, 2), scaling=1)
    sc_bp = sc_bp*55
    
    #      n = which(row.names(sc_sp) %in% c('Streptococcus_mitis','Streptococcus_infantis',
    #                                        'Escherichia_coli'))
    #      sc_sp[n,] = sc_sp[n,]/2
    
    rownames(sc_sp) = str_replace_all(rownames(sc_sp),'_'," ")
    
    write.csv(sc_sp, 'dbRDA_sc_sp.csv')
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
    
    {
      bg_color <- adjustcolor("#f2bd33", alpha.f = 0.5) # 0.5 makes it 50% transparent
      
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

      # add points for species scores
      points(sc_sp, 
             pch = 22, # set shape (here, square with a fill colour)
             col = "black",
             bg = bg_color, 
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
setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/case_match/')
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

  IPW = data.frame(Comparison = NA, IPW = NA); n_IPW = 0
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
          
          # Inverse Probability Weighting (IPW)
          {
            alpha = cbind(alpha, metadata_x)
            alpha = alpha[,-c(5,6,7,10,11)]
            alpha = alpha[,-c(2,3)]
            
            alpha$Outcome <- as.factor(alpha$Outcome)
            
            # Estimate Propensity Scores
            ps_model <- glm(Outcome ~ Age, 
                            family = binomial, data = alpha)
            
            # Obtain predicted probabilities (propensity scores)
            alpha$ps <- predict(ps_model, type = "response")
            
            # Compute Inverse Probability Weights (IPW)
            alpha$ipw <- ifelse(alpha$Outcome == "PTB", 1 / alpha$ps, 1 / (1 - alpha$ps))
            
            # Perform Weighted Regression
            weighted_data <- svydesign(ids = ~1, weights = ~ipw, data = alpha)

            model <- svyglm(Shannon ~ Outcome, design = weighted_data)
            x = summary(model)
            
            n_IPW = n_IPW+1
            IPW[n_IPW, 1] = paste0("Race:",aa, "; Age:",Age_quantile_list[bb])
            IPW[n_IPW, 2] = x$coefficients[2, 4]
          }
        }
        
      }
    }
  }
  write.csv(IPW, 'IPW.csv')
  
}

# analysis_all_dif_race_dif_ages_dif_cohort_correct
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
  
  cohorts = unique(metadata_all_3$BioProject)
  for (cohort in cohorts) {
    keep = metadata_all_3$BioProject == cohort; sum(keep)
    if (sum(keep) < 10) {next}
    metadata_all_22 = metadata_all_3[keep,]; reads_table_all_22 = reads_table_all_3[,keep]
    
    for (bb in 1:3) {
      keep = metadata_all_22$Age > Age_quantile[bb] & metadata_all_22$Age <= Age_quantile[bb+1]; sum(keep)
      if (sum(keep) < 10) {next}
      metadata_all_2 = metadata_all_22[keep,]; reads_table_all_2 = reads_table_all_22[,keep]
      
      for (aa in race) {
        keep = metadata_all_2$Race == aa; sum(keep)
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
              ggtitle(paste0("Cohort: ", cohort, '\n',"Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            
            ggsave(paste0(cohort,"_Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Age.pdf'), height = 2,width = 2)
            
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
              ggtitle(paste0("Cohort: ", cohort, '\n',"Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            ggsave(paste0(cohort,"_Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
            
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
        if (nrow(metadata_x) < 10) {next}
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
            write.csv(x.all,paste0(cohort,"_Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Differential_abundance.csv'))
            
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
              ggtitle(paste0("Cohort: ", cohort, '\n',"Race: ", aa, '\n', "Age: ", Age_quantile_list[bb])) +
              scale_color_brewer(palette="Set1")+
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(hjust = 0.5, size = 7))
            ggsave(paste0(cohort,"_Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Differential_abundance.pdf'), width=4, height=1+nrow(das)*0.1)
            
            
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
              ggtitle(paste0("Cohort: ", cohort, '\n',"Race: ", aa, '\n', "Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
              theme(
                axis.title.x = element_text( size=12),
                axis.title.y = element_text( size=12),
                legend.text = element_text(size=12),
                legend.title = element_text(size=12),
                plot.title = element_text(hjust = 0.5, size = 12)
              ) 
            
            ggsave(paste0(cohort,"_Race_",Age_quantile_list[bb],"_", str_replace(aa,'\\/','_'),'tsne_Outcome.pdf'), width=4, height=6)
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
              ggtitle(paste0("Cohort: ", cohort, '\n',"Race: ", aa, '\n',"Age: ", Age_quantile_list[bb],'\n',"n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            ggsave(paste0(cohort,"_Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Shannon.pdf'),width=2.5, height=2)
          }
          
        }
      }
    }
  }
  

  
}

# analysis_all_dif_ages_dif_cohort_correct
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
  
  cohorts = unique(metadata_all_3$BioProject)
  for (cohort in cohorts) {
    keep = metadata_all_3$BioProject == cohort; sum(keep)
    if (sum(keep) < 10) {next}
    metadata_all_22 = metadata_all_3[keep,]; reads_table_all_22 = reads_table_all_3[,keep]
    
    for (bb in 1:3) {
      keep = metadata_all_22$Age > Age_quantile[bb] & metadata_all_22$Age <= Age_quantile[bb+1]; sum(keep)
      if (sum(keep) < 10) {next}
      metadata_all_2 = metadata_all_22[keep,]; reads_table_all_2 = reads_table_all_22[,keep]
        
        metadata_x = metadata_all_2
        reads_table_x = reads_table_all_2
        
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
              ggtitle(paste0("Cohort: ", cohort, '\n',"Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            
            ggsave(paste0(cohort,"_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Age.pdf'), height = 2,width = 2)
            
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
              ggtitle(paste0("Cohort: ", cohort, '\n',"Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            ggsave(paste0(cohort,Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
            
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
        if (nrow(metadata_x) < 10) {next}
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
            write.csv(x.all,paste0(cohort,Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Differential_abundance.csv'))
            
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
              ggtitle(paste0("Cohort: ", cohort, '\n', "Age: ", Age_quantile_list[bb])) +
              scale_color_brewer(palette="Set1")+
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(hjust = 0.5, size = 7))
            ggsave(paste0(cohort,Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Differential_abundance.pdf'), width=4, height=1+nrow(das)*0.1)
            
            
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
              ggtitle(paste0("Cohort: ", cohort, '\n',"Age: ", Age_quantile_list[bb],'\n','P = ',P)) +
              theme(
                axis.title.x = element_text( size=12),
                axis.title.y = element_text( size=12),
                legend.text = element_text(size=12),
                legend.title = element_text(size=12),
                plot.title = element_text(hjust = 0.5, size = 12)
              ) 
            
            ggsave(paste0(cohort,Age_quantile_list[bb],"_", str_replace(aa,'\\/','_'),'tsne_Outcome.pdf'), width=4, height=6)
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
              ggtitle(paste0("Cohort: ", cohort, '\n',"Age: ", Age_quantile_list[bb],'\n',"n = ", nrow(alpha)/2, " × 2", '\n','P = ',P)) +
              theme(axis.title = element_text(size = 7), 
                    axis.text = element_text(size = 7), 
                    legend.text = element_text(size = 7), 
                    legend.title = element_text(size = 7),
                    plot.title = element_text(size = 7))
            ggsave(paste0(cohort,Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Shannon.pdf'),width=2.5, height=2)
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
  
  IPW = data.frame(Comparison = NA, IPW = NA); n_IPW = 0
  
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
            
            alpha = cbind(alpha, metadata_x)
            alpha = alpha[,-c(5,6,7,10,11)]
            alpha = alpha[,-c(2,3)]
            
            # Inverse Probability Weighting (IPW)
            {
              alpha$Outcome <- as.factor(alpha$Outcome)
              
              # Estimate Propensity Scores
              ps_model <- glm(Outcome ~ gest_day_collection + Age, 
                              family = binomial, data = alpha)
              
              # Obtain predicted probabilities (propensity scores)
              alpha$ps <- predict(ps_model, type = "response")
              
              # Compute Inverse Probability Weights (IPW)
              alpha$ipw <- ifelse(alpha$Outcome == "PTB", 1 / alpha$ps, 1 / (1 - alpha$ps))
              
              # Perform Weighted Regression
              weighted_data <- svydesign(ids = ~1, weights = ~ipw, data = alpha)
              
              model <- svyglm(Shannon ~ Outcome, design = weighted_data)
              x = summary(model)
              
              n_IPW = n_IPW+1
              IPW[n_IPW, 1] = paste0("Ga_",gest_day_collection_quantile_list[cc],"_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Shannon')
              IPW[n_IPW, 2] = x$coefficients[2, 4]
            }
          }
          
        }
      }
    }
  }
  write.csv(IPW, 'IPW.csv')
  
}




##### age #####
setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/age/')
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

# mixed-effects model # different cohort
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/age/cohort')
  cohorts = unique(metadata_all$BioProject)
  table(metadata_all$BioProject)
  
  for (cohort in cohorts) {
    
    keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
      !is.na(metadata_all$Race) & metadata_all$Race %in% c('B/AA', 'White', 'Asian') &
      !is.na(metadata_all$gest_day_collection) & metadata_all$BioProject == cohort; sum(keep)
    
    if (sum(keep) < 20) {next}
    
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
        
        
        data = data.frame(Gestational_age_in_sample_collection = c(metadata_PTB$gest_day_collection,metadata_Term$gest_day_collection), 
                          Outcome = c(rep('PTB',length(metadata_PTB$gest_day_collection)), rep('Term',length(metadata_Term$gest_day_collection))))
        
        P = kruskal.test(Gestational_age_in_sample_collection ~ Outcome,data = data)
        if (P$p.value < 0.001) {
          P = formatC(P$p.value, format = "e", digits = 3)
        } else {
          P = formatC(P$p.value, format = "f", digits = 3)
        }
        
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
    
    
    if (nrow(metadata_x) < 20) {next}
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
        FDR = results$Pval_Interaction[which(results$Bacterium == bacterium)]; FDR = formatC(FDR, format = "e", digits = 3)
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
                           'P-value for interaction = ', FDR),
            x = "Age",
            y = "CLR value of a bacterial taxon"
          ) +
          theme_minimal()
        ggsave(paste0(cohort,'_Interaction_',bacterium, '.pdf'), width = 4, height = 3)
      }
    }
    write.csv(results,paste0(cohort,'_Age_preterm_interaction.csv'))
  }

}

# bar plot # different cohort
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/age/cohort')
  cohorts = unique(metadata_all$BioProject)
  table(metadata_all$BioProject)
  
  for (cohort in cohorts) {
    keep = !is.na(metadata_all$Outcome) & !is.na(metadata_all$Age) &
      !is.na(metadata_all$Race) & metadata_all$Race %in% c('B/AA', 'White', 'Asian') &
      !is.na(metadata_all$gest_day_collection) & metadata_all$BioProject == cohort; sum(keep)
    
    if (sum(keep) < 20) {next}
    
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
      if (sum(keep) < 20) {next}
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
            
            ggsave(paste0(cohort, "_Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Age.pdf'), height = 2,width = 2)
            
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
            ggsave(paste0(cohort, "_Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'Gestational_age_in_sample_collection.pdf'), height = 2,width = 2)
            
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
        
        if (nrow(metadata_x) < 20) {next}
        
        table(metadata_x$Race,metadata_x$Outcome)
        table(metadata_x$BioProject,metadata_x$Outcome)
        
        x = data.frame(Outcome = metadata_x$Outcome)
        data = barplot(reads_table_x, x, type_th = 0.1, taxa_num = 9, VMB = T)
        data$p1
        ggsave(paste0(cohort, "_Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'_bar_1.pdf'), height = 3,width = 8)
        
        data$p2
        ggsave(paste0(cohort, "_Race_",Age_quantile_list[bb],"_",str_replace(aa,'\\/','_'),'_bar_2.pdf'), height = 6,width = 6)
        
      }
    }
  }

  
}

##### ML #####
library(pROC)
library(doSNOW)
library(parallel)
library(matrixStats)
library(erer) # export list
library(ggplot2)
library(randomForest)
library(forestError)

library(caret)
library(tidyr) # for 'gather'

case_number = vector()
# all
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race); sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[10] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )
        
        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
      }
      
      resamps_all <- resamples(list(RF = model_rf,
                                  GLM = model_glm,
                                  GLMB = model_glmboost,
                                  SVM = model_svm,
                                  KNN = model_knn,
                                  XGB_tree = model_xgbtree,
                                  C5.0Cost = model_C5.0Cost))
      resamps_all
      summary(resamps_all)
      
      pdf('Kappa_ML_methods_all.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_all, metric = "Kappa")
      dev.off()
      
      pdf('Accuracy_ML_methods_all.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_all, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}

# black
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "B/AA"; sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
        
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[1] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )

        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
      }
      
      resamps_b <- resamples(list(RF = model_rf,
                                GLM = model_glm,
                                GLMB = model_glmboost,
                                SVM = model_svm,
                                KNN = model_knn,
                                XGB_tree = model_xgbtree,
                                C5.0Cost = model_C5.0Cost))
      resamps_b
      summary(resamps_b)
      
      pdf('Kappa_ML_methods_Black.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_b, metric = "Kappa")
      dev.off()
      
      pdf('Accuracy_ML_methods_Black.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_b, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}

# white
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "White"; sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[2] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )
        
        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
        model_cforest <- train(
          output_all ~ .,
          data = data,
          method = 'cforest',
          trControl = fitControl
        )
        
      }
      
      resamps_w <- resamples(list(RF = model_rf,
                                GLM = model_glm,
                                GLMB = model_glmboost,
                                SVM = model_svm,
                                KNN = model_knn,
                                XGB_tree = model_xgbtree,
                                CFOREST = model_cforest,
                                C5.0Cost = model_C5.0Cost))
      resamps_w
      summary(resamps_w)
      
      pdf('Kappa_ML_methods_White.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_w, metric = "Kappa")
      dev.off()
      
      
      pdf('Accuracy_ML_methods_White.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_w, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}

# Asian
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "Asian"; sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[3] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )
        
        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
        model_cforest <- train(
          output_all ~ .,
          data = data,
          method = 'cforest',
          trControl = fitControl
        )
        
      }
      
      resamps_a <- resamples(list(RF = model_rf,
                                GLM = model_glm,
                                GLMB = model_glmboost,
                                SVM = model_svm,
                                KNN = model_knn,
                                XGB_tree = model_xgbtree,
                                CFOREST = model_cforest,
                                C5.0Cost = model_C5.0Cost))
      resamps_a
      summary(resamps_a)
      
      pdf('Kappa_ML_methods_Asian.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_a, metric = "Kappa")
      dev.off()
      
      
      pdf('Accuracy_ML_methods_Asian.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_a, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}

# black 25+
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "B/AA" &
        metadata_all$Age >25 & !is.na(metadata_all$Age); sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[4] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )
        
        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
      }
      
      resamps_b_25 <- resamples(list(RF = model_rf,
                                GLM = model_glm,
                                GLMB = model_glmboost,
                                SVM = model_svm,
                                KNN = model_knn,
                                XGB_tree = model_xgbtree,
                                C5.0Cost = model_C5.0Cost))
      resamps_b_25
      summary(resamps_b_25)
      
      pdf('Kappa_ML_methods_Black_25.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_b_25, metric = "Kappa")
      dev.off()
      
      
      pdf('Accuracy_ML_methods_Black_25.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_b_25, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}

# white 25+
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "White"&
        metadata_all$Age >25 & !is.na(metadata_all$Age); sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[5] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )
        
        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
        model_cforest <- train(
          output_all ~ .,
          data = data,
          method = 'cforest',
          trControl = fitControl
        )
        
      }
      
      resamps_w_25 <- resamples(list(RF = model_rf,
                                GLM = model_glm,
                                GLMB = model_glmboost,
                                SVM = model_svm,
                                KNN = model_knn,
                                XGB_tree = model_xgbtree,
                                CFOREST = model_cforest,
                                C5.0Cost = model_C5.0Cost))
      resamps_w_25
      summary(resamps_w_25)
      
      pdf('Kappa_ML_methods_White_25.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_w_25, metric = "Kappa")
      dev.off()
      
      
      pdf('Accuracy_ML_methods_White_25.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_w_25, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}

# Asian 25+
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "Asian"&
        metadata_all$Age >25 & !is.na(metadata_all$Age); sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[6] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )
        
        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
        model_cforest <- train(
          output_all ~ .,
          data = data,
          method = 'cforest',
          trControl = fitControl
        )
        
      }
      
      resamps_a_25 <- resamples(list(RF = model_rf,
                                GLM = model_glm,
                                GLMB = model_glmboost,
                                SVM = model_svm,
                                KNN = model_knn,
                                XGB_tree = model_xgbtree,
                                CFOREST = model_cforest,
                                C5.0Cost = model_C5.0Cost))
      resamps_a_25
      summary(resamps_a_25)
      
      pdf('Kappa_ML_methods_asian_25.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_a_25, metric = "Kappa")
      dev.off()
      
      
      pdf('Accuracy_ML_methods_asian_25.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_a_25, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}

# black 25+ ga_135-
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "B/AA" &
        metadata_all$Age >25 & !is.na(metadata_all$Age) &
        metadata_all$gest_day_collection <=135 & !is.na(metadata_all$gest_day_collection); sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[7] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )
        
        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
      }
      
      resamps_b_25_135 <- resamples(list(RF = model_rf,
                                     GLM = model_glm,
                                     GLMB = model_glmboost,
                                     SVM = model_svm,
                                     KNN = model_knn,
                                     XGB_tree = model_xgbtree,
                                     C5.0Cost = model_C5.0Cost))
      resamps_b_25_135
      summary(resamps_b_25)
      
      pdf('Kappa_ML_methods_Black_25_135.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_b_25_135, metric = "Kappa")
      dev.off()
      
      
      pdf('Accuracy_ML_methods_Black_25_135.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_b_25_135, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}

# white 25+ ga_135-
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "White"&
        metadata_all$Age >25 & !is.na(metadata_all$Age)&
        metadata_all$gest_day_collection <=135 & !is.na(metadata_all$gest_day_collection); sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[8] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )
        
        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
        model_cforest <- train(
          output_all ~ .,
          data = data,
          method = 'cforest',
          trControl = fitControl
        )
        
      }
      
      resamps_w_25_135 <- resamples(list(RF = model_rf,
                                     GLM = model_glm,
                                     GLMB = model_glmboost,
                                     SVM = model_svm,
                                     KNN = model_knn,
                                     XGB_tree = model_xgbtree,
                                     CFOREST = model_cforest,
                                     C5.0Cost = model_C5.0Cost))
      resamps_w_25_135
      summary(resamps_w_25_135)
      
      pdf('Kappa_ML_methods_White_25_135.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_w_25_135, metric = "Kappa")
      dev.off()
      
      
      pdf('Accuracy_ML_methods_White_25_135.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_w_25_135, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}

# Asian 25+ ga_135-
{
  setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
  {
    
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "Asian"&
        metadata_all$Age >25 & !is.na(metadata_all$Age)&
        metadata_all$gest_day_collection <=135 & !is.na(metadata_all$gest_day_collection); sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      case_number[9] = length(output_all)
    }
    
    # Between-Models
    {
      all_model_info <- getModelInfo()
      all_model_names <- names(all_model_info)
      
      data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
      
      set.seed(2024)
      
      fitControl <- trainControl(## 5-fold CV
        method = "repeatedcv",
        number = 5,
        ## repeated ten times
        repeats = 10)
      
      {
        cl <- makeSOCKcluster(detectCores())
        registerDoSNOW(cl)
        
        model_rf <- train(
          output_all ~ .,
          data = data,
          method = 'rf',
          keep.inbag = TRUE,  # Pass additional parameters to randomForest
          trControl = fitControl
        )
        
        model_glmboost <- train(
          output_all ~ .,
          data = data,
          method = 'glmboost',
          trControl = fitControl
        )
        
        model_svm <- train(output_all ~ .,
                           data = data,
                           method = "svmLinear",
                           trControl = fitControl
        )
        
        model_knn <- train(
          output_all ~ .,
          data = data,
          method = 'knn',
          trControl = fitControl
        )
        
        model_xgbtree <- train(
          output_all ~ .,
          data = data,
          method = 'xgbTree',
          trControl = fitControl
        )
        
        model_C5.0Cost <- train(
          output_all ~ .,
          data = data,
          method = 'C5.0Cost',
          trControl = fitControl
        )
        
        model_glm <- train(
          output_all ~ .,
          data = data,
          method = 'glm',
          trControl = fitControl
        )
        
        model_cforest <- train(
          output_all ~ .,
          data = data,
          method = 'cforest',
          trControl = fitControl
        )
        
      }
      
      resamps_a_25_135 <- resamples(list(RF = model_rf,
                                     GLM = model_glm,
                                     GLMB = model_glmboost,
                                     SVM = model_svm,
                                     KNN = model_knn,
                                     XGB_tree = model_xgbtree,
                                     CFOREST = model_cforest,
                                     C5.0Cost = model_C5.0Cost))
      resamps_a_25_135
      summary(resamps_a_25_135)
      
      pdf('Kappa_ML_methods_asian_25_135.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_a_25_135, metric = "Kappa")
      dev.off()
      
      
      pdf('Accuracy_ML_methods_asian_25_135.pdf',width = 3.5, height = 2.5)
      trellis.par.set(caretTheme())
      dotplot(resamps_a_25_135, metric = "Accuracy")
      dev.off()
    }
    
  }
  
}
write.csv(case_number,'ML_case_number.csv')



##### ML RFE #####
# run 'ML_RFE_server.R' and then continue
setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_race/')
library(stringr) 
library(ggplot2)
library(parallel)
library(compositions) # for 'clr' function
library(vegan) # for diversity

library(pROC)
library(doSNOW)
library(parallel)
library(matrixStats)
library(erer) # export list
library(ggplot2)
library(randomForest)
library(forestError)

library(caret)
library(tidyr) # for 'gather'
# Black
{
  # summary
  {
    # selected features
    data = read.csv('selected_features_list_Black.csv')
    data = data$z..k..[data$Result != '' & data$Result != 'Result']
    
    data = as.data.frame(table(data))
    data = data[order(data$Freq, decreasing = T),]
    write.csv(data,'feature_freq_Black.csv')
    
    colnames(data) = c('Taxa', 'Frequency')
    data = data[data$Frequency>=50,]
    data$Taxa = factor(data$Taxa, levels = data$Taxa)
    
    ggplot(data, aes(Taxa, Frequency)) + 
      geom_bar(stat="identity", position="stack", width=1, color = 'grey') +theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1),
            axis.ticks=element_blank(),legend.position="none") + labs(y = 'Frequency / 100 selections')
    ggsave('feature_freq_Black.pdf',width = 5, height = 5)

    
  }
  
  # Between-Models
  {
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "B/AA"; sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      #input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      input_all = input_all[,colnames(input_all) %in% data$Taxa[data$Frequency >=50]]
    }
    
    all_model_info <- getModelInfo()
    all_model_names <- names(all_model_info)
    
    data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
    
    set.seed(2024)
    
    fitControl <- trainControl(## 5-fold CV
      method = "repeatedcv",
      number = 5,
      ## repeated ten times
      repeats = 10)
    
    {
      cl <- makeSOCKcluster(detectCores())
      registerDoSNOW(cl)
      
      model_rf <- train(
        output_all ~ .,
        data = data,
        method = 'rf',
        keep.inbag = TRUE,  # Pass additional parameters to randomForest
        trControl = fitControl
      )
      
      model_glmboost <- train(
        output_all ~ .,
        data = data,
        method = 'glmboost',
        trControl = fitControl
      )
      
      model_svm <- train(output_all ~ .,
                         data = data,
                         method = "svmLinear",
                         trControl = fitControl
      )
      
      model_knn <- train(
        output_all ~ .,
        data = data,
        method = 'knn',
        trControl = fitControl
      )
      
      model_xgbtree <- train(
        output_all ~ .,
        data = data,
        method = 'xgbTree',
        trControl = fitControl
      )
      
      model_C5.0Cost <- train(
        output_all ~ .,
        data = data,
        method = 'C5.0Cost',
        trControl = fitControl
      )
      
      model_glm <- train(
        output_all ~ .,
        data = data,
        method = 'glm',
        trControl = fitControl
      )
      
      model_cforest <- train(
        output_all ~ .,
        data = data,
        method = 'cforest',
        trControl = fitControl
      )
      
    }
    
    resamps_b_RFE <- resamples(list(RF = model_rf,
                                GLM = model_glm,
                                GLMB = model_glmboost,
                                SVM = model_svm,
                                KNN = model_knn,
                                XGB_tree = model_xgbtree,
                                CFOREST = model_cforest,
                                C5.0Cost = model_C5.0Cost))
    resamps_b_RFE
    summary(resamps_b_RFE)
    
    pdf('Kappa_ML_methods_resamps_b_RFE.pdf',width = 3.5, height = 2.5)
    trellis.par.set(caretTheme())
    dotplot(resamps_b_RFE, metric = "Kappa")
    dev.off()
    
    
    pdf('Accuracy_ML_methods_resamps_b_RFE.pdf',width = 3.5, height = 2.5)
    trellis.par.set(caretTheme())
    dotplot(resamps_b_RFE, metric = "Accuracy")
    dev.off()
  }
}

# White
{
  # summary
  {
    # selected features
    data = read.csv('selected_features_list_White.csv')
    data = data$z..k..[data$Result != '' & data$Result != 'Result']
    
    data = as.data.frame(table(data))
    data = data[order(data$Freq, decreasing = T),]
    write.csv(data,'feature_freq_White.csv')
    
    colnames(data) = c('Taxa', 'Frequency')
    data = data[data$Frequency>=50,]
    data$Taxa = factor(data$Taxa, levels = data$Taxa)
    
    ggplot(data, aes(Taxa, Frequency)) + 
      geom_bar(stat="identity", position="stack", width=1, color = 'grey') +theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1),
            axis.ticks=element_blank(),legend.position="none") + labs(y = 'Frequency / 100 selections')
    ggsave('feature_freq_White.pdf',width = 3, height = 4)
    
    
  }
  
  # Between-Models
  {
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "White"; sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      #input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      input_all = input_all[,colnames(input_all) %in% data$Taxa[data$Frequency >=50]]
      
      table(output_all)
      length(output_all)
      dim(input_all)
      
    }
    
    all_model_info <- getModelInfo()
    all_model_names <- names(all_model_info)
    
    data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
    
    set.seed(2024)
    
    fitControl <- trainControl(## 5-fold CV
      method = "repeatedcv",
      number = 5,
      ## repeated ten times
      repeats = 10)
    
    {
      cl <- makeSOCKcluster(detectCores())
      registerDoSNOW(cl)
      
      model_rf <- train(
        output_all ~ .,
        data = data,
        method = 'rf',
        keep.inbag = TRUE,  # Pass additional parameters to randomForest
        trControl = fitControl
      )
      
      model_glmboost <- train(
        output_all ~ .,
        data = data,
        method = 'glmboost',
        trControl = fitControl
      )
      
      model_svm <- train(output_all ~ .,
                         data = data,
                         method = "svmLinear",
                         trControl = fitControl
      )
      
      model_knn <- train(
        output_all ~ .,
        data = data,
        method = 'knn',
        trControl = fitControl
      )
      
      model_xgbtree <- train(
        output_all ~ .,
        data = data,
        method = 'xgbTree',
        trControl = fitControl
      )
      
      model_C5.0Cost <- train(
        output_all ~ .,
        data = data,
        method = 'C5.0Cost',
        trControl = fitControl
      )
      
      model_glm <- train(
        output_all ~ .,
        data = data,
        method = 'glm',
        trControl = fitControl
      )
      
      model_cforest <- train(
        output_all ~ .,
        data = data,
        method = 'cforest',
        trControl = fitControl
      )
      
    }
    
    resamps_w_RFE <- resamples(list(RF = model_rf,
                                    GLM = model_glm,
                                    GLMB = model_glmboost,
                                    SVM = model_svm,
                                    KNN = model_knn,
                                    XGB_tree = model_xgbtree,
                                    CFOREST = model_cforest,
                                    C5.0Cost = model_C5.0Cost))
    resamps_w_RFE
    summary(resamps_w_RFE)
    
    pdf('Kappa_ML_methods_resamps_w_RFE.pdf',width = 3.5, height = 2.5)
    trellis.par.set(caretTheme())
    dotplot(resamps_w_RFE, metric = "Kappa")
    dev.off()
    
    
    pdf('Accuracy_ML_methods_resamps_w_RFE.pdf',width = 3.5, height = 2.5)
    trellis.par.set(caretTheme())
    dotplot(resamps_w_RFE, metric = "Accuracy")
    dev.off()
  }
}

# asian
{
  # summary
  {
    # selected features
    data = read.csv('selected_features_list_asian.csv')
    data = data$z..k..[data$Result != '' & data$Result != 'Result']
    
    data = as.data.frame(table(data))
    data = data[order(data$Freq, decreasing = T),]
    write.csv(data,'feature_freq_asian.csv')
    
    colnames(data) = c('Taxa', 'Frequency')
    data = data[data$Frequency>=10,]
    data$Taxa = factor(data$Taxa, levels = data$Taxa)
    
    ggplot(data, aes(Taxa, Frequency)) + 
      geom_bar(stat="identity", position="stack", width=1, color = 'grey') +theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            axis.ticks=element_blank(),legend.position="none") + labs(y = 'Frequency / 100 selections')
    ggsave('feature_freq_asian.pdf',width = 2, height = 5)
    
    
  }
  
  # Between-Models
  {
    cl <- makeSOCKcluster(detectCores())
    registerDoSNOW(cl)
    
    # preparing input and output
    {
      keep = !is.na(metadata_all$Race) & metadata_all$Race == "Asian"; sum(keep)
      metadata_x = metadata_all[keep,]; 
      reads_table_x = reads_table_all_batch_rm[,colnames(reads_table_all_batch_rm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      
      #input_all = FilterFeatures(input_all, output_all, p.value_th = 0.05)
      
      table(output_all)
      length(output_all)
      dim(input_all)
      input_all = input_all[,colnames(input_all) %in% data$Taxa[data$Frequency >=50]]
    }
    
    all_model_info <- getModelInfo()
    all_model_names <- names(all_model_info)
    
    data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
    
    set.seed(2024)
    
    fitControl <- trainControl(## 5-fold CV
      method = "repeatedcv",
      number = 5,
      ## repeated ten times
      repeats = 10)
    
    {
      cl <- makeSOCKcluster(detectCores())
      registerDoSNOW(cl)
      
      model_rf <- train(
        output_all ~ .,
        data = data,
        method = 'rf',
        keep.inbag = TRUE,  # Pass additional parameters to randomForest
        trControl = fitControl
      )
      
      model_glmboost <- train(
        output_all ~ .,
        data = data,
        method = 'glmboost',
        trControl = fitControl
      )
      
      model_svm <- train(output_all ~ .,
                         data = data,
                         method = "svmLinear",
                         trControl = fitControl
      )
      
      model_knn <- train(
        output_all ~ .,
        data = data,
        method = 'knn',
        trControl = fitControl
      )
      
      model_xgbtree <- train(
        output_all ~ .,
        data = data,
        method = 'xgbTree',
        trControl = fitControl
      )
      
      model_C5.0Cost <- train(
        output_all ~ .,
        data = data,
        method = 'C5.0Cost',
        trControl = fitControl
      )
      
      model_glm <- train(
        output_all ~ .,
        data = data,
        method = 'glm',
        trControl = fitControl
      )
      
      model_cforest <- train(
        output_all ~ .,
        data = data,
        method = 'cforest',
        trControl = fitControl
      )
      
    }
    
    resamps_a_RFE <- resamples(list(RF = model_rf,
                                    GLM = model_glm,
                                    GLMB = model_glmboost,
                                    SVM = model_svm,
                                    KNN = model_knn,
                                    XGB_tree = model_xgbtree,
                                    CFOREST = model_cforest,
                                    C5.0Cost = model_C5.0Cost))
    resamps_a_RFE
    summary(resamps_a_RFE)
    
    pdf('Kappa_ML_methods_resamps_a_RFE.pdf',width = 3.5, height = 2.5)
    trellis.par.set(caretTheme())
    dotplot(resamps_a_RFE, metric = "Kappa")
    dev.off()
    
    
    pdf('Accuracy_ML_methods_resamps_a_RFE.pdf',width = 3.5, height = 2.5)
    trellis.par.set(caretTheme())
    dotplot(resamps_a_RFE, metric = "Accuracy")
    dev.off()
  }
}

##### all factors #####
setwd('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/ML_with_all_variables/')
load("~/Desktop/2_KH/VMB_PTB/others/Untitled.RData")
# association network
{
  metadata = metadata_all[,c(7,1,8,10,11)]
  colnames(metadata)[colnames(metadata) == 'disease'] = 'BV'
  
  keep = sapply(1:ncol(metadata), function(j) (is.numeric(metadata[,j]))); sum(keep)
  metadata_1 = metadata[,keep]; metadata_2 = metadata[,!keep];  # metadata_1 is interval
  
  keep = sapply(1:ncol(metadata_2), function(j) (length(unique(metadata_2[!is.na(metadata_2[,j]),j])) ==2)); sum(keep)
  metadata_3 = metadata_2[,keep]; metadata_2 = metadata_2[,!keep];  # metadata_2 is ordinal  # metadata_3 is binary
  metadata_3 = data.frame(Outcome = metadata_3)
  
  # The spearman's correlation for interval and ordinal
  {
    data =  t(cbind(metadata_1, metadata_2)); 
    
    output = newwork_rcorr(data, normalization_method = NA, type = 'spearman', pvalue = 0.05, 
                           cor_parameter= 0, style = 1, bar_max = 1, bar_min = -1, 
                           pheatmap_fontsize = 5, treeheight = 20, alpha = 0.05,FDR = T)
    
    correlated_variables_network_1 = data.frame(p_value = output$pvalue, adj_pvalue = output$adj_pvalue, R_value = output$Rvalue)
    correlated_variables_network_1 = correlated_variables_network_1[,c(3,5,4,1,2)]
    colnames(correlated_variables_network_1) = c('Variable 1','Variable 2', 'R-value','P-value','Adj-P-value')
    write.csv(correlated_variables_network_1,'correlation_interval_and_ordinal.csv',row.names = F)
  }
  
  
  # The Point-Biserial Correlation between binary and interval/ordinal
  {
    correlated_variables_network_2 = as.data.frame(matrix(data = NA, ncol = 5, nrow = 100000))
    colnames(correlated_variables_network_2) = c('Variable 1','Variable 2', 'R-value','P-value','Adj-P-value')
    
    data =  as.data.frame(cbind(metadata_1, metadata_2))
    n=0
    for (a in 1: ncol(data)) {
      for (b in 1: ncol(metadata_3)) {
        
        data_2 = data.frame(V1 = xtfrm(data[,a]), V2 = xtfrm(metadata_3[,b]))
        data_2 = data_2[!is.na(data_2$V1) & !is.na(data_2$V2),]
        
        if (nrow(data_2) < nrow(data) /10) {next}
        
        n=n+1
        correlated_variables_network_2$`Variable 1`[n] = colnames(data)[a]
        correlated_variables_network_2$`Variable 2`[n] = colnames(metadata_3)[b]
        
        x = cor.test(data_2$V1, data_2$V2)
        correlated_variables_network_2$`R-value`[n] = as.numeric(as.character(x$estimate))
        correlated_variables_network_2$`P-value`[n] = as.numeric(as.character(x$p.value))
      }
    }
    correlated_variables_network_2 = correlated_variables_network_2[!is.na(correlated_variables_network_2$`Variable 1`),]
    x = adjust.p(correlated_variables_network_2$`P-value`)
    correlated_variables_network_2$`Adj-P-value` = x$adjp$adjusted.p
    write.csv(correlated_variables_network_2,'correlation_binary_and_interval_or_ordinal.csv',row.names = F)
  }
  
  # The Phi Coefficient among binary
#  {
#    output = newwork_chi_squared(t(metadata_3), style = 1,pheatmap_fontsize = 5, treeheight = 50, alpha = 0.05,FDR = T)
    
#    correlated_variables_network_3 = output$pvalue
#    colnames(correlated_variables_network_3) = c('Variable 1','Variable 2', 'R-value','P-value','Adj-P-value')
#    write.csv(correlated_variables_network_3,'correlation_binary.csv',row.names = F)
#  }
  
#  colnames(correlated_variables_network_3) = colnames(correlated_variables_network_1)
  correlated_variables_network = rbind(correlated_variables_network_1,correlated_variables_network_2)
  correlated_variables_network = correlated_variables_network[correlated_variables_network$`Adj-P-value`<=0.05 & !is.na(correlated_variables_network$`Adj-P-value`),]
  correlated_variables_network = data.frame(Source = correlated_variables_network$`Variable 1`,
                                            Weight = abs(correlated_variables_network$`R-value`),
                                            Target = correlated_variables_network$`Variable 2`,
                                            Correlation= correlated_variables_network$`R-value`)
  
  correlated_variables_network$V1 = paste0(correlated_variables_network$Source,correlated_variables_network$Target)
  correlated_variables_network$V2 = paste0(correlated_variables_network$Target,correlated_variables_network$Source)
  for (a in 1:nrow(correlated_variables_network)) {
    if (correlated_variables_network$V1[a] > correlated_variables_network$V2[a]) {
      x = correlated_variables_network$V1[a]; correlated_variables_network$V1[a] = correlated_variables_network$V2[a]
      correlated_variables_network$V2[a] = x
    }
  }
  correlated_variables_network = correlated_variables_network[!duplicated(correlated_variables_network$V1),]
  
  write.csv(correlated_variables_network,'correlated_variables_network.csv', row.names = F)
}

# association
{
  metadata = metadata_all[,c(7,1,8,10,11)]
  data = metadata[,c(1,3)]
  data = data[!is.na(data$gest_day_collection),]
  ggplot(data=data, aes(x=Outcome, y=gest_day_collection)) +geom_violin(trim=T, aes(fill=Outcome))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1, alpha= 0.1)+theme_bw()+
    theme_classic()+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  ggsave('Correlation_gest_day_collection.pdf',width=2.5, height=2.5)
  
  data = metadata[,c(1,4)]
  data = data[!is.na(data$Age),]
  ggplot(data=data, aes(x=Outcome, y=Age)) +geom_violin(trim=T, aes(fill=Outcome))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1, alpha= 0.1)+theme_bw()+
    theme_classic()+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  ggsave('Correlation_Age.pdf',width=2.5, height=2.5)
  
  data = metadata[,c(1,5)]
  data = data[!is.na(data$Race),]
  data = as.data.frame(table(data))
  data = data[data$Race != 'American Indian' & data$Race != 'Mixed',]
  
  ggplot(data, aes(x = interaction(Outcome, Race), y = Freq, aes(fill=Race))) +
    geom_bar(stat = "identity", aes(fill=Race)) +
    labs( y = "Case number") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave('Correlation_Race.pdf',width=4, height=2.5)
}

# establish the model
{
  # preparing input and output
  {
    keep = !is.na(metadata_cm$gest_day_collection) & metadata_cm$gest_day_collection <= 91 &
      !is.na(metadata_cm$Race) & metadata_cm$Race %in% c("Asian", "B/AA","White") &
      !is.na(metadata_cm$Age); sum(keep)
    metadata_x = metadata_cm[keep,]; 
    reads_table_x = reads_table_cm[,colnames(reads_table_cm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
    y = metadata_x$Run
    reads_table = reads_table_x
    reads_table = as.data.frame(t(reads_table))
    reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
    input_all = reads_table
    
    output_all = as.character(metadata_x$Outcome)
    keep = output_all == 'PTB'
    output_all[keep] = 1; output_all[!keep] = 0; 
    data = as.data.frame(cbind(output_all, as.data.frame(reads_table)))
    
    table(output_all)
    
    data = read.csv('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/case_match_ML/feature_freq_RFE_0.68up.csv')
    data = data$data[data$Freq >= 5]
    input_all = input_all[,colnames(input_all) %in% data]
    input_all = cbind(input_all, metadata_x[,c(4:6)])
    dim(input_all)
  }

  # Between-Models
  {
    # preparing input and output
    {
      keep = !is.na(metadata_cm$gest_day_collection) & metadata_cm$gest_day_collection <= 91 &
        !is.na(metadata_cm$Race) & metadata_cm$Race %in% c("Asian", "B/AA","White") &
        !is.na(metadata_cm$Age); sum(keep)
      metadata_x = metadata_cm[keep,]; 
      reads_table_x = reads_table_cm[,colnames(reads_table_cm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      y = metadata_x$Run
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      data = as.data.frame(cbind(output_all, as.data.frame(reads_table)))
      
      table(output_all)
      
      data = read.csv('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/case_match_ML/feature_freq_RFE_0.68up.csv')
      data = data$data[data$Freq >= 5]
      input_all = input_all[,colnames(input_all) %in% data]
      input_all = cbind(input_all, metadata_x[,c(4:6)])
      dim(input_all)
    }
    
    all_model_info <- getModelInfo()
    all_model_names <- names(all_model_info)
    
    data = as.data.frame(cbind(output_all, as.data.frame(input_all)))
    
    set.seed(2024)
    
    fitControl <- trainControl(## 5-fold CV
      method = "repeatedcv",
      number = 5,
      ## repeated ten times
      repeats = 10)
    
    model_rf <- train(
      output_all ~ .,
      data = data,
      method = 'rf',
      keep.inbag = TRUE,  # Pass additional parameters to randomForest
      trControl = fitControl
    )
    
    model_glmboost <- train(
      output_all ~ .,
      data = data,
      method = 'glmboost',
      trControl = fitControl
    )
    
    model_svm <- train(output_all ~ .,
                       data = data,
                       method = "svmLinear",
                       trControl = fitControl
    )
    
    model_knn <- train(
      output_all ~ .,
      data = data,
      method = 'knn',
      trControl = fitControl
    )
    
    model_xgbLinear <- train(
      output_all ~ .,
      data = data,
      method = 'xgbLinear',
      trControl = fitControl
    )
    
    model_xgbtree <- train(
      output_all ~ .,
      data = data,
      method = 'xgbTree',
      trControl = fitControl
    )
    
    model_xgbDART <- train(
      output_all ~ .,
      data = data,
      method = 'xgbDART',
      trControl = fitControl
    )
    
    model_bayesglm <- train(
      output_all ~ .,
      data = data,
      method = 'bayesglm',
      trControl = fitControl
    )
    
    model_C5.0 <- train(
      output_all ~ .,
      data = data,
      method = 'C5.0',
      trControl = fitControl
    )
    
    model_cforest <- train(
      output_all ~ .,
      data = data,
      method = 'cforest',
      trControl = fitControl
    )
    
    model_ctree <- train(
      output_all ~ .,
      data = data,
      method = 'ctree',
      trControl = fitControl
    )
    
    model_ctree2 <- train(
      output_all ~ .,
      data = data,
      method = 'ctree2',
      trControl = fitControl
    )
    
    model_C5.0Cost <- train(
      output_all ~ .,
      data = data,
      method = 'C5.0Cost',
      trControl = fitControl
    )
    
    model_glm <- train(
      output_all ~ .,
      data = data,
      method = 'glm',
      trControl = fitControl
    )
    
    model_Rborist <- train(
      output_all ~ .,
      data = data,
      method = 'Rborist',
      trControl = fitControl
    )
    
    model_cforest <- train(
      output_all ~ .,
      data = data,
      method = 'cforest',
      trControl = fitControl
    )
    
    model_parRF <- train(
      output_all ~ .,
      data = data,
      method = 'parRF',
      trControl = fitControl
    )
    
    model_rFerns <- train(
      output_all ~ .,
      data = data,
      method = 'rFerns',
      trControl = fitControl
    )
    
    model_ordinalRF <- train(
      output_all ~ .,
      data = data,
      method = 'ordinalRF',
      trControl = fitControl
    )
    
    resamps <- resamples(list(RF = model_rf,
                              Rborist = model_Rborist,
                              PARRF = model_parRF,
                              rFerns = model_rFerns,
                              CFOREST = model_cforest,
                              GLM = model_glm,
                              GLMB = model_glmboost,
                              SVM = model_svm,
                              KNN = model_knn,
                              XGB_linear = model_xgbLinear,
                              XGB_tree = model_xgbtree,
                              XGB_DART = model_xgbDART,
                              CFOREST = model_cforest,
                              CTREE = model_ctree,
                              CTREE2 = model_ctree2,
                              C5.0Cost = model_C5.0Cost))
    resamps
    summary(resamps)
    
    pdf('Kappa_ML_methods_RFE_metadata.pdf',width = 4, height = 4)
    trellis.par.set(caretTheme())
    dotplot(resamps, metric = "Kappa")
    dev.off()
    
    
    pdf('Accuracy_ML_methods_RFE_metadata.pdf',width = 4, height = 4)
    trellis.par.set(caretTheme())
    dotplot(resamps, metric = "Accuracy")
    dev.off()
  }
  
  # testing >85
  {
    {
      keep = !is.na(metadata_cm$gest_day_collection) & metadata_cm$gest_day_collection <= 91 &
        !is.na(metadata_cm$Race) & metadata_cm$Race %in% c("Asian", "B/AA","White") &
        !is.na(metadata_cm$Age); sum(keep)
      metadata_x = metadata_cm[keep,]; 
      reads_table_x = reads_table_cm[,colnames(reads_table_cm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      y = metadata_x$Run
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      data = as.data.frame(cbind(output_all, as.data.frame(reads_table)))
      
      table(output_all)
      
      data = read.csv('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/case_match_ML/feature_freq_RFE_0.68up.csv')
      data = data$data[data$Freq >= 5]
      input_all = input_all[,colnames(input_all) %in% data]
      input_all = cbind(input_all, metadata_x[,c(4:6)])
      dim(input_all)
    }
    set.seed(2028)
    training = as.data.frame(cbind(output_all, as.data.frame((input_all))))
    
    model_1 <- train(
      output_all ~ .,
      data = training,
      method = 'rf',
      keep.inbag = TRUE,  # Pass additional parameters to randomForest
      trControl = trainControl(method = "cv", number = 5)  # 5-fold cross-validation control
    )
    
    # preparing input and output
    {
      keep = !is.na(metadata_cm$gest_day_collection) & metadata_cm$gest_day_collection > 85 &
        !is.na(metadata_cm$Race) & metadata_cm$Race %in% c("B/AA","White") &
        !is.na(metadata_cm$Age); sum(keep)
      metadata_x = metadata_cm[keep,]; 
      reads_table_x = reads_table_cm[,colnames(reads_table_cm) %in% metadata_x$Run]; reads_table_x = reads_table_x[metadata_x$Run]
      y = metadata_x$Run
      reads_table = reads_table_x
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      input_all = reads_table
      
      output_all = as.character(metadata_x$Outcome)
      keep = output_all == 'PTB'
      output_all[keep] = 1; output_all[!keep] = 0; 
      data = as.data.frame(cbind(output_all, as.data.frame(reads_table)))
      
      table(output_all)
      
      data = read.csv('/Users/binzhu/Desktop/1_KH/VMB_PTB/results/case_match_ML/feature_freq_RFE_0.68up.csv')
      data = data$data[data$Freq >= 5]
      input_all = input_all[,colnames(input_all) %in% data]
      input_all = cbind(input_all, metadata_x[,c(4:6)])
      dim(input_all)
    }
    
    testing = as.data.frame(cbind(output_all, as.data.frame((input_all))))
    predictions <- predict(model_1, newdata = testing, type = "prob")  
    auc_1 = roc(testing[,1], predictions$`1`)$auc      # It builds a ROC curve and returns a “roc” object, a list of class “roc”.
    #        auc_2
    auc_1
  }

}

##### Figures #####
# fig 1. Profiles of VMB. bar plot.   No case match compared with case match. dbRDA and marginal Adonis to test impact factors.
# fig 2. Different time point; longitudinal changes of microbes. Significant change in different races; whether transition associated with PTB.
# fig 3. Machine learning XGB_DART. Case match all, only one race.   No case match, only one race;  ML using transitional ratio between first tri and second tri. and add vagitype
# fig 4. In first trimester, case match linked by lines, show vagitype.  find details like PCA which PCA components are significantly different, and what is the linear formular for the PCA.



