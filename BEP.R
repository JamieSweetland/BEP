setwd("D:/TUDelft/3rd Year Nanobiology/BEP/R directory/Main table analysis")

library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(zoo)

#Preparing the main table
{
binded_table <- read.csv(file="binded_table.csv",header=TRUE ,sep=",")

tmp<-binded_table[binded_table$transcript_biotype=="protein_coding" & 
                    binded_table$sp_transcript_biotype=="protein_coding" &
                    binded_table$chromosome_name!='X' &
                    binded_table$chromosome_name!='Y' ,]


main_t <- tmp %>% 
  mutate(len_intron = len_maxmin - transcript_length) %>% 
  mutate(sp_len_intron = sp_len_maxmin - sp_transcript_length) %>%
  mutate(ratio_intron_transcript = len_intron/len_maxmin) %>%
  mutate(sp_ratio_intron_transcript = sp_len_intron/sp_len_maxmin) %>%
  mutate(ratio_exon_transcript = transcript_length/len_maxmin) %>%
  mutate(sp_ratio_exon_transcript = sp_transcript_length/sp_len_maxmin) %>%
  mutate(sp_cds_maxmin_ratio = sp_cds_length/sp_len_maxmin) %>%
  mutate(cds_maxmin_ratio = cds_length/len_maxmin) %>%
  #main based on the following
  mutate(dif_len_intron = abs(len_intron - sp_len_intron)) %>%
  mutate(dif_len_exon = abs(transcript_length - sp_transcript_length)) %>%
  mutate(dif_num_exon = abs(sp_exon_num - exon_num)) %>%
  mutate(dif_intron = dif_len_intron/len_intron) %>%
  mutate(dif_rank = abs(rank_len_maxmin - rank_sp_len_maxmin)) %>%
  mutate(dif_ratio_intron = abs(ratio_intron_transcript - sp_ratio_intron_transcript)) %>%
  mutate(dif_gc = abs(sp_percentage_gene_gc_content - percentage_gene_gc_content)) %>%
  mutate(dif_ratio_exon = abs(ratio_exon_transcript - sp_ratio_exon_transcript)) %>%
  mutate(dif_cds_ratio = abs(sp_cds_maxmin_ratio - cds_maxmin_ratio))
}

#Length difference filter
{
  main_t <- main_t %>%
    group_by(ensembl_transcript_id) %>%
    mutate(num_sp = length(ensembl_transcript_id[sp_homolog_orthology_confidence==1 &
                                                   sp_homolog_perc_id > 50 &
                                                   sp_homolog_perc_id_r1 > 50])) %>%
    
    mutate(num_similar_sp1 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.5 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50])) %>%
        
    mutate(num_similar_sp2 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.25 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50]))
  
  data_t <- unique(main_t[,c(22:25,30:36,38:41,76:78)])
  
  num_long_genes_0.5 = unique(main_t[main_t$real_transcript_len>50000 &
                                            main_t$sp_real_transcript_len>50000 &
                                            main_t$num_similar_sp1 >= 9
                                          ,"ensembl_gene_id"])
  
  num_long_genes_0.25 = unique(main_t[main_t$real_transcript_len>50000 &
                                            main_t$sp_real_transcript_len>50000 &
                                            main_t$num_similar_sp2 >= 9
                                          ,"ensembl_gene_id"])
  
  #Plotting length difference filters
  {
    p1 <- ggplot(data_t) +
      geom_point(aes(x=log10(real_transcript_len),y=num_sp)) +
      geom_point(data=data_t[data_t$external_gene_name=='IGF1',],
                 aes(y=num_sp, x=log10(real_transcript_len)), alpha=1,color='blue',size=5) +
      ylab("number of species") +
      xlab(expression(paste(log[10], " transcript length (bp)"))) +
      theme_minimal()
    
    p2 <- ggplot(data_t) +
      geom_point(aes(x=log10(real_transcript_len),y=num_similar_sp1)) +
      geom_point(data=data_t[data_t$external_gene_name=='IGF1',],
                 aes(y=num_similar_sp1, x=log10(real_transcript_len)), alpha=1,color='blue',size=5) +
      ylab("number of similar species") +
      xlab(expression(paste(log[10], " transcript length (bp)"))) +
      theme_minimal()
    
    p3 <- ggplot(data_t) +
      geom_point(aes(x=log10(real_transcript_len),y=num_similar_sp2)) +
      geom_point(data=data_t[data_t$external_gene_name=='IGF1',],
                 aes(y=num_similar_sp2, x=log10(real_transcript_len)), alpha=1,color='blue',size=5) +
      ylab("number of similar species") +
      xlab(expression(paste(log[10], " transcript length (bp)"))) +
      theme_minimal()
    
    p4 <- ggplot(data_t, mapping=aes(x=log10(real_transcript_len), y=num_sp)) +
      stat_summary_bin(fun = "mean", geom="bar", bins=20) +
      ylab("mean number of species") +
      xlab(expression(paste(log[10], " transcript length (bp)"))) +
      theme_minimal()
    
    p5 <- ggplot(data_t, mapping=aes(x=log10(real_transcript_len), y=num_similar_sp1)) +
      stat_summary_bin(fun = "mean", geom="bar", bins=20) +
      ylab("mean number of similar species") +
      xlab(expression(paste(log[10], " transcript length (bp)"))) +
      theme_minimal()
    
    p6 <- ggplot(data_t, mapping=aes(x=log10(real_transcript_len), y=num_similar_sp2)) +
      stat_summary_bin(fun = "mean", geom="bar", bins=20) +
      ylab("mean number of similar species") +
      xlab(expression(paste(log[10], " transcript length (bp)"))) +
      theme_minimal()
    
    grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3,nrow=2)
    
  }
  
  
}

#CDS distribution
{
  data_t <- unique(main_t[,c(22:25,30:36,38:41,76:78)])
  
  ggplot(data_t) +
    geom_density(data=data_t[data_t$real_transcript_len>50000 ,],
                 aes(x=100*cds_length/real_transcript_len), fill='yellow') +
    xlim(0, 15) +
    xlab('% of transcript which is cds') +
    geom_vline(xintercept=0.5, color='blue', linetype='dashed', size=1.2) +
    geom_text(aes(0, 0.29, label='IGF1', vjust=-1), size=4, color='blue') +
    theme_minimal()
}

#CDS and length difference filter
{
  main_t <- main_t %>%
    group_by(ensembl_transcript_id) %>%
    mutate(num_sp = length(ensembl_transcript_id[sp_homolog_orthology_confidence==1 &
                                                   sp_homolog_perc_id > 50 &
                                                   sp_homolog_perc_id_r1 > 50])) %>%
    
    mutate(num_similar_sp1 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.5 &
        cds_maxmin_ratio < 0.007 &
        sp_cds_maxmin_ratio < 0.014 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50])) %>%
    
    mutate(num_similar_sp2 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.5 &
        cds_maxmin_ratio < 0.006 &
        sp_cds_maxmin_ratio < 0.013 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50]))
  
  data_t <- unique(main_t[,c(22:25,30:36,38:41,76:78)]) 
  
  goi_cds <- data_t[data_t$real_transcript_len>50000 & data_t$num_similar_sp1>=9 ,]
  goi_cds_strict <- data_t[data_t$real_transcript_len>50000 & data_t$num_similar_sp2>=9 ,]
  
  #Plotting CDS and length difference filters
  {
    p1 <- ggplot(data_t) +
      geom_point(aes(x=log10(real_transcript_len),y=num_similar_sp1)) +
      geom_point(data=data_t[data_t$external_gene_name=='IGF1',],
                 aes(y=num_sp, x=log10(real_transcript_len)), alpha=1,color='blue',size=5) +
      ylab("number of species") +
      xlab(expression(paste(log[10], " transcript length (bp)"))) +
      theme_minimal()
    
    p2 <- ggplot(data_t) +
      geom_point(aes(x=log10(real_transcript_len),y=num_similar_sp2)) +
      geom_point(data=data_t[data_t$external_gene_name=='IGF1',],
                 aes(y=num_similar_sp1, x=log10(real_transcript_len)), alpha=1,color='blue',size=5) +
      ylab("number of similar species") +
      xlab(expression(paste(log[10], " transcript length (bp)"))) +
      theme_minimal()
    
    grid.arrange(p1,p2,ncol=2,nrow=1)
  }
}

#Longest intron length distribution
{
  data_t <- unique(main_t[,c(22:25,30:36,38:41,76:78)])
  
  ggplot(data_t) +
    geom_density(data=data_t[data_t$real_transcript_len>50000 ,],
                 aes(x=log10(longest_intron)),fill='yellow') +
    geom_vline(xintercept = log10(data_t$longest_intron[data_t$external_gene_name=='IGF1']),
               color='blue', linetype='dashed', size=1.2) +
    geom_text(aes(4.85, 0.6, label='IGF1', vjust=-1), size=4, color='blue') +
    xlab(expression(paste(log[10], " longest intron length (bp)"))) +
    theme_minimal()
}

#CDS plotted against longest intron
{
  main_t <- main_t %>%
    group_by(ensembl_transcript_id) %>%
    mutate(num_sp = length(ensembl_transcript_id[sp_homolog_orthology_confidence==1 &
                                                   sp_homolog_perc_id > 50 &
                                                   sp_homolog_perc_id_r1 > 50])) %>%
    
    mutate(num_similar_sp1 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.5 &
        cds_maxmin_ratio < 0.007 &
        sp_cds_maxmin_ratio < 0.014 &
        longest_intron > 56000 &
        sp_longest_intron > 44000 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50])) %>%
    
    mutate(num_similar_sp2 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.5 &
        cds_maxmin_ratio < 0.006 &
        sp_cds_maxmin_ratio < 0.013 &
        longest_intron > 56000 &
        sp_longest_intron > 44000 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50]))
  
  data_t <- unique(main_t[,c(22:25,30:36,38:41,76:78)])
  
  ggplot(data_t) +
    geom_point(data=data_t[data_t$real_transcript_len>50000 ,],
               aes(y=100*cds_length/real_transcript_len, x=log10(longest_intron)), color='red') +
    
    geom_point(data=data_t[data_t$real_transcript_len>50000 & data_t$num_similar_sp2>=9 ,],
               aes(y=100*cds_length/real_transcript_len, x=log10(longest_intron)), color='black') +
    
    geom_point(data=data_t[data_t$external_gene_name=='IGF1',],
               aes(y=100*cds_length/real_transcript_len, x=log10(longest_intron)),
               color='blue', size=5) +
    
    geom_vline(xintercept = log10(data_t$longest_intron[data_t$external_gene_name=='IGF1']),
               color='black', linetype='dashed', size=1.2) +
    
    geom_hline(yintercept = 0.6, color='black', linetype='dashed', size=1.2) +
    xlab(expression(paste(log[10], " longest intron length (bp)"))) +
    ylab('% of transcript which is CDS') +
    ylim(0, 0.8) +
    theme_minimal()
}

#CDS plotted against longest intron proportion
{
  main_t <- main_t %>%
    group_by(ensembl_transcript_id) %>%
    mutate(num_sp = length(ensembl_transcript_id[sp_homolog_orthology_confidence==1 &
                                                   sp_homolog_perc_id > 50 &
                                                   sp_homolog_perc_id_r1 > 50])) %>%
    
    mutate(num_similar_sp1 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.5 &
        cds_maxmin_ratio < 0.007 &
        sp_cds_maxmin_ratio < 0.014 &
        longest_intron/real_transcript_len > 0.65 &
        sp_longest_intron/sp_real_transcript_len > 0.65 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50])) %>%
    
    mutate(num_similar_sp2 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.5 &
        cds_maxmin_ratio < 0.006 &
        sp_cds_maxmin_ratio < 0.013 &
        longest_intron/real_transcript_len > 0.65 &
        sp_longest_intron/sp_real_transcript_len > 0.65 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50]))
  
  data_t <- unique(main_t[,c(22:25,30:36,38:41,76:78)])
  
  ggplot(data_t) +
    geom_point(data=data_t[data_t$real_transcript_len>50000 ,],
               aes(y=100*cds_length/real_transcript_len, x=100*longest_intron/real_transcript_len), 
               color='red') +
    
    geom_point(data=data_t[data_t$real_transcript_len>50000 & data_t$num_similar_sp2>=9 &
                             data_t$longest_intron/data_t$real_transcript_len > 0.665,],
               aes(y=100*cds_length/real_transcript_len, x=100*longest_intron/real_transcript_len),
               color='black') +
    
    geom_point(data=data_t[data_t$external_gene_name=='IGF1',],
               aes(y=100*cds_length/real_transcript_len, x=100*longest_intron/real_transcript_len),
               color='blue', size=5) +
    
    geom_vline(xintercept = 100*data_t$longest_intron[data_t$external_gene_name=='IGF1']/
                 data_t$real_transcript_len[data_t$external_gene_name=='IGF1'],
               color='black', linetype='dashed', size=1.2) +
    
    geom_hline(yintercept = 0.6, color='black', linetype='dashed', size=1.2) +
    xlab('% of transcript which is longest intron') +
    ylab('% of transcript which is CDS') +
    ylim(0, 0.8) +
    theme_minimal()
  
  goi_cds_lip <- data_t[data_t$real_transcript_len>50000 & data_t$num_similar_sp1>=9 ,]
  goi_cds_lip_strict <- data_t[data_t$real_transcript_len>50000 & data_t$num_similar_sp2>=9 ,]
}

#GC content against the transcript length
{
  main_t <- main_t %>%
    group_by(ensembl_transcript_id) %>%
    mutate(num_sp = length(ensembl_transcript_id[sp_homolog_orthology_confidence==1 &
                                                   sp_homolog_perc_id > 50 &
                                                   sp_homolog_perc_id_r1 > 50])) %>%
    
    mutate(num_similar_sp1 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.5 &
        cds_maxmin_ratio < 0.007 &
        sp_cds_maxmin_ratio < 0.014 &
        longest_intron > 56000 &
        sp_longest_intron > 44000 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50])) %>%
    
    mutate(num_similar_sp2 = length(ensembl_transcript_id[
      abs(real_transcript_len - sp_real_transcript_len)/real_transcript_len < 0.5 &
        cds_maxmin_ratio < 0.006 &
        sp_cds_maxmin_ratio < 0.013 &
        longest_intron > 56000 &
        sp_longest_intron > 44000 &
        sp_homolog_orthology_confidence==1 &
        sp_homolog_perc_id > 50 &
        sp_homolog_perc_id_r1 > 50]))
  
  data_t <- unique(main_t[,c(22:25,30:36,38:41,76:78)])
  
  ggplot(data_t) +
    geom_point(data=data_t[data_t$real_transcript_len>50000 ,],
               aes(y=percentage_gene_gc_content, x=log10(real_transcript_len)), color='red') +
    
    geom_point(data=data_t[data_t$real_transcript_len>50000 & data_t$num_similar_sp1>=9 ,],
               aes(y=percentage_gene_gc_content, x=log10(real_transcript_len)), color='black') +
    
    geom_point(data=data_t[data_t$external_gene_name=='IGF1',],
               aes(y=percentage_gene_gc_content, x=log10(real_transcript_len)),
               color='blue', size=5) +
    
    geom_hline(yintercept = data_t$percentage_gene_gc_content[data_t$external_gene_name=='IGF1'], 
               color='blue', linetype='dashed', size=1.2) +
    
    xlab(expression(paste(log[10], " transcript length (bp)"))) +
    ylab('% GC content') +
    theme_minimal()
}