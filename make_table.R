#Build species features table for all species
#Build species structures table for all species
#Merge them
#Build species homolog table for human vs species
#Mutate to get the desired data
#Merge by the canonical transcript 
#This should give final table to be used elsewhere

setwd("D:/TUDelft/3rd Year Nanobiology/BEP/R directory/Make table")

library(biomaRt)
library(dplyr)
library(stringr)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
la = listAttributes(human)

structure_attr <- c("ensembl_transcript_id","cds_length","ensembl_exon_id","exon_chrom_start",
                    "exon_chrom_end")

feature_attr <- c("ensembl_gene_id","external_gene_name","ensembl_transcript_id","ensembl_peptide_id",
                  "start_position","end_position","transcript_start","transcript_end","transcript_length","strand","chromosome_name",
                  "percentage_gene_gc_content","gene_biotype","transcript_biotype")


structure_t <- getBM(attributes = structure_attr, values=T, mart=human, uniqueRows=T)
feature_t <- getBM(attributes = feature_attr, values=T, mart=human, uniqueRows=T)
hsapiens_structure_feature_t <- merge(feature_t,structure_t,by="ensembl_transcript_id")

hsapiens_structure_feature_t <- hsapiens_structure_feature_t %>% 
  group_by(ensembl_transcript_id) %>% 
  mutate(min_start = min(exon_chrom_start)) %>% 
  mutate(max_end = max(exon_chrom_end)) %>% 
  mutate(exon_num = length(ensembl_transcript_id)) %>%
  mutate(exon_len = exon_chrom_end - exon_chrom_start) %>%
  mutate(intron_len = abs(lead(exon_chrom_start, n=1) - exon_chrom_end)) %>%
  mutate(longest_intron = max(intron_len, na.rm =T)) %>%
  mutate(real_transcript_len = transcript_end - transcript_start) %>%
  mutate(longest_exon = max(exon_len))

hsapiens_structure_feature_t$ensembl_exon_id = NULL
hsapiens_structure_feature_t$exon_chrom_start = NULL
hsapiens_structure_feature_t$exon_chrom_end = NULL
hsapiens_structure_feature_t$exon_len = NULL
hsapiens_structure_feature_t$intron_len = NULL

hsapiens_structure_feature_t <- unique(hsapiens_structure_feature_t)

rm(structure_t)
rm(feature_t)

#All species that need features and structure table info.
splist = c("rnorvegicus","mmusculus","ptroglodytes","clfamiliaris","sscrofa","btaurus","oarambouillet",
           "cporcellus","fcatus","pabelii","lafricana") 

binded_t <- data.frame()

for (sp in splist){
  mart = useMart("ensembl", dataset = paste0(sp,"_gene_ensembl"))
  
  homolog_attr <-c("ensembl_transcript_id","ensembl_peptide_id",paste0(sp,"_homolog_ensembl_peptide"),
                   paste0(sp,"_homolog_canonical_transcript_protein"),
                   paste0(sp,"_homolog_perc_id"),paste0(sp,"_homolog_perc_id_r1"),paste0(sp,"_homolog_goc_score"),
                   paste0(sp,"_homolog_wga_coverage"),paste0(sp,"_homolog_orthology_confidence"),
                   paste0(sp,"_homolog_orthology_type"))
  
  sp_structure_t <- getBM(attributes = structure_attr, values=T, mart=mart, uniqueRows=T)
  sp_feature_t <- getBM(attributes = feature_attr, values=T, mart=mart, uniqueRows=T)
  sp_structure_feature_t <- merge(sp_feature_t,sp_structure_t,by="ensembl_transcript_id")
  
  sp_structure_feature_t <- sp_structure_feature_t %>% 
    group_by(ensembl_transcript_id) %>% 
    mutate(min_start = min(exon_chrom_start)) %>% 
    mutate(max_end = max(exon_chrom_end)) %>% 
    mutate(exon_num = length(ensembl_transcript_id)) %>%
    mutate(exon_len = exon_chrom_end - exon_chrom_start) %>%
    mutate(intron_len = abs(lead(exon_chrom_start, n=1) - exon_chrom_end)) %>%
    mutate(longest_intron = max(intron_len, na.rm =T)) %>%
    mutate(real_transcript_len = transcript_end - transcript_start) %>%
    mutate(longest_exon = max(exon_len))
  
  sp_structure_feature_t$ensembl_exon_id = NULL
  sp_structure_feature_t$exon_chrom_start = NULL
  sp_structure_feature_t$exon_chrom_end = NULL
  sp_structure_feature_t$exon_len = NULL
  sp_structure_feature_t$intron_len = NULL
  
  sp_structure_feature_t <- unique(sp_structure_feature_t)
  
  colnames(sp_structure_feature_t) <- paste("sp", colnames(sp_structure_feature_t), sep = "_")
  
  rm(sp_structure_t)
  rm(sp_feature_t)
  
  homolog_t <- getBM(attributes = homolog_attr, values=T, mart=human, uniqueRows=T)
  
  cn<-colnames(homolog_t)
  cn2<-str_replace(cn,paste0(sp,"_"),"sp_")
  colnames(homolog_t)<-cn2
  
  canonical_homolog_t <- homolog_t[homolog_t$ensembl_peptide_id==homolog_t$sp_homolog_canonical_transcript_protein ,]
  canonical_homolog_t <- canonical_homolog_t[!(is.na(canonical_homolog_t$sp_homolog_orthology_confidence)),]
  canonical_homolog_t$ensembl_peptide_id=NULL
  
  cn<-colnames(sp_structure_feature_t)
  cn[3]<-"sp_homolog_associated_gene_name"
  cn[4]<-"sp_homolog_ensembl_peptide"
  colnames(sp_structure_feature_t)<-cn
  
  #Merge the 3 tables
  hsapiens_homolog_t <- merge(hsapiens_structure_feature_t, canonical_homolog_t, by="ensembl_transcript_id")
  
  sp_structure_feature_t <- sp_structure_feature_t[sp_structure_feature_t$sp_gene_biotype=="protein_coding",]
  sp_structure_feature_t <- sp_structure_feature_t[!(is.na(sp_structure_feature_t$sp_homolog_ensembl_peptide) | 
                             sp_structure_feature_t$sp_homolog_ensembl_peptide==""), ]
  hsapiens_homolog_t <- hsapiens_homolog_t[hsapiens_homolog_t$gene_biotype=="protein_coding",]
  
  #Filter first before getting rank data (orthology type TBD)
  final_t <- merge(sp_structure_feature_t,hsapiens_homolog_t, by="sp_homolog_ensembl_peptide")
  final_t <- final_t[final_t$sp_homolog_orthology_type=="ortholog_one2one" &
                           final_t$exon_num>1 & final_t$sp_exon_num>1 ,]
  
  final_t$species <- sp
  final_t$len_maxmin <- final_t$max_end - final_t$min_start
  final_t$len_startend <- final_t$end_position - final_t$start_position
  final_t$sp_len_maxmin <- final_t$sp_max_end - final_t$sp_min_start
  final_t$sp_len_startend <- final_t$sp_end_position - final_t$sp_start_position
  
  final_t$rank_len_maxmin <- rank(final_t$len_maxmin,ties.method="min")
  final_t$rank_len_startend <- rank(final_t$len_startend,ties.method="min")
  final_t$rank_sp_len_maxmin <- rank(final_t$sp_len_maxmin,ties.method="min")
  final_t$rank_sp_len_startend <- rank(final_t$sp_len_startend,ties.method="min")
  
  write.csv(final_t,file=paste0(sp,"_main_table"))

  binded_t <- rbind.data.frame(binded_t,final_t)
}

write.csv(binded_t, file="binded_table.csv",row.names = F)









