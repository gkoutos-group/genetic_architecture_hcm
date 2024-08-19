#This script joins the information from VEP with that of genebe
#A consequence level is also created in order to stratify the risk of a particular mutation
#this is used to filter the variants later

args = commandArgs(trailingOnly = T)

library(data.table)
library(tidyverse)

df1 = fread(args[1])
df2 = fread(args[2])

outp = gsub('genebe.tsv.gz','vepiv.tsv',args[1])

colnames(df2)[1] = 'Uploaded_variation'

ref = df1 %>% mutate(Location = paste0(CHROM,':',POS)) %>% dplyr::select(Location, REF)

#classification based on https://doi.org/10.1038/s41593-019-0530-0

df_out = df2 %>% 
	separate(Uploaded_variation, into = c('CHROM','POS','ALT','REF'), remove = F) %>%
    mutate(ID = Uploaded_variation,
           p2 = POS,
            HGVSc = gsub('.*:','',HGVSc),
            HGVSp = gsub('.*:','',HGVSp),
	consq_level = ifelse(str_detect(Consequence, 'feature_truncation|frameshift_truncation|inframe_deletion|conservative_inframe_deletion|disruptive_inframe_deletion|stop_gained|stop_gained_NMD_escaping|stop_gained_NMD_triggering|splice_donor_variant|splice_acceptor_variant|frameshift_variant|inframe_insertion|inframe_deletion'),1,3),
	consq_level = ifelse(consq_level == 3 & str_detect(Consequence, 'inframe_deletion|inframe_insertion|missense_variant|stop_lost|start_lost|protein_altering_variant'), 2, consq_level)
	) %>%
            group_by(ID) %>%
	slice(which.min(consq_level)) %>%
	unique %>%
            mutate(Consequence = paste(unique(Consequence), collapse = ','),
            HGVSc = paste(unique(HGVSc), collapse = ','),
            HGVSp = paste(unique(HGVSp), collapse = ','),
        Uploaded_variation = paste(unique(Uploaded_variation), collapse = ','),
            LoF = gsub('-,','', LoF),
            Consequence = gsub(',-','', Consequence),
            Consequence = gsub('-,','', Consequence),
            HGVSc = gsub(',-','', HGVSc),
            HGVSc = gsub('-,','', HGVSc),
            HGVSp = gsub(',-','', HGVSp),
            HGVSp = gsub('-,','', HGVSp)) %>%
        unique %>%
ungroup %>%
        select(-REF, -CHROM, -POS, -ALT) %>%
     left_join(mutate(df1,ID = paste0(CHROM,'_',POS,'_',REF,'_',ALT)), by = 'ID') %>%
     dplyr::select(`#CHROM` = CHROM, POS, p2, REF, ALT, ID, Consequence, Ref.Gene, vep_snp = Uploaded_variation, ACMG, HGVSc, HGVSp, consq_level)  %>%
     arrange(as.numeric(POS)) %>% 
     arrange(as.numeric(`#CHROM`)) %>% 
     drop_na %>% 
     unique

colnames(df_out)[1] = 'CHROM'

write_tsv(df_out, outp)
