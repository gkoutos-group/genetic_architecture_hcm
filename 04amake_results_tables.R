#This script generates tables of mutations for an individual population, cases and controls
#Several tables are generated as a result, summarising the findings

library(data.table)
library(tidyverse)
library(R.utils)
set.seed(1234)

#select population
dat = 'ukb'
#dat = 'chinese'

#read in genetic data
df_case = fread(paste0(dat,'_cases.tsv.gz')) 
df_control = fread(paste0('filtered_inputs/',dat,'_controls.tsv.gz'))

#update column names
colnames(df_case) = gsub('# ','', colnames(df_case))
colnames(df_case) = sub('-','_', colnames(df_case))
colnames(df_case) = gsub('\\[[0-9]*\\]','', colnames(df_case))
colnames(df_case) = gsub(':GT','', colnames(df_case))
colnames(df_case)[6:ncol(df_case)] = paste0('x',colnames(df_case)[6:ncol(df_case)])

colnames(df_control) = gsub('# ','', colnames(df_control))
colnames(df_control) = sub('-','_', colnames(df_control))
colnames(df_control) = gsub('\\[[0-9]*\\]','', colnames(df_control))
colnames(df_control) = gsub(':GT','', colnames(df_control))
colnames(df_control)[6:ncol(df_control)] = paste0('x',colnames(df_control)[6:ncol(df_control)])

df = df_case %>%
    full_join(df_control, by = 'ID')


#read in case ids, others are controls
cases = colnames(df_case)[6:ncol(df_case)]

#controls are selected then, in the case of UK bb, sampled for now to match numbers
controls = colnames(df_control)[6:ncol(df_control)]

if (dat == 'ukb'){
  controls = sample(controls, 1000, replace = FALSE)
}

#remove unnecessary columns
df = df %>% dplyr::select(ID, all_of(cases), all_of(controls))

#loop to find how many alt alleles exist in data, this helps remove any that are not there
sm = c()

for (i in 1:nrow(df)){ 
    sm = c(sm,
           sum(as.numeric(unlist(df[i,2:ncol(df)])), na.rm = T)) 
}

df$sum = sm

write_tsv(df, paste0(dat,'gt.tsv'))
write_tsv(data.frame(id = c(controls,cases), status = c(rep('control',length(controls)),rep('case',length(cases)))),
          paste0(dat, 'status.tsv'))


#read in annotations
df2_cases = fread(paste0("annos/",dat,"_cases_annotations_vepiv.tsv"))
if (dat == 'ukb'){
  df2_controls = fread(paste0("annos/",dat,"_controls2_annotations_vepiv.tsv"))
} else {
  df2_controls = fread(paste0("annos/",dat,"_controls_annotations_vepiv.tsv"))
}

df2 = rbind(df2_cases, df2_controls)

#remove any variants that have 0 alt alleles present in samples and add annotations
df = df %>% 
    left_join(df2, by = 'ID') %>% 
    filter(consq_level < 3)

#update gnomaAD field to numeric and remove '-' as NA then filter for rare mutations
#change to long table for easier analysis
#define variant type as truncating or non-truncating, this is based on Consequence
#define zygosity and case/control status
#classification of variants according to - https://doi.org/10.1038/s41593-019-0530-0
#gnomad files are filtered using tabix by the genetic ranges used in this research, first is exomes, second is from genomes
frq = fread('gnomad_v4.txt') %>% 
  select(ID, ex_gnomAD_faf95 = fafmax_faf95_max)

frq_genomes = fread('gen_gnomad_v4.txt') %>% 
  select(ID, gen_gnomAD_faf95 = fafmax_faf95_max)

frq = frq %>%
  full_join(frq_genomes, by = 'ID') %>%
  mutate(ex_gnomAD_faf95 = ifelse(is.na(ex_gnomAD_faf95)|ex_gnomAD_faf95 == '.', 0, ex_gnomAD_faf95),
         gen_gnomAD_faf95 = ifelse(is.na(gen_gnomAD_faf95)|gen_gnomAD_faf95 == '.', 0, gen_gnomAD_faf95),
         ex_gnomAD_faf95 = as.numeric(ex_gnomAD_faf95),
         gen_gnomAD_faf95 = as.numeric(gen_gnomAD_faf95),
         gnomAD_faf95 = ifelse(ex_gnomAD_faf95 >= gen_gnomAD_faf95, ex_gnomAD_faf95, gen_gnomAD_faf95))

df = df %>%
    mutate(across(everything(), as.character)) %>%
    select(-gnomAD_faf95) %>%
    left_join(frq) %>%
    mutate(gnomAD_faf95 = ifelse(gnomAD_faf95 == '.', '0', gnomAD_faf95)) %>%
    mutate(gnomAD_faf95 = as.numeric(gnomAD_faf95)) %>% 
    filter(gnomAD_faf95 < 4e-5 | is.na(gnomAD_faf95)) %>% 
    pivot_longer(cols = starts_with("x"), names_to = 'sid', values_to = 'gt', values_drop_na = T) %>%
    filter(gt > 0) %>%
    mutate(variant_type = ifelse(str_detect(Consequence, 'feature_truncation|frameshift_truncation|inframe_deletion|conservative_inframe_deletion|disruptive_inframe_deletion|stop_gained|stop_gained_NMD_escaping|stop_gained_NMD_triggering|splice_donor_variant|splice_acceptor_variant|frameshift_variant|inframe_insertion|inframe_deletion'),'Truncating','Non-truncating'), # %in% c('feature_truncation','frameshift_truncation','inframe_deletion','conservative_inframe_deletion','disruptive_inframe_deletion','stop_gained','stop_gained_NMD_escaping','stop_gained_NMD_triggering')
	       zygosity = ifelse(gt == 2, 'Hom',''),
		   zygosity = ifelse(gt == 1, 'Het', zygosity),
		   status = ifelse(sid %in% cases, 'case', 'control')
        ) %>%
  unique 

#count case/control numbers per mutation
status_count = df %>% 
    group_by(ID, status) %>% 
    count() 

#add totals and calculate frequency
df = df %>% 
    left_join(status_count, by = c('ID', 'status')) %>%
	filter((gnomAD_faf95 < 4e-5 | is.na(gnomAD_faf95)) & zygosity != '') %>%
	mutate(freq = ifelse(status == 'case', n/length(cases), n/length(controls)))


#create supplementary table 4
s4 = df %>% 
    filter(status == 'case') %>% 
    select(sid, consq_level, ID, Gene = Ref.Gene, cDNA_variant = HGVSc, protein_variant = HGVSp, variant_consequence = Consequence, gnomAD_faf95, hcm_count = n, hcm_freq = freq, ACMG)

write_tsv(s4, paste0(dat,'supp_4_ids.tsv'))

#create supplementary table 5
s5 = df %>% 
    filter(status == 'control') %>% 
    select(sid, consq_level, ID, Gene = Ref.Gene, cDNA_variant = HGVSc, protein_variant = HGVSp, variant_consequence = Consequence, gnomAD_faf95, control_count = n, control_freq = freq, ACMG)

write_tsv(s5, paste0(dat,'supp_5_ids.tsv'))

#create supplementary table 6
s6 = df %>% 
    select(Gene = Ref.Gene, variant_type, status, sid) %>% 
    unique() %>% 
    select(-sid) %>% 
    group_by(Gene, variant_type, status) %>% 
    count()

s6 = rbind(
  s6,
  data.frame(
    Gene = c('Thin Filament','Thin Filament','Thin Filament','Thin Filament','Myosin light chain (MLC)','Myosin light chain (MLC)','Myosin light chain (MLC)','Myosin light chain (MLC)'),
    variant_type = c('Non-truncating','Truncating','Non-truncating','Truncating','Non-truncating','Truncating','Non-truncating','Truncating'),
    status = c('case','case','control','control','case','case','control','control'),
    n = c(
      sum(s6$n[s6$variant_type == 'Non-truncating' & s6$status == 'case' & s6$Gene %in% c('TNNI3', 'TNNT2', 'TNNC1', 'TPM1', 'ACTC1')]),
      sum(s6$n[s6$variant_type == 'Truncating' & s6$status == 'case' & s6$Gene %in% c('TNNI3', 'TNNT2', 'TNNC1', 'TPM1', 'ACTC1')]),
      sum(s6$n[s6$variant_type == 'Non-truncating' & s6$status == 'control' & s6$Gene %in% c('TNNI3', 'TNNT2', 'TNNC1', 'TPM1', 'ACTC1')]),
      sum(s6$n[s6$variant_type == 'Truncating' & s6$status == 'control' & s6$Gene %in% c('TNNI3', 'TNNT2', 'TNNC1', 'TPM1', 'ACTC1')]),
      sum(s6$n[s6$variant_type == 'Non-truncating' & s6$status == 'case' & s6$Gene %in% c('MYL2','MYL3')]),
      sum(s6$n[s6$variant_type == 'Truncating' & s6$status == 'case' & s6$Gene %in% c('MYL2','MYL3')]),
      sum(s6$n[s6$variant_type == 'Non-truncating' & s6$status == 'control' & s6$Gene %in% c('MYL2','MYL3')]),
      sum(s6$n[s6$variant_type == 'Truncating' & s6$status == 'control' & s6$Gene %in% c('MYL2','MYL3')])
    )
  )
)

#create template to join to 
s6dummy = data.frame(Gene = rep(c('MYH7','MYBPC3','Myosin light chain (MLC)',
                                  'Thin Filament', 'CSRP3', 'JPH2', 'ACTN2', 
                                  'PLN', 'TRIM63', 'ALPK3', 'FHOD3'), 4), 
                     variant_type = c(rep('Truncating', 22),
                                      rep('Non-truncating', 22 )), 
                     status = rep(
                       c(
                         rep('case', 11),
                         rep('control', 11)
                         ), 
                       2)
                     )

s6 = s6dummy %>% 
    left_join(s6, by = c('Gene', 'variant_type', 'status')) %>% 
    mutate(n = ifelse(is.na(n), 0, n)) %>% 
    pivot_wider(names_from = status, values_from = n) %>% 
    filter(!is.na(case) & !is.na(control)) %>% 
    arrange(Gene) %>% 
    mutate(case_without = length(cases) - case, control_without = length(controls) - control)

#iteratively carry out fisher's test
ft = c()

for (i in 1:nrow(s6)){

ft = c(ft,
         fisher.test(x = matrix(c(unlist(s6[i,3]), unlist(s6[i,5]), unlist(s6[i,4]), unlist(s6[i,6])), nrow = 2), alternative = 'greater')$p.value
	)

}

s6$fishers_test = ft

write_tsv(s6, paste0(dat,'supp_6.tsv'))

#create supplementary table 8
s8 = df %>% 
    select(sid, ID, Gene = Ref.Gene, cDNA_variant = HGVSc, protein_variant = HGVSp, zygosity, variant_consequence = Consequence, ACMG, gnomAD_faf95) %>% 
    unique %>% 
    group_by(sid) %>% 
    add_count() %>% 
    filter(n > 1) %>% 
    select(-n) %>% 
    arrange(sid)

#using anonymised patient ids
for (i in 1:length(unique(s8$sid))){

	s8$sid[s8$sid == unique(s8$sid)[i]]  = i

}

write_tsv(s8, paste0(dat,'supp_8.tsv'))







