#This script was used to visualize less stringent minor allele thresholds

library(data.table)
library(tidyverse)

chinese_gt = fread('chinesegt.tsv') %>% select(-sum) %>% as.data.frame

uk_gt = fread('ukbgt.tsv')  %>% select(-sum) %>% as.data.frame

status = fread('ukbstatus.tsv') %>% mutate(cohort = 'uk') %>%
  rbind(
    fread('chinesestatus.tsv') %>% mutate(cohort = 'chinese')
  )

all_genes = c('MYH7','MYBPC3','CSRP3','JPH2','ACTN2','PLN','TRIM63','ALPK3','FHOD3','MYL2','MYL3','TNNI3','TNNT2','TNNC1','TPM1','ACTC1')

chinese_gt = chinese_gt %>%
  pivot_longer(!ID, names_to = 'id', values_to = 'gt_bin') %>%
  filter(!is.na(gt_bin)) %>%
  mutate(gt_bin = ifelse(gt_bin > 0, 1, 0)) %>%
  filter(gt_bin > 0) %>%
  select(-gt_bin) %>% 
  unique %>%
  left_join(status, by = 'id') 

uk_gt = uk_gt %>%
  pivot_longer(!ID, names_to = 'id', values_to = 'gt_bin') %>%
  filter(!is.na(gt_bin)) %>%
  mutate(gt_bin = ifelse(gt_bin > 0, 1, 0)) %>%
  filter(gt_bin > 0) %>%
  select(-gt_bin) %>% 
  unique %>%
  left_join(status, by = 'id') 

all_gt = rbind(chinese_gt, uk_gt)

for (f in list.files(pattern = '*.tsv.gz')){
  
  print(f)
  
  if(f != list.files(pattern = '*.tsv.gz')[1]){
    
    df = df %>%
      rbind(
        fread(f, header = T, skip = 'chr') %>%
          mutate(group = 
                   gsub('_fulltables_genebe.tsv.gz','',f)
                 )
      )
    
  } else {
  
    df = fread(f, header = T, skip = 'chr') %>%
      mutate(group = 
               gsub('.tsv.gz','',f)
      )
      
  }
  
}

df = df %>% 
  mutate(ID = paste(chr,pos,ref,alt, sep = "_")) %>%
  select(ID, gene_symbol, pos, consequences, gnomad_exomes_af, gnomad_genomes_af, clinvar_classification) %>%
  unique() %>%
  mutate(consq_clin = ifelse(
    str_detect(consequences, 'conservative_inframe_deletion|conservative_inframe_insertion|disruptive_inframe_deletion|disruptive_inframe_insertion|exon_loss_variant|missense_variant|frameshift_variant|exon_loss_variant|splice_acceptor_variant|splice_donor_variant|start_lost|stop_gained|stop_lost')
  | (str_detect(clinvar_classification,'athogenic') & clinvar_classification != 'Conflicting classifications of pathogenicity' & !is.na(clinvar_classification))
  ,1,0
  )
  ) %>% 
  mutate(max_maf = pmax(gnomad_exomes_af,gnomad_genomes_af,na.rm = T)) %>% 
  filter(consq_clin == 1 & max_maf < 0.0005)



fig_cases = all_gt %>%
  left_join(df, by = c('ID' = 'ID')) %>% 
  filter(!is.na(pos) & status == 'case' & gene_symbol %in% all_genes) %>% 
  select(id, cohort, Gene = gene_symbol) 

templ = data.frame(
  Gene = rep(all_genes,2),
  cohort = c(rep('chinese', length(all_genes)),rep('uk', length(all_genes)))
)

fig_cases = templ %>%
  left_join(fig_cases, by = c('Gene' = 'Gene', 'cohort' = 'cohort')) %>%
  unique %>% 
  group_by(Gene,cohort) %>%
  mutate(num = 1) %>%
  summarise(num = sum(num)) %>% 
  ungroup %>%
  mutate(num = ifelse(cohort == 'chinese', num / 643 * 100, num),
         num = ifelse(cohort == 'uk', num / 406 * 100, num)) %>%
  group_by(Gene) %>%
  mutate(sum = sum(num))

fig_cases %>%
  ungroup %>%
  mutate(Gene = as.factor(Gene) %>% fct_reorder(desc(sum)),
         cohort = factor(cohort, levels = c('chinese','uk'))) %>%
  ggplot(aes(y = num, x = Gene, fill = cohort)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  ggtitle('Rare Damaging Mutations in HCM Cases (MAF<5e-4)') +
  ylab('Percent of Cohort') +
  scale_fill_manual(values = c("firebrick3","#0072B2"),labels = c('chinese' = "Chinese", 'uk' =  "UK"), name = 'Cohort') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -45, size=7))

ggsave('rare_genes_cases.png', dpi = 600)


