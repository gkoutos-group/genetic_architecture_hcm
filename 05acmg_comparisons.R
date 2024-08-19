#this script is for clarity of the analysis. Please be aware that cardioclassifier data was fairly crudely extracted from the website
library(tidyverse)
library(data.table)

#import cardioclassifier 
for (i in c('ukcases','ukcontrols','chinesecases', 'chinesecontrols')){
  
  if(i == 'ukcases'){
    cdf = fread(paste0('../small_cc_',i,'.txt'), header = F, col.names = c('cc_acmg', 'gene','consequence','hgvs','prot','pos','q','ab')) %>%
      select(pos, cc_acmg, gene, pos)
    cdf$cohort = 'uk_cases'
  } else {
    df2 = fread(paste0('../small_cc_',i,'.txt'), header = F, col.names = c('cc_acmg', 'gene','consequence','hgvs','prot','pos','q','ab')) %>%
      select(pos, cc_acmg, gene, pos)
    df2$cohort = gsub('uk','uk_',i)
    df2$cohort = gsub('chinese','chinese_',df2$cohort)
    
    cdf = rbind(cdf, df2)
  }
  
  
}

#hg38 coords and alleles
for (i in c('ukb_cases','ukb_controls2','chinese_cases', 'chinese_controls')){
  
  if(i == 'ukb_cases'){
    hg38_lu = fread(paste0('../hg19/',i,'_hg19.tsv'), header = F, col.names = c('chr','pos_hg19','chr_hg38','pos_hg38','ref','alt')) 
    hg38_lu$cohort = 'uk_cases'
  } else {
    hg38_lu2 = fread(paste0('../hg19/',i,'_hg19.tsv'), header = F, col.names = c('chr','pos_hg19','chr_hg38','pos_hg38','ref','alt')) 
    hg38_lu2$cohort = gsub('ukb','uk',i)
    hg38_lu2$cohort = gsub('2','',hg38_lu2$cohort)
    hg38_lu2$cohort = gsub('chinese','chinese',hg38_lu2$cohort)
    
    hg38_lu = rbind(hg38_lu, hg38_lu2)
  }
  
  
}

hg38_lu = hg38_lu %>% select(-chr_hg38)

cdf = cdf %>%
  separate(pos, sep = ":", into = c("chr","pos_hg19")) %>%
  mutate(pos_hg19 = as.integer(pos_hg19)) %>%
  left_join(hg38_lu, by = c('chr' = 'chr', 'pos_hg19' = 'pos_hg19', 'cohort' = 'cohort'))

df = cdf

#get intervar results

int_files = list.files('../basic_intervar/')

for (i in int_files){
  
  if(i == int_files[1]){
    
    iv = fread(paste0('../basic_intervar/',i))
    colnames(iv) = gsub('[# :]','',colnames(iv))
    iv = iv %>%
      select(chr = Chr, pos = Start, ref = Ref, alt = Alt, acmg_iv = InterVarInterVarandEvidence) %>%
      mutate(acmg_iv = gsub(' PVS.*','',acmg_iv),
             acmg_iv = gsub('InterVar: ','',acmg_iv),
             cohort = gsub('_ivin.*','',i),
             cohort = gsub('ukb','uk',cohort),
             cohort = gsub('2','',cohort))
    
  } else {
    
    iv2 = fread(paste0('../basic_intervar/',i))
    colnames(iv2) = gsub('[# :]','',colnames(iv2))
    iv2 = iv2 %>%
      select(chr = Chr, pos = Start, ref = Ref, alt = Alt, acmg_iv = InterVarInterVarandEvidence) %>%
      mutate(acmg_iv = gsub(' PVS.*','',acmg_iv),
             acmg_iv = gsub('InterVar: ','',acmg_iv),
             cohort = gsub('_ivin.*','',i),
             cohort = gsub('ukb','uk',cohort),
             cohort = gsub('2','',cohort))
    
    iv = rbind(iv, iv2)
  }
  
}

iv$chr = paste0('chr',iv$chr)

df = df %>%
  full_join(iv, by = c('chr' = 'chr', 'pos_hg38' = 'pos', 'ref' = 'ref', 'alt' = 'alt', 'cohort' = 'cohort'))


#bring in genebe

gb_files = list.files('../annos/')

for (i in gb_files){
  
  if(i == gb_files[1]){
    
    gb = fread(paste0('../annos/',i))

    gb = gb %>%
      select(chr = CHROM, pos = POS, ref = REF, alt = ALT, acmg_gb = ACMG, gene = Ref.Gene) %>%
      mutate(cohort = gsub('_anno.*','',i),
             cohort = gsub('ukb','uk',cohort),
             cohort = gsub('2','',cohort))
    
  } else {
    
    gb2 = fread(paste0('../annos/',i))
   
    gb2 = gb2 %>%
      select(chr = CHROM, pos = POS, ref = REF, alt = ALT, acmg_gb = ACMG, gene = Ref.Gene) %>%
      mutate(cohort = gsub('_anno.*','',i),
             cohort = gsub('ukb','uk',cohort),
             cohort = gsub('2','',cohort))
    
    gb = rbind(gb, gb2)
  }
  
}

gb$chr = paste0('chr',gb$chr)

df = df %>%
  select(-gene) %>%
  full_join(gb, by = c('chr' = 'chr', 'pos_hg38' = 'pos', 'ref' = 'ref', 'alt' = 'alt', 'cohort' = 'cohort'))

df_all = df %>%
  select(chr, pos_hg38,cohort,ref,alt, gene, acmg_gb,cc_acmg,acmg_iv) %>%
  mutate(acmg_gb = gsub('_',' ',tolower(acmg_gb)),
         acmg_iv = tolower(acmg_iv),
         cc_acmg = tolower(cc_acmg)) %>%
  filter(alt != '*')

df_nona = df_all %>% drop_na() %>% select(chr, pos_hg38, ref,alt, acmg_gb,cc_acmg,acmg_iv) %>% unique
df_nocc = df_all %>% select(chr, pos_hg38, ref,alt, acmg_gb,acmg_iv) %>% unique() 

per_cl = rbind(
  df_nona %>% select(cc_acmg) %>% mutate(num = 1) %>% group_by(cc_acmg) %>% summarise(n = sum(num)) %>% rename('acmg' = 1) %>% mutate(classifier = 'CardioClassifier'),
  df_nona %>% select(acmg_iv) %>% mutate(num = 1) %>% group_by(acmg_iv) %>% summarise(n = sum(num)) %>% rename('acmg' = 1) %>% mutate(classifier = 'InterVar'),
  df_nona %>% select(acmg_gb) %>% mutate(num = 1) %>% group_by(acmg_gb) %>% summarise(n = sum(num)) %>% rename('acmg' = 1) %>% mutate(classifier = 'genebe')
)

per_cl2 = rbind(
  df_nocc %>% select(acmg_iv) %>% mutate(num = 1) %>% group_by(acmg_iv) %>% summarise(n = sum(num)) %>% rename('acmg' = 1) %>% mutate(classifier = 'InterVar'),
  df_nocc %>% select(acmg_gb) %>% mutate(num = 1) %>% group_by(acmg_gb) %>% summarise(n = sum(num)) %>% rename('acmg' = 1) %>% mutate(classifier = 'genebe')
)




