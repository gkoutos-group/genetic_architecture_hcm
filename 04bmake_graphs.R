library(data.table)
library(tidyverse)
library(webr)
library(ggsunburst)
library(ggrepel)


ch_cases = fread('chinesesupp_4_ids.tsv') %>%
  mutate(grp = 'chinese_cases')
colnames(ch_cases)[c(9,10)] = gsub('hcm_','', colnames(ch_cases)[c(9,10)] )

ch_controls = fread('chinesesupp_5_ids.tsv') %>%
  mutate(grp = 'chinese_controls')
colnames(ch_controls)[c(9,10)] = gsub('control_','', colnames(ch_controls)[c(9,10)] )

uk_cases = fread('ukbsupp_4_ids.tsv') %>%
  mutate(grp = 'uk_cases')
colnames(uk_cases)[c(9,10)] = gsub('hcm_','', colnames(uk_cases)[c(9,10)] )

uk_controls = fread('ukbsupp_5_ids.tsv') %>%
  mutate(grp = 'uk_controls')
colnames(uk_controls)[c(9,10)] = gsub('control_','', colnames(uk_controls)[c(9,10)] )

df = rbind(
  ch_cases,
  ch_controls,
  uk_cases,
  uk_controls
)

n = c(chinese_cases = 643, chinese_controls = 500, uk_cases = 406, uk_controls = 1000)

fig1a_df = data.frame(group = factor(c("Chinese HCM", "Chinese Control", "UK HCM", "UK Control"),
                                     levels = c("Chinese HCM", "Chinese Control", "UK HCM", "UK Control")),
                      proportion = c(length(unique(ch_cases$sid))/n['chinese_cases'], length(unique(ch_controls$sid))/n['chinese_controls'],
                                     length(unique(uk_cases$sid))/n['uk_cases'], length(unique(uk_controls$sid))/n['uk_controls']) * 100)



ggplot(fig1a_df, aes(x = group, y = proportion)) +
  geom_bar(stat = 'identity', fill = 'slateblue') +
  ggtitle('Yield of Rare Variants Across Validated HCM Genes') +
  xlab('Group') +
  ylab('Percentage of Individuals with Rare Variant') +
  theme_light()
ggsave('fig1a.png', dpi = 600)

prop_per_gene = function(gen, typ, cohort, status){
  num = df %>%
    filter(Gene %in% gen, consq_level == typ, grp == paste0(cohort,'_',status)) %>%
    select(sid) %>%
    unique %>%
    nrow() / n[paste0(cohort,'_',status)]
  
  return(num)
}

MLC = c('MYL2','MYL3')
TF = c('ACTC1','TNNC1','TNNI3','TNNT2','TPM1')

gene_list = c()
prop_list = c()
fill_list = c()
pop_list = c()

for (g in c("ALPK3","MYBPC3","MYH7","ACTN2","JPH2","CSRP3","PLN","FHOD3","TRIM63","MLC","Thin Filament")){
  for (s in c('cases', 'controls')){
    for (c in c('chinese','uk')){
      for (t in c(1,2)){
        
        gene_list = c(gene_list, g)
        fill_list = c(fill_list, paste0(c,t))
        if(c == 'uk'){
          c_cap = 'UK'
        } else {
          c_cap = str_to_title(c)
        }
        pop_list = c(pop_list, paste0(c_cap,' ',str_to_title(s)))
        
        if (g == 'MLC'){
          
          gene = MLC
          
        } else if (g == 'Thin Filament'){
          
          gene = TF
        
        } else {
          gene = g
        }
        
        prop = prop_per_gene(gene, t, c, s)
        
        prop_list = c(prop_list, prop)
        
      }
    }
  }
  
}

fig1b_df = data.frame(
                      gene = gene_list,
                      pop = pop_list,
                      type = fill_list,
                      prop = prop_list
                      ) 
gene_order = fig1b_df %>%
  group_by(gene) %>% 
  summarise(n = sum(prop)) %>%
  arrange(desc(n)) %>%
  filter(n > 0) %>%
  select(gene) %>% 
  unique %>%
  unlist %>%
  as.character

fig1b_df = fig1b_df %>%
  filter(gene %in% gene_order)

fig1b_df$gene = factor(fig1b_df$gene, levels = gene_order)

pal = c("#CC79A7","firebrick3","#56B4E9", "#0072B2")

ggplot(fig1b_df, aes(x = pop, y = prop * 100, fill = type )) +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_grid(~ gene) +
  geom_text_repel(aes(label = ifelse(prop > 0, signif(100 * prop, 2),'')), position = position_stack(vjust = 0.5),
                  direction="y", hjust=0.5, size = 3, box.padding=0.1, fontface = 'bold', color = 'black') +
  ggtitle('Proportion of Individuals with Potentially Damaging Mutations in HCM Genes') +
  xlab('Group') +
  ylab('Percentage of Individuals with a Rare Variant') +
  theme_light() +
  scale_fill_manual(values=pal, labels = c('Truncating','Non-Truncating','Truncating','Non-Truncating')) +
  theme(axis.text.x = element_text(angle = -75, size=7)) +
  guides(fill=guide_legend(title="Mutation Type"))

ggsave('fig1b.png', width = 12, height = 6, dpi = 600)                                  

ts = c('Chinese Cases', 'Chinese Controls', 'UK Cases', 'UK Controls')
names(ts) = c("chinese_cases","chinese_controls","uk_cases","uk_controls")

for (i in c("chinese_cases","chinese_controls","uk_cases","uk_controls")){

ttl = as.character(ts[i])

fig2 = df %>% 
  select(sid, grp, Gene, ACMG) %>% 
  unique %>% 
  group_by(ACMG, Gene) %>%
  mutate(tal = 1) %>%
  mutate(Gene = ifelse(Gene %in% TF, "TF", Gene),
         Gene = ifelse(Gene %in% MLC, "MLC", Gene))

fig2 = fig2 %>%
  filter(grp == i) %>%
  summarise(n = sum(tal))
                    
                    
#https://stackoverflow.com/questions/50004058/multiple-dependent-level-sunburst-doughnut-chart-using-ggplot2    
fig2 = fig2 %>%
  mutate(ACMG = as.factor(ACMG) %>% fct_reorder(n, sum)) %>%
  arrange(ACMG, n) %>%
  mutate(type = as.factor(Gene) %>% fct_reorder2(Gene, n),
         ACMG = gsub('_',' ', ACMG))

lvl0 = tibble(name = "Parent", value = 0, level = 0, fill = '')

lvl1 = fig2 %>%
  select(name = ACMG, Gene, n, type) %>%
  group_by(name) %>%
  summarise(value = sum(n)) %>%
  ungroup() %>%
  mutate(level = 1) %>%
  mutate(fill = name)

lvl2 = fig2 %>%
  select(name = Gene, value = n, fill = ACMG) %>%
  mutate(level = 2)

colour_list = c('#66A61E','#1B9E77','#666666','#D95F02','#E7298A')
names(colour_list) = c('Benign','Likely benign','Uncertain significance','Likely pathogenic','Pathogenic')


Data = bind_rows(lvl0, lvl1, lvl2) %>%
  mutate(name = as.factor(name) %>% fct_reorder2(fill, value)) %>%
  mutate(fill = ifelse(name == 'Parent', NA, fill)) %>%
  arrange(fill, name) %>% 
  mutate(level = as.factor(level),
         name = ifelse(level == 2 & fill %in% c("Benign","Likely benign"),'', as.character(name)),
         name = as.factor(name)) 
  
  
ggplot(Data, aes(x = level, y = value, fill = fill, alpha = level)) +
  geom_col(width = 1, color = "gray90", size = 0.25, position = position_stack()) +
  geom_text_repel(aes(label = name), size = 3, position = position_stack(vjust = 0.5), box.padding = 0., fontface = 'bold', max.overlaps = 10) +
  coord_polar(theta = "y") +
  scale_alpha_manual(values = c("0" = 0, "1" = 1, "2" = 0.7), guide = F) +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_manual(values = colour_list, guide = F) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  ggtitle(paste0('Types of Mutation Found in ',ttl))

ggsave(paste0('fig2_',i,'.png'), dpi = 600)

}




