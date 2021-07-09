# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Boxplot of homozygous and heterozygous proportions
# # 
# Author: Linnéa Smeds 
# Date: 6 July 2021

#Clear the workspace
rm(list=ls())

# # # # 
# Load packages
require(vcfR)
require(tidyverse)
require(tictoc)
require(viridis)
require(gtable)
require(grid)

## Reading in the data
vcf_file<-"data/100S95F15R.vepfinal.withMaleFounders.chr1-38.vcf"
pop_file<-"data/5TempClasses11I95F15R.pop.list"
founder_file<-"data/founders.list"
gen_file<-"data/generation.list"
vep_file<-"data/100S95F15R.filt1.vepTypes.chr1-38.txt"
vcf <- read.vcfR(vcf_file)
pop_tib<-pop_file %>% read.table() %>% as_tibble() %>% rename(Indiv=V1, Pop=V2)
gen_tib<-gen_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Generation=generation)
fou_tib <- founder_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Founder=Ftype)
vep_tib <-vep_file %>% read.table(header=TRUE) %>% as_tibble()

## Convert to tidy format
tidy_vcf <- vcfR2tidy(vcf, 
                      info_fields=c("AA"), 
                      format_fields=c("GT"), 
                      dot_is_NA=TRUE)


## Extracting relevant data and inverting sites with alternative polarization
tidy_gt3<-tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>% 
  inner_join(tidy_vcf$fix)  %>% 
  mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
  group_by(Indiv, new_gt) %>% summarize(count=n()) %>% ungroup() %>%
  mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
  pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count)

# Add population names and some extra columns to facilitate plotting!
final_gt<- tidy_gt3 %>% mutate(sum=gt00+gt01+gt11) %>% inner_join(pop_tib) %>% mutate(hetfrq=gt01/sum, homderfrq=gt11/sum, deralfreq=(gt11*2+gt01)/(sum*2),Grp=case_when(Pop=='Finland' | Pop=='Russia' ~ Pop, TRUE ~ 'Scandinavia'))  %>% left_join(fou_tib) %>% left_join(gen_tib) 
# Rearrange the population in a custom order
final_gt$Pop <- factor(final_gt$Pop, levels = c("1983-1990", "1991-1998", "1999-2006", "2007-2014S", "2007-2014I", "Immigrant", "Finland", "Russia"))


######################################################################################################
## Boxplots with founder and generation info

# HETEROZYGOUS FREQUENCY
p1<- ggplot(final_gt,aes(Pop, hetfrq))+
    geom_violin(aes(fill=Pop))+
    scale_fill_viridis(discrete=TRUE, guide=FALSE) +
    facet_grid(~factor(Grp, level=c("Scandinavia", "Finland", "Russia")), scales="free_x", space="free_x") +
    labs(x=NULL,y="proportion heterozygous sites")+
    theme(panel.grid.major = element_line(colour = 'white'),
          panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
          axis.text.x = element_text(angle = 45, hjust=1),
          legend.box = "vertical") +
    geom_point(data=filter(final_gt, !is.na(Founder)),aes(shape=Founder), size=3, col='darkgrey', fill='lightyellow')+
    scale_shape_manual(values=c(22,24))+
    theme(legend.key=element_rect(fill='#3acf73', colour='NA')) +
    geom_point(data=filter(final_gt, !is.na(Generation)), aes(color=Generation)) + 
    scale_colour_continuous(low='#fffee3', high='#0214a3', breaks = c(1,2,3,4,5,6), labels = c(1,2,3,4,5,6))
plot_tab1 <- ggplotGrob(p1)
plot_filtered1 <- gtable_filter(plot_tab, 
                                 "(background|panel|strip-t|axis-l|xlab|ylab|guide-box|title|subtitle|caption|tag|axis-b-1)",
                                 trim=FALSE)
grid.newpage()
grid.draw(plot_filtered1)

# HOMOZYGOUS DERIVED FREQUENCY
p2<- ggplot(final_gt,aes(Pop, homderfrq))+
    geom_violin(aes(fill=Pop))+
    scale_fill_viridis(discrete=TRUE, guide=FALSE) +
    facet_grid(~factor(Grp, level=c("Scandinavia", "Finland", "Russia")), scales="free_x", space="free_x") +
    labs(x=NULL,y="proportion homozygous derived")+
    theme(panel.grid.major = element_line(colour = 'white'),
          panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
          axis.text.x = element_text(angle = 45, hjust=1),
          legend.box = "vertical") +
    geom_point(data=filter(final_gt, !is.na(Founder)),aes(shape=Founder), size=3, col='darkgrey', fill='lightyellow')+
    scale_shape_manual(values=c(22,24))+
    theme(legend.key=element_rect(fill='#3acf73', colour='NA')) +
    geom_point(data=filter(final_gt, !is.na(Generation)), aes(color=Generation)) + 
    scale_colour_continuous(low='#fffee3', high='#0214a3', breaks = c(1,2,3,4,5,6), labels = c(1,2,3,4,5,6))
plot_tab2 <- ggplotGrob(p2)
plot_filtered2 <- gtable_filter(plot_tab2, 
                                 "(background|panel|strip-t|axis-l|xlab|ylab|guide-box|title|subtitle|caption|tag|axis-b-1)",
                                 trim=FALSE)
grid.newpage()
grid.draw(plot_filtered2)
  
# FREQ OF DERIVED ALLELES 
p3<- ggplot(final_gt,aes(Pop, deralfreq))+
  geom_violin(aes(fill=Pop))+
  scale_fill_viridis(discrete=TRUE, guide=FALSE) +
  facet_grid(~factor(Grp, level=c("Scandinavia", "Finland", "Russia")), scales="free_x", space="free_x") +
  labs(x=NULL,y="derived allele frequency")+
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.box = "vertical") +
  geom_point(data=filter(final_gt, !is.na(Founder)),aes(shape=Founder), size=3, col='darkgrey', fill='lightyellow')+
  scale_shape_manual(values=c(22,24))+
  theme(legend.key=element_rect(fill='#3acf73', colour='NA')) +
  geom_point(data=filter(final_gt, !is.na(Generation)), aes(color=Generation)) + 
  scale_colour_continuous(low='#fffee3', high='#0214a3', breaks = c(1,2,3,4,5,6), labels = c(1,2,3,4,5,6))
plot_tab3 <- ggplotGrob(p3)
plot_filtered3 <- gtable_filter(plot_tab3, 
                                "(background|panel|strip-t|axis-l|xlab|ylab|guide-box|title|subtitle|caption|tag|axis-b-1)",
                                trim=FALSE)
#grid.newpage()
grid.draw(plot_filtered3)
