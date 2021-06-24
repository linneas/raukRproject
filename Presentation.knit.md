---
title: "Genomic data in R"
author: "Linnéa Smeds"
date: "23/06/2021"
output: ioslides_presentation
logo: "images/raukR_logo.svg"

---



# Part 1  | Data and Aims

## Data

### Whole genome SNP data from 210 wolves:
    - *Scandinavian=Recently founded*
    - *Finnish/Russian=Source population*
   
### Detailed information about Scandinavian population: 
  - *sample year*
  - *founding individuals*
  - *generation to founder*
  
> (in this project, I only use ~5000 deleterious SNPs)
  
## Aim

To rewrite bash&python code that extract frequencies of homozygous derived/heterozygous sites for each individual from the vcf file, and plot population differences along with all info given for the focal population

### My personal goal

Learn how to use tidyverse (especially pipes) and ggplot

# Part 2 | Results

## Packages used

```r
require(vcfR)
require(tictoc)
require(viridis)
require(tidyverse)
```

## Read in the vcf file and convert it to tidy 

```r
vcf_file="data/Wolf.deleterious.missense.polarized.vcf"
vcf <- vcfR::read.vcfR(vcf_file)
tidy_vcf <- vcfR::vcfR2tidy(vcf, 
                      info_fields=c("AA"), 
                      format_fields=c("GT"), 
                      dot_is_NA=TRUE)
```
The population, founder and generation info from plain textfiles are saved as tibbles directly

```r
# example
pop_tib<- "data/populations.txt" %>% 
      read.table() %>% 
      as_tibble() %>% 
      rename(Indiv=V1, Pop=V2)
```



## Inspecting the tables

```r
str(tidy_vcf)
```

```
## List of 3
##  $ fix : tibble[,9] [4,922 × 9] (S3: tbl_df/tbl/data.frame)
##   ..$ ChromKey: int [1:4922] 1 1 1 1 1 1 1 1 1 1 ...
##   ..$ CHROM   : chr [1:4922] "1" "1" "1" "1" ...
##   ..$ POS     : int [1:4922] 2997440 4310746 4311229 4838054 5191478 8427892 8535137 8590825 8590856 9112132 ...
##   ..$ ID      : chr [1:4922] NA NA NA NA ...
##   ..$ REF     : chr [1:4922] "G" "C" "T" "G" ...
##   ..$ ALT     : chr [1:4922] "C" "T" "G" "T" ...
##   ..$ QUAL    : num [1:4922] 2998 2732 2720 38729 14608 ...
##   ..$ FILTER  : chr [1:4922] "PASS" "PASS" "PASS" "PASS" ...
##   ..$ AA      : chr [1:4922] "G" "C" "T" "G" ...
##  $ gt  : tibble[,5] [1,033,620 × 5] (S3: tbl_df/tbl/data.frame)
##   ..$ ChromKey     : int [1:1033620] 1 1 1 1 1 1 1 1 1 1 ...
##   ..$ POS          : int [1:1033620] 2997440 4310746 4311229 4838054 5191478 8427892 8535137 8590825 8590856 9112132 ...
##   ..$ Indiv        : chr [1:1033620] "100-G47-11" "100-G47-11" "100-G47-11" "100-G47-11" ...
##   ..$ gt_GT        : chr [1:1033620] "0/0" "0/0" "0/0" "1/1" ...
##   ..$ gt_GT_alleles: chr [1:1033620] "G/G" "C/C" "T/T" "T/T" ...
##  $ meta: tibble[,5] [2 × 5] (S3: tbl_df/tbl/data.frame)
##   ..$ Tag        : chr [1:2] "INFO" "FORMAT"
##   ..$ ID         : chr [1:2] "AA" "gt_GT"
##   ..$ Number     : chr [1:2] "1" "1"
##   ..$ Type       : chr [1:2] "String" "String"
##   ..$ Description: chr [1:2] "Ancestral Allele" "Genotype"
```

## Inspecting the tables (cont.)

```r
head(tidy_vcf$gt, 3)
```

```
## # A tibble: 3 x 5
##   ChromKey     POS Indiv      gt_GT gt_GT_alleles
##      <int>   <int> <chr>      <chr> <chr>        
## 1        1 2997440 100-G47-11 0/0   G/G          
## 2        1 4310746 100-G47-11 0/0   C/C          
## 3        1 4311229 100-G47-11 0/0   T/T
```

```r
head(tidy_vcf$fix, 3)
```

```
## # A tibble: 3 x 9
##   ChromKey CHROM     POS ID    REF   ALT    QUAL FILTER AA   
##      <int> <chr>   <int> <chr> <chr> <chr> <dbl> <chr>  <chr>
## 1        1 1     2997440 <NA>  G     C     2998. PASS   G    
## 2        1 1     4310746 <NA>  C     T     2732. PASS   C    
## 3        1 1     4311229 <NA>  T     G     2720. PASS   T
```



## Extract relevant info - Trial #1

```r
tic()
okpos<- tidy_vcf$fix %>%filter(REF==AA)
invertpos<-tidy_vcf$fix %>%filter(ALT==AA)
homref_aa_ref<-tidy_vcf$gt %>% inner_join(okpos) %>% filter(gt_GT=="0/0") %>% group_by(Indiv) %>%summarize(num_homref_aa=n())
homalt_aa_ref<-tidy_vcf$gt %>% inner_join(okpos) %>% filter(gt_GT=="1/1") %>% group_by(Indiv) %>%summarize(num_homalt_der=n())
homref_aa_alt<-tidy_vcf$gt %>% inner_join(invertpos) %>% filter(gt_GT=="0/0") %>% group_by(Indiv) %>%summarize(num_homref_der=n())
homalt_aa_alt<-tidy_vcf$gt %>% inner_join(invertpos) %>% filter(gt_GT=="1/1") %>% group_by(Indiv) %>%summarize(num_homalt_aa=n())
aa_tib<- homref_aa_ref %>% inner_join(homalt_aa_alt) %>% mutate(num_aa=num_homref_aa+num_homalt_aa) %>% select(Indiv, num_aa)
der_tib<-homalt_aa_ref %>% inner_join(homref_aa_alt) %>% mutate(num_der=num_homref_der+num_homalt_der) %>% select(Indiv, num_der)
het<-tidy_vcf$gt %>% filter(gt_GT=="0/1") %>% group_by(Indiv) %>%summarize(num_het=n())
tidy_gt<-aa_tib %>% inner_join(het) %>% inner_join(der_tib) 
toc()
```

```
## 3.712 sec elapsed
```

## Extract relevant info - Trial #2

```r
tic()
refAA_tib<-tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>% inner_join(okpos) %>% group_by(Indiv, gt_GT) %>% summarize(count=n()) %>% ungroup() %>% mutate(gt = case_when(gt_GT=='0/0' ~ 'refAA_homanc', gt_GT=='0/1' ~'refAA_het', gt_GT=='1/1' ~'refAA_homder')) %>% select(-gt_GT) %>% pivot_wider(names_from=gt, values_from=count)
altAA_tib<-tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>% inner_join(invertpos) %>% group_by(Indiv, gt_GT) %>% summarize(count=n()) %>% ungroup() %>% mutate(gt = case_when(gt_GT=='0/0' ~ 'altAA_homder', gt_GT=='0/1' ~'altAA_het', gt_GT=='1/1' ~'altAA_homanc')) %>% select(-gt_GT) %>% pivot_wider(names_from=gt, values_from=count)
tidy_gt2 <- refAA_tib %>% inner_join(altAA_tib) %>% transmute(Indiv, num_aa=altAA_homanc+refAA_homanc, num_het=refAA_het+altAA_het, num_der=refAA_homder+altAA_homder) 
toc()
```

```
## 1.777 sec elapsed
```

## Extract relevant info - Trial #3

```r
tic()
tidy_gt3<-tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>% 
  inner_join(tidy_vcf$fix)  %>% 
  mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                          (ALT==AA & gt_GT=='1/1') ~'0/0',
                          TRUE ~ gt_GT)) %>%
  group_by(Indiv, new_gt) %>% summarize(count=n()) %>% 
  ungroup() %>% mutate_if(is.character, str_replace_all, 
                          pattern = "/", replacement = "") %>%
  pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count)
toc()
```

```
## 1.089 sec elapsed
```

## Finally!

```r
head(tidy_gt3, 3)
```

```
## # A tibble: 3 x 4
##   Indiv        gt00  gt01  gt11
##   <chr>       <int> <int> <int>
## 1 100-G47-11   4062   646   208
## 2 102-G34-10   4162   276   288
## 3 103-D-11-17  4280   398   225
```
Merge with other info tibbles

```r
final_gt<- tidy_gt3 %>% mutate(sum=gt00+gt01+gt11) %>% 
  inner_join(pop_tib) %>% 
  mutate(hetfrq=gt01/sum, homderfrq=gt11/sum, 
         deralfreq=(gt11*2+gt01)/(sum*2),
         Grp=case_when(Pop=='Finland' | Pop=='Russia' ~ Pop, TRUE ~ 'Scandinavia'))  %>% 
  left_join(fou_tib) %>% left_join(gen_tib)
```

## Now to the plotting

```r
ggplot(final_gt,aes(Pop,hetfrq))+
  geom_boxplot(aes(fill=Pop))+
    geom_point() +
    labs(x=NULL,y="proportion heterozygous sites")+
    theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
```

<img src="Presentation_files/figure-html/unnamed-chunk-13-1.png" width="720" />

