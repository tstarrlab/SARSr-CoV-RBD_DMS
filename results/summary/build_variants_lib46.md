build final barcode-variant lookup table for lib46 (SARSr wildtypes and
mini-SSM pool v2)
================
Tyler Starr
12/28/2022

-   <a href="#setup" id="toc-setup">Setup</a>
-   <a href="#data-input" id="toc-data-input">Data input</a>
-   <a href="#process-pacbio-sequencing-variant-parsing-library-coverage"
    id="toc-process-pacbio-sequencing-variant-parsing-library-coverage">Process
    PacBio sequencing: variant parsing, library coverage</a>
-   <a href="#save-codon-variant-table"
    id="toc-save-codon-variant-table">Save codon-variant table</a>

## Setup

This notebook reads in the nt-variant-lookup table containing both the
wildtypes pools of SARSr-CoVs as well as the SSM mutants and parses
nt-variant lookup table

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","seqinr")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("config.yaml")

#read in file giving concordance between RBD numbering and SARS-CoV-2 Spike numbering
RBD_sites <- read.csv(file=config$RBD_annotation_file,stringsAsFactors=F)

#make output directory
if(!file.exists(config$variants_dir)){
 dir.create(file.path(config$variants_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_haswellp-r0.3.7.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] seqinr_3.6-1      gridExtra_2.3     forcats_0.4.0     stringr_1.4.0    
    ##  [5] dplyr_0.8.3       purrr_0.3.3       readr_1.3.1       tidyr_1.0.0      
    ##  [9] tibble_3.0.2      ggplot2_3.3.0     tidyverse_1.3.0   data.table_1.12.8
    ## [13] yaml_2.2.0        knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0 xfun_0.11        haven_2.2.0      colorspace_1.4-1
    ##  [5] vctrs_0.3.1      generics_0.0.2   htmltools_0.4.0  rlang_0.4.7     
    ##  [9] pillar_1.4.5     glue_1.3.1       withr_2.1.2      DBI_1.1.0       
    ## [13] dbplyr_1.4.2     modelr_0.1.5     readxl_1.3.1     lifecycle_0.2.0 
    ## [17] munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0 rvest_0.3.5     
    ## [21] evaluate_0.14    fansi_0.4.0      broom_0.7.0      Rcpp_1.0.3      
    ## [25] scales_1.1.0     backports_1.1.5  jsonlite_1.6     fs_1.3.1        
    ## [29] hms_0.5.2        digest_0.6.23    stringi_1.4.3    ade4_1.7-13     
    ## [33] grid_3.6.2       cli_2.0.0        tools_3.6.2      magrittr_1.5    
    ## [37] crayon_1.3.4     pkgconfig_2.0.3  MASS_7.3-51.4    ellipsis_0.3.0  
    ## [41] xml2_1.3.3       reprex_0.3.0     lubridate_1.7.4  assertthat_0.2.1
    ## [45] rmarkdown_2.0    httr_1.4.1       rstudioapi_0.10  R6_2.4.1        
    ## [49] compiler_3.6.2

## Data input

Read in tables from processed sequencing data. The `dt` table gives the
barcode, RBD background, and any nucleotide mutations found for each
variant in each background. Remove any clashing barcodes.

## Process PacBio sequencing: variant parsing, library coverage

First, let’s assess how well we got all of the desired site-saturation
mutagenesis positions covered in our variants coming off of the PacBio
sequencing. We need to parse nucleotide mutations to their amino acid
mutations, and flag as invalid double mutants or mutants at unintended
positions/backgrounds. To do this, we create a function that checks
whether a background should be mutated, and if so, indexes the positions
that should be mutated as targeted positiosn differ in index across
backgrounds due to variations in length over evolutionary time. We then
output the variant_class as wildtype, mutant, invalid (unintended
mutation), indel (also unintended), synonymous, or stop mutant.

``` r
#load a table giving the indexing of nucleotide numbers for sites targeted in each mutated background
index <- read.csv(file=config$mutant_indexing_file,stringsAsFactors=F)

#set empty columns for filling with mutant information
dt[,variant_class:=as.character(NA)];dt[,wildtype:=as.character(NA)];dt[,position:=as.numeric(NA)];dt[,mutant:=as.character(NA)]

#set a function that returns variant_class, wildtype aa, SARS2 indexed position, and mutant aa for each single nt barcode
parse_aamut <- function(nt_substitutions,background){
  variant_class_return <- as.character(NA)
  wildtype_return <- as.character(NA)
  position_return <- as.numeric(NA)
  mutant_return <- as.character(NA)
  subs <- strsplit(as.character(nt_substitutions),split=" ")[[1]]
  if(length(subs)==0){ #if wildtype
    variant_class_return <- "wildtype"
  }else{ #if mutations
    if(background %in% config$mutated_targets){ #if background with mutation intended
      positions <- vector(length=length(subs), mode="numeric")
      for(k in 1:length(subs)){
        positions[k] <- as.numeric(paste(strsplit(subs[k],split="")[[1]][2:(length(strsplit(subs[k],split="")[[1]])-1)],collapse=""))
      }
      aa_pos <- unique(ceiling(positions/3))
      if(length(aa_pos)>1){ #if multiple codon mutations
        variant_class_return <- "invalid"
      }
      if(length(aa_pos)==1){ #if single codon mutation, assign site in SARS-CoV-2 numbering, wildtype AA, and mutant AA
        if(aa_pos==index[index$target==background,"index_455"]){ #if mutant at site 455
          variant_class_return <- "mutant"
          wildtype_return <- index[index$target==background,"aa_455"]
          position_return <- 455
          codon <- strsplit(index[index$target==background,"codon_455"],split="")[[1]]
          for(j in 1:length(positions)){ #iterate through positions in codon and mutate nt in codon string if needed
            if(positions[j]==index[index$target==background,"nt_455"]){ #if first position of codon is mutated
              codon[1] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_455"]+1){ #if second position of codon is mutated
              codon[2] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_455"]+2){ #if third position of codon is mutated
              codon[3] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
          }
          mutant_return <- translate(tolower(codon))
        }else if(aa_pos==index[index$target==background,"index_486"]){ #if mutant at site 486
          variant_class_return <- "mutant"
          wildtype_return <- index[index$target==background,"aa_486"]
          position_return <- 486
          codon <- strsplit(index[index$target==background,"codon_486"],split="")[[1]]
          for(j in 1:length(positions)){ #iterate through positions in codon and mutate nt in codon string if needed
            if(positions[j]==index[index$target==background,"nt_486"]){ #if first position of codon is mutated
              codon[1] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_486"]+1){ #if second position of codon is mutated
              codon[2] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_486"]+2){ #if third position of codon is mutated
              codon[3] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
          }
          mutant_return <- translate(tolower(codon))
        }else if(aa_pos==index[index$target==background,"index_493"]){ #if mutant at site 493
          variant_class_return <- "mutant"
          wildtype_return <- index[index$target==background,"aa_493"]
          position_return <- 493
          codon <- strsplit(index[index$target==background,"codon_493"],split="")[[1]]
          for(j in 1:length(positions)){ #iterate through positions in codon and mutate nt in codon string if needed
            if(positions[j]==index[index$target==background,"nt_493"]){ #if first position of codon is mutated
              codon[1] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_493"]+1){ #if second position of codon is mutated
              codon[2] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_493"]+2){ #if third position of codon is mutated
              codon[3] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
          }
          mutant_return <- translate(tolower(codon))
        }else if(aa_pos==index[index$target==background,"index_494"]){ #if mutant at site 494
          variant_class_return <- "mutant"
          wildtype_return <- index[index$target==background,"aa_494"]
          position_return <- 494
          codon <- strsplit(index[index$target==background,"codon_494"],split="")[[1]]
          for(j in 1:length(positions)){ #iterate through positions in codon and mutate nt in codon string if needed
            if(positions[j]==index[index$target==background,"nt_494"]){ #if first position of codon is mutated
              codon[1] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_494"]+1){ #if second position of codon is mutated
              codon[2] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_494"]+2){ #if third position of codon is mutated
              codon[3] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
          }
          mutant_return <- translate(tolower(codon))
        }else if(aa_pos==index[index$target==background,"index_498"]){ #if mutant at site 498
          variant_class_return <- "mutant"
          wildtype_return <- index[index$target==background,"aa_498"]
          position_return <- 498
          codon <- strsplit(index[index$target==background,"codon_498"],split="")[[1]]
          for(j in 1:length(positions)){ #iterate through positions in codon and mutate nt in codon string if needed
            if(positions[j]==index[index$target==background,"nt_498"]){ #if first position of codon is mutated
              codon[1] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_498"]+1){ #if second position of codon is mutated
              codon[2] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_498"]+2){ #if third position of codon is mutated
              codon[3] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
          }
          mutant_return <- translate(tolower(codon))
        }else if(aa_pos==index[index$target==background,"index_501"]){ #if mutant at site 501
          variant_class_return <- "mutant"
          wildtype_return <- index[index$target==background,"aa_501"]
          position_return <- 501
          codon <- strsplit(index[index$target==background,"codon_501"],split="")[[1]]
          for(j in 1:length(positions)){ #iterate through positions in codon and mutate nt in codon string if needed
            if(positions[j]==index[index$target==background,"nt_501"]){ #if first position of codon is mutated
              codon[1] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_501"]+1){ #if second position of codon is mutated
              codon[2] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
            if(positions[j]==index[index$target==background,"nt_501"]+2){ #if third position of codon is mutated
              codon[3] <- strsplit(subs[j],split="")[[1]][length(strsplit(subs[j],split="")[[1]])]
            }
          }
          mutant_return <- translate(tolower(codon))
        }else{ #if mutant at other (unintended) position
          variant_class_return <- "invalid"
        }
      }
    }else{ #if background with mutation unintended
      variant_class_return <- "invalid"
    }
  }
  return(list(variant_class_return, wildtype_return, position_return, mutant_return))
}

dt[,c("variant_class", "wildtype", "position", "mutant") := parse_aamut(nt_substitutions=substitutions,background=target),by=c("library","barcode")]

dt[variant_class=="mutant" & mutant==wildtype,variant_class:="synonymous"]
dt[mutant=="*",variant_class:="stop"]
dt[number_of_indels>0,variant_class:="indel"]

#side thing -- for position 455, we introduced a secondary pool of variants ("pool 6") replacing missed mutations at positions 455 in GD-Pangolin and RaTG13
#the NNS mutagenesis should reintroduce wildtype for comparison within the pool-6-only experiments (because we need to normalize mutants to wildtype within the experiment), but it appears the wildtype that is reintroduced is a synonymous mut.
#(This is not a problem for BtKY72, becasue apparently between the 6 positiosn in pool 6 in this background, there are sufficient true-wildtype representatives, even if we are discarding some additional synonymous variants)
#because we filter out synonymous mutations below for simplicity, let's recode specifically 455 synonymous muts in these two backgrounds so the proper wt comparitors in the pool-6-only experiments make it through to the phenotyping script where they are needed.
dt[variant_class=="synonymous" & target %in% c("RaTG13","GD-Pangolin") & position==455,variant_class:="wildtype"]
```

For the mutated backgrounds, let’s look at the fraction of variants in
each class. In this PacBio sequencing, we have two backgrounds currently
present for SARS-CoV-1_Urbani and SARS-CoV-2. For SARS-CoV-1, we created
our own “in-house” SSM libraries in the 2693 background coding sequence
due to delay in shipment of the product from Genscript. The Genscript
product did arrive in time for PacBio sequencing, so we included it for
barcode attribution in case there were problems with our “in-house”
assembly, but because the “in-house” assembly looks good, we will move
forward just that background (\_2693, but we will rename it to the full
\_Urbani_HP03L name further below). On the SARS-CoV-2 side, the tube
from Genscript corresponding to the SSM library for position S494 had
zero volume, and it does not appear that DNA was deposited in any of the
other tubes inadvertently since I do not see these mutated positions
present in this background. In expectation that this was the case, we
performed SSM at position S494 in our 2649 background SARS-CoV-2
sequence, so further below we will pool the mutants across these two
backgrounds to get our final set of mini-mutant scanning mutants for
this background.

We can see there is perhaps a slightly higher proportion of invalid
mutants in the GD-Pangolin and RaTG13 backgrounds (which we see below is
attributable to Genscript targeting the wrong mutation with their SSM!).
Can see our SARS-CoV-2_2649 assembly has virtually only mutants,
consistent with the fact that we introduced mutations only at a single
position here and didn’t really pool in extra wildtype of this sequence
(since we have the Genscript SARS-CoV-2 wildtype).

<img src="build_variants_lib46_files/figure-gfm/variant_fractions_mutated-1.png" style="display: block; margin: auto;" />

For backgrounds not targeted with mutations, do any stand out as having
many “invalid” (i.e. mutated) and indel variants? Most of these were
cloned in bulk directly from Twist-synthesized oligos – so it’s
encouraging to see most assemblies are correct.

<img src="build_variants_lib46_files/figure-gfm/variant_fractions_unmutated-1.png" style="display: block; margin: auto;" />

Let’s collapse the SARS-CoV-1_Urbani and SARS-CoV-2 information, as
described above. For SARS-CoV-1_Urbani, we will remove the
“SARS-CoV-1_Urbani_HP03L” mutant variants – these were pooled into the
PacBio but were not part of the titrations. (The wildtype of this
background is, though, as it was pooled separately from the delayed
mutants shipment.) For SARS-CoV-2, we can simply pool the \_2649 mutants
into the remaining SARS-CoV-2_WH1 mutants from Genscript.

``` r
dt[target=="SARS-CoV-2_2649",target:="SARS-CoV-2_WH1"]
dt[target=="SARS-CoV-1_Urbani_HP03L" & variant_class=="wildtype",target:="SARS-CoV-1_2693"]
dt <- dt[target!="SARS-CoV-1_Urbani_HP03L"]
```

Let’s look at coverage of mutants at each position across each
background.

    ## Warning: Removed 84 rows containing missing values (position_stack).

<img src="build_variants_lib46_files/figure-gfm/variant_fractions-1.png" style="display: block; margin: auto;" />

Check out if we are missing any mutants at our intended positions. Apart
from the missing 455 mutations in RaTG13 and GD-Pangolin due to
Genscript targeting L452 mistakenly, we are only missing one mutation!
N501K in AncSARS1a.

``` r
kable(aa_coverage[aa_coverage$count==0 & !is.na(aa_coverage$count),])
```

|     | target        | site | mutant | count |
|-----|:--------------|-----:|:-------|------:|
| 746 | AncSARS1a_MAP |  501 | K      |     0 |

Check out some other coverage statistics: average (median) \# barcodes
for a mutant across the two libraries is 27, minimum is 2. The
distribution of bc number per mutant is shown in the histogram below.

<img src="build_variants_lib46_files/figure-gfm/coverage_stats-1.png" style="display: block; margin: auto;" />

## Save codon-variant table

Save the codon-variant lookup table. This table should have target,
library, barcode, variant_call_support, codon_substitutions,
aa_substitutions, n_codon_substitutions, n_aa_substitutions in that
order. We make these formatting changes below

``` r
dt[,codon_substitutions:=NA]
dt[variant_class=="mutant", aa_substitutions:=paste(wildtype,position,mutant,sep=""),by=c("library","barcode")]
dt[,n_codon_substitutions:=as.numeric(NA)]
dt[variant_class %in% c("wildtype"),n_codon_substitutions:=0]
dt[,n_aa_substitutions:=as.numeric(NA)]
dt[variant_class %in% c("synonymous","wildtype"),n_aa_substitutions:=0]
dt[variant_class=="mutant", n_aa_substitutions:=1]
```

``` r
write.csv(dt[variant_class != "invalid",.(target, library, barcode, variant_call_support, codon_substitutions, aa_substitutions, n_codon_substitutions, n_aa_substitutions)], file=config$codon_variant_table_file_lib46, quote=F, row.names=F)
```
