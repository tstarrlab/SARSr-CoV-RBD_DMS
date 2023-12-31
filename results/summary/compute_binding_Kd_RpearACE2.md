Compute per-barcode RpearACE2 binding affinity
================
Tyler Starr
01/02/2023

This notebook reads in per-barcode counts from `count_variants.ipynb`
for RpearACE2-binding Sort-seq experiments, computes functional scores
for RBD ACE2-binding affiniity, and does some basic QC on variant
binding functional scores.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","fitdistrplus","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$Titeseq_Kds_dir)){
  dir.create(file.path(config$Titeseq_Kds_dir))
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
    ##  [1] gridExtra_2.3       fitdistrplus_1.0-14 npsurv_0.4-0       
    ##  [4] lsei_1.2-0          survival_3.1-8      MASS_7.3-51.4      
    ##  [7] forcats_0.4.0       stringr_1.4.0       dplyr_0.8.3        
    ## [10] purrr_0.3.3         readr_1.3.1         tidyr_1.0.0        
    ## [13] tibble_3.0.2        ggplot2_3.3.0       tidyverse_1.3.0    
    ## [16] data.table_1.12.8   yaml_2.2.0          knitr_1.26         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0 xfun_0.11        lattice_0.20-38  splines_3.6.2   
    ##  [5] haven_2.2.0      colorspace_1.4-1 vctrs_0.3.1      generics_0.0.2  
    ##  [9] htmltools_0.4.0  rlang_0.4.7      pillar_1.4.5     glue_1.3.1      
    ## [13] withr_2.1.2      DBI_1.1.0        dbplyr_1.4.2     modelr_0.1.5    
    ## [17] readxl_1.3.1     lifecycle_0.2.0  munsell_0.5.0    gtable_0.3.0    
    ## [21] cellranger_1.1.0 rvest_0.3.5      evaluate_0.14    fansi_0.4.0     
    ## [25] broom_0.7.0      Rcpp_1.0.3       scales_1.1.0     backports_1.1.5 
    ## [29] jsonlite_1.6     fs_1.3.1         hms_0.5.2        digest_0.6.23   
    ## [33] stringi_1.4.3    grid_3.6.2       cli_2.0.0        tools_3.6.2     
    ## [37] magrittr_1.5     crayon_1.3.4     pkgconfig_2.0.3  Matrix_1.2-18   
    ## [41] ellipsis_0.3.0   xml2_1.3.3       reprex_0.3.0     lubridate_1.7.4 
    ## [45] assertthat_0.2.1 rmarkdown_2.0    httr_1.4.1       rstudioapi_0.10 
    ## [49] R6_2.4.1         compiler_3.6.2

## Setup

First, we will read in metadata on our sort samples, the table giving
number of reads of each barcode in each of the sort bins, and the
barcode-variant lookup tables, and merge these tables together.

``` r
#read dataframe with list of barcode runs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

#eliminate rows from barcode_runs that are not from a titeseq expt.
barcode_runs <- barcode_runs[barcode_runs$sample_type == "TiteSeq_RpearACE2",]

#read file giving count of each barcode in each sort partition
counts <- data.table(read.csv(file=config$variant_counts_file,stringsAsFactors=F))

#eliminate rows from counts that are not part of an titration bin sample. Also remove the ACE2 monomer data.
counts <- subset(counts, sample %in% barcode_runs[barcode_runs$sample_type=="TiteSeq_RpearACE2","sample"])

#read in barcode-variant lookup tables
dt_Wh1 <- data.table(read.csv(file=config$WH1_bc_variant_lookup,stringsAsFactors=F))
dt_Wh1[library=="pool1",library:="lib40"]
dt_Wh1[library=="pool2",library:="lib41"]
dt_RshSTT182 <- data.table(read.csv(file=config$codon_variant_table_file_RshSTT182,stringsAsFactors=F))
dt_PRD0038 <- data.table(read.csv(file=config$`codon_variant_table_file_PRD-0038`,stringsAsFactors=F))
dt_SARS1 <- data.table(read.csv(file=config$`codon_variant_table_file_SARS-CoV-1_2693`,stringsAsFactors=F))
dt_RsYN04 <- data.table(read.csv(file=config$codon_variant_table_file_RsYN04,stringsAsFactors=F))
dt_RmYN02 <- data.table(read.csv(file=config$codon_variant_table_file_RmYN02,stringsAsFactors=F))
dt_lib46 <- data.table(read.csv(file=config$codon_variant_table_file_lib46,stringsAsFactors=F))

#merge, eliminate barcodes duplicated within a library
dt <- rbind(dt_Wh1,dt_RshSTT182,dt_PRD0038,dt_SARS1,dt_RsYN04,dt_RmYN02,dt_lib46); setkey(dt,barcode,library)
duplicates <- dt[duplicated(dt,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplciates, and then remove
dt[,duplicate:=FALSE]
for(i in 1:nrow(duplicates)){
  dt[library==duplicates[i,library] & barcode==duplicates[i,barcode],duplicate:=TRUE]
}
dt <- dt[duplicate==FALSE,]; dt[,duplicate:=NULL]

dt <- merge(counts, dt, by=c("library","barcode")); rm(dt_Wh1);rm(dt_RshSTT182);rm(dt_PRD0038);rm(dt_SARS1);rm(dt_RsYN04);rm(dt_RmYN02); rm(dt_lib46);rm(counts); rm(duplicates)
```

Convert from Illumina read counts to estimates of the number of cells
that were sorted into a bin, and add some other useful information to
our data frame.

``` r
#for each bin, normalize the read counts to the observed ratio of cell recovery among bins
for(i in 1:nrow(barcode_runs)){
  lib <- as.character(barcode_runs$library[i])
  bin <- as.character(barcode_runs$sample[i])
  ratio <- sum(dt[library==lib & sample==bin,"count"])/barcode_runs$number_cells[i]
  if(ratio<1){ #if there are fewer reads from a FACS bin than cells sorted
    dt[library==lib & sample==bin, count.norm := as.numeric(count)] #don't normalize cell counts, make count.norm the same as count
    print(paste("reads < cells for",lib,bin,", un-normalized (ratio",ratio,")")) #print to console to inform of undersampled bins
  }else{
    dt[library==lib & sample==bin, count.norm := as.numeric(count/ratio)] #normalize read counts by the average read:cell ratio, report in new "count.norm" column
    print(paste("read:cell ratio for",lib,bin,"is",ratio))
  }
}
```

    ## [1] "read:cell ratio for lib40 TiteSeq_RpearACE2_01_bin1 is 1.27695606975348"
    ## [1] "read:cell ratio for lib40 TiteSeq_RpearACE2_01_bin2 is 1.48524836496688"
    ## [1] "reads < cells for lib40 TiteSeq_RpearACE2_01_bin3 , un-normalized (ratio 0.653892922599746 )"
    ## [1] "reads < cells for lib40 TiteSeq_RpearACE2_01_bin4 , un-normalized (ratio 0.0928662420382166 )"
    ## [1] "read:cell ratio for lib46 TiteSeq_RpearACE2_01_bin1 is 1.24565627650844"
    ## [1] "read:cell ratio for lib46 TiteSeq_RpearACE2_01_bin2 is 1.45684779419004"
    ## [1] "reads < cells for lib46 TiteSeq_RpearACE2_01_bin3 , un-normalized (ratio 0.164391043145822 )"
    ## [1] "reads < cells for lib46 TiteSeq_RpearACE2_01_bin4 , un-normalized (ratio 0.13677130044843 )"

``` r
#annotate each barcode as to whether it's a homolog variant, SARS-CoV-2 wildtype, synonymous muts only, stop, nonsynonymous, >1 nonsynonymous mutations
dt[,variant_class:=as.character(NA)]
dt[n_codon_substitutions==0, variant_class := "wildtype"]
dt[n_codon_substitutions > 0 & n_aa_substitutions==0, variant_class := "synonymous"]
dt[n_aa_substitutions>0 & grepl("*",aa_substitutions,fixed=T), variant_class := "stop"]
dt[n_aa_substitutions == 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := "1 nonsynonymous"]
dt[n_aa_substitutions > 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := ">1 nonsynonymous"]

#cast the data frame into wide format
dt <- dcast(dt, library + barcode + target + variant_class + aa_substitutions + n_aa_substitutions ~ sample, value.var="count.norm")

#add total count corresponding to count across the four bins for each barcode
dt[,RpearACE2_count := sum(TiteSeq_RpearACE2_01_bin1,TiteSeq_RpearACE2_01_bin2,TiteSeq_RpearACE2_01_bin3,TiteSeq_RpearACE2_01_bin4),by=c("library","barcode")]

#add indicator if count>1 in >1 bin
dt[,total_bins_w_count := sum(.(TiteSeq_RpearACE2_01_bin1,TiteSeq_RpearACE2_01_bin2,TiteSeq_RpearACE2_01_bin3,TiteSeq_RpearACE2_01_bin4)>0),by=c("library","barcode")]
```

## Calculating mean fluorescence

Next, for each barcode, calculate its mean fluorescence as an indicator
of RpearACE2 binding level. We will use a maximum likelihood approach to
determine the mean and standard deviation of fluorescence for a barcode,
given its distribution of cell counts across sort bins, and the known
fluorescence boundaries of those sort bins from the sorting log. The
package `fitdistcens` enables this ML estimation for these type of
*censored* observations, where we know we observed a cell within some
fluorescence interval but do not know the exact fluorescence value
attributed to that observation. The counts are multiplied by 20 so that
there is not a large rounding effect when they are rounded to integers.

Fluor boundaries:

    (-288, 172), (173, 1476), (1477, 12889), (12890, 262143)

``` r
#define function to calculate ML meanF
calc.MLmean <- function(b1,b2,b3,b4,min.b1,min.b2,min.b3,min.b4,max.b4,min.count=1){ #b1-4 gives observed cell counts in bins 1-4; remaining arguments give fluorescence boundaries of the respective bins; min.count gives minimum number of total observations needed across bins in order to calculate meanF (default 1)
  data <- data.frame(left=c(rep(min.b1,round(b1)),rep(min.b2,round(b2)),rep(min.b3,round(b3)),rep(min.b4,round(b4))),
                     right=c(rep(min.b2,round(b1)),rep(min.b3,round(b2)),rep(min.b4,round(b3)),rep(max.b4,round(b4)))) #define data input in format required for fitdistcens
  if(nrow(unique(data))>1 & nrow(data)>min.count){ #only fits if above user-specified min.count, and if the data satisfies the fitdistcens requirement that cells are observed in at least two of the censored partitions to enable ML estimation of identifiable parameters
    fit <- fitdistcens(data,"norm")
    return(list(as.numeric(summary(fit)$estimate["mean"]),as.numeric(summary(fit)$estimate["sd"])))
  } else {
    return(list(as.numeric(NA),as.numeric(NA)))
  }
}

#fit ML mean and sd fluorescence for each barcode, and calculate total cell count as the sum across the four bins. Multiply cell counts by a factor of 20 to minimize rounding errors since fitdistcens requires rounding to integer inputs
invisible(dt[,c("MFI_RpearACE2","sd_RpearACE2") := tryCatch(calc.MLmean(b1=TiteSeq_RpearACE2_01_bin1*20,b2=TiteSeq_RpearACE2_01_bin2*20,
                                                               b3=TiteSeq_RpearACE2_01_bin3*20,b4=TiteSeq_RpearACE2_01_bin4*20,
                                                               min.b1=log(20),min.b2=log(172.5),min.b3=log(1477.5),
                                                               min.b4=log(12889.5),max.b4=log(229000)),
                                                   error=function(e){return(list(as.numeric(NA),as.numeric(NA)))}),by=c("library","barcode")])

#save temp data file for downstream troubleshooting since the ML meanF took >1hr to calculate -- don't use these for final anlaysis though for reproducibility!
save(dt,file=paste(config$Titeseq_Kds_dir,"/dt.temp.RpearACE2.Rda",sep=""))
```

Next let’s look at the distributon of cell counts across the four bins
for each barcode, for those for which estimates were generated on the
bottom.

``` r
#histograms
par(mfrow=c(2,2))
hist(log10(dt[library=="lib40",RpearACE2_count]+0.1),xlab="cell count (log10, plus 0.1 pseudocount)",main="lib40, all bc",col="gray50")
hist(log10(dt[library=="lib46",RpearACE2_count]+0.1),xlab="cell count (log10, plus 0.1 pseudocount)",main="lib46, all bc",col="gray50")
hist(log10(dt[library=="lib40" & !is.na(MFI_RpearACE2),RpearACE2_count]+0.1),xlab="cell count (log10, plus 0.1 pseudocount)",main="lib40, determined",col="gray50")
hist(log10(dt[library=="lib46" & !is.na(MFI_RpearACE2),RpearACE2_count]+0.1),xlab="cell count (log10, plus 0.1 pseudocount)",main="lib46, determined",col="gray50")
```

<img src="compute_binding_Kd_RpearACE2_files/figure-gfm/cell_count_coverage-1.png" style="display: block; margin: auto;" />

Filter out expression measurements determined from \<5 estimated cells

``` r
min_count <- 5
dt[RpearACE2_count<min_count, c("MFI_RpearACE2","sd_RpearACE2","RpearACE2_count") := NA]
```

Distribution of MFI estimates for lib40 backgrounds, with wt/syn
barcodes in purple:

``` r
par(mfrow=c(3,2))
hist(dt[library=="lib40" & target=="Wuhan_Hu_1",MFI_RpearACE2],col="gray40",breaks=60,xlab="MFI, RpearACE2",main="Wuhan_Hu_1, lib40",xlim=c(4,8))
hist(dt[library=="lib40" & target=="Wuhan_Hu_1" & variant_class %in% (c("synonymous","wildtype")),MFI_RpearACE2],col="#92278F",add=T,breaks=60)

hist(dt[library=="lib40" & target=="RshSTT182",MFI_RpearACE2],col="gray40",breaks=20,xlab="MFI, RpearACE2",main="RshSTT182, lib40",xlim=c(4,8))
hist(dt[library=="lib40" & target=="RshSTT182" & variant_class %in% (c("synonymous","wildtype")),MFI_RpearACE2],col="#92278F",add=T,breaks=20)

hist(dt[library=="lib40" & target=="PRD-0038",MFI_RpearACE2],col="gray40",breaks=20,xlab="MFI, RpearACE2",main="PRD-0038, lib40",xlim=c(4,8))
hist(dt[library=="lib40" & target=="PRD-0038" & variant_class %in% (c("synonymous","wildtype")),MFI_RpearACE2],col="#92278F",add=T,breaks=20)

hist(dt[library=="lib40" & target=="SARS-CoV-1_2693",MFI_RpearACE2],col="gray40",breaks=20,xlab="MFI, RpearACE2",main="SARS-CoV-1_2693, lib40",xlim=c(4,8))
hist(dt[library=="lib40" & target=="SARS-CoV-1_2693" & variant_class %in% (c("synonymous","wildtype")),MFI_RpearACE2],col="#92278F",add=T,breaks=20)

hist(dt[library=="lib40" & target=="RsYN04",MFI_RpearACE2],col="gray40",breaks=20,xlab="MFI, RpearACE2",main="RsYN04, lib40",xlim=c(4,8))
hist(dt[library=="lib40" & target=="RsYN04" & variant_class %in% (c("synonymous","wildtype")),MFI_RpearACE2],col="#92278F",add=T,breaks=20)

hist(dt[library=="lib40" & target=="RmYN02",MFI_RpearACE2],col="gray40",breaks=20,xlab="MFI, RpearACE2",main="RmYN02, lib40",xlim=c(4,8))
hist(dt[library=="lib40" & target=="RmYN02" & variant_class %in% (c("synonymous","wildtype")),MFI_RpearACE2],col="#92278F",add=T,breaks=20)
```

<img src="compute_binding_Kd_RpearACE2_files/figure-gfm/MFI_distribution-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$Titeseq_Kds_dir,"/hist_MFI-per-barcode_lib40_RpearACE2.pdf",sep="")))
```

Distribution of MFI estimates for lib46 backgrounds:

``` r
par(mfrow=c(1,2))
hist(dt[library=="lib46" & variant_class=="wildtype",MFI_RpearACE2],col="gray40",breaks=60,xlab="MFI, RpearACE2",main="wildtypes, lib46",xlim=c(4,8))

hist(dt[library=="lib46" & variant_class=="1 nonsynonymous",MFI_RpearACE2],col="gray40",breaks=20,xlab="MFI, RpearACE2",main="mutants, lib46",xlim=c(4,8))
```

<img src="compute_binding_Kd_RpearACE2_files/figure-gfm/Kd_distribution_lib46-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$Titeseq_Kds_dir,"/hist_MFI-per-barcode_lib46_RpearACE2.pdf",sep="")))
```

Remove stops

``` r
#remove stop variants, which even if they eke through, either a) still have low counts and give poor fits as a result, or b) seem to be either dubious PacBio calls (lower variant_call_support) or have late stop codons which perhaps don't totally ablate funciton. Either way, the vast majority were purged before this step and we don't want to deal with the remaining ones!
dt[variant_class == "stop",c("MFI_RpearACE2","sd_RpearACE2","RpearACE2_count") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA))]
```

Let’s visualize the final binding measurements as violin plots for the
different wildtype targets. In next notebook, we’ll evaluate count depth
and possibly apply further filtering to remove low-count expression
estimates

``` r
p1 <- ggplot(dt[library=="lib40" & !is.na(MFI_RpearACE2),],aes(x=variant_class,y=MFI_RpearACE2))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("RpearACE2, MFI")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~target+library,nrow=6)

grid.arrange(p1,ncol=1)
```

<img src="compute_binding_Kd_RpearACE2_files/figure-gfm/binding_distribution_vioplot_lib40-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$Titeseq_Kds_dir,"/violin-plot_MFI-by-target_lib40_RpearACE2.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[library=="lib46" & !is.na(MFI_RpearACE2) & variant_class=="wildtype",],aes(x=target,y=MFI_RpearACE2))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("RpearACE2, MFI (wildtypes)")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))

grid.arrange(p1,ncol=1)
```

<img src="compute_binding_Kd_RpearACE2_files/figure-gfm/binding_distribution_vioplot_lib46-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$Titeseq_Kds_dir,"/violin-plot_MFI-by-target_lib46_RpearACE2.pdf",sep="")))
```

We have generated binding measurements for 80.92% of the barcodes in our
libraries.

## Data Output

Finally, let’s output our measurements for downstream analyses.

``` r
dt[,.(library,barcode,target,variant_class,aa_substitutions,n_aa_substitutions,
     RpearACE2_count,MFI_RpearACE2)] %>%
  mutate_if(is.numeric, round, digits=6) %>%
  write.csv(file=config$Titeseq_Kds_file_RpearACE2, row.names=F)
```
