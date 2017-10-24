Genes mapped per species
================

Purpose:
--------

To determine the number of genes with at least 1x coverage for each non-human primate when using the human genome as the reference sequence for alignment.

``` r
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(stringr)
library(genefilter)
library(data.table)
library(genefilter)
library(dplyr)
library(tibble)
```

Upload file of transcript lengths and create function to upload and prep gene counts.

``` r
data_dir <- "Compiled_by_Species"

##Downloading the transcript length for ENSEMBL IDs, GRCh38 - the download can take a while, so ultimtely saved as a csv for easier access later.
##library(biomaRt)
##library(dplyr)
##txdb <- makeTxDbFromBiomart(dataset="hsapiens_gene_ensembl")
##exons.list.per.gene <- exonsBy(txdb, by = "gene")
##gene_sizes <- sum(width(reduce(exons.list.per.gene)))
##write.csv(gene_sizes, file = "transcript length.csv")

gene_sizes <- read.table("transcript length.csv", sep = ",", header = TRUE)
colnames(gene_sizes) <- c("rowname", "length")
gene_sizes_df <- data.frame(gene_sizes) 
gene_sizes_df[sapply(gene_sizes_df, is.factor)] <- lapply(gene_sizes_df[sapply(gene_sizes_df, is.factor)], as.character)

species_count <- function(species_dir) {
  ##File list
  sampleFiles <- basename(Sys.glob(file.path(species_dir, "*Gene_counts.tabular")))
  sampleFiles
  ##Creating a sample table for further analysis
  sampleNames <- sampleFiles %>% str_replace("_*htseq-count_Gene_counts.tabular", "") %>% str_replace("^[0-9]*_", "")
  sampleNames
  sampleDonor <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[1]))
  sampleTreatment <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[2]))
  sampleReplicate <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[3]))
  sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, 
                            donor = sampleDonor, treatment = sampleTreatment, 
                            replicate = sampleReplicate)
  sampleTable
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = species_dir, design = ~ donor + treatment)
  contrast <- c("treatment", "treated", "mock")
  dds <- estimateSizeFactors(dds)
  dds@colData
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  counts <- data.frame(counts(dds, normalized = TRUE), stringsAsFactors = FALSE)
  counts_df <- as_tibble(rownames_to_column(counts))
}
```

Now make a new function to assess the count coverage, taking into account each donor so that if at least one donor of a species has 1x coverage, it will be included. Then evaluate each species using this function and export the final number of genes per species with at least 1x coverage for at least one donor as a csv.

``` r
count_coverage <- function(w, donor1, donor2, donor3) {
joined <- left_join(w, gene_sizes_df, by = "rowname") %>%
  mutate(x = rowMeans(select(., contains(donor1)))) %>%
  mutate(y = rowMeans(select(., contains(donor2)))) %>%
  mutate(z = rowMeans(select(., contains(donor3)))) %>%
  mutate(RPK_1 = (((x)/(length))*1000), 
         RPK_2 = (((y)/(length))*1000), 
         RPK_3 = (((z)/(length))*1000)) %>%
  filter(RPK_1 >= 13.3 | RPK_2 >= 13.3 | RPK_3 >= 13.3)
joined
}

bonobo_list <- species_count("Compiled_by_Species/Bonobo/data") %>%
  count_coverage("PR00248", "PR111", "PR235")
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
human_list <- species_count("Compiled_by_Species/Human/data") %>%
  count_coverage("AF", "NHDF", "SR")
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
chimp_list <- species_count("Compiled_by_Species/Chimpanzee/data") %>%
  count_coverage("S003611", "S003649", "S004933")
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
gorilla_list <- species_count("Compiled_by_Species/Gorilla/data") %>%
  count_coverage("PR00107", "PR0230", "PR00573")
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
baboon_list <- species_count("Compiled_by_Species/Olive baboon/data") %>%
  count_coverage("PR00033", "PR00036", "PR00039")
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
orangutan_list <- species_count("Compiled_by_Species/Orangutan/data") %>%
  count_coverage("AG06105", "PR00054", "PR01109")
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
ptmac_list <- species_count("Compiled_by_Species/Pigtailed macaque/data") %>%
  count_coverage("AG08490", "AG07923", "PR0058")
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
rhesus_list <- species_count("Compiled_by_Species/Rhesus macaque/data") %>%
  count_coverage("AG08305", "AG08308", "AG08312")
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
sqmonk_list <- species_count("Compiled_by_Species/Squirrel monkey/data") %>%
  count_coverage("AG05311", "SQMA", "SQMB")
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
NHP_list <- list(sqmonk_list, baboon_list, chimp_list, gorilla_list, rhesus_list, 
                 ptmac_list, human_list, orangutan_list, bonobo_list)
names(NHP_list) <- c("Squirrel monkey", "Baboon", "Chimpanzee", "Gorilla", 
                     "Rhesus", "Pigtailed macaque", "Human", "Orangutan", "Bonobo")
l <- lapply(NHP_list, nrow)
l
```

    ## $`Squirrel monkey`
    ## [1] 8329
    ## 
    ## $Baboon
    ## [1] 11199
    ## 
    ## $Chimpanzee
    ## [1] 11860
    ## 
    ## $Gorilla
    ## [1] 11471
    ## 
    ## $Rhesus
    ## [1] 10533
    ## 
    ## $`Pigtailed macaque`
    ## [1] 10609
    ## 
    ## $Human
    ## [1] 12041
    ## 
    ## $Orangutan
    ## [1] 11436
    ## 
    ## $Bonobo
    ## [1] 11780

``` r
write.table(l, file = "NHP 1x coverage.txt", row.names = FALSE, sep="\t",quote=FALSE)
```

Session info

``` r
sessionInfo()
```

    ## R version 3.3.3 (2017-03-06)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2               tibble_1.3.3              
    ##  [3] dplyr_0.7.3                data.table_1.10.0         
    ##  [5] genefilter_1.54.2          stringr_1.2.0             
    ##  [7] RColorBrewer_1.1-2         gplots_3.0.1              
    ##  [9] DESeq2_1.12.4              SummarizedExperiment_1.2.3
    ## [11] Biobase_2.32.0             GenomicRanges_1.24.3      
    ## [13] GenomeInfoDb_1.8.7         IRanges_2.6.1             
    ## [15] S4Vectors_0.10.3           BiocGenerics_0.18.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.10         locfit_1.5-9.1       lattice_0.20-35     
    ##  [4] gtools_3.5.0         assertthat_0.2.0     rprojroot_1.2       
    ##  [7] digest_0.6.12        R6_2.2.0             plyr_1.8.4          
    ## [10] backports_1.0.5      acepack_1.4.1        RSQLite_1.1-2       
    ## [13] evaluate_0.10        ggplot2_2.2.1        zlibbioc_1.18.0     
    ## [16] rlang_0.1.2          lazyeval_0.2.0       annotate_1.50.1     
    ## [19] gdata_2.17.0         rpart_4.1-10         Matrix_1.2-8        
    ## [22] checkmate_1.8.2      rmarkdown_1.4        splines_3.3.3       
    ## [25] BiocParallel_1.6.6   geneplotter_1.50.0   foreign_0.8-67      
    ## [28] htmlwidgets_0.9      RCurl_1.95-4.8       munsell_0.4.3       
    ## [31] pkgconfig_2.0.1      base64enc_0.1-3      htmltools_0.3.5     
    ## [34] nnet_7.3-12          gridExtra_2.2.1      htmlTable_1.9       
    ## [37] Hmisc_4.0-2          XML_3.98-1.9         bitops_1.0-6        
    ## [40] grid_3.3.3           xtable_1.8-2         gtable_0.2.0        
    ## [43] DBI_0.6-1            magrittr_1.5         scales_0.4.1        
    ## [46] KernSmooth_2.23-15   stringi_1.1.5        XVector_0.12.1      
    ## [49] latticeExtra_0.6-28  Formula_1.2-1        tools_3.3.3         
    ## [52] glue_1.1.1           survival_2.41-3      yaml_2.1.14         
    ## [55] AnnotationDbi_1.34.4 colorspace_1.3-2     cluster_2.0.6       
    ## [58] caTools_1.17.1       memoise_1.0.0        bindr_0.1           
    ## [61] knitr_1.16
