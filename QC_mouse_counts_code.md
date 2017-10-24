QC of mouse counts
================

Purpose:
--------

To see how all the mouse read counts cluster.

``` r
library(DESeq2)
library(gplots)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(genefilter)
library(data.table)
```

Upload gene count files and generate a table of the samples broken up by donor ID, replicate designation, and treatment (polyI:C treated or mock).

``` r
data_dir <- "Mouse counts"
output_dir <- "QC_mouse_counts_output"

## File list
sampleFiles <- basename(Sys.glob(file.path(data_dir, "*Gene_counts.tabular")))
sampleFiles
```

    ##  [1] "11_C57A_mock_A__htseq-count_Gene_counts.tabular"   
    ##  [2] "11_C57A_mock_B__htseq-count_Gene_counts.tabular"   
    ##  [3] "11_C57A_mock_C__htseq-count_Gene_counts.tabular"   
    ##  [4] "11_C57A_treated_A__htseq-count_Gene_counts.tabular"
    ##  [5] "11_C57A_treated_B__htseq-count_Gene_counts.tabular"
    ##  [6] "11_C57A_treated_C__htseq-count_Gene_counts.tabular"
    ##  [7] "11_C57B_mock_A__htseq-count_Gene_counts.tabular"   
    ##  [8] "11_C57B_mock_B__htseq-count_Gene_counts.tabular"   
    ##  [9] "11_C57B_mock_C__htseq-count_Gene_counts.tabular"   
    ## [10] "11_C57B_treated_A__htseq-count_Gene_counts.tabular"
    ## [11] "11_C57B_treated_B__htseq-count_Gene_counts.tabular"
    ## [12] "11_C57B_treated_C__htseq-count_Gene_counts.tabular"
    ## [13] "11_C57C_mock_A__htseq-count_Gene_counts.tabular"   
    ## [14] "11_C57C_mock_B__htseq-count_Gene_counts.tabular"   
    ## [15] "11_C57C_mock_C__htseq-count_Gene_counts.tabular"   
    ## [16] "11_C57C_treated_A__htseq-count_Gene_counts.tabular"
    ## [17] "11_C57C_treated_B__htseq-count_Gene_counts.tabular"
    ## [18] "11_C57C_treated_C__htseq-count_Gene_counts.tabular"

``` r
## Create table of samples with donor, treatment, and replicate information
sampleNames <- sampleFiles %>% str_replace("_*htseq-count_Gene_counts.tabular", "") %>% str_replace("^[0-9]*_", "")
sampleNames
```

    ##  [1] "C57A_mock_A"    "C57A_mock_B"    "C57A_mock_C"    "C57A_treated_A"
    ##  [5] "C57A_treated_B" "C57A_treated_C" "C57B_mock_A"    "C57B_mock_B"   
    ##  [9] "C57B_mock_C"    "C57B_treated_A" "C57B_treated_B" "C57B_treated_C"
    ## [13] "C57C_mock_A"    "C57C_mock_B"    "C57C_mock_C"    "C57C_treated_A"
    ## [17] "C57C_treated_B" "C57C_treated_C"

``` r
sampleDonor <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[1]))
sampleTreatment <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[2]))
sampleReplicate <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[3]))

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, donor = sampleDonor, treatment = sampleTreatment, replicate = sampleReplicate)
sampleTable
```

    ##        sampleName                                           fileName donor
    ## 1     C57A_mock_A    11_C57A_mock_A__htseq-count_Gene_counts.tabular  C57A
    ## 2     C57A_mock_B    11_C57A_mock_B__htseq-count_Gene_counts.tabular  C57A
    ## 3     C57A_mock_C    11_C57A_mock_C__htseq-count_Gene_counts.tabular  C57A
    ## 4  C57A_treated_A 11_C57A_treated_A__htseq-count_Gene_counts.tabular  C57A
    ## 5  C57A_treated_B 11_C57A_treated_B__htseq-count_Gene_counts.tabular  C57A
    ## 6  C57A_treated_C 11_C57A_treated_C__htseq-count_Gene_counts.tabular  C57A
    ## 7     C57B_mock_A    11_C57B_mock_A__htseq-count_Gene_counts.tabular  C57B
    ## 8     C57B_mock_B    11_C57B_mock_B__htseq-count_Gene_counts.tabular  C57B
    ## 9     C57B_mock_C    11_C57B_mock_C__htseq-count_Gene_counts.tabular  C57B
    ## 10 C57B_treated_A 11_C57B_treated_A__htseq-count_Gene_counts.tabular  C57B
    ## 11 C57B_treated_B 11_C57B_treated_B__htseq-count_Gene_counts.tabular  C57B
    ## 12 C57B_treated_C 11_C57B_treated_C__htseq-count_Gene_counts.tabular  C57B
    ## 13    C57C_mock_A    11_C57C_mock_A__htseq-count_Gene_counts.tabular  C57C
    ## 14    C57C_mock_B    11_C57C_mock_B__htseq-count_Gene_counts.tabular  C57C
    ## 15    C57C_mock_C    11_C57C_mock_C__htseq-count_Gene_counts.tabular  C57C
    ## 16 C57C_treated_A 11_C57C_treated_A__htseq-count_Gene_counts.tabular  C57C
    ## 17 C57C_treated_B 11_C57C_treated_B__htseq-count_Gene_counts.tabular  C57C
    ## 18 C57C_treated_C 11_C57C_treated_C__htseq-count_Gene_counts.tabular  C57C
    ##    treatment replicate
    ## 1       mock         A
    ## 2       mock         B
    ## 3       mock         C
    ## 4    treated         A
    ## 5    treated         B
    ## 6    treated         C
    ## 7       mock         A
    ## 8       mock         B
    ## 9       mock         C
    ## 10   treated         A
    ## 11   treated         B
    ## 12   treated         C
    ## 13      mock         A
    ## 14      mock         B
    ## 15      mock         C
    ## 16   treated         A
    ## 17   treated         B
    ## 18   treated         C

``` r
## Use DESeq2 to examine count data
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = data_dir, design = ~ donor + treatment)
dds
```

    ## class: DESeqDataSet 
    ## dim: 46078 18 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(46078): ENSMUSG00000000001 ENSMUSG00000000003 ...
    ##   ENSMUSG00000107391 ENSMUSG00000107392
    ## rowData names(0):
    ## colnames(18): C57A_mock_A C57A_mock_B ... C57C_treated_B
    ##   C57C_treated_C
    ## colData names(3): donor treatment replicate

``` r
##Determining size factor for each column
dds <- estimateSizeFactors(dds)
dds@colData
```

    ## DataFrame with 18 rows and 4 columns
    ##                   donor treatment replicate sizeFactor
    ##                <factor>  <factor>  <factor>  <numeric>
    ## C57A_mock_A        C57A      mock         A  1.1767032
    ## C57A_mock_B        C57A      mock         B  1.1619211
    ## C57A_mock_C        C57A      mock         C  0.8878449
    ## C57A_treated_A     C57A   treated         A  1.0687889
    ## C57A_treated_B     C57A   treated         B  1.2091169
    ## ...                 ...       ...       ...        ...
    ## C57C_mock_B        C57C      mock         B  1.1657031
    ## C57C_mock_C        C57C      mock         C  0.9479510
    ## C57C_treated_A     C57C   treated         A  0.7878581
    ## C57C_treated_B     C57C   treated         B  1.0291789
    ## C57C_treated_C     C57C   treated         C  0.7502351

``` r
rld <- rlogTransformation(dds, blind = TRUE)

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
paste(donor, treatment, replicate, sep=" : "))
```

Making a heat map of the sample to sample distances matrix, adding color bars to indicate treatment as well as species.

``` r
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
```

    ## SHA-1 hash of file is 015fc0457e61e3e93a903e69a24d96d2dac7b9fb

``` r
colcoloring = function(donor) {
  ifelse(donor == "C57A", "deeppink4",
                       ifelse(donor == "C57B", "deeppink2",
                              ifelse(donor == "C57C", "lightpink", "grey")))
                                
}

colcolors <- unlist(lapply(sampleTable$donor, colcoloring))
      
treatcoloring <- function(treatment) { if (treatment=="mock") "red" else "orange" }
treatcolors <- unlist(lapply(sampleTable$treatment, treatcoloring))

clab = cbind(treatcolors, colcolors)
colnames(clab) <- c("", "")
png(file = file.path(output_dir, paste(Sys.Date(), "mouse_distance clustering_labeled.png")), units = 'in', height = 16, width = 19, res = 300)
heatmap.3(mat, trace="none", keysize = 0.5, col = viridis(100), srtCol = 45, cexCol = 1.5, cexRow = 1.5, dendrogram = "column", density.info = "none", margin = c(11, 12),  Rowv = TRUE, Colv = TRUE, ColSideColors = clab, lhei = c(1.5, 7), lwid = c(1.5, 6), ColSideColorsSize=3) 
```

    ## Warning in plot.window(...): "srtCol" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "srtCol" is not a graphical parameter

    ## Warning in title(...): "srtCol" is not a graphical parameter

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png(file = file.path(output_dir, paste(Sys.Date(), "mouse_distance clustering_unlabeled.png")), units = 'in', height = 16, width = 18, res = 300)
heatmap.3(mat, trace="none", keysize = 0.8, col = viridis(100), dendrogram = "column", density.info = "none", margin = c(5, 5),  Rowv = TRUE, Colv = TRUE, ColSideColors = clab, lhei = c(1.5, 11), lwid = c(1.5,8), ColSideColorsSize=3.5, labCol = FALSE, labRow = FALSE) 
dev.off()
```

    ## quartz_off_screen 
    ##                 2

PCA plot with donor and treatment indicated by color and shape, respectively.

``` r
a <- plotPCA(rld, intgroup = c("donor", "treatment", "replicate"), returnData = TRUE)
donors_PCA_colors <- unlist(lapply(a$donor, colcoloring))
names(donors_PCA_colors) <- a$donor


b <- plotPCA(rld, intgroup = c("donor", "treatment", "replicate")) + aes(colour = a$donor, shape = a$treatment) +
  scale_colour_manual(values = donors_PCA_colors) +
  geom_point(aes(size = 5)) + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.x = element_text(face = "bold", size = 22), 
        axis.text = element_text(size = 20), 
        panel.grid.major = element_line(size = 0.65, 
        color = "gray69"), 
        panel.grid.minor = element_line(size = 0.3, color = "gray69"), 
        axis.line = element_line(size = 2), 
        axis.title.y = element_text(face = "bold", size = 18))
b
```

![](QC_mouse_counts_code_files/figure-markdown_github/PCA%20plot-1.png)

``` r
ggsave(file = file.path(output_dir, paste(Sys.Date(), "mouse_PCA.png")), plot = b, units = 'in', height = 5, width = 10, dpi = 300, device = "png")
```

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
    ##  [1] devtools_1.12.0            viridis_0.4.0             
    ##  [3] viridisLite_0.2.0          data.table_1.10.0         
    ##  [5] genefilter_1.54.2          stringr_1.2.0             
    ##  [7] RColorBrewer_1.1-2         dplyr_0.7.3               
    ##  [9] ggplot2_2.2.1              gplots_3.0.1              
    ## [11] DESeq2_1.12.4              SummarizedExperiment_1.2.3
    ## [13] Biobase_2.32.0             GenomicRanges_1.24.3      
    ## [15] GenomeInfoDb_1.8.7         IRanges_2.6.1             
    ## [17] S4Vectors_0.10.3           BiocGenerics_0.18.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.2.1           splines_3.3.3        gtools_3.5.0        
    ##  [4] Formula_1.2-1        assertthat_0.2.0     latticeExtra_0.6-28 
    ##  [7] yaml_2.1.14          RSQLite_1.1-2        backports_1.0.5     
    ## [10] lattice_0.20-35      glue_1.1.1           digest_0.6.12       
    ## [13] XVector_0.12.1       checkmate_1.8.2      colorspace_1.3-2    
    ## [16] htmltools_0.3.5      Matrix_1.2-8         plyr_1.8.4          
    ## [19] XML_3.98-1.9         pkgconfig_2.0.1      zlibbioc_1.18.0     
    ## [22] xtable_1.8-2         scales_0.4.1         gdata_2.17.0        
    ## [25] BiocParallel_1.6.6   htmlTable_1.9        tibble_1.3.3        
    ## [28] annotate_1.50.1      withr_1.0.2          nnet_7.3-12         
    ## [31] lazyeval_0.2.0       survival_2.41-3      magrittr_1.5        
    ## [34] memoise_1.0.0        evaluate_0.10        foreign_0.8-67      
    ## [37] tools_3.3.3          munsell_0.4.3        locfit_1.5-9.1      
    ## [40] cluster_2.0.6        AnnotationDbi_1.34.4 bindrcpp_0.2        
    ## [43] caTools_1.17.1       rlang_0.1.2          grid_3.3.3          
    ## [46] RCurl_1.95-4.8       htmlwidgets_0.9      labeling_0.3        
    ## [49] bitops_1.0-6         base64enc_0.1-3      rmarkdown_1.4       
    ## [52] gtable_0.2.0         curl_2.3             DBI_0.6-1           
    ## [55] R6_2.2.0             gridExtra_2.2.1      knitr_1.16          
    ## [58] bindr_0.1            Hmisc_4.0-2          rprojroot_1.2       
    ## [61] KernSmooth_2.23-15   stringi_1.1.5        Rcpp_0.12.10        
    ## [64] geneplotter_1.50.0   rpart_4.1-10         acepack_1.4.1
