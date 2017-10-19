Rhesus\_macaque\_DGE
================

Purpose:
--------

To determine the differential gene expression in primary rhesus macaque dermal fibroblasts transfected with 0.5 ug polyI:C for 24 h versus those mock transfected.

``` r
library(DESeq2)
library(gplots)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(genefilter)
library(data.table)
library(genefilter)
```

Upload gene count files and generate a table of the samples broken up by donor ID, replicate designation, and treatment (polyI:C treated or mock).

``` r
data_dir <- "Compiled_by_species/Rhesus macaque/data"
output_dir <- "DGE_per_species"
## File list
sampleFiles <- basename(Sys.glob(file.path(data_dir, "*Gene_counts.tabular")))
sampleFiles
```

    ##  [1] "20_AG08305_mock_B__htseq-count_Gene_counts.tabular"    
    ##  [2] "20_AG08308_treated_B__htseq-count_Gene_counts.tabular" 
    ##  [3] "20_AG08312_treated_B__htseq-count_Gene_counts.tabular" 
    ##  [4] "32_AG08305_mock_C__htseq-count_Gene_counts.tabular"    
    ##  [5] "32_AG08308_treated_C__htseq-count_Gene_counts.tabular" 
    ##  [6] "32_AG08312_treated_C__htseq-count_Gene_counts.tabular" 
    ##  [7] "44_AG08305_treated_A___htseq-count_Gene_counts.tabular"
    ##  [8] "44_AG08308_mock_A__htseq-count_Gene_counts.tabular"    
    ##  [9] "44_AG08312_mock_A__htseq-count_Gene_counts.tabular"    
    ## [10] "56_AG08305_treated_B__htseq-count_Gene_counts.tabular" 
    ## [11] "56_AG08308_mock_B__htseq-count_Gene_counts.tabular"    
    ## [12] "56_AG08312_mock_B__htseq-count_Gene_counts.tabular"    
    ## [13] "68_AG08305_treated_C__htseq-count_Gene_counts.tabular" 
    ## [14] "68_AG08308_mock_C__htseq-count_Gene_counts.tabular"    
    ## [15] "68_AG08312_mock_C__htseq-count_Gene_counts.tabular"    
    ## [16] "8_AG08305_mock_A___htseq-count_Gene_counts.tabular"    
    ## [17] "8_AG08308_treated_A__htseq-count_Gene_counts.tabular"  
    ## [18] "8_AG08312_treated_A__htseq-count_Gene_counts.tabular"

``` r
## Create table of samples with donor, treatment, and replicate information
sampleNames <- sampleFiles %>% str_replace("_Gene.counts.tabluar", "") %>% str_replace("__htseq-count", "") %>% str_replace("^[0-9]*_", "")
sampleNames
```

    ##  [1] "AG08305_mock_B_Gene_counts.tabular"    
    ##  [2] "AG08308_treated_B_Gene_counts.tabular" 
    ##  [3] "AG08312_treated_B_Gene_counts.tabular" 
    ##  [4] "AG08305_mock_C_Gene_counts.tabular"    
    ##  [5] "AG08308_treated_C_Gene_counts.tabular" 
    ##  [6] "AG08312_treated_C_Gene_counts.tabular" 
    ##  [7] "AG08305_treated_A__Gene_counts.tabular"
    ##  [8] "AG08308_mock_A_Gene_counts.tabular"    
    ##  [9] "AG08312_mock_A_Gene_counts.tabular"    
    ## [10] "AG08305_treated_B_Gene_counts.tabular" 
    ## [11] "AG08308_mock_B_Gene_counts.tabular"    
    ## [12] "AG08312_mock_B_Gene_counts.tabular"    
    ## [13] "AG08305_treated_C_Gene_counts.tabular" 
    ## [14] "AG08308_mock_C_Gene_counts.tabular"    
    ## [15] "AG08312_mock_C_Gene_counts.tabular"    
    ## [16] "AG08305_mock_A__Gene_counts.tabular"   
    ## [17] "AG08308_treated_A_Gene_counts.tabular" 
    ## [18] "AG08312_treated_A_Gene_counts.tabular"

``` r
sampleDonor <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[1]))
sampleTreatment <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[2]))
sampleReplicate <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[3]))

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, donor = sampleDonor, treatment = sampleTreatment, replicate = sampleReplicate)
sampleTable
```

    ##                                sampleName
    ## 1      AG08305_mock_B_Gene_counts.tabular
    ## 2   AG08308_treated_B_Gene_counts.tabular
    ## 3   AG08312_treated_B_Gene_counts.tabular
    ## 4      AG08305_mock_C_Gene_counts.tabular
    ## 5   AG08308_treated_C_Gene_counts.tabular
    ## 6   AG08312_treated_C_Gene_counts.tabular
    ## 7  AG08305_treated_A__Gene_counts.tabular
    ## 8      AG08308_mock_A_Gene_counts.tabular
    ## 9      AG08312_mock_A_Gene_counts.tabular
    ## 10  AG08305_treated_B_Gene_counts.tabular
    ## 11     AG08308_mock_B_Gene_counts.tabular
    ## 12     AG08312_mock_B_Gene_counts.tabular
    ## 13  AG08305_treated_C_Gene_counts.tabular
    ## 14     AG08308_mock_C_Gene_counts.tabular
    ## 15     AG08312_mock_C_Gene_counts.tabular
    ## 16    AG08305_mock_A__Gene_counts.tabular
    ## 17  AG08308_treated_A_Gene_counts.tabular
    ## 18  AG08312_treated_A_Gene_counts.tabular
    ##                                                  fileName   donor
    ## 1      20_AG08305_mock_B__htseq-count_Gene_counts.tabular AG08305
    ## 2   20_AG08308_treated_B__htseq-count_Gene_counts.tabular AG08308
    ## 3   20_AG08312_treated_B__htseq-count_Gene_counts.tabular AG08312
    ## 4      32_AG08305_mock_C__htseq-count_Gene_counts.tabular AG08305
    ## 5   32_AG08308_treated_C__htseq-count_Gene_counts.tabular AG08308
    ## 6   32_AG08312_treated_C__htseq-count_Gene_counts.tabular AG08312
    ## 7  44_AG08305_treated_A___htseq-count_Gene_counts.tabular AG08305
    ## 8      44_AG08308_mock_A__htseq-count_Gene_counts.tabular AG08308
    ## 9      44_AG08312_mock_A__htseq-count_Gene_counts.tabular AG08312
    ## 10  56_AG08305_treated_B__htseq-count_Gene_counts.tabular AG08305
    ## 11     56_AG08308_mock_B__htseq-count_Gene_counts.tabular AG08308
    ## 12     56_AG08312_mock_B__htseq-count_Gene_counts.tabular AG08312
    ## 13  68_AG08305_treated_C__htseq-count_Gene_counts.tabular AG08305
    ## 14     68_AG08308_mock_C__htseq-count_Gene_counts.tabular AG08308
    ## 15     68_AG08312_mock_C__htseq-count_Gene_counts.tabular AG08312
    ## 16     8_AG08305_mock_A___htseq-count_Gene_counts.tabular AG08305
    ## 17   8_AG08308_treated_A__htseq-count_Gene_counts.tabular AG08308
    ## 18   8_AG08312_treated_A__htseq-count_Gene_counts.tabular AG08312
    ##    treatment replicate
    ## 1       mock         B
    ## 2    treated         B
    ## 3    treated         B
    ## 4       mock         C
    ## 5    treated         C
    ## 6    treated         C
    ## 7    treated         A
    ## 8       mock         A
    ## 9       mock         A
    ## 10   treated         B
    ## 11      mock         B
    ## 12      mock         B
    ## 13   treated         C
    ## 14      mock         C
    ## 15      mock         C
    ## 16      mock         A
    ## 17   treated         A
    ## 18   treated         A

``` r
## Use DESeq2 to examine count data
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = data_dir, design = ~ donor + treatment)
dds
```

    ## class: DESeqDataSet 
    ## dim: 65217 18 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(65217): ENSG00000000003 ENSG00000000005 ...
    ##   ENSG00000281921 ENSG00000281922
    ## rowData names(0):
    ## colnames(18): AG08305_mock_B_Gene_counts.tabular
    ##   AG08308_treated_B_Gene_counts.tabular ...
    ##   AG08308_treated_A_Gene_counts.tabular
    ##   AG08312_treated_A_Gene_counts.tabular
    ## colData names(3): donor treatment replicate

``` r
##Interested in comparing treated versus mock
contrast <- c("treatment", "treated", "mock")

##How we want the output file named
output_basename <- sprintf("%s-%s_vs_%s_analysis", "Rhesus macaque", contrast[2], contrast[3])

##Determining size factor for each column
dds <- estimateSizeFactors(dds)
dds@colData
```

    ## DataFrame with 18 rows and 4 columns
    ##                                          donor treatment replicate
    ##                                       <factor>  <factor>  <factor>
    ## AG08305_mock_B_Gene_counts.tabular     AG08305      mock         B
    ## AG08308_treated_B_Gene_counts.tabular  AG08308   treated         B
    ## AG08312_treated_B_Gene_counts.tabular  AG08312   treated         B
    ## AG08305_mock_C_Gene_counts.tabular     AG08305      mock         C
    ## AG08308_treated_C_Gene_counts.tabular  AG08308   treated         C
    ## ...                                        ...       ...       ...
    ## AG08308_mock_C_Gene_counts.tabular     AG08308      mock         C
    ## AG08312_mock_C_Gene_counts.tabular     AG08312      mock         C
    ## AG08305_mock_A__Gene_counts.tabular    AG08305      mock         A
    ## AG08308_treated_A_Gene_counts.tabular  AG08308   treated         A
    ## AG08312_treated_A_Gene_counts.tabular  AG08312   treated         A
    ##                                       sizeFactor
    ##                                        <numeric>
    ## AG08305_mock_B_Gene_counts.tabular     0.8892133
    ## AG08308_treated_B_Gene_counts.tabular  1.1937314
    ## AG08312_treated_B_Gene_counts.tabular  1.0950224
    ## AG08305_mock_C_Gene_counts.tabular     1.0586643
    ## AG08308_treated_C_Gene_counts.tabular  0.9788578
    ## ...                                          ...
    ## AG08308_mock_C_Gene_counts.tabular     0.6504444
    ## AG08312_mock_C_Gene_counts.tabular     1.1856845
    ## AG08305_mock_A__Gene_counts.tabular    0.9324556
    ## AG08308_treated_A_Gene_counts.tabular  1.1934385
    ## AG08312_treated_A_Gene_counts.tabular  0.7486315

``` r
dds <- estimateDispersions(dds)
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
plotDispEsts(dds, main=sprintf("%s Dispersion Estimates", output_basename))
```

![](Rhesus_macaque_DGE_code_files/figure-markdown_github/Dispersion%20Estimates-1.png)

``` r
dds <- nbinomWaldTest(dds)
res <- results(dds, contrast=contrast)

sum(res$padj <= 0.05, na.rm=TRUE)
```

    ## [1] 4873

``` r
sum(res$padj <= 0.1, na.rm=TRUE)
```

    ## [1] 5662

``` r
res <- res[order(res$padj, -abs(res$log2FoldChange)),]
head(res)
```

    ## log2 fold change (MAP): treatment treated vs mock 
    ## Wald test p-value: treatment treated vs mock 
    ## DataFrame with 6 rows and 6 columns
    ##                  baseMean log2FoldChange     lfcSE      stat    pvalue
    ##                 <numeric>      <numeric> <numeric> <numeric> <numeric>
    ## ENSG00000157601 5856.0416       8.754289 0.1166586  75.04196         0
    ## ENSG00000111335 1881.5329       7.447484 0.1349374  55.19215         0
    ## ENSG00000137628  856.4702       6.758637 0.1519572  44.47726         0
    ## ENSG00000119917 2807.8443       6.248956 0.1658003  37.68967         0
    ## ENSG00000119922 4445.2762       5.994279 0.1512383  39.63467         0
    ## ENSG00000115267 1560.3931       5.768585 0.1138661  50.66111         0
    ##                      padj
    ##                 <numeric>
    ## ENSG00000157601         0
    ## ENSG00000111335         0
    ## ENSG00000137628         0
    ## ENSG00000119917         0
    ## ENSG00000119922         0
    ## ENSG00000115267         0

``` r
mcols(res, use.names=TRUE)
```

    ## DataFrame with 6 rows and 2 columns
    ##                        type
    ##                 <character>
    ## baseMean       intermediate
    ## log2FoldChange      results
    ## lfcSE               results
    ## stat                results
    ## pvalue              results
    ## padj                results
    ##                                                      description
    ##                                                      <character>
    ## baseMean               mean of normalized counts for all samples
    ## log2FoldChange log2 fold change (MAP): treatment treated vs mock
    ## lfcSE                  standard error: treatment treated vs mock
    ## stat                   Wald statistic: treatment treated vs mock
    ## pvalue              Wald test p-value: treatment treated vs mock
    ## padj                                        BH adjusted p-values

``` r
##Log-intensity ratios = M values, log-intensity averages = A values
##Red points indicate padj < 0.1. 
plotMA(res, alpha=0.1, main=sprintf(output_basename))
```

![](Rhesus_macaque_DGE_code_files/figure-markdown_github/MA%20Plot-1.png)

``` r
##FYI 
attr(res, "filterThreshold")
```

    ## NULL

``` r
metadata(res)$alpha
```

    ## [1] 0.1

``` r
metadata(res)$filterThreshold
```

    ## 76.94775% 
    ##  1.951283

``` r
plot(metadata(res)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
```

![](Rhesus_macaque_DGE_code_files/figure-markdown_github/FYI%20demonstration%20of%20independent%20filtering-1.png)

``` r
##source("http://bioconductor.org/biocLite.R")
##biocLite("org.Hs.eg.db")
require(org.Hs.eg.db)
```

    ## Loading required package: org.Hs.eg.db

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 

``` r
key = "ENSEMBL"
cols = c("SYMBOL", "ENTREZID", "GENENAME", "ALIAS", "REFSEQ", "ACCNUM")
for (col in cols) {
  # Get annotation data for column
  annotation_data <- AnnotationDbi::select(org.Hs.eg.db, rownames(res), col, keytype=key)
  # Collapse one-to-many relationships
  tmp <- aggregate(annotation_data[col], by=annotation_data[key],
                   # to a list
                   FUN=function(x)list(x))
  # Match on key and append to results
  idx <- match(rownames(res), tmp[[key]])
  res[[col]] <- tmp[idx,col]
}
```

    ## 'select()' returned 1:many mapping between keys and columns

    ## 'select()' returned 1:many mapping between keys and columns
    ## 'select()' returned 1:many mapping between keys and columns
    ## 'select()' returned 1:many mapping between keys and columns
    ## 'select()' returned 1:many mapping between keys and columns
    ## 'select()' returned 1:many mapping between keys and columns

``` r
head(res)
```

    ## log2 fold change (MAP): treatment treated vs mock 
    ## Wald test p-value: treatment treated vs mock 
    ## DataFrame with 6 rows and 12 columns
    ##                  baseMean log2FoldChange     lfcSE      stat    pvalue
    ##                 <numeric>      <numeric> <numeric> <numeric> <numeric>
    ## ENSG00000157601 5856.0416       8.754289 0.1166586  75.04196         0
    ## ENSG00000111335 1881.5329       7.447484 0.1349374  55.19215         0
    ## ENSG00000137628  856.4702       6.758637 0.1519572  44.47726         0
    ## ENSG00000119917 2807.8443       6.248956 0.1658003  37.68967         0
    ## ENSG00000119922 4445.2762       5.994279 0.1512383  39.63467         0
    ## ENSG00000115267 1560.3931       5.768585 0.1138661  50.66111         0
    ##                      padj   SYMBOL ENTREZID GENENAME    ALIAS   REFSEQ
    ##                 <numeric>   <list>   <list>   <list>   <list>   <list>
    ## ENSG00000157601         0 ######## ######## ######## ######## ########
    ## ENSG00000111335         0 ######## ######## ######## ######## ########
    ## ENSG00000137628         0 ######## ######## ######## ######## ########
    ## ENSG00000119917         0 ######## ######## ######## ######## ########
    ## ENSG00000119922         0 ######## ######## ######## ######## ########
    ## ENSG00000115267         0 ######## ######## ######## ######## ########
    ##                   ACCNUM
    ##                   <list>
    ## ENSG00000157601 ########
    ## ENSG00000111335 ########
    ## ENSG00000137628 ########
    ## ENSG00000119917 ########
    ## ENSG00000119922 ########
    ## ENSG00000115267 ########

``` r
output_data <- as.data.frame(res)
LIST_COLS <- sapply(output_data, is.list)
for (COL in colnames(output_data)[LIST_COLS]) {
  output_data[COL] <-
    sapply(output_data[COL],
           function(x)sapply(x, function(y) paste(unlist(y),
                                                  collapse=", ") ) )
}

# Save data frame above as tab-separated file
write.table(output_data,
            file=file.path(output_dir, paste(output_basename,
                                             "_results.txt", sep='')),
            quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
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
    ##  [1] org.Hs.eg.db_3.3.0         AnnotationDbi_1.34.4      
    ##  [3] data.table_1.10.0          genefilter_1.54.2         
    ##  [5] stringr_1.2.0              RColorBrewer_1.1-2        
    ##  [7] dplyr_0.7.3                gplots_3.0.1              
    ##  [9] DESeq2_1.12.4              SummarizedExperiment_1.2.3
    ## [11] Biobase_2.32.0             GenomicRanges_1.24.3      
    ## [13] GenomeInfoDb_1.8.7         IRanges_2.6.1             
    ## [15] S4Vectors_0.10.3           BiocGenerics_0.18.0       
    ## [17] knitr_1.16                
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.10        locfit_1.5-9.1      lattice_0.20-35    
    ##  [4] gtools_3.5.0        assertthat_0.2.0    rprojroot_1.2      
    ##  [7] digest_0.6.12       R6_2.2.0            plyr_1.8.4         
    ## [10] backports_1.0.5     acepack_1.4.1       RSQLite_1.1-2      
    ## [13] evaluate_0.10       ggplot2_2.2.1       zlibbioc_1.18.0    
    ## [16] rlang_0.1.2         lazyeval_0.2.0      annotate_1.50.1    
    ## [19] gdata_2.17.0        rpart_4.1-10        Matrix_1.2-8       
    ## [22] checkmate_1.8.2     rmarkdown_1.4       splines_3.3.3      
    ## [25] BiocParallel_1.6.6  geneplotter_1.50.0  foreign_0.8-67     
    ## [28] htmlwidgets_0.9     RCurl_1.95-4.8      munsell_0.4.3      
    ## [31] pkgconfig_2.0.1     base64enc_0.1-3     htmltools_0.3.5    
    ## [34] nnet_7.3-12         tibble_1.3.3        gridExtra_2.2.1    
    ## [37] htmlTable_1.9       Hmisc_4.0-2         XML_3.98-1.9       
    ## [40] bitops_1.0-6        grid_3.3.3          xtable_1.8-2       
    ## [43] gtable_0.2.0        DBI_0.6-1           magrittr_1.5       
    ## [46] scales_0.4.1        KernSmooth_2.23-15  stringi_1.1.5      
    ## [49] XVector_0.12.1      bindrcpp_0.2        latticeExtra_0.6-28
    ## [52] Formula_1.2-1       tools_3.3.3         glue_1.1.1         
    ## [55] survival_2.41-3     yaml_2.1.14         colorspace_1.3-2   
    ## [58] cluster_2.0.6       caTools_1.17.1      memoise_1.0.0      
    ## [61] bindr_0.1
