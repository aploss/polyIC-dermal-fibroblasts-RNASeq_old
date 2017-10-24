QC of NHP counts
================

Purpose:
--------

To see how all the NHP RNASeq read counts cluster.

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

Upload gene count files and generate a table of the samples broken up by donor ID, replicate designation, and treatment (polyI:C treated or mock). Note this is only for the NHP species which were all mapped to the human genome - mouse samples are run separately since the reads were mapped to the mouse genome.

``` r
data_dir <- "NHP counts"
output_dir <- "QC_NHP_counts_output"

## File list
sampleFiles <- basename(Sys.glob(file.path(data_dir, "*Gene_counts.tabular")))
sampleFiles
```

    ##   [1] "11_AG07923_mock_A__htseq-count_Gene_counts.tabular"      
    ##   [2] "11_AG07923_mock_B__htseq-count_Gene_counts.tabular"      
    ##   [3] "11_AG07923_mock_C__htseq-count_Gene_counts.tabular"      
    ##   [4] "11_AG07923_treated_A__htseq-count_Gene_counts.tabular"   
    ##   [5] "11_AG07923_treated_B__htseq-count_Gene_counts.tabular"   
    ##   [6] "11_AG07923_treated_C__htseq-count_Gene_counts.tabular"   
    ##   [7] "11_AG08490_mock_A__htseq-count_Gene_counts.tabular"      
    ##   [8] "11_AG08490_mock_B__htseq-count_Gene_counts.tabular"      
    ##   [9] "11_AG08490_treated_A__htseq-count_Gene_counts.tabular"   
    ##  [10] "11_AG08490_treated_B__htseq-count_Gene_counts.tabular"   
    ##  [11] "11_AG08490_treated_C__htseq-count_Gene_counts.tabular"   
    ##  [12] "11_PR0058_mock_A__htseq-count_Gene_counts.tabular"       
    ##  [13] "11_PR0058_mock_B__htseq-count_Gene_counts.tabular"       
    ##  [14] "11_PR0058_mock_C__htseq-count_Gene_counts.tabular"       
    ##  [15] "11_PR0058_treated_A___htseq-count_Gene_counts.tabular"   
    ##  [16] "11_PR0058_treated_B__htseq-count_Gene_counts.tabular"    
    ##  [17] "11_PR0058_treated_C__htseq-count_Gene_counts.tabular"    
    ##  [18] "20_AG05311_mock_B__htseq-count_Gene_counts.tabular"      
    ##  [19] "20_AG08305_mock_B__htseq-count_Gene_counts.tabular"      
    ##  [20] "20_AG08308_treated_B__htseq-count_Gene_counts.tabular"   
    ##  [21] "20_AG08312_treated_B__htseq-count_Gene_counts.tabular"   
    ##  [22] "20_PR00248_mock_B__htseq-count_Gene_counts.tabular"      
    ##  [23] "20_PR111_mock_B__htseq-count_Gene_counts.tabular"        
    ##  [24] "20_PR235_mock_B__htseq-count_Gene_counts.tabular"        
    ##  [25] "20_SQMA_mock_A__htseq-count_Gene_counts.tabular"         
    ##  [26] "20_SQMB_mock_B__htseq-count_Gene_counts.tabular"         
    ##  [27] "32_AG05311_mock_C__htseq-count_Gene_counts.tabular"      
    ##  [28] "32_AG08305_mock_C__htseq-count_Gene_counts.tabular"      
    ##  [29] "32_AG08308_treated_C__htseq-count_Gene_counts.tabular"   
    ##  [30] "32_AG08312_treated_C__htseq-count_Gene_counts.tabular"   
    ##  [31] "32_PR00248_mock_C__htseq-count_Gene_counts.tabular"      
    ##  [32] "32_PR111_mock_C__htseq-count_Gene_counts.tabular"        
    ##  [33] "32_PR235_mock_C__htseq-count_Gene_counts.tabular"        
    ##  [34] "32_SQMA_mock_C__htseq-count_Gene_counts.tabular"         
    ##  [35] "32_SQMB_mock_C__htseq-count_Gene_counts.tabular"         
    ##  [36] "44_AG05311_treated_A__htseq-count_Gene_counts.tabular"   
    ##  [37] "44_AG08305_treated_A___htseq-count_Gene_counts.tabular"  
    ##  [38] "44_AG08308_mock_A__htseq-count_Gene_counts.tabular"      
    ##  [39] "44_AG08312_mock_A__htseq-count_Gene_counts.tabular"      
    ##  [40] "44_PR00248_treated_A__htseq-count_Gene_counts.tabular"   
    ##  [41] "44_PR111_treated_A__htseq-count_Gene_counts.tabular"     
    ##  [42] "44_PR235_treated_A__htseq-count_Gene_counts.tabular"     
    ##  [43] "44_SQMA_treated_A__htseq-count_Gene_counts.tabular"      
    ##  [44] "44_SQMB_treated_B__htseq-count_Gene_counts.tabular"      
    ##  [45] "56_AG05311_treated_B__htseq-count_Gene_counts.tabular"   
    ##  [46] "56_AG08305_treated_B__htseq-count_Gene_counts.tabular"   
    ##  [47] "56_AG08308_mock_B__htseq-count_Gene_counts.tabular"      
    ##  [48] "56_AG08312_mock_B__htseq-count_Gene_counts.tabular"      
    ##  [49] "56_PR00248_treated_B__htseq-count_Gene_counts.tabular"   
    ##  [50] "56_PR111_treated_B__htseq-count_Gene_counts.tabular"     
    ##  [51] "56_PR235_treated_B__htseq-count_Gene_counts.tabular"     
    ##  [52] "56_SQMA_treated_B__htseq-count_Gene_counts.tabular"      
    ##  [53] "56_SQMB_treated_C__htseq-count_Gene_counts.tabular"      
    ##  [54] "68_AG05311_treated_C__htseq-count_Gene_counts.tabular"   
    ##  [55] "68_AG08305_treated_C__htseq-count_Gene_counts.tabular"   
    ##  [56] "68_AG08308_mock_C__htseq-count_Gene_counts.tabular"      
    ##  [57] "68_AG08312_mock_C__htseq-count_Gene_counts.tabular"      
    ##  [58] "68_PR00248_treated_C__htseq-count_Gene_counts.tabular"   
    ##  [59] "68_PR111_treated_C__htseq-count_Gene_counts.tabular"     
    ##  [60] "68_PR235_treated_C__htseq-count_Gene_counts.tabular"     
    ##  [61] "68_SQMA_treated_C__htseq-count_Gene_counts.tabular"      
    ##  [62] "68_SQMB_treated_A__htseq-count_Gene_counts.tabular"      
    ##  [63] "8_AF_mock_A__htseq-count_Gene_counts.tabular"            
    ##  [64] "8_AF_mock_B__htseq-count_Gene_counts.tabular"            
    ##  [65] "8_AF_mock_C__htseq-count_Gene_counts.tabular"            
    ##  [66] "8_AF_treated_A__htseq-count_Gene_counts.tabular"         
    ##  [67] "8_AF_treated_B__htseq-count_Gene_counts.tabular"         
    ##  [68] "8_AF_treated_C__htseq-count_Gene_counts.tabular"         
    ##  [69] "8_AG05311_mock_A__htseq-count_Gene_counts.tabular"       
    ##  [70] "8_AG06105_mock_A__htseq-count_Gene_counts.tabular"       
    ##  [71] "8_AG06105_mock_B__htseq-count_Gene_counts.tabular"       
    ##  [72] "8_AG06105_mock_C___htseq-count_Gene_counts.tabular"      
    ##  [73] "8_AG06105_treated_A__htseq-count_Gene_counts.tabular"    
    ##  [74] "8_AG06105_treated_B__htseq-count_Gene_counts.tabular"    
    ##  [75] "8_AG06105_treated_C__htseq-count_Gene_counts.tabular"    
    ##  [76] "8_AG08305_mock_A___htseq-count_Gene_counts.tabular"      
    ##  [77] "8_AG08308_treated_A__htseq-count_Gene_counts.tabular"    
    ##  [78] "8_AG08312_treated_A__htseq-count_Gene_counts.tabular"    
    ##  [79] "8_AG08490_mock_C_Sunday__htseq-count_Gene_counts.tabular"
    ##  [80] "8_NHDF_mock_A___htseq-count_Gene_counts.tabular"         
    ##  [81] "8_NHDF_mock_B_htseq-count_Gene_counts.tabular"           
    ##  [82] "8_NHDF_mock_C___htseq-count_Gene_counts.tabular"         
    ##  [83] "8_NHDF_treated_A___htseq-count_Gene_counts.tabular"      
    ##  [84] "8_NHDF_treated_B__htseq-count_Gene_counts.tabular"       
    ##  [85] "8_NHDF_treated_C___htseq-count_Gene_counts.tabular"      
    ##  [86] "8_PR00033_mock_A__htseq-count_Gene_counts.tabular"       
    ##  [87] "8_PR00033_mock_B__htseq-count_Gene_counts.tabular"       
    ##  [88] "8_PR00033_mock_C___htseq-count_Gene_counts.tabular"      
    ##  [89] "8_PR00033_treated_A__htseq-count_Gene_counts.tabular"    
    ##  [90] "8_PR00033_treated_B__htseq-count_Gene_counts.tabular"    
    ##  [91] "8_PR00033_treated_C___htseq-count_Gene_counts.tabular"   
    ##  [92] "8_PR00036_mock_A__htseq-count_Gene_counts.tabular"       
    ##  [93] "8_PR00036_mock_B__htseq-count_Gene_counts.tabular"       
    ##  [94] "8_PR00036_mock_C__htseq-count_Gene_counts.tabular"       
    ##  [95] "8_PR00036_treated_A__htseq-count_Gene_counts.tabular"    
    ##  [96] "8_PR00036_treated_B__htseq-count_Gene_counts.tabular"    
    ##  [97] "8_PR00036_treated_C___htseq-count_Gene_counts.tabular"   
    ##  [98] "8_PR00039_mock_A__htseq-count_Gene_counts.tabular"       
    ##  [99] "8_PR00039_mock_B__htseq-count_Gene_counts.tabular"       
    ## [100] "8_PR00039_mock_C__htseq-count_Gene_counts.tabular"       
    ## [101] "8_PR00039_treated_A__htseq-count_Gene_counts.tabular"    
    ## [102] "8_PR00039_treated_B__htseq-count_Gene_counts.tabular"    
    ## [103] "8_PR00039_treated_C___htseq-count_Gene_counts.tabular"   
    ## [104] "8_PR00054_mock_A__htseq-count_Gene_counts.tabular"       
    ## [105] "8_PR00054_mock_B__htseq-count_Gene_counts.tabular"       
    ## [106] "8_PR00054_mock_C__htseq-count_Gene_counts.tabular"       
    ## [107] "8_PR00054_treated_A___htseq-count_Gene_counts.tabular"   
    ## [108] "8_PR00054_treated_B__htseq-count_Gene_counts.tabular"    
    ## [109] "8_PR00054_treated_C__htseq-count_Gene_counts.tabular"    
    ## [110] "8_PR00107_mock_A__htseq-count_Gene_counts.tabular"       
    ## [111] "8_PR00107_mock_B__htseq-count_Gene_counts.tabular"       
    ## [112] "8_PR00107_mock_C__htseq-count_Gene_counts.tabular"       
    ## [113] "8_PR00107_treated_A__htseq-count_Gene_counts.tabular"    
    ## [114] "8_PR00107_treated_B__htseq-count_Gene_counts.tabular"    
    ## [115] "8_PR00107_treated_C__htseq-count_Gene_counts.tabular"    
    ## [116] "8_PR00248_mock_A__htseq-count_Gene_counts.tabular"       
    ## [117] "8_PR00573_mock_A__htseq-count_Gene_counts.tabular"       
    ## [118] "8_PR00573_mock_B__htseq-count_Gene_counts.tabular"       
    ## [119] "8_PR00573_mock_C__htseq-count_Gene_counts.tabular"       
    ## [120] "8_PR00573_treated_A__htseq-count_Gene_counts.tabular"    
    ## [121] "8_PR00573_treated_B__htseq-count_Gene_counts.tabular"    
    ## [122] "8_PR00573_treated_C___htseq-count_Gene_counts.tabular"   
    ## [123] "8_PR01109_mock_A__htseq-count_Gene_counts.tabular"       
    ## [124] "8_PR01109_mock_B__htseq-count_Gene_counts.tabular"       
    ## [125] "8_PR01109_mock_C__htseq-count_Gene_counts.tabular"       
    ## [126] "8_PR01109_treated_A__htseq-count_Gene_counts.tabular"    
    ## [127] "8_PR01109_treated_B__htseq-count_Gene_counts.tabular"    
    ## [128] "8_PR01109_treated_C__htseq-count_Gene_counts.tabular"    
    ## [129] "8_PR0230_mock_A__htseq-count_Gene_counts.tabular"        
    ## [130] "8_PR0230_mock_B__htseq-count_Gene_counts.tabular"        
    ## [131] "8_PR0230_mock_C___htseq-count_Gene_counts.tabular"       
    ## [132] "8_PR0230_treated_A__htseq-count_Gene_counts.tabular"     
    ## [133] "8_PR0230_treated_B__htseq-count_Gene_counts.tabular"     
    ## [134] "8_PR0230_treated_C___htseq-count_Gene_counts.tabular"    
    ## [135] "8_PR111_mock_A__htseq-count_Gene_counts.tabular"         
    ## [136] "8_PR235_mock_A__htseq-count_Gene_counts.tabular"         
    ## [137] "8_S003611_mock_A__htseq-count_Gene_counts.tabular"       
    ## [138] "8_S003611_mock_B__htseq-count_Gene_counts.tabular"       
    ## [139] "8_S003611_mock_C__htseq-count_Gene_counts.tabular"       
    ## [140] "8_S003611_treated_A___htseq-count_Gene_counts.tabular"   
    ## [141] "8_S003611_treated_B__htseq-count_Gene_counts.tabular"    
    ## [142] "8_S003611_treated_C__htseq-count_Gene_counts.tabular"    
    ## [143] "8_S003649_mock_A__htseq-count_Gene_counts.tabular"       
    ## [144] "8_S003649_mock_B__htseq-count_Gene_counts.tabular"       
    ## [145] "8_S003649_mock_C__htseq-count_Gene_counts.tabular"       
    ## [146] "8_S003649_treated_A__htseq-count_Gene_counts.tabular"    
    ## [147] "8_S003649_treated_B__htseq-count_Gene_counts.tabular"    
    ## [148] "8_S003649_treated_C___htseq-count_Gene_counts.tabular"   
    ## [149] "8_S004933_mock_A__htseq-count_Gene_counts.tabular"       
    ## [150] "8_S004933_mock_B__htseq-count_Gene_counts.tabular"       
    ## [151] "8_S004933_mock_C__htseq-count_Gene_counts.tabular"       
    ## [152] "8_S004933_treated_A__htseq-count_Gene_counts.tabular"    
    ## [153] "8_S004933_treated_B__htseq-count_Gene_counts.tabular"    
    ## [154] "8_S004933_treated_C___htseq-count_Gene_counts.tabular"   
    ## [155] "8_SQMA_mock_B__htseq-count_Gene_counts.tabular"          
    ## [156] "8_SQMB_mock_A__htseq-count_Gene_counts.tabular"          
    ## [157] "8_SR_mock_A__htseq-count_Gene_counts.tabular"            
    ## [158] "8_SR_mock_B__htseq-count_Gene_counts.tabular"            
    ## [159] "8_SR_mock_C__htseq-count_Gene_counts.tabular"            
    ## [160] "8_SR_treated_A__htseq-count_Gene_counts.tabular"         
    ## [161] "8_SR_treated_B__htseq-count_Gene_counts.tabular"         
    ## [162] "8_SR_treated_C__htseq-count_Gene_counts.tabular"

``` r
## Create table of samples with donor, treatment, and replicate information
sampleNames <- sampleFiles %>% str_replace("_*htseq-count_Gene_counts.tabular", "") %>% str_replace("^[0-9]*_", "")
sampleNames
```

    ##   [1] "AG07923_mock_A"        "AG07923_mock_B"       
    ##   [3] "AG07923_mock_C"        "AG07923_treated_A"    
    ##   [5] "AG07923_treated_B"     "AG07923_treated_C"    
    ##   [7] "AG08490_mock_A"        "AG08490_mock_B"       
    ##   [9] "AG08490_treated_A"     "AG08490_treated_B"    
    ##  [11] "AG08490_treated_C"     "PR0058_mock_A"        
    ##  [13] "PR0058_mock_B"         "PR0058_mock_C"        
    ##  [15] "PR0058_treated_A"      "PR0058_treated_B"     
    ##  [17] "PR0058_treated_C"      "AG05311_mock_B"       
    ##  [19] "AG08305_mock_B"        "AG08308_treated_B"    
    ##  [21] "AG08312_treated_B"     "PR00248_mock_B"       
    ##  [23] "PR111_mock_B"          "PR235_mock_B"         
    ##  [25] "SQMA_mock_A"           "SQMB_mock_B"          
    ##  [27] "AG05311_mock_C"        "AG08305_mock_C"       
    ##  [29] "AG08308_treated_C"     "AG08312_treated_C"    
    ##  [31] "PR00248_mock_C"        "PR111_mock_C"         
    ##  [33] "PR235_mock_C"          "SQMA_mock_C"          
    ##  [35] "SQMB_mock_C"           "AG05311_treated_A"    
    ##  [37] "AG08305_treated_A"     "AG08308_mock_A"       
    ##  [39] "AG08312_mock_A"        "PR00248_treated_A"    
    ##  [41] "PR111_treated_A"       "PR235_treated_A"      
    ##  [43] "SQMA_treated_A"        "SQMB_treated_B"       
    ##  [45] "AG05311_treated_B"     "AG08305_treated_B"    
    ##  [47] "AG08308_mock_B"        "AG08312_mock_B"       
    ##  [49] "PR00248_treated_B"     "PR111_treated_B"      
    ##  [51] "PR235_treated_B"       "SQMA_treated_B"       
    ##  [53] "SQMB_treated_C"        "AG05311_treated_C"    
    ##  [55] "AG08305_treated_C"     "AG08308_mock_C"       
    ##  [57] "AG08312_mock_C"        "PR00248_treated_C"    
    ##  [59] "PR111_treated_C"       "PR235_treated_C"      
    ##  [61] "SQMA_treated_C"        "SQMB_treated_A"       
    ##  [63] "AF_mock_A"             "AF_mock_B"            
    ##  [65] "AF_mock_C"             "AF_treated_A"         
    ##  [67] "AF_treated_B"          "AF_treated_C"         
    ##  [69] "AG05311_mock_A"        "AG06105_mock_A"       
    ##  [71] "AG06105_mock_B"        "AG06105_mock_C"       
    ##  [73] "AG06105_treated_A"     "AG06105_treated_B"    
    ##  [75] "AG06105_treated_C"     "AG08305_mock_A"       
    ##  [77] "AG08308_treated_A"     "AG08312_treated_A"    
    ##  [79] "AG08490_mock_C_Sunday" "NHDF_mock_A"          
    ##  [81] "NHDF_mock_B"           "NHDF_mock_C"          
    ##  [83] "NHDF_treated_A"        "NHDF_treated_B"       
    ##  [85] "NHDF_treated_C"        "PR00033_mock_A"       
    ##  [87] "PR00033_mock_B"        "PR00033_mock_C"       
    ##  [89] "PR00033_treated_A"     "PR00033_treated_B"    
    ##  [91] "PR00033_treated_C"     "PR00036_mock_A"       
    ##  [93] "PR00036_mock_B"        "PR00036_mock_C"       
    ##  [95] "PR00036_treated_A"     "PR00036_treated_B"    
    ##  [97] "PR00036_treated_C"     "PR00039_mock_A"       
    ##  [99] "PR00039_mock_B"        "PR00039_mock_C"       
    ## [101] "PR00039_treated_A"     "PR00039_treated_B"    
    ## [103] "PR00039_treated_C"     "PR00054_mock_A"       
    ## [105] "PR00054_mock_B"        "PR00054_mock_C"       
    ## [107] "PR00054_treated_A"     "PR00054_treated_B"    
    ## [109] "PR00054_treated_C"     "PR00107_mock_A"       
    ## [111] "PR00107_mock_B"        "PR00107_mock_C"       
    ## [113] "PR00107_treated_A"     "PR00107_treated_B"    
    ## [115] "PR00107_treated_C"     "PR00248_mock_A"       
    ## [117] "PR00573_mock_A"        "PR00573_mock_B"       
    ## [119] "PR00573_mock_C"        "PR00573_treated_A"    
    ## [121] "PR00573_treated_B"     "PR00573_treated_C"    
    ## [123] "PR01109_mock_A"        "PR01109_mock_B"       
    ## [125] "PR01109_mock_C"        "PR01109_treated_A"    
    ## [127] "PR01109_treated_B"     "PR01109_treated_C"    
    ## [129] "PR0230_mock_A"         "PR0230_mock_B"        
    ## [131] "PR0230_mock_C"         "PR0230_treated_A"     
    ## [133] "PR0230_treated_B"      "PR0230_treated_C"     
    ## [135] "PR111_mock_A"          "PR235_mock_A"         
    ## [137] "S003611_mock_A"        "S003611_mock_B"       
    ## [139] "S003611_mock_C"        "S003611_treated_A"    
    ## [141] "S003611_treated_B"     "S003611_treated_C"    
    ## [143] "S003649_mock_A"        "S003649_mock_B"       
    ## [145] "S003649_mock_C"        "S003649_treated_A"    
    ## [147] "S003649_treated_B"     "S003649_treated_C"    
    ## [149] "S004933_mock_A"        "S004933_mock_B"       
    ## [151] "S004933_mock_C"        "S004933_treated_A"    
    ## [153] "S004933_treated_B"     "S004933_treated_C"    
    ## [155] "SQMA_mock_B"           "SQMB_mock_A"          
    ## [157] "SR_mock_A"             "SR_mock_B"            
    ## [159] "SR_mock_C"             "SR_treated_A"         
    ## [161] "SR_treated_B"          "SR_treated_C"

``` r
sampleDonor <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[1]))
sampleTreatment <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[2]))
sampleReplicate <- sapply(strsplit(sampleNames, split = '_', fixed = TRUE), function(x) (x[3]))

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, donor = sampleDonor, treatment = sampleTreatment, replicate = sampleReplicate)
sampleTable
```

    ##                sampleName
    ## 1          AG07923_mock_A
    ## 2          AG07923_mock_B
    ## 3          AG07923_mock_C
    ## 4       AG07923_treated_A
    ## 5       AG07923_treated_B
    ## 6       AG07923_treated_C
    ## 7          AG08490_mock_A
    ## 8          AG08490_mock_B
    ## 9       AG08490_treated_A
    ## 10      AG08490_treated_B
    ## 11      AG08490_treated_C
    ## 12          PR0058_mock_A
    ## 13          PR0058_mock_B
    ## 14          PR0058_mock_C
    ## 15       PR0058_treated_A
    ## 16       PR0058_treated_B
    ## 17       PR0058_treated_C
    ## 18         AG05311_mock_B
    ## 19         AG08305_mock_B
    ## 20      AG08308_treated_B
    ## 21      AG08312_treated_B
    ## 22         PR00248_mock_B
    ## 23           PR111_mock_B
    ## 24           PR235_mock_B
    ## 25            SQMA_mock_A
    ## 26            SQMB_mock_B
    ## 27         AG05311_mock_C
    ## 28         AG08305_mock_C
    ## 29      AG08308_treated_C
    ## 30      AG08312_treated_C
    ## 31         PR00248_mock_C
    ## 32           PR111_mock_C
    ## 33           PR235_mock_C
    ## 34            SQMA_mock_C
    ## 35            SQMB_mock_C
    ## 36      AG05311_treated_A
    ## 37      AG08305_treated_A
    ## 38         AG08308_mock_A
    ## 39         AG08312_mock_A
    ## 40      PR00248_treated_A
    ## 41        PR111_treated_A
    ## 42        PR235_treated_A
    ## 43         SQMA_treated_A
    ## 44         SQMB_treated_B
    ## 45      AG05311_treated_B
    ## 46      AG08305_treated_B
    ## 47         AG08308_mock_B
    ## 48         AG08312_mock_B
    ## 49      PR00248_treated_B
    ## 50        PR111_treated_B
    ## 51        PR235_treated_B
    ## 52         SQMA_treated_B
    ## 53         SQMB_treated_C
    ## 54      AG05311_treated_C
    ## 55      AG08305_treated_C
    ## 56         AG08308_mock_C
    ## 57         AG08312_mock_C
    ## 58      PR00248_treated_C
    ## 59        PR111_treated_C
    ## 60        PR235_treated_C
    ## 61         SQMA_treated_C
    ## 62         SQMB_treated_A
    ## 63              AF_mock_A
    ## 64              AF_mock_B
    ## 65              AF_mock_C
    ## 66           AF_treated_A
    ## 67           AF_treated_B
    ## 68           AF_treated_C
    ## 69         AG05311_mock_A
    ## 70         AG06105_mock_A
    ## 71         AG06105_mock_B
    ## 72         AG06105_mock_C
    ## 73      AG06105_treated_A
    ## 74      AG06105_treated_B
    ## 75      AG06105_treated_C
    ## 76         AG08305_mock_A
    ## 77      AG08308_treated_A
    ## 78      AG08312_treated_A
    ## 79  AG08490_mock_C_Sunday
    ## 80            NHDF_mock_A
    ## 81            NHDF_mock_B
    ## 82            NHDF_mock_C
    ## 83         NHDF_treated_A
    ## 84         NHDF_treated_B
    ## 85         NHDF_treated_C
    ## 86         PR00033_mock_A
    ## 87         PR00033_mock_B
    ## 88         PR00033_mock_C
    ## 89      PR00033_treated_A
    ## 90      PR00033_treated_B
    ## 91      PR00033_treated_C
    ## 92         PR00036_mock_A
    ## 93         PR00036_mock_B
    ## 94         PR00036_mock_C
    ## 95      PR00036_treated_A
    ## 96      PR00036_treated_B
    ## 97      PR00036_treated_C
    ## 98         PR00039_mock_A
    ## 99         PR00039_mock_B
    ## 100        PR00039_mock_C
    ## 101     PR00039_treated_A
    ## 102     PR00039_treated_B
    ## 103     PR00039_treated_C
    ## 104        PR00054_mock_A
    ## 105        PR00054_mock_B
    ## 106        PR00054_mock_C
    ## 107     PR00054_treated_A
    ## 108     PR00054_treated_B
    ## 109     PR00054_treated_C
    ## 110        PR00107_mock_A
    ## 111        PR00107_mock_B
    ## 112        PR00107_mock_C
    ## 113     PR00107_treated_A
    ## 114     PR00107_treated_B
    ## 115     PR00107_treated_C
    ## 116        PR00248_mock_A
    ## 117        PR00573_mock_A
    ## 118        PR00573_mock_B
    ## 119        PR00573_mock_C
    ## 120     PR00573_treated_A
    ## 121     PR00573_treated_B
    ## 122     PR00573_treated_C
    ## 123        PR01109_mock_A
    ## 124        PR01109_mock_B
    ## 125        PR01109_mock_C
    ## 126     PR01109_treated_A
    ## 127     PR01109_treated_B
    ## 128     PR01109_treated_C
    ## 129         PR0230_mock_A
    ## 130         PR0230_mock_B
    ## 131         PR0230_mock_C
    ## 132      PR0230_treated_A
    ## 133      PR0230_treated_B
    ## 134      PR0230_treated_C
    ## 135          PR111_mock_A
    ## 136          PR235_mock_A
    ## 137        S003611_mock_A
    ## 138        S003611_mock_B
    ## 139        S003611_mock_C
    ## 140     S003611_treated_A
    ## 141     S003611_treated_B
    ## 142     S003611_treated_C
    ## 143        S003649_mock_A
    ## 144        S003649_mock_B
    ## 145        S003649_mock_C
    ## 146     S003649_treated_A
    ## 147     S003649_treated_B
    ## 148     S003649_treated_C
    ## 149        S004933_mock_A
    ## 150        S004933_mock_B
    ## 151        S004933_mock_C
    ## 152     S004933_treated_A
    ## 153     S004933_treated_B
    ## 154     S004933_treated_C
    ## 155           SQMA_mock_B
    ## 156           SQMB_mock_A
    ## 157             SR_mock_A
    ## 158             SR_mock_B
    ## 159             SR_mock_C
    ## 160          SR_treated_A
    ## 161          SR_treated_B
    ## 162          SR_treated_C
    ##                                                     fileName   donor
    ## 1         11_AG07923_mock_A__htseq-count_Gene_counts.tabular AG07923
    ## 2         11_AG07923_mock_B__htseq-count_Gene_counts.tabular AG07923
    ## 3         11_AG07923_mock_C__htseq-count_Gene_counts.tabular AG07923
    ## 4      11_AG07923_treated_A__htseq-count_Gene_counts.tabular AG07923
    ## 5      11_AG07923_treated_B__htseq-count_Gene_counts.tabular AG07923
    ## 6      11_AG07923_treated_C__htseq-count_Gene_counts.tabular AG07923
    ## 7         11_AG08490_mock_A__htseq-count_Gene_counts.tabular AG08490
    ## 8         11_AG08490_mock_B__htseq-count_Gene_counts.tabular AG08490
    ## 9      11_AG08490_treated_A__htseq-count_Gene_counts.tabular AG08490
    ## 10     11_AG08490_treated_B__htseq-count_Gene_counts.tabular AG08490
    ## 11     11_AG08490_treated_C__htseq-count_Gene_counts.tabular AG08490
    ## 12         11_PR0058_mock_A__htseq-count_Gene_counts.tabular  PR0058
    ## 13         11_PR0058_mock_B__htseq-count_Gene_counts.tabular  PR0058
    ## 14         11_PR0058_mock_C__htseq-count_Gene_counts.tabular  PR0058
    ## 15     11_PR0058_treated_A___htseq-count_Gene_counts.tabular  PR0058
    ## 16      11_PR0058_treated_B__htseq-count_Gene_counts.tabular  PR0058
    ## 17      11_PR0058_treated_C__htseq-count_Gene_counts.tabular  PR0058
    ## 18        20_AG05311_mock_B__htseq-count_Gene_counts.tabular AG05311
    ## 19        20_AG08305_mock_B__htseq-count_Gene_counts.tabular AG08305
    ## 20     20_AG08308_treated_B__htseq-count_Gene_counts.tabular AG08308
    ## 21     20_AG08312_treated_B__htseq-count_Gene_counts.tabular AG08312
    ## 22        20_PR00248_mock_B__htseq-count_Gene_counts.tabular PR00248
    ## 23          20_PR111_mock_B__htseq-count_Gene_counts.tabular   PR111
    ## 24          20_PR235_mock_B__htseq-count_Gene_counts.tabular   PR235
    ## 25           20_SQMA_mock_A__htseq-count_Gene_counts.tabular    SQMA
    ## 26           20_SQMB_mock_B__htseq-count_Gene_counts.tabular    SQMB
    ## 27        32_AG05311_mock_C__htseq-count_Gene_counts.tabular AG05311
    ## 28        32_AG08305_mock_C__htseq-count_Gene_counts.tabular AG08305
    ## 29     32_AG08308_treated_C__htseq-count_Gene_counts.tabular AG08308
    ## 30     32_AG08312_treated_C__htseq-count_Gene_counts.tabular AG08312
    ## 31        32_PR00248_mock_C__htseq-count_Gene_counts.tabular PR00248
    ## 32          32_PR111_mock_C__htseq-count_Gene_counts.tabular   PR111
    ## 33          32_PR235_mock_C__htseq-count_Gene_counts.tabular   PR235
    ## 34           32_SQMA_mock_C__htseq-count_Gene_counts.tabular    SQMA
    ## 35           32_SQMB_mock_C__htseq-count_Gene_counts.tabular    SQMB
    ## 36     44_AG05311_treated_A__htseq-count_Gene_counts.tabular AG05311
    ## 37    44_AG08305_treated_A___htseq-count_Gene_counts.tabular AG08305
    ## 38        44_AG08308_mock_A__htseq-count_Gene_counts.tabular AG08308
    ## 39        44_AG08312_mock_A__htseq-count_Gene_counts.tabular AG08312
    ## 40     44_PR00248_treated_A__htseq-count_Gene_counts.tabular PR00248
    ## 41       44_PR111_treated_A__htseq-count_Gene_counts.tabular   PR111
    ## 42       44_PR235_treated_A__htseq-count_Gene_counts.tabular   PR235
    ## 43        44_SQMA_treated_A__htseq-count_Gene_counts.tabular    SQMA
    ## 44        44_SQMB_treated_B__htseq-count_Gene_counts.tabular    SQMB
    ## 45     56_AG05311_treated_B__htseq-count_Gene_counts.tabular AG05311
    ## 46     56_AG08305_treated_B__htseq-count_Gene_counts.tabular AG08305
    ## 47        56_AG08308_mock_B__htseq-count_Gene_counts.tabular AG08308
    ## 48        56_AG08312_mock_B__htseq-count_Gene_counts.tabular AG08312
    ## 49     56_PR00248_treated_B__htseq-count_Gene_counts.tabular PR00248
    ## 50       56_PR111_treated_B__htseq-count_Gene_counts.tabular   PR111
    ## 51       56_PR235_treated_B__htseq-count_Gene_counts.tabular   PR235
    ## 52        56_SQMA_treated_B__htseq-count_Gene_counts.tabular    SQMA
    ## 53        56_SQMB_treated_C__htseq-count_Gene_counts.tabular    SQMB
    ## 54     68_AG05311_treated_C__htseq-count_Gene_counts.tabular AG05311
    ## 55     68_AG08305_treated_C__htseq-count_Gene_counts.tabular AG08305
    ## 56        68_AG08308_mock_C__htseq-count_Gene_counts.tabular AG08308
    ## 57        68_AG08312_mock_C__htseq-count_Gene_counts.tabular AG08312
    ## 58     68_PR00248_treated_C__htseq-count_Gene_counts.tabular PR00248
    ## 59       68_PR111_treated_C__htseq-count_Gene_counts.tabular   PR111
    ## 60       68_PR235_treated_C__htseq-count_Gene_counts.tabular   PR235
    ## 61        68_SQMA_treated_C__htseq-count_Gene_counts.tabular    SQMA
    ## 62        68_SQMB_treated_A__htseq-count_Gene_counts.tabular    SQMB
    ## 63              8_AF_mock_A__htseq-count_Gene_counts.tabular      AF
    ## 64              8_AF_mock_B__htseq-count_Gene_counts.tabular      AF
    ## 65              8_AF_mock_C__htseq-count_Gene_counts.tabular      AF
    ## 66           8_AF_treated_A__htseq-count_Gene_counts.tabular      AF
    ## 67           8_AF_treated_B__htseq-count_Gene_counts.tabular      AF
    ## 68           8_AF_treated_C__htseq-count_Gene_counts.tabular      AF
    ## 69         8_AG05311_mock_A__htseq-count_Gene_counts.tabular AG05311
    ## 70         8_AG06105_mock_A__htseq-count_Gene_counts.tabular AG06105
    ## 71         8_AG06105_mock_B__htseq-count_Gene_counts.tabular AG06105
    ## 72        8_AG06105_mock_C___htseq-count_Gene_counts.tabular AG06105
    ## 73      8_AG06105_treated_A__htseq-count_Gene_counts.tabular AG06105
    ## 74      8_AG06105_treated_B__htseq-count_Gene_counts.tabular AG06105
    ## 75      8_AG06105_treated_C__htseq-count_Gene_counts.tabular AG06105
    ## 76        8_AG08305_mock_A___htseq-count_Gene_counts.tabular AG08305
    ## 77      8_AG08308_treated_A__htseq-count_Gene_counts.tabular AG08308
    ## 78      8_AG08312_treated_A__htseq-count_Gene_counts.tabular AG08312
    ## 79  8_AG08490_mock_C_Sunday__htseq-count_Gene_counts.tabular AG08490
    ## 80           8_NHDF_mock_A___htseq-count_Gene_counts.tabular    NHDF
    ## 81             8_NHDF_mock_B_htseq-count_Gene_counts.tabular    NHDF
    ## 82           8_NHDF_mock_C___htseq-count_Gene_counts.tabular    NHDF
    ## 83        8_NHDF_treated_A___htseq-count_Gene_counts.tabular    NHDF
    ## 84         8_NHDF_treated_B__htseq-count_Gene_counts.tabular    NHDF
    ## 85        8_NHDF_treated_C___htseq-count_Gene_counts.tabular    NHDF
    ## 86         8_PR00033_mock_A__htseq-count_Gene_counts.tabular PR00033
    ## 87         8_PR00033_mock_B__htseq-count_Gene_counts.tabular PR00033
    ## 88        8_PR00033_mock_C___htseq-count_Gene_counts.tabular PR00033
    ## 89      8_PR00033_treated_A__htseq-count_Gene_counts.tabular PR00033
    ## 90      8_PR00033_treated_B__htseq-count_Gene_counts.tabular PR00033
    ## 91     8_PR00033_treated_C___htseq-count_Gene_counts.tabular PR00033
    ## 92         8_PR00036_mock_A__htseq-count_Gene_counts.tabular PR00036
    ## 93         8_PR00036_mock_B__htseq-count_Gene_counts.tabular PR00036
    ## 94         8_PR00036_mock_C__htseq-count_Gene_counts.tabular PR00036
    ## 95      8_PR00036_treated_A__htseq-count_Gene_counts.tabular PR00036
    ## 96      8_PR00036_treated_B__htseq-count_Gene_counts.tabular PR00036
    ## 97     8_PR00036_treated_C___htseq-count_Gene_counts.tabular PR00036
    ## 98         8_PR00039_mock_A__htseq-count_Gene_counts.tabular PR00039
    ## 99         8_PR00039_mock_B__htseq-count_Gene_counts.tabular PR00039
    ## 100        8_PR00039_mock_C__htseq-count_Gene_counts.tabular PR00039
    ## 101     8_PR00039_treated_A__htseq-count_Gene_counts.tabular PR00039
    ## 102     8_PR00039_treated_B__htseq-count_Gene_counts.tabular PR00039
    ## 103    8_PR00039_treated_C___htseq-count_Gene_counts.tabular PR00039
    ## 104        8_PR00054_mock_A__htseq-count_Gene_counts.tabular PR00054
    ## 105        8_PR00054_mock_B__htseq-count_Gene_counts.tabular PR00054
    ## 106        8_PR00054_mock_C__htseq-count_Gene_counts.tabular PR00054
    ## 107    8_PR00054_treated_A___htseq-count_Gene_counts.tabular PR00054
    ## 108     8_PR00054_treated_B__htseq-count_Gene_counts.tabular PR00054
    ## 109     8_PR00054_treated_C__htseq-count_Gene_counts.tabular PR00054
    ## 110        8_PR00107_mock_A__htseq-count_Gene_counts.tabular PR00107
    ## 111        8_PR00107_mock_B__htseq-count_Gene_counts.tabular PR00107
    ## 112        8_PR00107_mock_C__htseq-count_Gene_counts.tabular PR00107
    ## 113     8_PR00107_treated_A__htseq-count_Gene_counts.tabular PR00107
    ## 114     8_PR00107_treated_B__htseq-count_Gene_counts.tabular PR00107
    ## 115     8_PR00107_treated_C__htseq-count_Gene_counts.tabular PR00107
    ## 116        8_PR00248_mock_A__htseq-count_Gene_counts.tabular PR00248
    ## 117        8_PR00573_mock_A__htseq-count_Gene_counts.tabular PR00573
    ## 118        8_PR00573_mock_B__htseq-count_Gene_counts.tabular PR00573
    ## 119        8_PR00573_mock_C__htseq-count_Gene_counts.tabular PR00573
    ## 120     8_PR00573_treated_A__htseq-count_Gene_counts.tabular PR00573
    ## 121     8_PR00573_treated_B__htseq-count_Gene_counts.tabular PR00573
    ## 122    8_PR00573_treated_C___htseq-count_Gene_counts.tabular PR00573
    ## 123        8_PR01109_mock_A__htseq-count_Gene_counts.tabular PR01109
    ## 124        8_PR01109_mock_B__htseq-count_Gene_counts.tabular PR01109
    ## 125        8_PR01109_mock_C__htseq-count_Gene_counts.tabular PR01109
    ## 126     8_PR01109_treated_A__htseq-count_Gene_counts.tabular PR01109
    ## 127     8_PR01109_treated_B__htseq-count_Gene_counts.tabular PR01109
    ## 128     8_PR01109_treated_C__htseq-count_Gene_counts.tabular PR01109
    ## 129         8_PR0230_mock_A__htseq-count_Gene_counts.tabular  PR0230
    ## 130         8_PR0230_mock_B__htseq-count_Gene_counts.tabular  PR0230
    ## 131        8_PR0230_mock_C___htseq-count_Gene_counts.tabular  PR0230
    ## 132      8_PR0230_treated_A__htseq-count_Gene_counts.tabular  PR0230
    ## 133      8_PR0230_treated_B__htseq-count_Gene_counts.tabular  PR0230
    ## 134     8_PR0230_treated_C___htseq-count_Gene_counts.tabular  PR0230
    ## 135          8_PR111_mock_A__htseq-count_Gene_counts.tabular   PR111
    ## 136          8_PR235_mock_A__htseq-count_Gene_counts.tabular   PR235
    ## 137        8_S003611_mock_A__htseq-count_Gene_counts.tabular S003611
    ## 138        8_S003611_mock_B__htseq-count_Gene_counts.tabular S003611
    ## 139        8_S003611_mock_C__htseq-count_Gene_counts.tabular S003611
    ## 140    8_S003611_treated_A___htseq-count_Gene_counts.tabular S003611
    ## 141     8_S003611_treated_B__htseq-count_Gene_counts.tabular S003611
    ## 142     8_S003611_treated_C__htseq-count_Gene_counts.tabular S003611
    ## 143        8_S003649_mock_A__htseq-count_Gene_counts.tabular S003649
    ## 144        8_S003649_mock_B__htseq-count_Gene_counts.tabular S003649
    ## 145        8_S003649_mock_C__htseq-count_Gene_counts.tabular S003649
    ## 146     8_S003649_treated_A__htseq-count_Gene_counts.tabular S003649
    ## 147     8_S003649_treated_B__htseq-count_Gene_counts.tabular S003649
    ## 148    8_S003649_treated_C___htseq-count_Gene_counts.tabular S003649
    ## 149        8_S004933_mock_A__htseq-count_Gene_counts.tabular S004933
    ## 150        8_S004933_mock_B__htseq-count_Gene_counts.tabular S004933
    ## 151        8_S004933_mock_C__htseq-count_Gene_counts.tabular S004933
    ## 152     8_S004933_treated_A__htseq-count_Gene_counts.tabular S004933
    ## 153     8_S004933_treated_B__htseq-count_Gene_counts.tabular S004933
    ## 154    8_S004933_treated_C___htseq-count_Gene_counts.tabular S004933
    ## 155           8_SQMA_mock_B__htseq-count_Gene_counts.tabular    SQMA
    ## 156           8_SQMB_mock_A__htseq-count_Gene_counts.tabular    SQMB
    ## 157             8_SR_mock_A__htseq-count_Gene_counts.tabular      SR
    ## 158             8_SR_mock_B__htseq-count_Gene_counts.tabular      SR
    ## 159             8_SR_mock_C__htseq-count_Gene_counts.tabular      SR
    ## 160          8_SR_treated_A__htseq-count_Gene_counts.tabular      SR
    ## 161          8_SR_treated_B__htseq-count_Gene_counts.tabular      SR
    ## 162          8_SR_treated_C__htseq-count_Gene_counts.tabular      SR
    ##     treatment replicate
    ## 1        mock         A
    ## 2        mock         B
    ## 3        mock         C
    ## 4     treated         A
    ## 5     treated         B
    ## 6     treated         C
    ## 7        mock         A
    ## 8        mock         B
    ## 9     treated         A
    ## 10    treated         B
    ## 11    treated         C
    ## 12       mock         A
    ## 13       mock         B
    ## 14       mock         C
    ## 15    treated         A
    ## 16    treated         B
    ## 17    treated         C
    ## 18       mock         B
    ## 19       mock         B
    ## 20    treated         B
    ## 21    treated         B
    ## 22       mock         B
    ## 23       mock         B
    ## 24       mock         B
    ## 25       mock         A
    ## 26       mock         B
    ## 27       mock         C
    ## 28       mock         C
    ## 29    treated         C
    ## 30    treated         C
    ## 31       mock         C
    ## 32       mock         C
    ## 33       mock         C
    ## 34       mock         C
    ## 35       mock         C
    ## 36    treated         A
    ## 37    treated         A
    ## 38       mock         A
    ## 39       mock         A
    ## 40    treated         A
    ## 41    treated         A
    ## 42    treated         A
    ## 43    treated         A
    ## 44    treated         B
    ## 45    treated         B
    ## 46    treated         B
    ## 47       mock         B
    ## 48       mock         B
    ## 49    treated         B
    ## 50    treated         B
    ## 51    treated         B
    ## 52    treated         B
    ## 53    treated         C
    ## 54    treated         C
    ## 55    treated         C
    ## 56       mock         C
    ## 57       mock         C
    ## 58    treated         C
    ## 59    treated         C
    ## 60    treated         C
    ## 61    treated         C
    ## 62    treated         A
    ## 63       mock         A
    ## 64       mock         B
    ## 65       mock         C
    ## 66    treated         A
    ## 67    treated         B
    ## 68    treated         C
    ## 69       mock         A
    ## 70       mock         A
    ## 71       mock         B
    ## 72       mock         C
    ## 73    treated         A
    ## 74    treated         B
    ## 75    treated         C
    ## 76       mock         A
    ## 77    treated         A
    ## 78    treated         A
    ## 79       mock         C
    ## 80       mock         A
    ## 81       mock         B
    ## 82       mock         C
    ## 83    treated         A
    ## 84    treated         B
    ## 85    treated         C
    ## 86       mock         A
    ## 87       mock         B
    ## 88       mock         C
    ## 89    treated         A
    ## 90    treated         B
    ## 91    treated         C
    ## 92       mock         A
    ## 93       mock         B
    ## 94       mock         C
    ## 95    treated         A
    ## 96    treated         B
    ## 97    treated         C
    ## 98       mock         A
    ## 99       mock         B
    ## 100      mock         C
    ## 101   treated         A
    ## 102   treated         B
    ## 103   treated         C
    ## 104      mock         A
    ## 105      mock         B
    ## 106      mock         C
    ## 107   treated         A
    ## 108   treated         B
    ## 109   treated         C
    ## 110      mock         A
    ## 111      mock         B
    ## 112      mock         C
    ## 113   treated         A
    ## 114   treated         B
    ## 115   treated         C
    ## 116      mock         A
    ## 117      mock         A
    ## 118      mock         B
    ## 119      mock         C
    ## 120   treated         A
    ## 121   treated         B
    ## 122   treated         C
    ## 123      mock         A
    ## 124      mock         B
    ## 125      mock         C
    ## 126   treated         A
    ## 127   treated         B
    ## 128   treated         C
    ## 129      mock         A
    ## 130      mock         B
    ## 131      mock         C
    ## 132   treated         A
    ## 133   treated         B
    ## 134   treated         C
    ## 135      mock         A
    ## 136      mock         A
    ## 137      mock         A
    ## 138      mock         B
    ## 139      mock         C
    ## 140   treated         A
    ## 141   treated         B
    ## 142   treated         C
    ## 143      mock         A
    ## 144      mock         B
    ## 145      mock         C
    ## 146   treated         A
    ## 147   treated         B
    ## 148   treated         C
    ## 149      mock         A
    ## 150      mock         B
    ## 151      mock         C
    ## 152   treated         A
    ## 153   treated         B
    ## 154   treated         C
    ## 155      mock         B
    ## 156      mock         A
    ## 157      mock         A
    ## 158      mock         B
    ## 159      mock         C
    ## 160   treated         A
    ## 161   treated         B
    ## 162   treated         C

``` r
## Use DESeq2 to examine count data
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = data_dir, design = ~ donor + treatment)
dds
```

    ## class: DESeqDataSet 
    ## dim: 65217 162 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(65217): ENSG00000000003 ENSG00000000005 ...
    ##   ENSG00000281921 ENSG00000281922
    ## rowData names(0):
    ## colnames(162): AG07923_mock_A AG07923_mock_B ... SR_treated_B
    ##   SR_treated_C
    ## colData names(3): donor treatment replicate

``` r
##Determining size factor for each column
dds <- estimateSizeFactors(dds)
dds@colData
```

    ## DataFrame with 162 rows and 4 columns
    ##                      donor treatment replicate sizeFactor
    ##                   <factor>  <factor>  <factor>  <numeric>
    ## AG07923_mock_A     AG07923      mock         A  0.9465769
    ## AG07923_mock_B     AG07923      mock         B  0.8103831
    ## AG07923_mock_C     AG07923      mock         C  0.9291446
    ## AG07923_treated_A  AG07923   treated         A  0.9207055
    ## AG07923_treated_B  AG07923   treated         B  1.0009211
    ## ...                    ...       ...       ...        ...
    ## SR_mock_B               SR      mock         B  1.7642677
    ## SR_mock_C               SR      mock         C  1.2152078
    ## SR_treated_A            SR   treated         A  0.9954227
    ## SR_treated_B            SR   treated         B  1.0614094
    ## SR_treated_C            SR   treated         C  1.0801023

``` r
##To save time in running the rlogTransformation code, I saved this 
##from when it was first run. 
##rld <- rlogTransformation(dds, blind = TRUE)
##save(rld, file = "rld2_NHPs.Rdata")
rld_loaded <- get(load("rld2_NHPs.Rdata"))

distsRL <- dist(t(assay(rld_loaded)))
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
  ifelse(donor == "AG07923" | donor ==  "AG08490" | donor == "PR0058", "darkolivegreen3",
                       ifelse(donor == "PR00033" | donor ==  "PR00036" | donor == "PR00039", "darkgreen",
                              ifelse(donor == "AG05311" | donor == "SQMA" | donor == "SQMB", "purple",
                                     ifelse(donor == "AG06105" | donor == "PR00054" | donor == "PR01109", "deepskyblue", 
                                            ifelse(donor == "PR0230" | donor == "PR00573" | donor == "PR00107", "dodgerblue1", 
                                                   ifelse(donor == "PR111" | donor == "PR235" | donor == "PR00248", "dodgerblue2",
                                                          ifelse(donor == "S004933" | donor == "S003611" | donor == "S003649", "dodgerblue3",
                                                                 ifelse(donor == "AG08308" | donor == "AG08312" | donor == "AG08305", "forestgreen",
                                                                        ifelse(donor == "NHDF" | donor == "AF" | donor == "SR", "dodgerblue4", "grey")))))))))
}

colcolors <- unlist(lapply(sampleTable$donor, colcoloring))
      
treatcoloring <- function(treatment) { if (treatment=="mock") "red" else "orange" }
treatcolors <- unlist(lapply(sampleTable$treatment, treatcoloring))

clab  = cbind(treatcolors, colcolors)
colnames(clab) = c("", "")

png(file = file.path(output_dir, paste(Sys.Date(), "all species_distance clustering_labeled.png")), units = 'in', height = 16, width = 18, res = 300)
heatmap.3(mat, trace="none", keysize = 0.5, col = viridis(100), srtCol = 45, cexCol = 0.6, cexRow = 0.6, dendrogram = "column", density.info = "none", margin = c(6, 8),  Rowv = TRUE, Colv = TRUE, ColSideColors = clab, lhei = c(1.5, 11), ColSideColorsSize=3) 
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
png(file = file.path(output_dir, paste(Sys.Date(), "all NHPs_distance clustering_unlabeled.png")), units = 'in', height = 16, width = 18, res = 300)
heatmap.3(mat, trace="none", keysize = 0.8, col = viridis(100), dendrogram = "column", density.info = "none", margin = c(5, 5),  Rowv = TRUE, Colv = TRUE, ColSideColors = clab, lhei = c(1.5, 11), lwid = c(1.5,8), ColSideColorsSize=3.5, labCol = FALSE, labRow = FALSE) 
dev.off()    
```

    ## quartz_off_screen 
    ##                 2

PCA plot with donor and treatment indicated by color and shape, respectively.

``` r
a <- plotPCA(rld_loaded, intgroup = c("donor", "treatment", "replicate"), returnData = TRUE)
donors_PCA_colors <- unlist(lapply(a$donor, colcoloring))
names(donors_PCA_colors) <- a$donor

b <- plotPCA(rld_loaded, intgroup = c("donor", "treatment", "replicate")) + aes(colour = a$donor, shape = a$treatment) +
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

![](QC_of_NHP_counts_files/figure-markdown_github/PCA%20plot-1.png)

``` r
ggsave(file = file.path(output_dir, paste(Sys.Date(), "allNHPS_PCA.png")), plot = b, units = 'in', height = 10, width = 10, dpi = 300, device = "png")
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
