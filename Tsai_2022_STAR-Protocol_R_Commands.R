#### These commands are developed and tested with: 
####   R version 4.2.0 (2022-04-22 ucrt)
####   Platform: x86_64-w64-mingw32/x64 (64-bit)
####   Running under: Windows 10 x64 (build 19044)
#### All the R commands used in the STAR Protocol developed by Tsai et al 2022 are included.

### "Before you begin" section
## Step 2a. Get version 3.15 of Bioconductor
if(!require("BiocManager", quietly = TRUE)) 
install.packages("BiocManager") 
BiocManager::install(version = "3.15")

## Step 2b. Install the required packages 
BiocManager::install(c("ggplot2", "gplots", "gridExtra", "RColorBrewer", "rtracklayer", "GEOquery", "biomaRt"))

## Step 3. Install MFSelector from source code
download.file("http://microarray.ym.edu.tw:8080/tools/module/MFSelector/content/support/multicore/MFSelector_1.0.tar.gz", "MFSelector_1.0.tar.gz")
install.packages("MFSelector_1.0.tar.gz", repos = NULL, type = "source")

### "Step-by-step method details" section
## Step 1a. Download the processed data (GSE140914_RAW.tar) from GEO with GEOquery package
library(GEOquery)
getGEOSuppFiles("GSE140914", makeDirectory = F)
untar("GSE140914_RAW.tar")

## Step 1b. Extract files
list.files(pattern = "gz$") |> lapply(gunzip)

## Step 1c. Combine these data
# Import data 
data_list <- list.files(pattern="^GSM") |> lapply(read.delim, header = TRUE)
data_mat <- sapply(data_list, \(x) x[, 3]) 

# Get gene identifiers
gene_id <- data_list[[1]][, 1] 
gene_names <- data_list[[1]][, 2] 
gene_labels <- paste(gene_id, gene_names, sep = ":") 

# Get sample names, batch information, group information 
sample_names <- sapply(data_list, \(x) colnames(x)[3]) 
sample_names <- gsub("X", "", sample_names) 
batch_id <- gsub("\\.\\S+", "", sample_names) 
sample_group <- gsub("\\S+\\.", "", sample_names) 
sample_code <- paste0("B:", batch_id, "_D:", sample_group) 
colnames(data_mat) <- sample_code 
rownames(data_mat) <- gene_id 

# Generate the input data for MFSelector 
data <- cbind(gene_labels, as.data.frame(data_mat))

## Step 1d. Visualize the overall relationship between samples with a multidimensional scaling (MDS) plot
# Define a function to generate MDS plot 
get_MDS_plot <- function(prefix, data_mat, sample_group){ 
    cmd <- t(data_mat) |> dist() |> cmdscale()
   
    df <- data.frame(Sample = colnames(data_mat), Dimension_1 = cmd[, 1],
      Dimension_2 = cmd[, 2], Group = sample_group) 

    p1 <- ggplot(df, aes(Dimension_1, Dimension_2)) +
    geom_point(aes(colour = Group), size = 4) + 
    scale_colour_brewer(palette="Paired") + 
    ggtitle(paste(prefix, "MDS Plot")) 

    print(p1) 
} 

# Load required package 
library(ggplot2) 

# Generate a MDS plot for this data set
get_MDS_plot("GSE140914", data_mat, sample_group)
ggsave("GSE140914_MDSplot.tiff", dpi = "retina") # Figure 1

## Step 2. Prepare the inputs with R. 
# Load required package
library(MFSelector)
library(parallel)

# Prepare input arguments
nsc <- table(sample_group)
stageord <- order(sample_group)
stagename <- sample_group |> unique() |> sort()

# Detect the number of total cores in the current machine
total_cores <- detectCores()
half_cores <- total_cores/2

## Problem 2 & 3 block ##
## Run MFSelector with multiple cores in parallel on a Windows machine
BiocManager::install(c("foreach", "doSNOW"))
library(foreach)
library(doSNOW)
detach("package:MFSelector", unload=TRUE) # To ensure the original MFSelector has been detached
source("MFSelector_doSNOW.r") # https://github.com/yushuen/MFSelector_STAR-Protocol/blob/main/MFSelector_doSNOW.r
## End of Problem 2 & 3 block ##

## Step 3. Identify candidate descending monotonic key genes with MFSelecto
mfselector(data, nsc, stageord = stageord, stagename = stagename, 
  type = 1, nline = T, dline = T, pdf = 1:100, cmp = 0, permut = 100, svdenoise = 0.03, svdetimes = 4, cores = half_cores) 

## Step 4. Rename the output files
mf_outputs <- list.files(pattern="^mfselector")
mf_outputs_type1 <- gsub(".*\\.", "MFSelector_Type1.", mf_outputs)
mapply(file.rename, mf_outputs, mf_outputs_type1)
rm(mf_outputs, mf_outputs_type1)

## Step 5. Identify candidate ascending monotonic key genes
mfselector(data, nsc, stageord = stageord, stagename = stagename, 
  type = 2, nline = T, dline = T, pdf = 1:100, cmp = 0, permut = 100, svdenoise = 0.03,  svdetimes = 4, cores = half_cores)

## Step 6. Rename the second output files
mf_outputs <- list.files(pattern="^mfselector")
mf_outputs_type2 <- gsub(".*\\.", "MFSelector_Type2.", mf_outputs)
mapply(file.rename, mf_outputs, mf_outputs_type2)
rm(mf_outputs, mf_outputs_type2)

## Step 7. Define an R function to identify genes which have fulfilled the given criteria
candidate.genes.parser <- function(input_file, DE = NULL, SVDE = NULL, nline = NULL, out_col = NULL, ...){
    input_tab <- read.delim(input_file)
    rii <- 1:nrow(input_tab)
    if(!is.null(DE)){
        sii <- which(as.numeric(input_tab[, "DE"]) <= DE)
        rii <- intersect(rii, sii)
    }
    if(!is.null(SVDE)){
        sii <- which(as.numeric(input_tab[, "SVDE"]) <= SVDE)
        rii <- intersect(rii, sii)
    }
    if(!is.null(nline)){
        sii <- which(input_tab[, "with.N.1.distinct.lines"] == TRUE)
        rii <- intersect(rii, sii)
    }
    if(is.null(out_col)){
        candidate_genes <- as.character(input_tab[rii, 1])  
    }else{
        candidate_genes <- as.character(input_tab[rii, out_col]) 
    }    
    return(candidate_genes)
}

## Step 8. Get candidate monotonic key genes
MF_type1_genes <- candidate.genes.parser("MFSelector_Type1.txt", DE = 4)
MF_type2_genes <- candidate.genes.parser("MFSelector_Type2.txt", DE = 4)

## Step 9. Find the indices of candidate genes
mii_type1 <- match(MF_type1_genes, gene_labels)
mii_type2 <- match(MF_type2_genes, gene_labels)

## Step 10. Generate MDS plots with these candidate genes
get_MDS_plot("GSE140914_Type1", data_mat[mii_type1, ], sample_group)
ggsave("GSE140914_Type1_MDSplot.tiff", dpi = "retina") # Figure 2

get_MDS_plot("GSE140914_Type2", data_mat[mii_type2, ], sample_group)
get_MDS_plot("GSE140914_Type1+2", data_mat[c(mii_type1, mii_type2), ], sample_group)

## Step 11. Generate heatmap plots with these candidate genes
# Define a function for Z-score transformation
z_transformation <- function(x){(x-mean(x))/sd(x)}
# Define a function to perform hierarchical clustering with Ward's method
ward_hclust <- function(d){hclust(d, method = "ward.D")}

# Define a function to generate a heatmap
get_heatmap<-function(prefix, data, sample_group, gene_symbols = NULL, ...){
    sample_col_fac <- as.factor(sample_group)    
    n_group <- sample_group |> unique() |> length()                    
    dataZ <- apply(data, 1, z_transformation) |> t()

    my_col_palette <- brewer.pal(11, "RdBu") |> rev()
    levels(sample_col_fac) <- brewer.pal(n_group, "Paired")
    col_colors <- as.character(sample_col_fac)
        
    if(!is.null(gene_symbols)){
        rownames(dataZ) <- gene_symbols
    }
          
    heatmap.2(dataZ, col = my_col_palette, ColSideColors = col_colors, 
      Colv = TRUE, dendrogram = "both", hclustfun = ward_hclust, 
      margins = c(5, 7), trace = "none", keysize = 1.2, 
      lhei = c(1.5, 9.5), lwid = c(1.5, 6.5))                          
}

# Load required packages for generating heatmap
library(gplots)
library(RColorBrewer)

# Generate heatmap plots with these candidate genes
# It needs a larger plot device.
tiff("GSE140914_Type1_Heatmap.tiff", width=10, height=12, unit="in", res = 320) # Figure 3
get_heatmap("GSE140914_Type1", data_mat[mii_type1, ], sample_group, gene_symbols = gene_labels[mii_type1])
dev.off() 

get_heatmap("GSE140914_Type2", data_mat[mii_type2, ], sample_group, gene_symbols = gene_labels[mii_type2])
get_heatmap("GSE140914_Type1+2", data_mat[c(mii_type1, mii_type2), ], sample_group, gene_symbols = gene_labels[c(mii_type1, mii_type2)])

## Step 12. Generate scatter plots with these candidate genes
# Import the output text files of MFSelector
mf_tab_1 <- read.delim("MFSelector_Type1.txt", header = TRUE)
mf_tab_2 <- read.delim("MFSelector_Type2.txt", header = TRUE)

# Define a function to generate scatter plot with ggplot2
get_scatter_plots <- function(my_gene_id, data, gene_id, sample_group, mf_tab){
    mii1 <- which(gene_id == my_gene_id)
    mii2 <- which(mf_tab[, 1] == my_gene_id)

    if(length(mii1) > 0 & length(mii2) >0){
        DE <- mf_tab[mii2, 2]
        PVAL <- round(mf_tab[mii2, 3], 2)
        QVAL <- round(mf_tab[mii2, 4], 2)
        SVDE <- mf_tab[mii2, 5]
        TITLE <- paste(my_gene_id, "\nDE =", DE, " p-value =", PVAL, 
                 " q-value =", QVAL, " SVDE =", SVDE)

        expression <- as.numeric(data[mii1, ])
        ordered_idx <- order(sample_group)        

        df <- data.frame(Samples = 1:length(expression), 
              Expression = expression[ordered_idx], 
              Group = as.factor(sample_group[ordered_idx]))

        group_min_df <- data.frame(
                        Min = tapply(df$Expression, df$Group, min), 
                        Group = levels(df$Group))

        p <- ggplot(data = df, aes(x = Samples, y = Expression)) + 
             geom_point(aes(colour = Group, shape = Group), size = 4) + 
             scale_colour_brewer(palette="Paired") +
             scale_x_continuous(breaks = 1:length(expression), 
             	labels = colnames(data)[ordered_idx]) + 
             geom_hline(aes(yintercept = Min, colour = Group), 
             group_min_df, lty = "dashed") + ggtitle(TITLE) + 
             theme(plot.title = element_text(hjust = 0.5)) + 
             theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

        print(p)        
    }
}

# Generate scatter plots with these candidate genes
get_scatter_plots(MF_type1_genes[1], data_mat, gene_labels, sample_group, mf_tab_1)
ggsave("Figure_4.tiff", dpi = "retina") # Figure 4

# If “parallel” package is available in your system, you can use “mclapply()” instead of “lapply()”
pdf("Scatter-Plots_Type1.pdf")
lapply(MF_type1_genes, get_scatter_plots, data_mat, gene_labels, sample_group, mf_tab_1)
#mclapply(MF_type1_genes, get_scatter_plots, data_mat, gene_labels, sample_group, mf_tab_1, cores = half_cores)
dev.off()
pdf("Scatter-Plots_Type2.pdf")
lapply(MF_type2_genes, get_scatter_plots, data_mat, gene_labels, sample_group, mf_tab_2)
#mclapply(MF_type2_genes, get_scatter_plots, data_mat, gene_labels, sample_group, mf_tab_2, cores = half_cores)
dev.off()

## Step 13. Prepare the input files for TO-GCN
# Compute the mean CPM values for each group
sample_group_fac <- as.factor(sample_group)

get_group_mean <- function(x, sample_group_fac){
    tapply(x, sample_group_fac, mean)
}

group_mean_cpm <- apply(data_mat, 1, get_group_mean, sample_group_fac) |> t()

# Define a function to set negative values as zero
neg_to_zero <- function(x){ (abs(x)+x)/2 }

# Filter out lowly expressed genes (CPM <= 1)
check_mat <- neg_to_zero(group_mean_cpm - 1)
mode(check_mat) <- "logical"
check_ans <- apply(check_mat, 1, any)
group_mean_cpm_subset <- group_mean_cpm[check_ans, ]

## Problem 4 block ## 
library(biomaRt) # Load package
ensembl <- useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', host = "uswest.ensembl.org") # Select the database and create the ensembl object
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl) # Select the data set and update the ensembl object
mart_export <- getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), filters = 'ensembl_gene_id', values = gene_id, mart = ensembl) # Get gene biotype
write.table(mart_export, "gene_biotype.txt", row.names = F, col.names = T, sep = "\t", quote = F) # Generate an output table
## End of Problem 4 block ##

# Import Ensembl gene biotype annotation
gene_biotype <- read.delim("gene_biotype.txt", header = T)

# Check whether there are any gene ids not included in gene_biotype
bii <- match(rownames(group_mean_cpm_subset), gene_biotype$ensembl_gene_id, nomatch = 0)
nii <- which(bii!=0)

# Select only protein coding genes
feature_type <- gene_biotype[bii[nii], 2]
pii <- which(feature_type == "protein_coding")
group_mean_cpm_pro <- group_mean_cpm_subset[nii[pii], ]

# Download the full list of human transcript factors form AnimalTFDB3.0
download.file("http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF", "Homo_sapiens_TF.txt")

# Import full list of human transcript factors
human_tf <- read.delim("Homo_sapiens_TF.txt", header = T)

# Identify the indices of genes annotated as transcription factors
tii <- match(rownames(group_mean_cpm_pro), human_tf$Ensembl, nomatch = 0)
tf_vec <- as.logical(tii)

# Split the mean TPM table into two tables for 
# transcription factor coding genes and non-transcription factor coding genes, respectively.
tf_group_mean_cpm <- group_mean_cpm_pro[tf_vec, ]
non_tf_group_mean_cpm <- group_mean_cpm_pro[!tf_vec, ]

# Generate the input files for TO-GCN
write.table(tf_group_mean_cpm, "TF_gene_matrix.tsv", row.names = T, col.names = F, sep = "\t", quote = F)
write.table(non_tf_group_mean_cpm, "Non-TF_gene_matrix.tsv", row.names = T, col.names = F, sep = "\t", quote = F)
