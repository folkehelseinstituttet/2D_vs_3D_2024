#Supplementary material for: https://pubmed.ncbi.nlm.nih.gov/39551125/
'Transcriptomic characterization of 2D and 3D human induced pluripotent stem cell-based in vitro models as New Approach Methodologies for developmental neurotoxicity testing' 

#Abstract
The safety and developmental neurotoxicity (DNT) potential of chemicals remain critically understudied due to limitations of current in vivo testing guidelines, which are low throughput, resource-intensive, and hindered by species differences that limit their relevance to human health. To address these issues, robust New Approach Methodologies (NAMs) using deeply characterized cell models are essential. This study presents the comprehensive transcriptomic characterization of two advanced human-induced pluripotent stem cell (hiPSC)-derived models: a 2D adherent and a 3D neurosphere model of human neural progenitor cells (hiNPCs) differentiated up to 21 days. Using high-throughput RNA sequencing, we compared gene expression profiles of 2D and 3D models at three developmental stages (3, 14, and 21 days of differentiation). Both models exhibit maturation towards post-mitotic neurons, with the 3D model maturing faster and showing a higher prevalence of GABAergic neurons, while the 2D model is enriched with glutamatergic neurons. Both models demonstrate broad applicability domains, including excitatory and inhibitory neurons, astrocytes, and key endocrine and especially the understudied cholinergic receptors. Comparison with human fetal brain samples confirms their physiological relevance. This study provides novel in-depth applicability insights into the temporal and dimensional aspects of hiPSC-derived neural models for DNT testing. The complementary use of these two models is highlighted: the 2D model excels in synaptogenesis assessment, while the 3D model is particularly suited for neural network formation as observed as well in previous functional studies with these models. This research marks a significant advancement in developing human-relevant, high-throughput DNT assays for regulatory purposes.

#Keywords: Applicability domain; DNT IVB; Developing human brain; Hazard assessment; Human neural progenitor cells.


Containes python scripts and results, version in Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11397464.svg)](https://zenodo.org/records/11397464)

# Simplified pipeline
![Figure 1](https://github.com/bazyliszek/2D_vs_3D_2024.jpg)

# Summary
[script.r](https://github.com/bazyliszek/2D_vs_3D_2024/script.r) Main R script  



# 2D_vs_3D_2024
2D_vs_3D_2024

# https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html

#1375 check me

# Setting working directory and libraries  --------------------------------

setwd("C:/Users/maww/Documents/2D_vs_3D")

#Vectors needed

ribosomal <- "ribosomal_not_removed"
threshold <-  10 #thereshold for gene counts
compared_to <- "2D_Day3" #contrast are made 14_0, 21_0  #Day 3 and dose 0 etc (...)
contrast = "2D_3"

#setwd("C:/Users/maww/Downloads")
data_dir <- "C:/Users/maww/Documents/2D_vs_3D"
saving_dir <- paste0("C:/Users/maww/Documents/2D_vs_3D/results/", ribosomal, "/Threshold_", threshold, "/Compared_to_", compared_to, "/logFC2_final/", sep="")
#saving_dir <- paste0("C:/Users/maww/Downloads/meth_results/", ribosomal, "/Threshold_", threshold, "/Compared_to_", compared_to, "/logFC1", sep="")

#saving_dir <- "C:/Users/maww/Downloads/meth_results/treatment_without_ribosomal_RNA_filtered_10/0_14"
#saving_dir <- "C:/Users/maww/Downloads/meth_results/treatment_without_ribosomal_RNA_filtered_10/0_21"
setwd(saving_dir)
setwd(data_dir)

#load library
library(DESeq2)
library(readxl)
library(ggrepel)
library(pcaExplorer)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limmaquickpca2go) ##missing
library(pca2go) ##missing                                                       
library(limma)
library(BiocStyle)

# Loading the dds object --------------------------------------------------
load("deseq2.dds.RData")
dds_meth <- dds
rm(dds)
dds_meth
#pcaExplorer(dds_meth)

colData(dds_meth)
summary(dds_meth$sizeFactor)


# Adding metadata to the dds object ---------------------------------------
samples_all_meth <- read_xlsx("220704_samples_run1_2_met_Marcin_Lea_2D_vs_3D.xlsx")  #manually curated so that it fits as.data.frame(colData(dds_meth)) 
head(samples_all_meth)

subsamples_mymeth <- samples_all_meth[,c(1,8:13)] #was 12
subsamples_mymeth$Treatment <-as.factor(c(paste(subsamples_mymeth$Time, "_", subsamples_mymeth$Dose, sep = "")))

dds_meth$Dose <- as.factor(subsamples_mymeth$Dose)
dds_meth$Time <- as.factor(subsamples_mymeth$Time)
dds_meth$Replicate <- as.factor(subsamples_mymeth$Replicate)
dds_meth$Run <- as.factor(subsamples_mymeth$Run)
dds_meth$Dimension <- as.factor(subsamples_mymeth$Dimension)
dds_meth$Passage <- as.factor(subsamples_mymeth$Pasage)

dds_meth$Treatment <- as.factor(c(paste(dds_meth$Dimension, "_", dds_meth$Time, sep = ""))) #if only one contast is interesting 
dds_meth$Treatment


# Relevel according to the contrasts --------------------------------------
dds_meth$Treatment <- relevel(dds_meth$Treatment, ref = "2D_3")  #change contrasts
dds_meth$Dimension <- relevel(dds_meth$Dimension, ref = "2D")
dds_meth$Time <- relevel(dds_meth$Time, ref = "3")
dds_meth$Dose <- relevel(dds_meth$Dose, ref = "0")
dds_meth$Run <- relevel(dds_meth$Run, ref = "1")
dds_meth$Passage <- relevel(dds_meth$Passage, ref = "9")

# Choosing the design (s) ------------------------------------------------------

dds_meth@design <- ~Treatment

#### NB!!
# Load the list of interesting genes (move this down!!!!!) -------------------------------------------

#genes_lists <- read_xlsx("RNASeq_Genes_Analysis_Marcin_list.xlsx") 
#head(genes_lists)
genes_lists <- read_xlsx("RNASeq_Genes_Analysis_05102022.xlsx") 
dim(genes_lists)


# NB! Change location for saving the files -----------------------------------------
setwd(saving_dir)

pdf("SizeFactor.pdf")
hist(dds_meth$sizeFactor)
dev.off()

png("SizeFactor.png")
hist(dds_meth$sizeFactor)
dev.off()

colData(dds_meth)
dim(dds_meth) #29744 33 #was 72

# Remove low expressed genes ----------------------------------------------
threshold <-  10
dds_meth <- dds_meth[rowSums(counts(dds_meth)) > threshold, ] 
dim(dds_meth) #19626 33


# Running DEseq2 and saving the files -------------------------------------

dds_meth_res <- DESeq(dds_meth) #!re-run if needed
dim(dds_meth_res)
dds_meth_res <- dds_meth_res[which(mcols(dds_meth_res)$betaConv),] # NB- consider to implement if needed.
dim(dds_meth_res) #19621
saveRDS(dds_meth_res, file="with_ribo_dds_meth_res_design_treatment_filtered_2D_3.rds")
dds_meth_res <- readRDS(file="with_ribo_dds_meth_res_design_treatment_filtered_2D_3.rds") #2D_3

#original counts
head(counts(dds_meth_res))
dim(counts(dds_meth_res)) #19621
names_genes <- as.vector(rownames(counts(dds_meth_res)))
head(names_genes, 20)
length(names_genes) #19621
write.csv(names_genes, file="names_genes.csv")

library("gprofiler2")

newnames = gconvert(names_genes, organism = "hsapiens", target = "ENSG")
dim(newnames) #16796
head(newnames) #16796)
str(unique(newnames))

merged_genes <- merge() 


dataframe <- cbind(my_input=as.vector(rownames(counts(dds_meth_res))), as.data.frame(counts(dds_meth_res, normalized=T), row.names=F))
head(dataframe)

merged <- unique(merge(x=dataframe, y=newnames, by.x='my_input', 
                       by.y='input', rownames=F, all.x=F, all.y=F, no.dups=T))
dim(merged) #16796
ensg <- merged$target
finalframe <- as.data.frame(cbind(ensg, merged[ ,2:34]), rownames=F)
head(finalframe)

any(duplicated(finalframe$ensg))
sum(duplicated(finalframe$ensg))
length(which(duplicated(finalframe$ensg)))

table(duplicated(finalframe$ensg))


duplicated(finalframe$ensg) ##

finalframe[876:880, ]

finalframe2 <- finalframe[duplicated(finalframe$ensg)==F, ]
dim(finalframe2) #16678

getwd()
write.csv(finalframe2, file="deconvolution11.csv", row.names=FALSE)



ttps://biit.cs.ut.ee/gprofiler/convert

dim(newnames)

gene_annotation <- read.csv2("gProfiler_hsapiens_14.8.2023_21-05-31.csv", sep=",")
dim(gene_annotation) #19981

test <- cbind(gene_annotation$converted_alias, counts(dds_meth_res))
head(test)


summary(dds_meth_res) #19621


# Transformations of counts, PCA plots -------------------------------------------

rld <- rlog(dds_meth_res, blind=TRUE) #!re-run if needed
vst <- vst(dds_meth_res, blind=TRUE) #!re-run if needed
saveRDS(rld, file="without_ribo_rld_10_design_treatment_2D_3.rds")
saveRDS(vst, file="without_ribo_vst_10_design_treatment_2D_3.rds")

rld <- readRDS(file="without_ribo_rld_10_design_treatment_2D_3.rds") #
vst <- readRDS(file="without_ribo_vst_10_design_treatment_2D_3.rds") #

rld_values <- assay(rld) # values
write.csv(rld_values, file="rld_values_2D_3.csv")

vst_values <- assay(vst) # values
write.csv(vst_values, file="vst_values_2D_3.csv")


# Plot PCA for rld and vst 500, 1000, 2000, 3000, 19621 top genes --------------------

ntop_genes = 19621
pdf(paste0("PCA_rld_meth_top_", ntop_genes, "_genes_with_threshold_", threshold, "_contrast_", contrast, "_design_treatment_plus_time.pdf"))
plotPCA(rld, intgroup="Treatment", ntop=ntop_genes)
plotPCA(rld, intgroup="Dimension", ntop=ntop_genes)
plotPCA(rld, intgroup="Time", ntop=ntop_genes)
plotPCA(rld, intgroup=c("Time", "Dimension"), ntop=ntop_genes)
dev.off()

pdf(paste0("PCA_vst_meth_top_", ntop_genes, "_genes_with_threshold_", threshold, "_contrast_", contrast, "_design_treatment_plus_time.pdf"))
plotPCA(vst, intgroup="Dimension", ntop=ntop_genes)
plotPCA(vst, intgroup="Treatment", ntop=ntop_genes)
plotPCA(vst, intgroup="Time", ntop=ntop_genes)
plotPCA(vst, intgroup=c("Time", "Dimension"), ntop=ntop_genes)
dev.off()


# Making PCA based on rld object again (visualize different PCs) ----------
# Input is a matrix of log transformed values
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
#https://hbctraining.github.io/DGE_workshop/lessons/03_DGE_QC_analysis.html
df <- cbind(subsamples_mymeth, pca$x)

pdf(paste0("PCA_additional_meth_top_", "all", "_genes_with_threshold_", threshold, "_contrast_", contrast, "_design_treatment.pdf"))
ggplot(df) + geom_point(aes(x=PC1, y=PC2, color = as.factor(Dimension))) + theme_bw() + ggtitle("Dimension - all features", contrast)
ggplot(df) + geom_point(aes(x=PC1, y=PC2, color = as.factor(Time))) + theme_bw() + ggtitle("Time - all features", contrast)
ggplot(df) + geom_point(aes(x=PC1, y=PC2, color = as.factor(Dimension), shape=as.factor(Time))) + theme_bw() + ggtitle("Dimension - all features", contrast)
dev.off()

#Vizualization for                                                        
pdf(paste0("PCA_additional_meth_top_", "all_in_one", "_genes_with_threshold_", threshold, "_contrast_", contrast, "_design_treatment.pdf"))
ggplot(df) + geom_point(aes(x=PC1, y=PC2, color = as.factor(Dimension), shape=as.factor(Time))) + theme_bw() + ggtitle("Dimension - all features", contrast)
dev.off()


# Hierachical clustering --------------------------------------------------
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

### Plot heatmap
library(pheatmap)
library(RColorBrewer)

heat.colors <- brewer.pal(6, "Reds")

pdf(paste0("Heatmap_all", "_genes_with_threshold_", threshold, "_contrast_", contrast, "_design_treatment_plus_pasage.pdf"))
#png(paste0("Heatmap_all", "_genes_with_threshold", threshold, "_contrast_", contrast, "_design_treatment_plus_pasage.png"), width = 20, height = 20, units = "cm", res=600, pointsize = 9)
#jpeg(paste0("Heatmap_all", "_genes_with_threshold", threshold, "_contrast_", contrast, "_design_treatment_plus_pasage.jpg"), width = 20, height = 20, units = "cm", res=600, pointsize = 9)
pheatmap(rld_cor, color= heat.colors, border_color=NA, fontsize = 6, height = 2, main=paste0("Rld - Heatmap all for contrast", contrast ))
dev.off()

# Print all the comparisons made ------------------------------------------
resultsNames(dds_meth_res)

#only for the treatment
contrasts_to_check = resultsNames(dds_meth_res)
contrasts_to_check = contrasts_to_check[-1]
#####contrasts_to_check = contrasts_to_check[1:14]

results_to_check = c()
for (i in contrasts_to_check) {
  results_to_check = c(results_to_check, results(dds_meth_res, name=i))
}
results_to_check
names(results_to_check) <- resultsNames(dds_meth_res)[-1]


# Fint the relevant genes based on the baseMean, padj and abs(log2) -----------------------------------------
baseMean_cutoff <- 0
padj_cutoff <- 0.01  #make 0.01 and 0.05 and 0.1
log2FoldChange_cutoff <- 2 # or 1

relevant_genes <- c()

for (i in results_to_check){
  tmp =  rownames(i[!is.na(i$padj) & (i$baseMean>baseMean_cutoff) & 
                      (i$padj<padj_cutoff) & (abs(i$log2FoldChange)>log2FoldChange_cutoff), ])
  print(length(tmp))
  relevant_genes <- union(relevant_genes, tmp)
}
#relevant_genes
length(relevant_genes)  #7979

#padj 0.01:
#1.5 = 1153 
#2.0 = 1560
#3.0 = 6247   
#4.0 = 6617
#5.0 = 5283

# Write down to CSV file the genes with given p-values  ------------------------------------------

filtered_results <- list()
#for (index in 1:length(results_to_check)) {
for (index in names(results_to_check)) {
  newElement = results_to_check[[index]][c(relevant_genes), ]
  filtered_results[index] <- newElement
}
head(filtered_results)

printableResults = data.frame()
printableResults <- cbind(c(relevant_genes, printableResults))

for (index in names(filtered_results)) {
  printableResults <- cbind(printableResults, filtered_results[[index]][["baseMean"]])
  printableResults <- cbind(printableResults, filtered_results[[index]][["log2FoldChange"]])
  printableResults <- cbind(printableResults, filtered_results[[index]][["lfcSE"]])
  printableResults <- cbind(printableResults, filtered_results[[index]][["stat"]])
  printableResults <- cbind(printableResults, filtered_results[[index]][["pvalue"]])
  printableResults <- cbind(printableResults, filtered_results[[index]][["padj"]])
  
}

newColNames = c("Genes")
for (index in names(filtered_results)){
  newColNames = c(newColNames, paste0(index,"_baseMean"), paste0(index,"_log2FoldChange"), paste0(index,"_lfcSE"), paste0(index,"_stat"), paste0(index,"_pvalue"), paste0(index,"_padj"))
  
}

colnames(printableResults) <- newColNames

head(printableResults)
dim(printableResults)


#double check why this is not working:
write.table(printableResults, file=paste0("A_Significant_genes_with_baseMean_cutoff_", baseMean_cutoff,
                                          "_padj_cutoff_", padj_cutoff,
                                          "_log2FoldChange_cutoff_", log2FoldChange_cutoff,
                                          "_with_threshold_", threshold,
                                          "_in_contrast_",contrast,"_design_treatment.csv"), row.names = F, col.names = T, sep = ";")
#working 
write.table(printableResults, file=paste0("A_Significant_genes_with_baseMean_cutoff_", baseMean_cutoff, "_design_treatment.csv"), 
            row.names = F, col.names = T, sep = ";")



###results_to_check[[1]]["GDF1",]

#Print lists of significant genes in any of the contrasts 
length(relevant_genes)

write.table(relevant_genes, file=paste0("B_Significant_genes_names_with_baseMean_cutoff_", baseMean_cutoff,
                                        "_padj_cutoff_", padj_cutoff, 
                                        "_log2FoldChange_cutoff_", log2FoldChange_cutoff,
                                        "_with_threshold_", threshold, 
                                        "_in_contrast_", contrast, "_design_treatment.csv"), row.names = F, col.names = F, sep = ";")


# Make a heatmaps ---------------------------------------------------------
library(ComplexHeatmap)
mat <- counts(dds_meth_res, normalized=T)[relevant_genes, ]
mat.z <- t(apply(mat, 1, scale))  #head(t(scale(t(mat)))) scaling on the rows
head(mat.z)
colnames(mat.z) <-  colnames(mat)


dimension_colors <- unlist(lapply(colnames(mat.z),function(x){
  if(grepl('2D',x)) '#fff7c0'
  else if(grepl('3D',x)) '#d1ffc0'
}
))

dimension_names <- unlist(lapply(colnames(mat.z),function(x){
  if(grepl('2D',x)) '2D'
  else if(grepl('3D',x)) '3D'
}
))


time_colors <- unlist(lapply(colnames(mat.z),function(x){
  if(grepl('DIV3',x)) '#fff7c0'
  else if(grepl('DIV14',x)) '#d1ffc0'
  else if(grepl('DIV21',x)) '#c0f0ff'
}
))

time_names <- unlist(lapply(colnames(mat.z),function(x){
  if(grepl('DIV3',x)) 'DIV3'
  else if(grepl('DIV14',x)) 'DIV14'
  else if(grepl('DIV21',x)) 'DIV21'
}
))


#Something is wrong with colours here!
#https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html


my_colors <- list(Time = time_colors, 
                  Dimension = dimension_colors)

set.seed(1234)


h_test <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, 
                  column_labels=colnames(mat.z), name="Z-score", 
                  row_names_gp = gpar(fontsize = 5), 
                  top_annotation = HeatmapAnnotation(Time = time_names, 
                                                     Dimension= dimension_names),
                  bottom_annotation = HeatmapAnnotation(Dimension = dimension_names),
                  column_names_gp = gpar(fontsize = 7)) #, col=treatment_colors


colors = SingleAnnotation(value = 1:10)


# Printing the colors into the heatmap  ----------------------------------

my_colour_time = c("DIV3" = "#97f7c8", "DIV14" = "#42f5a7","DIV21" = "#1e784c")
my_colour_dimension = c("2D" = "#4287f5", "3D" = "#f59942")

#my_col_dose = c("0" = "#ffdfff", "0.07" = "#ffc0e0","2.2" = "#ff90b0",  "6.7" = "#ff5f77", "20" = "#ff0000")
#my_col_passage = c("9" = "#ffffff", "10" = "#e0ffe0","11" = "#b0ffb0",  "13" = "#77ff77")
#my_colour_time = c("DIV3" = "#00ff44", "DIV14" = "#00dd33","DIV21" = "#009922")
#my_col_passage = c("9" = "#BED0F0", "10" = "#91B4F1","11" = "#538CF0",  "13" = "#1564EF")

#col_Zscore = colorRamp2(c(-5, 0, 5), c("#F5B7B1", "#FCF3CF", "#3498DB"))

hh <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T,
              column_title =paste0("base_cf", baseMean_cutoff,
                     "_padj_cf_", padj_cutoff, 
                     "_log2Fold_cf_", log2FoldChange_cutoff,
                     "_thr_", threshold, 
                     "_in_contr_", contrast),
              #col =col_Zscore, 
              column_labels=colnames(mat.z), name="Z-score", 
              row_names_gp = gpar(fontsize = 3), 
              top_annotation = HeatmapAnnotation(Time = time_names, 
                                                 Dimension = dimension_names,
                                                 col = list(Time = my_colour_time, 
                                                            Dimension = my_colour_dimension)),
              #bottom_annotation = HeatmapAnnotation(Dose = dose_names, col = list(Dose = my_col_dose)),
              column_names_gp = gpar(fontsize = 3)) #, col=treatment_colors

#Print png and pdf (working nicely)

png(file=paste0("Heatmap_of_sig_genes_with_baseMean_cutoff_", baseMean_cutoff,
                "_padj_cf_", padj_cutoff, 
                "_log2FoldChange_cf_", log2FoldChange_cutoff,
                "_with_threshold_", threshold, 
                "_in_contrast_", contrast, "_des_treat.png")) #, res=400, width= 3000, height=3000)
plot(hh) 
dev.off()


png(file=paste0("Resolution_Heatmap_of_sig_genes_with_baseMean_cutoff_", baseMean_cutoff,
                "_padj_cf_", padj_cutoff, 
                "_log2FoldChange_cf_", log2FoldChange_cutoff,
                "_with_threshold_", threshold,
                "_in_contrast_", contrast, "_des_treat.png"), res=400, width= 3000, height=3000)
plot(hh) 
dev.off()


pdf(file=paste0("B_Heatmap_of_sig_genes_with_baseMean_cutoff_", baseMean_cutoff,
                "_padj_cf_", padj_cutoff, 
                "_log2FoldChange_cf_", log2FoldChange_cutoff,
                "_with_threshold_", threshold, 
                "_in_contrast_", contrast, "_des_treat.pdf"))
plot(hh) 
dev.off()



# Print also the subset of the genes --------------------------------------

mysubset <- 200
hh_sub <- hh[1:mysubset, ]

png(file=paste0("Heatmap_of_sig_genes_with_baseMean_cf_", baseMean_cutoff,
                "_padj_cf_", padj_cutoff, 
                "_log2FoldChange_cf_", log2FoldChange_cutoff,
                "_with_threshold_", threshold, 
                "_in_contrast_", contrast, "_subset", 
                mysubset, "_design.png"), 
    res=400, width= 3000, height=3000)
plot(hh_sub) 
dev.off()


pdf(file=paste0("Heatmap_of_sig_genes_with_baseMean_cf_", baseMean_cutoff,
                "_padj_cf_", padj_cutoff, 
                "_log2FoldChange_cf_", log2FoldChange_cutoff,
                "_with_threshold_", threshold, 
                "_in_contrast_", contrast, "_subset", 
               mysubset, "_design.pdf"))
plot(hh_sub) 
dev.off()


### INDIVIDUAL LISTS FROM JULIA AND MALENE HERE!
############################################

subsetname <- "Oddvar" 
Oddvar_genes <- c('ESR1', 'ESR2', 'AR', 'THRB', 'GLCCI1', 'RARA', 'RARB', 'RARG', 'PGR', 'PPARA', 'PPARG', 'PPARD', 'VDR', 
                  'NR1H3', 'PTGER2', 'PTGER1', 'PTGER3', 'PTGER4', 'PTGER4P1', 'PTGER4P2', 'PTGER4P3', 'PTGER4P2-CDK2AP2P2')

#other subsets from Julia
#genes <- read_xlsx("Genes.xlsx")  #list from Julia



setwd("C:/Users/maww/Documents/2D_vs_3D")
#genes <- read_xlsx("Gene_list_specified_30092023.xlsx") #list taken from Malene
#genes <- read_xlsx("Gene_list_specified_14112023.xlsx") #list taken from Malene
#genes <- read_xlsx("Gene_list_specified_20112023.xlsx") #with removed pseudogenes and shh correctly writen
#genes <- read_xlsx("Gene_list_specified_04122023.xlsx") #with removed pseudogenes and shh correctly writen and added cholinergic after email 412
#genes <- read_xlsx("Gene_list_specified_04122023.xlsx") #with removed pseudogenes and shh correctly writen and added cholinergic after email 412 and from Elisak 15.12
genes <- read_xlsx("Gene_list_specified_17122023.xlsx") #with removed pseudogenes and shh correctly writen and added cholinergic after email 412 and from Elisak 15.12
genes <- read_xlsx("Gene_list_specified_10012024.xlsx") #with removed pseudogenes and shh correctly writen and added cholinergic after email 412 and from Elisak 15.12
genes <- read_xlsx("Gene_list_specified_14012024.xlsx") #with removed pseudogenes and shh correctly writen and added cholinergic after email 412 and from Elisak 15.12



#Targeted_analysis
#saving_dir <- paste0("C:/Users/maww/Documents/2D_vs_3D/results/", ribosomal, "/Threshold_", threshold, "/Compared_to_", compared_to, "/logFC2_final/Targeted_analysis_20112023", sep="")
#saving_dir <- paste0("C:/Users/maww/Documents/2D_vs_3D/results/", ribosomal, "/Threshold_", threshold, "/Compared_to_", compared_to, "/logFC2_final/Targeted_analysis_04122023", sep="")
saving_dir <- paste0("C:/Users/maww/Documents/2D_vs_3D/results/", ribosomal, "/Threshold_", threshold, "/Compared_to_", compared_to, "/logFC2_final/Targeted_analysis_17122023", sep="")
saving_dir <- paste0("C:/Users/maww/Documents/2D_vs_3D/results/", ribosomal, "/Threshold_", threshold, "/Compared_to_", compared_to, "/logFC2_final/Targeted_analysis_10012024", sep="")
saving_dir <- paste0("C:/Users/maww/Documents/2D_vs_3D/results/", ribosomal, "/Threshold_", threshold, "/Compared_to_", compared_to, "/logFC2_final/Targeted_analysis_14012024", sep="")



setwd(saving_dir)

unique(genes$`Topic`)[1] #was GO_C Term


#sub = genes$`GO_C Term Desc`=='membrane' #make here for look and implement MAlene files also



#Chose all genes or only relevant_genes (that are significant)
mat <- counts(dds_meth_res, normalized=T)[relevant_genes, ]
dim(mat) #7983 33
mat <- counts(dds_meth_res, normalized=T)
dim(mat) #19621 33
mat.z <- t(apply(mat, 1, scale))  #head(t(scale(t(mat)))) scaling on the rows
head(mat.z)
colnames(mat.z) <-  colnames(mat)
dim(mat.z)

# Targeted lists from DEGs ----------------------------------------------------------
# STOP here -- developed 12.11 Category or Topic has to be done individually and manually, working on this when Malene noticed missing lists 12.11.2023 (fix error for Categories)


#for (i in 1:length(unique(genes$`Topic`))){ # Category or Topic
#for (i in c(1,2,3,4,5,6,7,8,9,11,13,14,15,16,17,18,19,20,21,22,23,24)){ # Category or Topic notice that for Categories significant you need to use this
for (i in c(1,2,3,4,5,6,7,8,9,11,13,14,15,16,17,18,19,20,21,22)){ # Category or Topic notice that for Categories significant you need to use this
  
    print(unique(genes$`Category`)[i]) # Category or Topic
  level <- "Category" # Category or Topic
  sub = genes$`Category`==unique(genes$`Category`)[i] #6 Topics or 20 Category Category or Topic
  subsetname <- unique(genes$`Category`)[i] #6 6 Topics or 20 Category # Category or Topic
  mysubset <- genes$`Gene symbol`[sub]
  
  
  mat.zz <- mat.z[rownames(mat.z)%in%unique(mysubset), ]
  if (length(mat.zz) >= 2) { 
  
    
    hhh <- Heatmap(mat.zz, cluster_rows = T, cluster_columns = T,
                  column_title =paste0(subsetname),
                  column_labels=colnames(mat.z), name="Z-score", 
                  row_names_gp = gpar(fontsize = 9), 
                  top_annotation = HeatmapAnnotation(Time = time_names, 
                                                     Dimension = dimension_names,
                                                     col = list(Time = my_colour_time, 
                                                                Dimension = my_colour_dimension)),
                  #bottom_annotation = HeatmapAnnotation(Dose = dose_names, col = list(Dose = my_col_dose)),
                  column_names_gp = gpar(fontsize = 0)) #, col=treatment_colors
    
    #Print png and pdf (working nicely)
    
    # png(file=paste0("Subset_", subsetname, "_sig_genes_with_baseMean_cutoff_", baseMean_cutoff,
    #                 "_padj_cf_", padj_cutoff, 
    #                 "_log2FoldChange_cf_", log2FoldChange_cutoff,
    #                 "_with_threshold_", threshold, 
    #                 "_in_contrast_", contrast, "_lr.png")) #, res=400, width= 3000, height=3000)
    # plot(hhh) 
    # dev.off()
    # 
    # 
    # png(file=paste0("Subset_", subsetname, "_sig_genes_with_baseMean_cutoff_", baseMean_cutoff,
    #                 "_padj_cf_", padj_cutoff, 
    #                 "_log2FoldChange_cf_", log2FoldChange_cutoff,
    #                 "_with_threshold_", threshold,
    #                 "_in_contrast_", contrast, "_hr.png"), res=400, width= 3000, height=3000)
    # plot(hhh) 
    # dev.off()
    
    
    pdf(file=paste0("Subset_", level, "_", i, "_with_unsignificant.pdf"))
    #png(file=paste0("Subset_", level, "_", i, "_only_significant.png")) #if you have significant genes
    plot(hhh) 
    dev.off()
    }
    else {
      
    }
}
  






##### Malene
genes <- read_xlsx("Gene_list_ML.xlsx")  #list from Malene 38
unique(genes$Functional_relation)[2]  ### was GO_C_Term Dewsc 

#sub = genes$`GO_C Term Desc`=='membrane' #make here for look and implement MAlene files also
sub = genes$Functional_relation==unique(genes$Functional_relation)[32]#32
subsetname = unique(genes$Functional_relation)[32] #32
mysubset <- genes$Gene_Identifer[sub]

mat.zz <- mat.z[rownames(mat.z)%in%unique(mysubset), ]
dim(mat.zz)
mat.zz
hhh <- Heatmap(mat.zz, cluster_rows = T, cluster_columns = T,
               column_title =paste0("base_cf", baseMean_cutoff,
                                    "_padj_cf_", padj_cutoff, 
                                    "_log2Fold_cf_", log2FoldChange_cutoff,
                                    "_thr_", threshold, 
                                    "_in_contr_", contrast),
               #col =col_Zscore, 
               column_labels=colnames(mat.z), name="Z-score", 
               row_names_gp = gpar(fontsize = 3), 
               top_annotation = HeatmapAnnotation(Time = time_names, 
                                                  Dimension = dimension_names,
                                                  col = list(Time = my_colour_time, 
                                                             Dimension = my_colour_dimension)),
               #bottom_annotation = HeatmapAnnotation(Dose = dose_names, col = list(Dose = my_col_dose)),
               column_names_gp = gpar(fontsize = 3)) #, col=treatment_colors

png(file=paste0("Sub_", substr(subsetname,1,14), "baseMean_", baseMean_cutoff,
                "_padj_cf_", padj_cutoff, 
                "_log2FC_cf_", log2FoldChange_cutoff,
                "_with_threshold_", threshold,
                "_in_", contrast, "_hr.png"), res=400, width= 3000, height=3000)
plot(hhh) 
dev.off()

pdf(file=paste0("Sub_", substr(subsetname,1,14), "baseMean_", baseMean_cutoff,
                "_padj_cf_", padj_cutoff, 
                "_log2FoldChange_cf_", log2FoldChange_cutoff,
                "_with_threshold_", threshold, 
                "_in_", contrast, ".pdf"))
plot(hhh) 
dev.off()









###########################################


#install.packages("gplots")
#install.packages("heatmap.plus")
#install.packages("RColorBrewer")

library("gplots")
library("heatmap.plus") #not working
library("RColorBrewer")

#Making the figures for major cell types markers 
setwd(saving_dir)
subfolder_name <- file.path(getwd(), "major_cell_type_marker")

# Create the subfolder if it doesn't exist
if (!file.exists(subfolder_name)) {
  dir.create(subfolder_name)
}

# Set the working directory to the subfolder
setwd(subfolder_name)

#####CHECK THIS!!!! ###add this ADRA2C

#Printing individual genes 
# GABRA1 that was GABAergic is not among significant genes, so I removed it to avoid the error 
# malene feedback change to GABBR1  and added 2 extra genes email from Malene 20.10.2023
major_cell_type_marker <- data.frame(Category=c("Neural Progenitor", "Neural Progenitor", "Neurons","Neurons", "Astrocytes", "Astrocytes", "Astrocytes",
                                                "Pre-synaptic", "Pre-synaptic", "Dopaminergic", "Glutamatergic", "Serotonergic", "GABAergic", "Noradrenergic ", "Glycerinergic", 
                                                "Cholinergic", "Cholinergic", "Cholinergic", "Cholinergic"), #added after malene email
                                     Gene=c("PAX6", "NES", "TUBB3", "MAP2", "AQP4", "GFAP", "S100B", "SYN1","DLG4","TH","SLC17A6","TPH1", "GABBR1", "ADRA2C", "GLRB", 
                                            "CHRNA4", #added after Malene email 04.12
                                            "ACHE",
                                            "CHAT",
                                            "SLC18A3"))


for (a in 1:length(unique(major_cell_type_marker$`Category`))){ # Category only
  print((major_cell_type_marker$`Category`)[a]) # Category only
  level <- "Category" # Category only
  sub = major_cell_type_marker$`Category`==unique(major_cell_type_marker$`Category`)[a] #11 Category Category only
  subsetname <- unique(major_cell_type_marker$`Category`)[a] #6 6 Topics or 20 Category # Category or Topic
  mysubset <- major_cell_type_marker$`Gene`[sub]
  
  for (i in 1:length(mysubset)){
    print(mysubset[i])
    data <- plotCounts(dds_meth_res, normalized=T, 
                       gene=mysubset[i],
                       intgroup=c("Dimension","Time"),
                       returnData=TRUE)
    #print(data)
    myplot <- ggplot(data, mapping = aes(x=Time:Dimension, y=log10(count))) + ggtitle(mysubset[i]) + theme_bw() + geom_boxplot() #+ scale_y_log10() 
    ggsave(myplot, file=paste0(mysubset[i], "-", subsetname, "-major_cell_type_marker-2112", level, ".png"))
    ggsave(myplot, file=paste0(mysubset[i], "-", subsetname, "-major_cell_type_marker-2112", level, ".pdf"))
    #plot(myplot)
  }
}

# Initialize an empty data frame to store the data for further vizualization
combined_data <- data.frame()

# Loop through each category
for (a in 1:length(unique(major_cell_type_marker$`Category`))) {
  category <- unique(major_cell_type_marker$`Category`)[a]
  level <- "Category"
  sub <- major_cell_type_marker$`Category` == unique(major_cell_type_marker$`Category`)[a]
  subsetname <- unique(major_cell_type_marker$`Category`)[a]
  mysubset <- major_cell_type_marker$`Gene`[sub]
  
  # Loop through genes in the current category
  for (i in 1:length(mysubset)) {
    gene <- mysubset[i]
    print(paste("Category:", category, "Gene:", gene))
    
    # Use plotCounts to get data for the current gene and category
    data <- plotCounts(dds_meth_res, normalized = TRUE, 
                       gene = gene, 
                       intgroup = c("Dimension", "Time"), 
                       returnData = TRUE)
    
    # Add columns for Category and Gene
    data$Category <- category
    data$Gene <- gene
    
    # Append the data to the combined_data data frame
    combined_data <- rbind(combined_data, data)
  }
}

# Save the combined data frame as a CSV file
write.csv(combined_data, file = "combined_data_major_cell_type_marker-2112.csv", row.names = FALSE)

# Show me statistics on log scale
# Load the required libraries
library(dplyr)

# Compute statistics for each group (category and gene)
summary_data <- combined_data %>%
  group_by(Category, Gene, Dimension, Time) %>%
  summarize(
    Mean_log10_count = mean(log10(count)),
    SD_log10_count = sd(log10(count)),
    Median_log10_count = median(log10(count)),
    Min_log10_count = min(log10(count)),
    Max_log10_count = max(log10(count))
  )

# Save the summary data as a CSV file
write.csv(summary_data, file = "summary_combined_data_major_cell_type_marker-2112.csv", row.names = FALSE)


# # Create a box plot using ggplot2
# ggplot(combined_data, aes(x = factor(Time), y = log10(count), fill = Dimension)) +
#   geom_boxplot() +
#   facet_grid(Category ~ Gene, scales = "free_y") +
#   ggtitle("Box Plot of log10(count) by Category") +
#   labs(x = "Time", y = "log10(count)") +
#   scale_fill_manual(values = c("2D" = "blue", "3D" = "red")) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   guides(fill = guide_legend(title = "Dimension"))


# Load the required libraries
library(ggplot2)
library(cowplot)

# Create a list to store individual plots
plot_list <- list()

# Loop through each category
for (category in unique(combined_data$Category)) {
  # Subset the data for the current category
  category_data <- combined_data[combined_data$Category == category, ]
  
  # Create a box plot for the current category
  plot <- ggplot(category_data, aes(x = factor(Time), y = log10(count), fill = Gene)) +
    #facet_grid(. ~ Dimension) +
    geom_boxplot() +
    ggtitle(paste("Box Plot of log10(count) for", category)) +
    labs(x = "Time", y = "log10(count)") +
    scale_fill_discrete(guide = guide_legend(title = "Genes")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Add the plot to the list
  plot_list[[category]] <- plot
}

# Arrange and display the individual plots with legends in 2 columns and 4 rows
plot_grid(plotlist = plot_list, ncol = 2)

### Statified by 2D and 3D  ####
# Load the required libraries


library(ggplot2)
library(cowplot)

# Create a list to store individual plots
plot_list <- list()

# Loop through each category
for (category in unique(combined_data$Category)) {
  # Subset the data for the current category
  category_data <- combined_data[combined_data$Category == category, ]
  
  # Create a box plot for the current category with stratification by Dimension
  plot <- ggplot(category_data, aes(x = factor(Time), y = log10(count), fill = Gene)) +
    geom_boxplot() +
    ggtitle(paste("Box Plot of log10(count) for", category)) +
    labs(x = "Time", y = "log10(count)") +
    scale_fill_discrete(guide = guide_legend(title = "Genes")) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_grid(. ~ Dimension)  # Stratify by Dimension
  
  # Add the plot to the list
  plot_list[[category]] <- plot
}

# Arrange and display the individual plots with legends in 2 columns and 4 rows
myplot_grid <- plot_grid(plotlist = plot_list, ncol = 3, nrow=4)
ggsave(myplot_grid, file=paste0("combined_data_major_cell_type_marker-2112", ".png"), width = 30, height = 30, units = "cm")
ggsave(myplot_grid, file=paste0("combined_data_major_cell_type_marker-2112", ".pdf"), width = 30, height = 30, units = "cm")







#########################################################################
########################################################################
###### CLEAN ME START 13.11.2023
#########################################################################



for (i in 1:length(unique(genes$`Topic`))){ # Category or Topic
  #for (i in c(1,2,3,4,5,6,7,8,9,11,13,14,15,16,17,18,19,20)){ # Category or Topic notice that for Categories significant you need to use this
  print(unique(genes$`Topic`)[i]) # Category or Topic
  level <- "Topic" # Category or Topic
  sub = genes$`Topic`==unique(genes$`Topic`)[i] #6 Topics or 20 Category Category or Topic
  subsetname <- unique(genes$`Topic`)[i] #6 6 Topics or 20 Category # Category or Topic
  mysubset <- genes$`Gene symbol`[sub]
  
  
  mat.zz <- mat.z[rownames(mat.z)%in%unique(mysubset), ]
  if (length(mat.zz) >= 2) { 
    
    
    
genes    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
Topic <- c(
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Cell type markers',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
#'Neurotransmitter receptors signaling', #
#'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
#'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Development markers', # introduced back after corection ssh to shh 
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Development markers',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Neural network formation',
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
#'Endocrine disruption', 
'Endocrine disruption', #
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', #
#'Endocrine disruption', #
#'Endocrine disruption', #
#'Endocrine disruption', #
#'Endocrine disruption', #
#'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption', 
'Endocrine disruption',
'Endocrine disruption', 
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Endocrine disruption',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling', #added after
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
#'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling',
'Neurotransmitter receptors signaling')

  
length(Topic)           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           

           
           
Category <- c(
'Neural Progenitor', 
'Neural Progenitor',
'Neural Progenitor', 
'Neural Progenitor', 
'Neural Progenitor', 
'Neural Progenitor', 
'Immature neuron',
'Immature neuron',
'Immature neuron',
'Astrocytes',
'Astrocytes',
'Astrocytes',
'Astrocytes',
'Astrocytes',
'Astrocytes',
'Mature neurons',
'Mature neurons',
'Mature neurons',
'Mature neurons',
'Pre-synaptic',
'Pre-synaptic',
'Pre-synaptic',
'Pre-synaptic',
'Pre-synaptic',
'Pre-synaptic',
'Pre-synaptic',
'Pre-synaptic',
'Post-synaptic', 
'Post-synaptic',
'Post-synaptic',
'Post-synaptic', 
'Post-synaptic', 
'Post-synaptic', 
'Post-synaptic', 
'Glutamatergic', 
'Glutamatergic', 
'Glutamatergic', 
'Glutamatergic', 
'Glutamatergic', 
'Glycerinergic', 
'Glycerinergic', 
'Glycerinergic', 
'Glycerinergic', 
'Glycerinergic', 
'Dopaminergic',
'Dopaminergic',
'Dopaminergic',
'Noradrenergic', 
'Noradrenergic',
'Noradrenergic', 
'Noradrenergic', #
#'GABAergic',
#'GABAergic',
'GABAergic',
'GABA switch',
'GABA switch',
'Serotonergic', 
#'Serotonergic',
'Serotonergic', 
'shh pathway', #introduced back after ssh correction to shh
'shh pathway', 
'shh pathway', 
'shh pathway', 
'shh pathway', 
'shh pathway', 
'Serotonergic', 
'Serotonergic', 
'Dopaminergic',
'Dopaminergic',
'Neural Development',
'Neural Development',
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural fate commitment', 
'Neural Development',
'Neural Development',
'Axon Guidance  Signaling',
'Axon Guidance  Signaling',
'Axon Guidance  Signaling',
'Axon Guidance  Signaling',
'Axon Guidance  Signaling',
'Axon Guidance  Signaling',
'Axon Guidance  Signaling',
'Axon Guidance  Signaling',
'Axon Guidance  Signaling',
'Axon Guidance  Signaling',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Cell Adhesion  Motility',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
#'Endocrine receptors',
'Endocrine receptors', #
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors', #
#'Endocrine receptors', #
#'Endocrine receptors', #
#'Endocrine receptors', #
#'Endocrine receptors', #
#'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Endocrine receptors',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Cell Cycle',
'Endocrine receptors',
'Cholinergic',
'Cholinergic',
'Cholinergic', 
'Cholinergic',
'Cholinergic receptor nicotinic alpha X subunit',
'Cholinergic receptor nicotinic alpha X subunit',
'Cholinergic receptor nicotinic alpha X subunit',
'Cholinergic receptor nicotinic alpha X subunit',
'Cholinergic receptor nicotinic alpha X subunit',
'Cholinergic receptor nicotinic alpha X subunit',
'Cholinergic receptor nicotinic alpha X subunit',
'Cholinergic receptor nicotinic alpha X subunit',
'Cholinergic receptor nicotinic alpha X subunit',
'Cholinergic receptor nicotinic beta X subunit',
'Cholinergic receptor nicotinic beta X subunit',
#'Cholinergic receptor nicotinic beta X subunit',
'Cholinergic receptor nicotinic beta X subunit',
'Cholinergic receptor muscarinic X',
'Cholinergic receptor muscarinic X',
'Cholinergic receptor muscarinic X',
'Cholinergic receptor muscarinic X',
'Cholinergic receptor muscarinic X')







length(Category)




















#Removed:
# 'SLC6A2', - I do not see these   
# 'GABRA1',
# 'PGR',
# 'SSH',
# 'PTGER4P1',
# 'PTGER4P2',
# 'PTGER4P3',
# 'PTGER4P2',
# 'CDK2AP2P2',



plotCounts(dds_meth_res, gene="PRPH", intgroup="Time", returnData=TRUE)
plotCounts(dds_meth_res, gene="CHRNB5", intgroup="Time", returnData=FALSE)


plotCounts(dds_meth_res, gene="CDK2AP2P2", intgroup="Time", returnData=TRUE)
plotCounts(dds_meth_res, gene="CDK2AP2P2", intgroup="Time", returnData=FALSE)



check me:
plotCounts(dds_meth_res, gene="CHRNB3", intgroup="Time", returnData=FALSE)

  
CHRNA1
CHRNA2
CHRNA3
CHRNA4
CHRNA5
CHRNA6
CHRNA7
CHRNA9
CHRNA10
CHRNB1
CHRNB2
CHRNB3 ##!!!
CHRNB4
CHRM1
CHRM2
CHRM3
CHRM4
CHRM5





Gene <- 
c('CDH2',
'FABP7',
'HES5',
'PAX6',
'NES',
'TNC',
'DCX',
'TBR1',
'TUBB3',
'AQP4',
'ALDH1L1',
'GFAP',
'SLC1A2',
'SLC1A3',
'S100B',
'MAP2',
'NEFL',
'NEFM',
'RBFOX3',
'BSN',
'PCLO',
'RIMBP2',
'RIMS2',
'SNAP25',
'SYP',
'SYN1',
'UNC13A',
'DLG4',
'HOMER1',
'HOMER2',
'HOMER3',
'SHANK1',
'SHANK3',
'NLGN1',
'GRIA1',
'GLS',
'GRIN1',
'GRIN2B',
'SLC17A6',
'GLRA1',
'GLRA2',
'GLRB',
'SLC6A5',
'SLC6A9',
'DRD4',
'SLC6A3',
'TH',
'ADRA2A',
'ADRA2B',
'ADRA2C',
#'SLC6A2',
#'GABRA1',
'GABBR1',
'GAD1',
'SLC12A4',
'SLC12A1',
#'HTR1A',
'SLC6A4',
'TPH1',
'SHH',  #CHANGE!!!
'PTCH1',
'SMO',
'GLI1',
'GLI2',
'GLI3',
'NKX2-2',
'GATA2',
'PHOX2A',
'LMX1A',
'EN1',
'EN2',
'WNT1',
'TGFB2',
'TGFBR1',
'SOX1',
'RBPJ',
'PAX3',
'ISL2',
'NRG1',
'OTX2',
'DLL4',
'TLX3',
'EYA1',
'NTRK3',
'HOXC10',
'FOXG1',
'HOXA2',
'SEMA6A',
'EPHA8',
'SLIT1',
'NTNG1',
'NTNG1',
'NTN1',
'DCC',
'NRP1',
'ROBO1',
'NRTN',
'TLN1',
'VCL',
'ITGB1',
'NCAM1',
'MYH10',
'TNC',
'CNTN1',
'L1CAM',
'SDK1',
'NPNT',
'BCAN',
'NRXN1',
'NLGN1',
'ESR1',
'ESR2',
'AR',
'THRB',
'GLCCI1',
'RARA',
'RARG',
#'PGR',
'PPARA',
'PPARG',
'PPARD',
'VDR',
'NR1H3',
'PTGER2',
'PTGER1',
'PTGER3',
'PTGER4',
#'PTGER4P1',
#'PTGER4P2',
#'PTGER4P3',
#'PTGER4P2',
#'CDK2AP2P2',
'NR3C1',
'NR1H2',
'RXRA',
'RXRB',
'RXRG',
'THRA',
'RARB',
'CDK1',
'CDK2',
'CDK4',
'CDC20',
'CDC25B',
'CDC7',
'CDC16',
'CCND1',
'CCNB3',
'CCNA1',
'FOXM1',
'BUB1',
'E2F1',
'AHR',
'CHRNA4',
'ACHE',
'CHAT',
'SLC18A3', #added
'CHRNA1',
'CHRNA2',
'CHRNA3',
'CHRNA4',
'CHRNA5',
'CHRNA6',
'CHRNA7',
'CHRNA9',
'CHRNA10',
'CHRNB1',
'CHRNB2',
#'CHRNB3',
'CHRNB4',
'CHRM1',
'CHRM2',
'CHRM3',
'CHRM4',
'CHRM5')




length(Gene)






































# Gene <- c(
# 'CDH2',
# 'FABP7',
# 'HES5',
# 'PAX6',
# 'NES',
# 'TNC',
# 'DCX',
# 'TBR1',
# 'TUBB3',
# 'AQP4',
# 'ALDH1L1',
# 'GFAP',
# 'SLC1A2',
# 'SLC1A3',
# 'S100B',
# 'MAP2',
# 'NEFL',
# 'NEFM',
# 'RBFOX3',
# 'BSN',
# 'PCLO',
# 'RIMBP2',
# 'RIMS2',
# 'SNAP25',
# 'SYP',
# 'SYN1',
# 'UNC13A',
# 'DLG4',
# 'HOMER1',
# 'HOMER2',
# 'HOMER3',
# 'SHANK1',
# 'SHANK3',
# 'NLGN1',
# 'GRIA1',
# 'GLS',
# 'GRIN1',
# 'GRIN2B',
# 'SLC17A6',
# 'GLRA1',
# 'GLRA2',
# 'GLRB',
# 'SLC6A5',
# 'SLC6A9',
# 'DRD4',
# 'SLC6A3',
# 'TH',
# 'ADRA2A',
# 'ADRA2B',
# 'ADRA2C',
# 'ADRA2C', #
# 'ADRA2C',
# 'GABBR1',
# 'GAD1',
# 'SLC12A4',
# 'SLC12A1',
# 'SLC6A4',
# 'SLC6A4',
# 'TPH1',
# 'PTCH1', #
# 'PTCH1',
# 'SMO',
# 'GLI1',
# 'GLI2',
# 'GLI3',
# 'NKX2-2',
# 'GATA2',
# 'PHOX2A',
# 'LMX1A',
# 'EN1',
# 'EN2',
# 'WNT1',
# 'TGFB2',
# 'TGFBR1',
# 'SOX1',
# 'RBPJ',
# 'PAX3',
# 'ISL2',
# 'NRG1',
# 'OTX2',
# 'DLL4',
# 'TLX3',
# 'EYA1',
# 'NTRK3',
# 'HOXC10',
# 'FOXG1',
# 'HOXA2',
# 'SEMA6A',
# 'EPHA8',
# 'SLIT1',
# 'NTNG1',
# 'LRRC4C',
# 'NTN1',
# 'DCC',
# 'NRP1',
# 'ROBO1',
# 'NRTN',
# 'TLN1',
# 'VCL',
# 'ITGB1',
# 'NCAM1',
# 'MYH10',
# 'TNC',
# 'CNTN1',
# 'L1CAM',
# 'SDK1',
# 'NPNT',
# 'BCAN',
# 'NRXN1',
# 'NLGN1',
# 'ESR1',
# 'ESR2',
# 'AR',
# 'THRB',
# 'GLCCI1',
# 'RARA',
# 'RARG',
# 'PPARA',
# 'PPARA', #
# 'PPARG',
# 'PPARD',
# 'VDR',
# 'NR1H3',
# 'PTGER2',
# 'PTGER1',
# 'PTGER3',
# 'PTGER4',
# 'PTGER4', #
# 'PTGER4', #
# 'PTGER4', #
# 'PTGER4', #
# 'PTGER4', #
# 'NR3C1',
# 'NR1H2',
# 'RXRA',
# 'RXRB',
# 'RXRG',
# 'THRA',
# 'RARB',
# 'CDK1',
# 'CDK2',
# 'CDK4',
# 'CDC20',
# 'CDC25B',
# 'CDC7',
# 'CDC16',
# 'CCND1',
# 'CCNB3',
# 'CCNA1',
# 'FOXM1',
# 'BUB1',
# 'E2F1'
# )
                      
length(Gene)
    

### Change level !!!!

major_cell_type_marker_B <- data.frame(Category=Category, Gene=Gene)
level = "Category"

major_cell_type_marker_B <- data.frame(Category=Topic, Gene=Gene)  # Be careful here with Category as I used it for Category first and Topic later
level = "Topic"


for (a in 1:length(unique(genes$`Topic`))){ # Category or Topic
  print((genes$`Topic`)[a]) # Category or Topic
  level <- "Topic" # Category or Topic
  #level <- "Topic" # Category or Topic
  sub = major_cell_type_marker_B$`Category`==unique(major_cell_type_marker_B$`Category`)[a] #6 Topics or 20 Category Category or Topic
  subsetname <- unique(major_cell_type_marker_B$`Category`)[a] #6 6 Topics or 20 Category # Category or Topic
  mysubset <- major_cell_type_marker_B$`Gene`[sub]
  
  for (i in 1:length(mysubset)){
    print(mysubset[i])
    data <- plotCounts(dds_meth_res, normalized=T, 
                       gene=mysubset[i],
                       intgroup=c("Dimension","Time"),
                       returnData=TRUE)
    #print(data)
    myplot <- ggplot(data, mapping = aes(x=Time:Dimension, y=log10(count))) + ggtitle(mysubset[i]) + theme_bw() + geom_boxplot() #+ scale_y_log10() 
    ggsave(myplot, file=paste0(mysubset[i], "-", subsetname, "-targeted-2112_", level, ".png"))
    ggsave(myplot, file=paste0(mysubset[i], "-", subsetname, "-tergeted-2112_", level, ".pdf"))
    #plot(myplot)
  }
}

# Initialize an empty data frame to store the data for further vizualization
combined_data <- data.frame()

# Loop through each category
for (a in 1:length(unique(major_cell_type_marker_B$`Category`))) { #Topic
  category <- unique(major_cell_type_marker_B$`Category`)[a]
  level <- "Category"
  sub <- major_cell_type_marker_B$`Category` == unique(major_cell_type_marker_B$`Category`)[a]
  subsetname <- unique(major_cell_type_marker_B$`Category`)[a]
  mysubset <- major_cell_type_marker_B$`Gene`[sub]
  
  # Loop through genes in the current category
  for (i in 1:length(mysubset)) {
    gene <- mysubset[i]
    print(paste("Category:", category, "Gene:", gene))
    
    # Use plotCounts to get data for the current gene and category
    data <- plotCounts(dds_meth_res, normalized = TRUE, 
                       gene = gene, 
                       intgroup = c("Dimension", "Time"), 
                       returnData = TRUE)
    
    # Add columns for Category and Gene
    data$Category <- category
    data$Gene <- gene
    
    # Append the data to the combined_data data frame
    combined_data <- rbind(combined_data, data)
  }
}

# Save the combined data frame as a CSV file
write.csv(combined_data, file = "combined_target_Category_marker-2112.csv", row.names = FALSE)
#write.csv(combined_data, file = "combined_target_Topic_marker-2112.csv", row.names = FALSE)

# Show me statistics on log scale
# Load the required libraries
library(dplyr)

# Compute statistics for each group (category and gene)
summary_data <- combined_data %>%
  group_by(Category, Gene, Dimension, Time) %>%
  summarize(
    Mean_log10_count = mean(log10(count)),
    SD_log10_count = sd(log10(count)),
    Median_log10_count = median(log10(count)),
    Min_log10_count = min(log10(count)),
    Max_log10_count = max(log10(count))
  )

# Save the summary data as a CSV file
write.csv(summary_data, file = "summary_combined_target_Category_marker-2112.csv", row.names = FALSE)
#write.csv(summary_data, file = "summary_combined_target_Topic_marker-2112.csv", row.names = FALSE)


#######

# # Create a box plot using ggplot2
# ggplot(combined_data, aes(x = factor(Time), y = log10(count), fill = Dimension)) +
#   geom_boxplot() +
#   facet_grid(Category ~ Gene, scales = "free_y") +
#   ggtitle("Box Plot of log10(count) by Category") +
#   labs(x = "Time", y = "log10(count)") +
#   scale_fill_manual(values = c("2D" = "blue", "3D" = "red")) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   guides(fill = guide_legend(title = "Dimension"))


# Load the required libraries
library(ggplot2)
library(cowplot)

# Create a list to store individual plots
plot_list <- list()

# Loop through each category
for (category in unique(combined_data$Category)) {
  # Subset the data for the current category
  category_data <- combined_data[combined_data$Category == category, ]
  
  # Create a box plot for the current category
  plot <- ggplot(category_data, aes(x = factor(Time), y = log10(count), fill = Gene)) +
    #facet_grid(. ~ Dimension) +
    geom_boxplot() +
    ggtitle(paste("Box Plot of log10(count) for", category)) +
    labs(x = "Time", y = "log10(count)") +
    scale_fill_discrete(guide = guide_legend(title = "Genes")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Add the plot to the list
  plot_list[[category]] <- plot
}

# Arrange and display the individual plots with legends in 2 columns and 4 rows
plot_grid(plotlist = plot_list, ncol = 2)

### Statified by 2D and 3D  ####
# Load the required libraries


library(ggplot2)
library(cowplot)

# Create a list to store individual plots
plot_list <- list()

# Loop through each category
for (category in unique(combined_data$Category)) {
  # Subset the data for the current category
  category_data <- combined_data[combined_data$Category == category, ]
  
  # Create a box plot for the current category with stratification by Dimension
  plot <- ggplot(category_data, aes(x = factor(Time), y = log10(count), fill = Gene)) +
    geom_boxplot() +
    ggtitle(paste("Box Plot of log10(count) for", category)) +
    labs(x = "Time", y = "log10(count)") +
    scale_fill_discrete(guide = guide_legend(title = "Genes")) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_grid(. ~ Dimension)  # Stratify by Dimension
  
  # Add the plot to the list
  plot_list[[category]] <- plot
}

# Arrange and display the individual plots with legends in 2 columns and 4 rows
myplot_grid <- plot_grid(plotlist = plot_list, ncol = 2)
ggsave(myplot_grid, file=paste0("combined_data_major_cell_type_marker-2410", ".png"))
ggsave(myplot_grid, file=paste0("combined_data_major_cell_type_marker-2410", ".pdf"))









##########################################################################
###### CLEAN ME STOP 13.11.2023
##########################################################################






# Make xy plot 
# Create a list of unique categories

library(DESeq2)
library(ggplot2)
library(patchwork)












### end


###https://biohpc.cornell.edu/doc/RNA-Seq-2020-exercise2.html

# NB: In case you need individual gene ----------------------------------------
#dds_meth_res <- readRDS(file="with_ribo_dds_meth_res_design_treatment_filtered_2D_3.rds")


for (i in 1:length(mysubset)){
  print(mysubset[i])
}


targeted_list <- list()



# Counst data normalized on 2D and 3D targeted genes ----------------------
#Print all genes into separate figures  #I mixed SLIT1 it was small letters and NGL-1 was replaced by LRRC4C
for (a in 1:length(unique(genes$`Category`))){ # Category or Topic
  print((genes$`Category`)[a]) # Category or Topic
  level <- "Category" # Category or Topic
  sub = genes$`Category`==unique(genes$`Category`)[a] #6 Topics or 20 Category Category or Topic
  subsetname <- unique(genes$`Category`)[a] #6 6 Topics or 20 Category # Category or Topic
  mysubset <- genes$`Gene symbol`[sub]
  
    for (i in 1:length(mysubset)){
      print(mysubset[i])
      data <- plotCounts(dds_meth_res, normalized=T, 
                         gene=mysubset[i],
                         intgroup=c("Dimension","Time"),
                         returnData=TRUE)
      #print(data)
      myplot <- ggplot(data, mapping = aes(x=Time:Dimension, y=count)) + ggtitle(mysubset[i]) + theme_bw() + geom_boxplot() #+ scale_y_log10() 
      ggsave(myplot, file=paste0(mysubset[i], "-", subsetname, "-", level, ".png"))

    }
}

# Trying to plot all categories on one plot!

library(gridExtra)

data_list <- list()

for (i in 1:length(unique(genes$`Category`))){ # Category or Topic
  level <- "Category" # Category or Topic
  sub = genes$`Category`==unique(genes$`Category`)[i] #6 Topics or 20 Category Category or Topic
  subsetname <- unique(genes$`Category`)[i] #6 6 Topics or 20 Category # Category or Topic
  mysubset <- unique(genes$`Gene symbol`[sub]) #gives me genes that comes for Category
  length(mysubset) #prints number of genes will be ploted 
  
  for (j in 1:length(mysubset)){
    #print(mysubset[j])
    data <- plotCounts(dds_meth_res, normalized=T, 
                       gene=mysubset[j],
                       intgroup=c("Dimension","Time"),
                       returnData=TRUE)
    data$geneID <- mysubset[j]
    data$Category <- subsetname
    
    data_list[[length(data_list) + 1]] <- data
    
    
    #datalists <- data[[j]]
    #data$geneID <- rbind(data, mysubset[j]) 
    
    #print(data)
    #myplots <- as.list(ggplot(data, mapping = aes(x=Time:Dimension, y=count)) + ggtitle(subsetname[j]) + theme_bw() + geom_boxplot()) #+ scale_y_log10() 
    #ggsave(myplots, file=paste0(mysubset, "-", subsetname, "-", level, ".png"))
  }
  #datalists[i] <- data
}


# Genes_by_Categories 01102023 -----------------------------------------------------

setwd("C:/Users/maww/Documents/2D_vs_3D/results/ribosomal_not_removed/Threshold_10/Compared_to_2D_Day3/logFC2_final/Targeted_analysis")
### WORKI here the 8th 
# Combine all data frames into one
combined_data <- do.call(rbind, data_list)
custom_colors <- c("lightblue", "orange")
indices <- list(a=1:6, b=7:9, c=10:15, d=16:19, e=20:27, f=28:34, g=35:39, 
                h=40:44, i=45:49, j=50:52, k=53:55, l=56:57, m=58:61, n=62:66, o=67:70,
                p=71:84, r=85:94, s=95:107, t=108:130, u=131:143)
names(indices) <- c("Neural Progenitor", "Immature neuron", "Astrocytes", "Mature neurons", "Pre-synaptic", "Post-synaptic", "Glutamatergic",
                    "Glycerinergic", "Dopaminergic", "Noradrenergic", "GABAergic", "GABA switch", "Serotonergic", "shh pathway", "Neural Development", 
                    "Neural fate commitment", "Axon Guidance / Signaling", "Cell Adhesion / Motility", "Endocrine receptors", "Cell Cycle")

# Create the subfolder path relative to the current working directory
subfolder_name <- file.path(getwd(), "Canva_Joint_08102023")

# Create the subfolder if it doesn't exist
if (!file.exists(subfolder_name)) {
  dir.create(subfolder_name)
}

# Set the working directory to the subfolder
setwd(subfolder_name)


for (z in indices) {
print(z)

subset_data_list <- data_list[z]
  
combined_data <- do.call(rbind, subset_data_list) 
#print("This will be", unique(combined_data$Category))

# Create a boxplot with facets
gg <- ggplot(combined_data, aes(x = Time:Dimension, y = count, fill = Dimension)) +
  geom_boxplot() +
  #scale_fill_discrete() +
  scale_fill_manual(values = setNames(custom_colors, levels(combined_data$Dimension))) +
  ggtitle(unique(combined_data$Category)) +
  theme_bw() +
  facet_wrap(~ geneID, ncol = 4, nrow = 6, scales = "free_y") + # Adjust ncol as needed
  xlab("Days") 

#gg <- gg + scale_x_discrete(labels = function(x) gsub("^\\d+:","",x))
gg <- gg + scale_x_discrete(labels = function(x) gsub(":.*$", "", x))


# Print or save the ggplot object
#print(gg)
ggsave(gg, file=paste0(unique(combined_data$Category),"nexttoeachother.png"))
ggsave(gg, file=paste0(unique(combined_data$Category),"nexttoeachother.pdf"))
}



# ### Make changes so that the 2D is desplayed on left and 3D on r --------

setwd("C:/Users/maww/Documents/2D_vs_3D/results/ribosomal_not_removed/Threshold_10/Compared_to_2D_Day3/logFC2_final/Targeted_analysis")

# Combine all data frames into one
combined_data <- do.call(rbind, data_list)

# Create custom colors
custom_colors <- c("lightblue", "orange")

indices <- list(a=1:6, b=7:9, c=10:15, d=16:19, e=20:27, f=28:34, g=35:39, 
                h=40:44, i=45:49, j=50:52, k=53:55, l=56:57, m=58:61, n=62:66, o=67:70,
                p=71:84, r=85:94, s=95:107, t=108:130, u=131:143)
names(indices) <- c("Neural Progenitor", "Immature neuron", "Astrocytes", "Mature neurons", "Pre-synaptic", "Post-synaptic", "Glutamatergic",
                    "Glycinergic", "Dopaminergic", "Noradrenergic", "GABAergic", "GABA switch", "Serotonergic", "shh pathway", "Neural Development", 
                    "Neural fate commitment", "Axon Guidance / Signaling", "Cell Adhesion / Motility", "Endocrine receptors", "Cell Cycle")


# Create the subfolder path relative to the current working directory
original_working_directory <- saving_dir
setwd(original_working_directory)
subfolder_name <- file.path(getwd(), "Canva_2Dleft_3Dright_8102023")

# Create the subfolder if it doesn't exist
if (!file.exists(subfolder_name)) {
  dir.create(subfolder_name)
}

# Set the working directory to the subfolder
setwd(subfolder_name)


for (z in indices) {
  print(z)
  
  subset_data_list <- data_list[z]
  
  combined_data <- do.call(rbind, subset_data_list)
  
  # Create a boxplot with facets and switch the order of Dimension within each subplot
  gg <- ggplot(combined_data, aes(x = Time, y = count, fill = Dimension)) +
    geom_boxplot() +
    scale_fill_manual(values = setNames(custom_colors, levels(combined_data$Dimension))) +
    ggtitle(unique(combined_data$Category)) +
    theme_bw() +
    facet_grid(geneID ~ Dimension, scales = "free_y", switch = "x") + # Switch the order of Dimension within each subplot
    xlab("Days")
  
  # Print or save the ggplot object
  ggsave(gg, file = paste0(unique(combined_data$Category), "_2Dleft_3Dright.png"))
  ggsave(gg, file = paste0(unique(combined_data$Category), "_2Dleft_3Dright.pdf"))
}











#### Notice that in the excel SLC6A2 was replaced to SLC6A5





data <- plotCounts(dds_meth_res, normalized=T, 
                 gene="PTGER4", #SLC6A2
                 intgroup=c("Dimension","Time"),
                 returnData=TRUE)




data <- plotCounts(dds_meth_res, normalized=T, 
                   gene="RARB",
                   intgroup=c("Dimension", "Time"),
                   returnData=TRUE)

ggplot(data, mapping = aes(x=Time, y=count, z=Dimension)) + theme_bw() + geom_boxplot() #scale_y_log10() 
ggsave("ADRB2.png")




### Volcano plots for different comparisons 
#https://www.biostars.org/p/282295/
https://www.biostars.org/p/282295/
  

# Volcano plots I red and black (based on the results to check) -------------

#Making the figures for major cell types markers 
setwd(saving_dir)
subfolder_name_vulcI <- file.path(getwd(), "Vulcanos_30102023_I")

# Create the subfolder if it doesn't exist
if (!file.exists(subfolder_name_vulcI)) {
  dir.create(subfolder_name_vulcI)
}

setwd(subfolder_name_vulcI)

for (i in 1:length(results_to_check)) {
  print(names(results_to_check[i]))
  tmp = results_to_check[[i]]
  #png(paste0("My_Vulcano_", names(results_to_check[i]), "_log2FC2.png"))
  pdf(paste0("My_Vulcano_", names(results_to_check[i]), "_log2FC2.pdf"))
  par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
  topT <- as.data.frame(tmp)
  #Adjusted P values (FDR Q values)
  with(topT, plot(log2FoldChange, -log10(padj), pch=20, main=names(results_to_check[i]), 
                  cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
  with(subset(topT, padj<0.01 & abs(log2FoldChange)>4), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5)) #was 2

  #with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))
  #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
  abline(v=0, col="black", lty=3, lwd=1.0)
  abline(v=-2, col="black", lty=4, lwd=2) #was lwd=2.0
  abline(v=2, col="black", lty=4, lwd=2) #was lwd=2.0
  abline(h=-log10(max(topT$pvalue[topT$padj<0.01], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
  dev.off()
}




# Volcano plots II with genes names and cuttoffs (based on the results to check) -------------
setwd(saving_dir)
subfolder_name_vulcII <- file.path(getwd(), "Vulcanos_30102023_II")

# Create the subfolder if it doesn't exist
if (!file.exists(subfolder_name_vulcII)) {
  dir.create(subfolder_name_vulcII)
}
setwd(subfolder_name_vulcII)


for (i in 1:length(results_to_check)) {
  print(names(results_to_check[i]))
  res = results_to_check[[i]]
  
  png(paste0("Genes_my_Vulcano_", names(results_to_check[i]), ".png"), width=15,height=20,units="cm",res=600)
  #pdf(paste0("Genes_my_Vulcano_", names(results_to_check[i]), ".pdf"), width=15,height=20)
  
  alpha <- 0.01 # Threshold on the adjusted p-value
  log2FoldChange_cutoff <- 2
  cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
  plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
       main=names(results_to_check[i]), xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
       pch=20, cex=0.6)
  abline(v=0)
  abline(v=c(-2,2), col="brown")
  abline(h=-log10(alpha), col="brown")
  
  gn.selected <- abs(res$log2FoldChange) > log2FoldChange_cutoff & res$padj < alpha 
  text(res$log2FoldChange[gn.selected],
       -log10(res$padj)[gn.selected],
       lab=rownames(res)[gn.selected ], cex=0.3)
  
  dev.off()
  
}



# Volcano plots III with genes names and cuttoffs (based on the results to check) with implemented ifstatment for cases with 0 genes!  -------------
setwd(saving_dir)
subfolder_name_vulcIII <- file.path(getwd(), "Vulcanos_30102023_III")

# Create the subfolder if it doesn't exist
if (!file.exists(subfolder_name_vulcIII)) {
  dir.create(subfolder_name_vulcIII)
}
setwd(subfolder_name_vulcIII)


for (i in 1:length(results_to_check)) {
  print(names(results_to_check[i]))
  res = results_to_check[[i]]
  
  #png(paste0("Style_II_Genes_my_Vulcano_", names(results_to_check[i]), ".png"), width=15,height=20,units="cm",res=600)
  pdf(paste0("Style_II_Genes_my_Vulcano_", names(results_to_check[i]), ".pdf"), width=15,height=20)
  
  alpha <- 0.01 # Threshold on the adjusted p-value
  log2FoldChange_cutoff <- 2 #was2
  
  cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
  plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
       main=names(results_to_check[i]), xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
       pch=20, cex=0.6)
  abline(v=0)
  abline(v=c(-2,2), col="brown")
  abline(h=-log10(alpha), col="brown")
  
  gn.selected <- abs(res$log2FoldChange) > log2FoldChange_cutoff & res$padj < alpha
    if (unique(gn.selected) %in% c('TRUE')) {
      
   
    text(res$log2FoldChange[gn.selected],
         -log10(res$padj)[gn.selected],
         lab=rownames(res)[gn.selected ], cex=0.5)
    #else i+1
    }
   else i+1
  dev.off()
  
}


#### Can I make a venn diagrams for the seperate 2D and 3D
argument 2D









# Print the lists of the genes for different contrasts in a given design ----------------
alpha_cf <- 0.01 # was 0.01## Threshold on the adjusted p-value
fold_cf <- 2 #1 #was 2

#https://stackoverflow.com/questions/14735078/how-to-handle-conditional-data-frame-subscripts-that-contain-nas

for (i in 1:length(results_to_check)) 
  {
  print(names(results_to_check[i]))
  res = results_to_check[[i]]
  res$padj[is.na(res$padj)] <- 1 #replaceing na values with 1
  gn.selected <- abs(res$log2FoldChange) > fold_cf & res$padj < alpha_cf
  newlist = res[gn.selected, ]
  print(dim(newlist))
  print(names(newlist))
}


####### 
# Write down to CSV file the genes with given p-values only for future venn diagrams  ------------------------------------------

filtered_results_XX <- list()
#for (index in 1:length(results_to_check)) {
for (index in names(results_to_check)) { #17contrasts
  newElement = results_to_check[[index]]
  newElement$padj[is.na(newElement$padj)] <- 1
  newElement <- newElement[abs(newElement$log2FoldChange) > fold_cf & newElement$padj < alpha_cf, ]
  filtered_results_XX[index] <- newElement
}
warnings()
head(filtered_results_XX)
length(filtered_results_XX) #list with 5 data frame 5*7=35 ####   was list with 17 data frames  17*7=119


#printableResults_XX = data.frame(NA, nrow = 5000, ncol = 120)
printableResults_XX = data.frame()

printableResults_XX <- cbind(c(relevant_genes, printableResults_XX))
dim(printableResults_XX) #1 column with 119 comming 

#Function for empty vectors or data frames 
isEmpty <- function(x) {
  return(length(x)==0)
}

#save the results 

#############################################################################
# Producing csv2 files with results for all contrasts ---------------------
#############################################################################

resultsNames(dds_meth_res)

alpha_cutoff <- 0.1 #this should be set to 0.1
FoldChange_cutoff <- 2

# 3D vs 2D at given time points  "3D_3 vs 2D_3"
res1 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_3", "2D_3")) #lfcThreshold=0, 
sum(res1$padj <= 0.01 & abs(res1$log2FoldChange >= FoldChange_cutoff), na.rm=TRUE)
resFilt1 <- res1[which(res1$padj < 0.01 & abs(res1$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt1) #5283
write.csv2(resFilt1, file="DE_Treatment_3D_vs_2D_day3_filtered.csv")

res2 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_14", "2D_14"))
sum(res2$padj < 0.01 & abs(res2$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #3574
resFilt2 <- res2[which(res2$padj < 0.01 & abs(res2$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt2) #5515
write.csv2(resFilt2, file="DE_Treatment_3D_vs_2D_day14_filtered.csv")

res3 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_21", "2D_21"))
sum(res3$padj < 0.01 & abs(res3$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #3710
resFilt3 <- res3[which(res3$padj < 0.01 & abs(res3$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt3) #5608
write.csv2(resFilt3, file="DE_Treatment_3D_vs_2D_day21_filtered.csv")


# Development at 2D
res4 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_3", "2D_14"))
sum(res4$padj < 0.01 & abs(res4$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #353
resFilt4 <- res4[which(res4$padj < 0.01 & abs(res4$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt4) #1152
write.csv2(resFilt4, file="DE_Treatment_day3_vs_day14_2D_filtered.csv")


res4bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_14", "2D_3")) #was opposite
sum(res4bis$padj < 0.01 & abs(res4bis$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #353
resFilt4bis <- res4bis[which(res4bis$padj < 0.01 & abs(res4bis$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt4bis) #1152
write.csv2(resFilt4bis, file="DE_Treatment_day14_vs_day3_2D_filtered.csv") #check this


res5 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_3", "2D_21"))
sum(res5$padj < 0.01 & abs(res5$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #427
resFilt5 <- res5[which(res5$padj < 0.01 & abs(res5$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt5) #1559
write.csv2(resFilt5, file="DE_Treatment_day3_vs_day21_2D_filtered.csv")

res5bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_21", "2D_3"))
sum(res5bis$padj < 0.01 & abs(res5bis$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #427
resFilt5bis <- res5bis[which(res5bis$padj < 0.01 & abs(res5bis$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt5bis) #1559
write.csv2(resFilt5bis, file="DE_Treatment_day21_vs_day3_2D_filtered.csv")

res6 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_14", "2D_21"))
sum(res6$padj < 0.01 & abs(res6$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #7
resFilt6 <- res6[which(res6$padj < 0.01 & abs(res6$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt6) #62
write.csv2(resFilt6, file="DE_Treatment_day14_vs_day21_2D_filtered.csv")

res6bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_21", "2D_14"))
sum(res6bis$padj < 0.01 & abs(res6bis$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #7
resFilt6bis <- res6bis[which(res6bis$padj < 0.01 & abs(res6bis$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt6bis) #62
write.csv2(resFilt6bis, file="DE_Treatment_day21_vs_day14_2D_filtered.csv")


#Development at 3D
res7 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_3", "3D_14"))
sum(res7$padj < 0.01 & abs(res7$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #51
resFilt7 <- res7[which(res7$padj < 0.01 & abs(res7$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt7) #522
write.csv2(resFilt7, file="DE_Treatment_day3_vs_day14_3D_filtered.csv")

res7bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_14", "3D_3"))
sum(res7bis$padj < 0.01 & abs(res7bis$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #51
resFilt7bis <- res7bis[which(res7bis$padj < 0.01 & abs(res7bis$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt7bis) #522
write.csv2(resFilt7bis, file="DE_Treatment_day14_vs_day3_3D_filtered.csv")


res8 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_3", "3D_21"))
sum(res8$padj < 0.01 & abs(res8$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #297
resFilt8 <- res8[which(res8$padj < 0.01 & abs(res8$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt8) #1313
write.csv2(resFilt8, file="DE_Treatment_day3_vs_day21_3D_filtered.csv")


res8bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_21", "3D_3"))
sum(res8bis$padj < 0.01 & abs(res8bis$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #297
resFilt8bis <- res8bis[which(res8bis$padj < 0.01 & abs(res8bis$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt8bis) #1313
write.csv2(resFilt8bis, file="DE_Treatment_day21_vs_day3_3D_filtered.csv")


res9 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_14", "3D_21"))
sum(res9$padj < 0.01 & abs(res9$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #16
resFilt9 <- res9[which(res9$padj < 0.01 & abs(res9$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt9) #30
write.csv2(resFilt9, file="DE_Treatment_day14_vs_day21_3D_filtered.csv")

res9bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_21", "3D_14"))
sum(res9bis$padj < 0.01 & abs(res9bis$log2FoldChange >FoldChange_cutoff), na.rm=TRUE) #16
resFilt9bis <- res9bis[which(res9bis$padj < 0.01 & abs(res9bis$log2FoldChange) >FoldChange_cutoff), ]
dim(resFilt9bis) #30
write.csv2(resFilt9bis, file="DE_Treatment_day21_vs_day14_3D_filtered.csv")



####Extra Vulcanos! (30.10.2023)

setwd(saving_dir)
subfolder_name_vulc_extra <- file.path(getwd(), "Vulcanos_30102023_extra")

# Create the subfolder if it doesn't exist
if (!file.exists(subfolder_name_vulc_extra)) {
  dir.create(subfolder_name_vulc_extra)
}
setwd(subfolder_name_vulc_extra)



extra_results_to_check <- list(Treatment_2D_3_2D_14=res4,
                              Treatment_2D_14_2D_3=res4bis,
                              Treatment_2D_3_2D_21=res5,
                              Treatment_2D_21_2D_3 = res5bis,
                              Treatment_2D_14_2D_21 = res6,
                              Treatment_2D_21_2D_14 = res6bis,
                              Treatment_3D_3_3D_14 = res7,
                              Treatment_3D_14_3D_3 = res7bis,
                              Treatment_3D_3_3D_21 = res8,
                              Treatment_3D_21_3D_3 = res8bis,
                              Treatment_3D_14_3D_21 = res9,
                              Treatment_3D_21_3D_14 = res9bis)

# # Development at 2D
# res4 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_3", "2D_14")) #ok
# res4bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_14", "2D_3")) #ok
# res5 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_3", "2D_21")) #ok
# res5bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_21", "2D_3")) #ok
# res6 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_14", "2D_21"))
# res6bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "2D_21", "2D_14"))
# 
# #Development at 3D
# res7 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_3", "3D_14"))
# res7bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_14", "3D_3"))
# res8 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_3", "3D_21"))
# res8bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_21", "3D_3"))
# res9 = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_14", "3D_21"))
# res9bis = results(dds_meth_res, alpha = alpha_cutoff, contrast=c("Treatment", "3D_21", "3D_14"))



for (i in 1:length(extra_results_to_check)) {
  print(names(extra_results_to_check[i]))
  res = extra_results_to_check[[i]]
  
  png(paste0("Genes_my_Vulcano_", names(extra_results_to_check[i]), "_extra.png"), width=15,height=20,units="cm",res=600)
  #pdf(paste0("Genes_my_Vulcano_", names(extra_results_to_check[i]), "_extra.pdf"), width=15,height=20)
  
  alpha <- 0.01 # Threshold on the adjusted p-value
  log2FoldChange_cutoff <- 2
  cols <- densCols(res$log2FoldChange, -log10(res$padj)) #padjust? res$pvalue
  plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
       main=names(extra_results_to_check[i]), xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
       pch=20, cex=0.6)
  abline(v=0)
  abline(v=c(-2,2), col="brown")
  abline(h=-log10(alpha), col="brown")
  
  gn.selected <- abs(res$log2FoldChange) > log2FoldChange_cutoff & res$padj < alpha 
  text(res$log2FoldChange[gn.selected],
       -log10(res$padj)[gn.selected],
       lab=rownames(res)[gn.selected ], cex=0.4)
  
  dev.off()
  
}








################################################################################
# Producing Venn diagrams based on DEGs -----------------------------------
################################################################################

#https://r-graph-gallery.com/14-venn-diagramm
# Load library
library(VennDiagram)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# 3D vs 2D 
venn.diagram(
  x = list(A=rownames(resFilt1), B=rownames(resFilt2), C=rownames(resFilt3)),
  category.names = c("3Dvs2D#Day3" , "3Dvs2D#Day14" , "3Dvs2D#Day21"),
  filename = '3D_vs_2D_different_days_venn_diagramm_controls.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1080 , 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

### 2D only development
venn.diagram(
  x = list(A=rownames(resFilt4), B=rownames(resFilt5), C=rownames(resFilt6)),
  category.names = c("3vs14#2D", "3vs21#2D" , "14vs21#2D"),
  filename = '2D_only_different_days_venn_diagramm_controls.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1080 , 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

### 3D only development
venn.diagram(
  x = list(A=rownames(resFilt7), B=rownames(resFilt8), C=rownames(resFilt9)),
  category.names = c("3vs14#3D", "3vs21#3D" , "14vs21#3D"),
  filename = '3D_only_different_days_venn_diagramm_controls.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1080 , 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


# Get a table for genes ---------------------------------------------------------------------

#rownames(filtered_results_XX[[index]])

for (index in names(filtered_results_XX)) {
  #print(index)
  if (isEmpty(rownames(filtered_results_XX[[index]]))==T) {
    printableResults_XX <- cbind(printableResults_XX, "0") #next # printableResults_XX <- cbind(printableResults_XX, "NA")
  } 
        
  else if (isEmpty(filtered_results_XX[[index]][["baseMean"]])==T) {
    printableResults_XX <- cbind(printableResults_XX, "1") #printableResults_XX <- "NA" #  next # printableResults_XX <- cbind(printableResults_XX, "NA")
    }
  
  else if (isEmpty(filtered_results_XX[[index]][["log2FoldChange"]])==T) {
        next # printableResults_XX <- cbind(printableResults_XX, "NA")
        }
  else if (isEmpty(filtered_results_XX[[index]][["lfcSE"]])==T) {
          next #printableResults_XX <- cbind(printableResults_XX, "NA")
        }
  else if (isEmpty(filtered_results_XX[[index]][["stat"]])==T) {
          next #printableResults_XX <- cbind(printableResults_XX, "NA")
        }
        
  else if (isEmpty(filtered_results_XX[[index]][["pvalue"]])==T) {
          next # printableResults_XX <- cbind(printableResults_XX, "NA")
        }
        
  else if (isEmpty(filtered_results_XX[[index]][["padj"]])==T) {
          next # printableResults_XX <- cbind(printableResults_XX, "NA")
        }
    
  else {
    printableResults_XX <- cbind(printableResults_XX, rownames(filtered_results_XX[[index]]))
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["baseMean"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["log2FoldChange"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["lfcSE"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["stat"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["pvalue"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["padj"]])
    
  }
}


dim(printableResults_XX) #176 and was 114
printableResults_XX




for (index in names(filtered_results_XX)) {
    print(index)
    printableResults_XX <- cbind(printableResults_XX, rownames(filtered_results_XX[[index]]))
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["baseMean"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["log2FoldChange"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["lfcSE"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["stat"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["pvalue"]])
    printableResults_XX <- cbind(printableResults_XX, filtered_results_XX[[index]][["padj"]])
}

warnings()
dim(printableResults_XX) #7979 and 211 and was 2378, 118

newColNames = c("Genes")
for (index in names(filtered_results_XX)){
  newColNames = c(newColNames, 
                  paste0(names(filtered_results_XX[index]),"_GeneName"), 
                  paste0(index,"_baseMean"), 
                  paste0(index,"_log2FoldChange"), 
                  paste0(index,"_lfcSE"), paste0(index,"_stat"), 
                  paste0(index,"_pvalue"), paste0(index,"_padj"))
  
}

length(newColNames) #36 and was (120 This is correct)

colnames(printableResults_XX) <- newColNames

head(printableResults_XX)



write.csv2(printableResults_XX, file=paste0("TEST_Significant_genes_with_baseMean_cutoff_", baseMean_cutoff,
                                             "_padj_cutoff_", padj_cutoff, 
                                             "_log2FoldChange_cutoff_", log2FoldChange_cutoff,
                                             "_with_threshold_", threshold, 
                                             "_in_contrast_", contrast, "_design_treatment_plus_pasage.csv"))



write.table(printableResults_XX, file=paste0("TEST_Significant_genes_with_baseMean_cutoff_", baseMean_cutoff,
                                          "_padj_cutoff_", padj_cutoff, 
                                          "_log2FoldChange_cutoff_", log2FoldChange_cutoff,
                                          "_with_threshold_", threshold, 
                                          "_in_contrast_", contrast, "_design_treatment_plus_pasage.csv"),
            row.names = F,col.names = T, sep = ";")


###############################################################################
# Enrichment analysis using enrichR ---------------------------------------
# Enrichment analysis using enrichR ---------------------------------------
###############################################################################


#install.packages("enrichR")
#https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html
library(enrichR)
  
listEnrichrSites()
websiteLive <- TRUE
dbs <- listEnrichrDbs()
  
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)
    

# Define lists with genes for the enrichment NB! Individually ------------------------------

#mygenes = rownames(mat.z)


setwd("C:/Users/maww/Documents/2D_vs_3D/Enrichments_corrected_17102023")

mygenes<- rownames(resFilt9) #from resFilt1 to resFilt9  DE_Treatment_3D_3 vs 2D_3


#Define more libraries here !!!
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "Panther_2016", 
         "KEGG_2021_Human", "Reactome_2016")
if (websiteLive) {
  #enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), dbs) # examples of data bases
  enriched <- enrichr(c(mygenes), dbs)
}

dbs


# Saving the csv files with enrichments -----------------------------------

#gene_vector <-  "DE_Treatment_3D_3_vs_2D_3" #1
#gene_vector <-  "DE_Treatment_3D_14_vs_2D_14" #2
#gene_vector <-  "DE_Treatment_3D_21_vs_2D_21" #3
#gene_vector <-  "DE_Treatment_2D_3_vs_2D_14" #4
#gene_vector <-  "DE_Treatment_2D_3_vs_2D_21" #5
#gene_vector <-  "DE_Treatment_2D_14_vs_2D_21" #6
#gene_vector <-  "DE_Treatment_3D_3_vs_3D_14" #7
#gene_vector <-  "DE_Treatment_3D_3_vs_3D_21" #8
gene_vector <-  "DE_Treatment_3D_14_vs_3D_21" #9



library("xlsx")
write.xlsx(x=enriched[[1]], file=paste0("Enrich", dbs[1], "_", gene_vector, "1710.xlsx"), row.names = F)
write.xlsx(x=enriched[[2]], file=paste0("Enrich", dbs[2], "_", gene_vector, "1710.xlsx"), row.names = F)
write.xlsx(x=enriched[[3]], file=paste0("Enrich", dbs[3], "_", gene_vector, "1710.xlsx"), row.names = F)
write.xlsx(x=enriched[[4]], file=paste0("Enrich", dbs[4], "_", gene_vector, "1710.xlsx"), row.names = F)
write.xlsx(x=enriched[[5]], file=paste0("Enrich", dbs[5], "_", gene_vector, "1710.xlsx"), row.names = F)
write.xlsx(x=enriched[[6]], file=paste0("Enrich", dbs[6], "_", gene_vector, "1710.xlsx"), row.names = F)


#if (websiteLive) enriched[[1]]

#Make a loop for this for all genes and specific contrast


# Plotting enrichments for all the databases chosen -----------------------

pdf(file=paste0("Enrichments_ordered_by_p_value_baseMean_", baseMean_cutoff, "_padj_", padj_cutoff, "_log2FoldChange_", log2FoldChange_cutoff, "_", gene_vector, ".pdf"))
if (websiteLive) plotEnrich(enriched[[1]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.value", title= paste0("Enrich: ", dbs[1]))
if (websiteLive) plotEnrich(enriched[[2]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.value", title=paste0("Enrich: ", dbs[2]))
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.value", title=paste0("Enrich: ", dbs[3]))
if (websiteLive) plotEnrich(enriched[[4]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.value", title=paste0("Enrichm: ", dbs[4]))
if (websiteLive) plotEnrich(enriched[[5]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.value", title=paste0("Enrichm: ", dbs[5]))
if (websiteLive) plotEnrich(enriched[[6]], showTerms = 15, numChar = 50, y = "Count", orderBy = "P.value", title=paste0("Enrichm: ", dbs[6]))
dev.off()


pdf(file=paste0("Enrichments_ordered_by_combined_score_baseMean", baseMean_cutoff, "_padj_", padj_cutoff, "_log2FoldChange_", log2FoldChange_cutoff, "_", gene_vector, ".pdf"))
if (websiteLive) plotEnrich(enriched[[1]], showTerms = 15, numChar = 50, y = "Count", orderBy = "Combined.Score", title= paste0("Enrich: ", dbs[1]))
if (websiteLive) plotEnrich(enriched[[2]], showTerms = 15, numChar = 50, y = "Count", orderBy = "Combined.Score", title=paste0("Enrich: ", dbs[2]))
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 15, numChar = 50, y = "Count", orderBy = "Combined.Score", title=paste0("Enrich: ", dbs[3]))
if (websiteLive) plotEnrich(enriched[[4]], showTerms = 15, numChar = 50, y = "Count", orderBy = "Combined.Score", title=paste0("Enrichm: ", dbs[4]))
if (websiteLive) plotEnrich(enriched[[5]], showTerms = 15, numChar = 50, y = "Count", orderBy = "Combined.Score", title=paste0("Enrichm: ", dbs[5]))
if (websiteLive) plotEnrich(enriched[[6]], showTerms = 15, numChar = 50, y = "Count", orderBy = "Combined.Score", title=paste0("Enrichm: ", dbs[6]))
dev.off()

##############################################################################
# Single vulcano plot ---------------------------------------------------
##############################################################################


res = results(dds_meth_res, contrast=c("Treatment", "3D_3", "2D_3"))

myvector <- "my contrast A"

png(paste0("Vulcano_", myvector, "_.png", sep=""), width = 2300, height = 2300, units = "px" )

alpha <- 0.001 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main=paste0("Volcano plot", myvector, sep=" "), xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res$log2FoldChange) > 3.5 & res$padj < alpha 
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.5)

dev.off()





# Export normalized counts ------------------------------------------------
normalized_counts <- counts(dds_meth_res, normalized=TRUE)
head(normalized_counts)
View(normalized_counts)

#####setwd("C:/Users/maww/Downloads/meth_results")

#write.table(normalized_counts, file="normalized_counts_treatment_withRNA_10cutoff.txt", sep="\t", quote=F, col.names=T)

write.table(normalized_counts, file="normalized_counts_meth_withoutRNA_10cutoff.txt", sep="\t", quote=F, col.names=T)
write.csv(normalized_counts, file="normalized_counts_meth_withoutRNA_10cutoff.csv")


#DEseq results
res <- results(dds_meth_res, alpha=0.01)
summary(res)
hist(res$padj)

#ordered by padjusted
resOrdered <- res[order(res$padj), ]
write.csv(resOrdered, file="normalized_counts_meth_withoutRNA_10cutoff_ordered.csv")
resOrdered
rownames(resOrdered[1:50, ])

table(resOrdered$padj<0.01)


####
plotMA(dds_meth_res)
plotMA(res, xlim=c(0, 50000), ylim=c(0, 50000))
plotMA(res)

#https://www.biostars.org/p/353618/#356111

#Building the results table
alpha_threshold = 0.1
##contrast_ts <- c("Time", "14", "3")
contrast_ts <- c("Dose", "0.07", "0")

contrast_ts <- c("Time", "14", "Dose", "0.07")

results(dds_meth_res, contrast=c("Time", "14", "0", "Dose", "0.07", "0"))

z = results(dds_meth_res, name="Time14.Dose0.07")
z = results(dds_meth_res, name="Time21.Dose0.07")
z = results(dds_meth_res, name="Time14.Dose2.2")
z = results(dds_meth_res, name="Time21.Dose2.2")
z = results(dds_meth_res, name="Time14.Dose6.7")
z = results(dds_meth_res, name="Time21.Dose6.7")
z = results(dds_meth_res, name="Time14.Dose20")
z = results(dds_meth_res, name="Time21.Dose20")


hist(log2(z$padj))

res_table <- results(dds_meth_res, contrast = contrast_ts, alpha = alpha_threshold)
dim(res_table)
class(res_table)

#look up how many genes are 
lookmeup = table(res_table$padj<0.01&abs(res_table$log2FoldChange)>2&res_table$baseMean>100)
lookmeup %>% data.frame() %>% View()

mcols(res_table, use.names=T)

library(dplyr)
res_table %>% data.frame() %>% View()




#Produce heatmap for interesting contrast

res_table_df <- as.data.frame(res_table)
res_table_df

res_table.df <- res_table_df[(res_table_df$baseMean > 60) & (abs(res_table_df$log2FoldChange) > 2), ]
dim(res_table.df)

mat <- counts(dds_meth, normalized=T)[rownames(res_table.df), ]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <-  colnames(mat)

library(ComplexHeatmap)

Heatmap(mat.z[1:10, ], cluster_rows = T, cluster_columns = T, 
        column_labels=colnames(mat.z), name="Z-score", row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 8))


Heatmap(mat.z, cluster_rows = diana(mat.z),
        cluster_columns = agnes(t(mat.z)),
        column_labels=colnames(mat.z), name="Z-score", row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 8))




h <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, 
             column_labels=colnames(mat.z), name="Z-score", 
             row_names_gp = gpar(fontsize = 2), 
             column_names_gp = gpar(fontsize = 5))

png('simple_heatmap_meth_Time_14_vs_3.png' , res=250, width= 3000, height=3000)
print(h)
dev.off()


####

#https://pubmed.ncbi.nlm.nih.gov/32124324/

# Plot gemes from Lauvas paper (left genes) -------------------------------
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Plot expression for single gene
# Nur genes

LL_genes <- c("CYP2E1", "GAP43", "MAPT", "NEUROD4", "NRG1", "NTRK2", "PLP1", "S100B",  "SOX9") # From Left panel
LL_genes <- c("CDK2", "CDKN1A", "CDKN1B", "CXCR4", "DDIT3", "E2F1", "NFKBIB", "PDK1",  "PFN1", "PLA2G12A", "RIT1", "TNFRSF1A")

#Genes from Lea
LL_genes <- genes_lists$Gene_Identifer[-c(1,4,10,25,28,31,43,50,51,54,56,60)]


# Markers for different genes Lauvas paper
neural_stem_cells -  c('NES','PAX6','SOX2')
astrocytes <- c('SLC1A2')
immature_neurones <- c('TBR1', 'GAP43')
committed_and_mature_neurons <- c('MAP2','MAPT','SYP')
GABAergic_neurons <-  c('SLC6A1','GABRA3','SLC32A1')
glutamatergic_neurons <- c('GRIN1', 'GRIN2b', 'GRIA1')
dopaminergic_neurons <- c('NR4A','TH') 
reference_genes <-  c('ACTB', 'GAPDH')


LL_genes_tab <- genes_lists[genes_lists$Functional_relation=="Synaptogenesis", ]
LL_genes_tab <- genes_lists[genes_lists$Functional_relation=="Neurite Outgrowth", ]
LL_genes_tab <- genes_lists[genes_lists$Functional_relation=="Metabolism Methadone", ]
LL_genes_tab <- genes_lists[genes_lists$Functional_relation=="Opioid receptor", ]

LL_genes_tab

#LL_genes <- LL_genes_tab$Gene_Identifer[-1] #synaptogenesis
LL_genes <- LL_genes_tab$Gene_Identifer #synaptogenesis




#setwd("C:/Users/maww/Downloads/meth_results/treatment_without_ribosomal_RNA_filtered_10/0_14/Nur_Right")

#Printing individual genes for loop up


LL_genes <- genes_lists$Gene_Identifer

for (i in 1:length(LL_genes)) {
  print(LL_genes[i])
  
  
  d <- plotCounts(dds_meth_res, gene=LL_genes[i], intgroup=c("Time", "Dose"), returnData=TRUE)
  
  
  g1 <- ggplot(d, aes(x = Time, y = log2(count), color = Time)) + 
    geom_point(position=position_jitter(w = 0.1,h = 0)) +
    geom_text_repel(aes(label = rownames(d))) + 
    theme_bw() +
    ggtitle(LL_genes[i]) +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  #png(file=paste0("Detailed_Expression_for_", LL_genes[i], ".png") , res=250, width= 3000, height=3000)
  pdf(file=paste0("Detailed_Expression_for_", LL_genes[i], ".pdf"))
  print(g1)
  dev.off()
  
  #box plot per given gene
  g2 <- ggplot(d, aes(Time, log2(count))) +
    geom_boxplot(aes(fill=factor(Dose))) + theme_bw() +
    theme(axis.text.x = element_text(angle=0, vjust=0.6)) +
    labs(title=LL_genes[i], 
         subtitle=" ",
         caption=" ",
         x="Development time",
         y="Log2 Normalized count")
  
  #png(file=paste0("Expression_for_", LL_genes[i], ".png") , res=250, width= 3000, height=3000)
  pdf(file=paste0("Expression_for_", LL_genes[i], ".pdf"))
  print(g2)
  dev.off()
  }
  
######IMPLEMENT MORE GENES

####Heatmap for chosen genes!!!

###start me


#OPEN HERE FOR annotation left
saveRDS(mat.z, file="mat.z.rds")
saveRDS(genes_lists, file="genes_lists.rds")


mat.z <- readRDS(file="mat.z.rds")
genes_lists <- readRDS(file="genes_lists.rds")


my_colour_time = c("DIV3" = "#B4ECC1", "DIV14" = "#62ED83","DIV21" = "#0BB835")
my_col_passage = c("9" = "#BED0F0", "10" = "#91B4F1","11" = "#538CF0",  "13" = "#1564EF")
my_col_dose = c("0" = "#FCFAFA", "0.07" = "#F7CDCD","2.2" = "#ED8B8B",  "6.7" = "#EF7272", "20" = "#EE3333")

#col_Zscore = colorRamp2(c(-5, 0, 5), c("#F5B7B1", "#FCF3CF", "#3498DB"))

genes_present = unique(rownames(mat.z[rownames(mat.z)%in%genes_lists$Gene_Identifer, ]))
genes_present
length(genes_present)


forheatmap <- unique(genes_lists[genes_lists$Gene_Identifer %in% genes_present, ]$Gene_Identifer)
length(unique(forheatmap)) #74 vs 83   I am finding now 56 genes

forleftannotation <- unique(genes_lists[genes_lists$Gene_Identifer %in% genes_present, ])

unique(forleftannotation$Functional_relation)


#https://htmlcolorcodes.com/

my_col_left_anno = c("Opioid receptor" = "#800080", #purple
                     "Shared LTP NMDA, Opioid receptor" = "#FF00FF", #fuchsia
                     "Reference Gene" = "#000080", #navy
                     "LTP NMDA" = '#0000FF', #blue
                     "Marker developmental status" = '#008080', #teal
                     "Marker differentiation process" = '#00FFFF',  #aqua
                     "Marker astrocytes" = '#008000', #green                
                     "Marker neural stem cells" = '#00FF00', #lime         
                     "Marker immature neurons" = '#808000', #olive
                     "Marker mature neurons"  = '#FFFF00', #yellow         
                     "Marker GABAergic neurons" = '#800000', #maroon     
                     "Marker glutamatergic neurons"  = '#FF0000', #red 
                     "Marker dopaminergic neurons" = '#000000', #black
                     "LTP Ca2+" = '#808080',  #gray                        
                     "Synaptogenesis" = "#C0C0C0", #silver                   
                     "Neurite Outgrowth" = "#FFFFFF") #white



arranged_colours = c("#800080", 
    '#008000',
    '#008000',
    '#808080',
    "#FF00FF",
    '#FFFF00',
    "#C0C0C0",
    "#FF00FF",
    "#C0C0C0",
    "#800080",
    "#FF00FF",
    '#00FF00',
    "#800080",
    '#808080',
    "#800080",
    '#00FF00',
    '#00FFFF',
    '#00FFFF',
    "#800080",
    "#C0C0C0",
    '#00FF00',
    '#808000',
    "#800080",
    '#008080',
    "#800080",
    "#FF00FF",
    '#FF0000',
    "#800080",
    "#800080",
    "#C0C0C0",
    '#808080',
    "#C0C0C0",
    "#800080",
    "#800080",
    "#C0C0C0",
    "#800080",
    '#0000FF',
    "#FF00FF",
    '#000000',
    '#800000',
    "#C0C0C0",
    "#C0C0C0",
    '#0000FF',
    "#FF00FF",
    "#C0C0C0",
    '#008000',
    '#808000',
    "#800080",
    "#FF00FF",
    '#0000FF',
    '#800000',
    "#C0C0C0",
    '#FF0000',
    "#800080",
    '#800000',
    "#C0C0C0")



forheatmap <- 

#Produce subset of matrix for genes that changes in time
hhh <- Heatmap(mat.z[rownames(mat.z)%in%forheatmap, ], cluster_rows = T, cluster_columns = T,
              column_labels=colnames(mat.z), name="Z-score", 
              row_names_gp = gpar(fontsize = 7, col = arranged_colours),
              top_annotation = HeatmapAnnotation(Time = condition_names, 
                                                 Dose = dose_names,
                                                 Passage= pasage_names,
                                                 col = list(Time = my_colour_time, 
                                                            Dose = my_col_dose, 
                                                            Passage=my_col_passage)),
              column_names_gp = gpar(fontsize = 3))

plot(hhh)
hhh@matrix  # use me 
rownames(hhh@matrix)

#Print png and pdf  ###make this more aoutomatized!   

png("20221021_subset.png", res=600, width= 7000, height=7000) #, res=400, width= 3000, height=3000)
plot(hhh) 
dev.off()

pdf("20221021_subset.pdf") #, res=400, width= 3000, height=3000)
plot(hhh) 
dev.off()

jpeg("20221021_subset.jpg", res=600, width= 7000, height=7000) #, res=400, width= 3000, height=3000)
plot(hhh) 
dev.off()



###Important genes for Lea
##CYP3A4, CYP2B6, UGT2B7

###Genes for Oddvar 

pdf(file=paste0("TEST2Heatmap_of_sig_genes_with_baseMean_cutoff_", baseMean_cutoff,
                "_padj_cf_", padj_cutoff, 
                "_log2FoldChange_cf_", log2FoldChange_cutoff,
                "_with_threshold_", threshold, 
                "_in_contrast_", contrast, "_des_treat_plus_passage.pdf"))


pdf(file=paste0("TEST2Heatmap_of_sig_genes_with_baseMean_cutoff_", baseMean_cutoff,
                "_padj_cf_", padj_cutoff, 
                "_log2FoldChange_cf_", log2FoldChange_cutoff,
                "_with_threshold_", threshold, 
                "_in_contrast_", contrast, "_des_treat_plus_passage.pdf"))
plot(hhh) 
dev.off()


###testing the annotation left

CYP3A4, CYP2B6, UGT2B7

hhh_left_test <- Heatmap(mat.z[rownames(mat.z)%in%forheatmap, ], cluster_rows = F, cluster_columns = T, row_km = 16,
               column_title =paste0("baseMean_cf", baseMean_cutoff,
                                    "_padj_cf_", padj_cutoff, 
                                    "_log2Fold_cf_", log2FoldChange_cutoff,
                                    "_thr_", threshold, 
                                    "_in_contr_", contrast),
               #col =col_Zscore, 
               column_labels=colnames(mat.z), name="Z-score", 
               #row_names_gp = gpar(fontsize = 5), #original
               #row_names_gp = gpar(fontsize = 5, col = as.factor(forleftannotation$Functional_relation)),  #wrong colours
               row_names_gp = gpar(fontsize = 5), 
               top_annotation = HeatmapAnnotation(Time = condition_names, 
                                                  Dose = dose_names,
                                                  Passage= pasage_names,
                                                  col = list(Time = my_colour_time, 
                                                             Dose = my_col_dose, 
                                                             Passage=my_col_passage)),
               #left_annotation = HeatmapAnnotation(Anno= as.matrix(forleftannotation[c(1,3)], 
               #                                   col = list(Anno = my_col_left_anno))),
               #bottom_annotation = HeatmapAnnotation(Dose = dose_names, col = list(Dose = my_col_dose)),
               column_names_gp = gpar(fontsize = 5)) #, col=treatment_colors


plot(hhh_left_test) 


###end me



# Plot expression for single gene (fix me)


Lea_genes = c("CYP3A4", "CYP2B6", "UGT2B7", "CYP1A1")

#dds_seq_subset <- dds_seq_subset[rownames(dds_seq)%in%genes_lists$Gene_Identifer, ]
dds_seq_subset <- dds_meth_res[rownames(dds_meth_res)%in%Lea_genes, ]



for (i in 1:3) {
  d <- plotCounts(dds_seq_subset, gene=rownames(dds_seq_subset)[i], intgroup=c("Treatment", "Time"), returnData=TRUE)
  d$Time <- ordered(d$Time, levels=c("3", "14", "21"))
  
  g <- ggplot(d, aes(x = Time, y = count, color = Time)) + 
    geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) +
    stat_summary(fun=mean, geom="point", shape=23, size=4) +
    theme_bw() +
    ggtitle(rownames(dds_seq_subset)[i]) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_brewer(palette="Dark2") +
    theme(legend.position=" ")
  pdf(file=paste0("gene_name_", rownames(dds_seq_subset)[i], ".pdf"))
  print(g)
  dev.off()
  png(file=paste0("gene_name_", rownames(dds_seq_subset)[i], ".png"))
  print(g)
  dev.off()
}





Oddvar_genes <- c('ESR1', 'ESR2', 'AR', 'THRB', 'GLCCI1', 'RARA', 'RARB', 'RARG', 'PGR', 'PPARA', 'PPARG', 'PPARD', 'VDR', 
            'NR1H3', 'PTGER2', 'PTGER1', 'PTGER3', 'PTGER4', 'PTGER4P1', 'PTGER4P2', 'PTGER4P3', 'PTGER4P2-CDK2AP2P2')




dds_meth_res


#Produce subset of matrix for genes that changes in time
hhhh <- Heatmap(mat.z[rownames(mat.z)%in%Oddvar_genes, ], cluster_rows = F, cluster_columns = T,
               column_labels=colnames(mat.z), name="Z-score", 
               row_names_gp = gpar(fontsize = 7),
               top_annotation = HeatmapAnnotation(Time = condition_names, 
                                                  Dose = dose_names,
                                                  Passage= pasage_names,
                                                  col = list(Time = my_colour_time, 
                                                             Dose = my_col_dose, 
                                                             Passage=my_col_passage)),
               column_names_gp = gpar(fontsize = 3))

plot(hhhh)




pdf(file="Oddvar.pdf")
print(hhhh)
dev.off()


png(file="Oddvar.png")
print(hhhh)
dev.off()


jpeg(file="Oddvar.jpeg")
print(hhhh)
dev.off()





library(ComplexHeatmap)
matOdd <- counts(dds_meth_res, normalized=T)
matOdd_tmp <- matOdd[rownames(matOdd)%in%Oddvar_genes, ]
dim(matOdd_tmp)

mat.zOdd.XX <- t(apply(matOdd_tmp, 1, scale))  #head(t(scale(t(mat)))) scaling on the rows
head(mat.zOdd.XX)
colnames(mat.zOdd.XX) <-  colnames(matOdd)


#Produce subset of matrix for genes that changes in time
hhhhh <- Heatmap(mat.zOdd.XX, cluster_rows = F, cluster_columns = T,
                column_labels=colnames(mat.z), name="Z-score", 
                row_names_gp = gpar(fontsize = 7),
                top_annotation = HeatmapAnnotation(Time = condition_names, 
                                                   Dose = dose_names,
                                                   Passage= pasage_names,
                                                   col = list(Time = my_colour_time, 
                                                              Dose = my_col_dose, 
                                                              Passage=my_col_passage)),
                column_names_gp = gpar(fontsize = 3))

plot(hhhhh)












###Next contrast


#Building the results table
alpha_threshold = 0.05
contrast_ts <- c("Dose", "20", "0")
res_table <- results(dds_meth_res, contrast = contrast_ts, alpha = alpha_threshold)
dim(res_table)
class(res_table)
mcols(res_table, use.names=T)

library(dplyr)
res_table %>% data.frame() %>% View()




###


res_table_df <- as.data.frame(res_table)
res_table_df

res_table.df <- res_table_df[(res_table_df$baseMean > 60) & (abs(res_table_df$log2FoldChange) > 1), ]
dim(res_table.df)

mat <- counts(dds_meth, normalized=T)[rownames(res_table.df), ]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <-  colnames(mat)

library(ComplexHeatmap)

Heatmap(mat.z[1:5, ], cluster_rows = T, cluster_columns = T, 
        column_labels=colnames(mat.z), name="Z-score", row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 8))

h <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, 
             column_labels=colnames(mat.z), name="Z-score", 
             row_names_gp = gpar(fontsize = 8), 
             column_names_gp = gpar(fontsize = 8))
png('simple_heatmap_meth_Time_14_vs_3.png' , res=250, width= 3000, height=3000)
print(h)
dev.off()


####
# Plot expression for single gene
d <- plotCounts(dds_meth_res, gene="PRPH", intgroup="Time", returnData=TRUE) 

# What is the data output of plotCounts()?
d %>% View()

# Plot the PRPH normalized counts, using the samplenames (rownames(d) as labels)
ggplot(d, aes(x = Time, y = count, color = Time)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("PRPH") +
  theme(plot.title = element_text(hjust = 0.5))

