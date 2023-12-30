install.packages("umap")
library(umap)
library(ggplot2)
library(ggpubr)  
library(dplyr)

#############################

exp <- read.delim("D:/BulkRNA/Data/data_mrna_seq_v2_rsem.txt" )
exp<- na.omit(exp)
exp<- data.frame(t(exp))

colnames(exp) <- exp[ 1,]
exp<- exp[-1:-2, ]


################PatientInfo 

PatientInfo <- read.delim("D:/BulkRNA/Data/tcga_clinical_data.tsv")
PatientInfo$Sample.ID <- gsub("-" ,"." , PatientInfo$Sample.ID)

## Defining a new column named Er.status
exp$ER.status <- PatientInfo$ER.Status.By.IHC[match(rownames(exp), PatientInfo$Sample.ID)]
exp<- na.omit(exp)
exp<- exp[exp$ER.status!="Indeterminate",]

 
log_exp <- as.data.frame(lapply(exp[,-20531], as.numeric))
#log-transform the expression data
# Convert to numeric and then apply log2 transformation
log_exp<-as.data.frame(log2(log_exp+1)) # Adding 1 to avoid log(0)
#dim(log_exp)
#log_exp <- t(scale(t(log_exp)))


 
######################################## UMAP

umap_result <- umap(log_exp, n_neighbors =120)
plot(umap_result$layout, main="UMAP Plot")


# Assuming 'labels' is a vector of sample labels (e.g., cancer types)
labels <- exp$ER.status
labels <- as.factor(labels)


plot(umap_result$layout, col = as.numeric(labels), pch = 16, main = "UMAP Plot")
# Add a legend if labels are factors
if (is.factor(labels)) {
  legend("topright", legend = levels(labels), col=1:length(levels(labels)), pch=16, cex=0.65)
}
