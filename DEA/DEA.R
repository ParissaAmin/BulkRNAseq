
library(pheatmap) 
library(ggplot2)
library(reshape2)  
library(ggpubr)  
library(dplyr)
library(DESeq2)
library(ggrepel)
#############################

exp <- read.delim("D:/BulkRNA/Data/data_mrna_seq_v2_rsem.txt")
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

remove(PatientInfo)

#############Creating a metadata table with information of our samples 

metadata<- as.data.frame(cbind(rownames(exp), exp[,"ER.status"]))
colnames(metadata)<- c("Sample.ID", "ER")

#####Deleting last column(the genes' names) of exp
expression_data<-  exp[,-20531]
expression_data <- data.frame(lapply(expression_data, as.numeric))
expression_data<-as.data.frame(log2(expression_data+1)) # Adding 1 to avoid log(0)

expression_data<- round(expression_data)
expression_data<- as.data.frame(t(expression_data))

remove(exp)
######################################Create dds object by DESeqDataSet function

dds<-DESeqDataSetFromMatrix(countData=expression_data,colData=metadata,design=~ER)

dds <- DESeq(dds)

My_Results<- results(dds)
My_Results_df<- as.data.frame(My_Results)

My_Results_df$GeneName <- rownames(My_Results_df)

#################Set the significance threshold
padj_threshold <- 0.05
logFC_threshold <- 2

#library(ggrepel)
pdf("D:/BulkRNA/DEA_ER_Stat-Pari.pdf", width = 12, height = 9)
ggplot(My_Results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size=2, aes(color = ifelse(abs(log2FoldChange) > logFC_threshold & padj < padj_threshold, "red", "grey")), alpha = 0.6) +
  geom_text_repel(
    aes(label = ifelse(abs(log2FoldChange) > logFC_threshold & padj < padj_threshold, GeneName, "")),
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.2,
    segment.color = "grey50",
    box.color = "grey50"
  ) +
  scale_color_identity() +
  theme_bw() +
  labs(
    title = "Volcano Plot for DEA ER Positive vs ER Negative in TCGA",
    x = "log2 Fold Change",
    y = "-log10(padj)",
    color = "Significant"
  )

dev.off()

















