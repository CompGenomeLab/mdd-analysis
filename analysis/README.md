# Downregulation in Glutamatergic Signaling Decreases Synaptic Plasticity through NPAS4 in Major Depression

## Differential Gene Expression (DEG) Analysis

- Import sample table as data.frame

```R
sampleTable <- data.frame(sampleName = sampleName,fileName = sampleFiles,condition = sampleCondition, sampleGender=sampleGender, sampleAge=sampleAge, PMI = PMI, brainRegion = brainRegion, sampleDataset = sampleDataset)
View(sampleTable)
```

| sampleName | fileName                         | condition | sampleGender | sampleAge | PMI  | brainRegion | sampleDataset |
| ---------- | -------------------------------- | --------- | ------------ | --------- | ---- | ----------- | ------------- |
| SRR5831944 | SRR5831944_count_control_1_m.txt | CTRL      | male         | 17        | 14   | dlpfc       | first         |
| SRR5831945 | SRR5831945_count_control_1_f.txt | CTRL      | female       | 30        | 13   | dlpfc       | first         |
| SRR5831946 | SRR5831946_count_control_1_m.txt | CTRL      | male         | 66        | 8    | dlpfc       | first         |
| SRR3438661 | SRR3438661_count.txt             | CTRL      | male         | 32        | 26   | nacc        | second        |
| SRR3438672 | SRR3438672_count.txt             | CTRL      | male         | 66        | 18   | nacc        | second        |
| SRR3438673 | SRR3438673_count.txt             | CTRL      | male         | 49        | 27.5 | nacc        | second        |

```R

sampleTable$condition <- factor(sampleTable$condition)
sampleTable$sampleDataset <- factor(sampleTable$sampleDataset)
sampleTable$sampleGender <- factor(sampleTable$sampleGender)
sampleTable$brainRegion <- factor(sampleTable$brainRegion)
library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory,design= ~ sampleDataset + sampleGender + sampleAge + brainRegion + PMI + condition)
dds <- DESeq(dds)
results <- results(dds)

contrast <- c("condition","MDD","CTRL")
results(dds, contrast = contrast, alpha = 0.05)
res_tableOE_unshrunken <- results(dds, contrast=contrast, alpha = 0.05)
res_tableOE <- lfcShrink(dds, contrast=contrast, type = "normal",  res=res_tableOE_unshrunken)
padj.cutoff <- 0.05
res_tableOE_df <- data.frame(res_tableOE)
res_tableOE_df$ensembl_gene_id <- rownames(res_tableOE_df)
rownames(res_tableOE_df) <- NULL
sig <- subset(res_tableOE_df, padj <= padj.cutoff)

FUN <- function(x){
x$oe = ifelse(x$padj <= 0.05, TRUE, FALSE)
return(x)
}
library("dplyr")
library("ggplot2")
res_tableOE_volcano <- res_tableOE_df %>% mutate(oe = padj < 0.05)
ggplot(res_tableOE_volcano) +
geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = oe)) +
ggtitle("Brain Regions") +
xlab("log2 fold change") +
ylab("-log10 adjusted p-value") +
theme(legend.position = "none",
plot.title = element_text(size = rel(1.5), hjust = 0.5),
axis.title = element_text(size = rel(1.25)))
write.csv(sig, "C:/Users/haksu/Desktop/Adebali Lab/Brain Regions/Brain Regions + PMI\\sig-batch(gender+age+PMI).csv")
```

## Co-expression Analysis

- Import sample annotation list

  | SampleName | Class |
  | ---------- | ----- |
  | SRR5831944 | CTRL  |
  | SRR5831945 | CTRL  |
  | SRR5831958 | MDD   |
  | SRR5831959 | MDD   |

```R
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts)
library("CEMiTool")
cem <- cemitool(normalized_counts, sample_annot)
generate_report(cem)
write_files(cem)
```

