library("DESeq2")
data <- read.delim("firre_tfea_30_vs_0/temp_files/count_file.header.bed", sep="	", header=TRUE)
countsTable <- subset(data, select=c(5, 6, 7, 8, 9, 10))

rownames(countsTable) <- data$region
cond_vector <- c("0min", "0min", "0min", "30min", "30min", "30min")
batch <- c()
if (length(batch) == 0) {
    conds <- data.frame(cond_vector)
    colnames(conds) <- c("treatment")
    ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countsTable, 
                                                colData = conds, 
                                                design = ~ treatment)
} else {
    conds <- data.frame(cond_vector, batch)
    colnames(conds) <- c("treatment", "batch")
    ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countsTable, 
                                                colData = conds, 
                                                design = ~ batch+treatment)
}

dds <- DESeq(ddsFullCountTable)
print("Size Factors")
print(sizeFactors(dds))
res <- results(dds, alpha = 0.05, contrast=c("treatment", "30min",
                                                            "0min"))
res$fc <- 2^(res$log2FoldChange)
res <- res[c(1:3,7,4:6)]

write.table(res, file = "firre_tfea_30_vs_0/temp_files/DESeq.res.txt", append = FALSE, sep= "	" )
sink()