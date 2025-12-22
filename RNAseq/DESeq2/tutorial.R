
# quick and basic
dds <- DESeqDataSetFromMatrix(countData = cts,                      # matrix: rows are GENES, coloumns are samples  #tximport can aggregate counts from transcript level to gene level.
                              colData = coldata,                    # dataframe: Rows of colData correspond to columns ofcountData # row.names should be sample ID. columns for design should be factors.
                                                                    # all(rownames(coldata) == colnames(cts)) should be True
                              design= ~ batch + condition)          # columns of coldata   
                                                                        # put the variable of interest at the end of the formula 
                                                                        # make sure the control level is the first level: coldata$condition <- relevel(coldata$condition, ref = "untreated")
                                                                                        
        
        #library(tximport)
        # counts.imported <- tximport(files = as.character(file.list), type = 'salmon', tx2gene = tx2gene)      # salmon output: Name   Length  EffectiveLength TPM NumReads
                    # tx2gene is required to aggregate transcripts into genes
                    # counts.imported is a list
                    # counts.imported$abundance -> tpm from salmon            counts.imported$length  -> EffectiveLength from salmon
                    # counts.imported$counts -> NumReads from salmon              counts.imported$countsFromAbundance  "No" a tag
                    # Aggregation from transcripts: 
                    #       --> abundacne: directly sum 
                    #       --> length: effectivelength weight by tpm
                    #       --> counts: directly sum by default, as "No" in countsFromAbundance. 
                    #       check countsFromAbundance options later
        # dds <- DESeqDataSetFromTximport(counts.imported, colData = samples, design = ~ condition)

        # other ways importting the data
        # tmimeta + DESeqDataSet()
        # DESeqDataSetFromHTSeq()
        # featureCounts?

# optional pre-filtering
    dds
    # 1) keep genes that have at least m reads among samples.
    gene_reads_thred <- 10
    keep <- rowSums(counts(dds)) > gene_reads_thred
    dds <- dds[keep, ]
    # 2) keep geens that have at least m reads for each of n samples
    gene_reads_thred <- 10
    smallestGroupSize <- 3
    keep <- rowSums(counts(dds) >= gene_reads_thred) >= smallestGroupSize
    dds <- dds[keep,]
    # see the dimension change
    dds

dds <- DESeq(dds)   # do the DEG analysis and produce data for plotting
    # dds[,] can be accessed by range, but can't be seen directly, use methods/fucntions below instead
    # mcols(dds)  metadata associated with genes, produced by DESeq2, can also be extended by user: mcols(dds) <- DataFrame(mcols(dds), featureData)
    # results(dds) # baseMean log2FoldChange     lfcSE      stat      pvalue        padj  #  mcols(dds)$description to get the description for each column 
        # use name='' or contrst=c("condition","treated","untreated"), to specify the coefficient, default is the last variable in the design formula and the last level of this variable over the reference level . 
        # contract is more advanced. 
        # other usage of results()
    # counts(dds)
    # colData(dds)
    # plotDispEsts(dds)
    # rld <- rlog(dds, blind = TRUE)
    # rld_mat <- assay(rld)
    # pca <- prcomp(t(rld_mat))
    # plot(pca$x[,1],pca$x[,2],col=dds$condition)
resultsNames(dds) # lists the coefficients                          # "Intercept", Comparison based on design
res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")  # change lfc but not p values
resOrdered <- res[order(res$pvalue),]
summary(res)  # break down as from sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05) # break down as from sum(res05$padj < 0.05, na.rm=TRUE)


plotMA(res, ylim=c(-2,2)) # Points will be colored blue if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotCounts(dds, gene=which.min(res$padj), intgroup="condition") # normalizes counts by the estimated size factors (or normalization factors if these were used) and adds a pseudocount of 1/2 to allow for log scale plotting


#exporting results
write.csv(as.data.frame(resOrdered), file="condition_treated_results.csv")
write.csv(as.data.frame(subset(resOrdered, padj < 0.1)), file="condition_treated_results.padj0.1.csv")
