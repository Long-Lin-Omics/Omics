library(DiffBind)
# setwd(system.file('extra',package='DiffBind'))

# tmpdir <- tempdir()
# url <- 'https://content.cruk.cam.ac.uk/bioinformatics/software/DiffBind/DiffBind_vignette_data.tar.gz'
# file <- basename(url)
# options(timeout=600)
# download.file(url, file.path(tmpdir,file))
# untar(file.path(tmpdir,file), exdir = tmpdir )
# setwd(file.path(tmpdir,"DiffBind_Vignette"))
setwd("DiffBind_Vignette")

#in a single step
tamoxifen <- dba.analyze("tamoxifen.csv")
tamoxifen.DB <- dba.report(tamoxifen)
        #header in tamoxifen.csv
        #SampleID Tissue Factor Condition Treatment Replicate bamReads ControlID bamControl Peaks PeakCaller

#seperately
tamoxifen <- dba(sampleSheet="tamoxifen.csv") %>%
    # insteadly
    # samples <- read.csv("tamoxifen.csv")
    # dba(sampleSheet=samples)
+ dba.blacklist() %>%
+ dba.count() %>%
+ dba.normalize() %>%
+ dba.contrast() %>%
+ dba.analyze()


#seperately in details
tamoxifen <- dba(sampleSheet="tamoxifen.csv") 
    # tamoxifen to show the metadata with this DBA object
    # * Samples, * overlapping peaks
    # header: ID Tissue Factor Condition Treatment Replicate Intervals
plot(tamoxifen) #or dba.plotHeatmap(tamoxifen) #heatmap with clusters
    #dba.plotHeatmap(tamoxifen,score=DBA_SCORE_RPKM_FOLD) # change the default normalisez read counts DBA_SCORE_NORMALIZED. 
# wait 
#tamoxifen = dba.blacklist()
tamoxifen <- dba.count(tamoxifen)
    # tamoxifen to show the metadata with this DBA object
    # * Samples, * overlapping peaks
    # header: ID Tissue Factor Condition Treatment Replicate Reads FRiP
plot(tamoxifen)
    #higher level of correlation using read counts matrix
tamoxifen <- dba.normalize(tamoxifen) 
    # by default, sequencing depth 
    # for details in normalization process
    # norm <- dba.normalize(tamoxifen, bRetrieve=TRUE)
    # The default library-size based methods results in all the library sizes being normalized to be the same (the mean library size): lib.size/norm.factors
tamoxifen <- dba.contrast(tamoxifen, reorderMeta=list(Condition="Responsive"))
    # add a design to the metadata
tamoxifen <- dba.analyze(tamoxifen)
    # by default: DESeq2
dba.show(tamoxifen, bContrasts=TRUE)
    # show total number of significantly differentially bound sites with FDR <= 5%
plot(tamoxifen, contrast=1)
    # heatmap with cluster using only differentially bound sites
    # Using only the differentially bound sites, we now can see the clusters to reflect the condition. 
    # This is because, all sites together, the number of random differentially bound sites, will dilute the difference between differentially bound sites (2845 vs 246).
readscores <- dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE, scale="row", colScheme = hmap)
    #heatmap between sites
tamoxifen.DB <- dba.report(tamoxifen)
    # differentially bound sites
    # header: seqnames ranges strand Conc Conc_Resistant Conc_Responsive Fold p-value FDR
    # Conc: mean read concentration over all the samples (the default calculation uses log2 normalized read counts)
dba.plotVenn(tamoxifen, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
    #sum(tamoxifen.DB$Fold>0)
    #sum(tamoxifen.DB$Fold<0)
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_CONDITION)
    # what is DBA_TISSUE??
dba.plotPCA(tamoxifen, contrast=1, label=DBA_TISSUE)
    # using only differentailly bound sites
dba.plotPCA(tamoxifen, attributes=c(DBA_TISSUE, DBA_CONDITION), label=DBA_REPLICATE)
dba.plotPCA(tamoxifen,contrast=1,b3D=TRUE)
dba.plotMA(tamoxifen)
dba.plotMA(tamoxifen, bXY=TRUE)
dba.plotVolcano(tamoxifen)

