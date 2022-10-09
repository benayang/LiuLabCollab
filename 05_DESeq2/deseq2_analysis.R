library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(tximport)
library(EnsDb.Hsapiens.v75)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(sva)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)

projdir <- "/nas/homes/benyang/LiuLabCollab"

# Prepare input table ---------------------
sampleTable <- data.frame(path = list.files(file.path(projdir, "03_Kallisto"), full.names=T))
sampleTable <- sampleTable %>% 
    mutate(sample = unlist(lapply(basename(path), function(x) unlist(strsplit(x, split="_"))[1])),
            pressure = factor(rep(c("0Pa","600Pa","1200Pa"), each=2), levels=c("0Pa","600Pa","1200Pa")),
            BNP = factor(rep(c("minus","plus"), times=3), levels=c("minus","plus")),
            group = factor(paste(BNP, pressure, sep="_")))
rownames(sampleTable) <- sampleTable$sample
saveRDS(sampleTable, file.path(projdir, "05_DESeq2", "sampleTable.RDS"))

# Create DESeq2 object ---------------------
# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto

edb = EnsDb.Hsapiens.v75
k <- keys(edb, keytype = "TXNAME")
tx2gene = genes(edb, filter = TxIdFilter(k), columns = c("tx_id","gene_id","gene_name"), return.type="data.frame")

files <- file.path(sampleTable$path, "abundance.h5")
names(files) <- sampleTable$sample
txi <- tximport(files, type="kallisto", ignoreTxVersion=T, tx2gene=tx2gene)
saveRDS(txi, file.path(projdir, "05_DESeq2", "txi_kallisto.RDS"))

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~pressure)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
dds <- DESeq(dds)

saveRDS(dds, file.path(projdir, "05_DESeq2", "pressure_dds.RDS"))

# Sample QC and comparisons ---------------------
for (i in c("spearman","pearson")) {
    dds_cor = cor(counts(dds, normalized=T), method=i)
    rownames(dds_cor) <- paste(dds$sample, sep="-")
    colnames(dds_cor) <- NULL
    png(file.path(projdir, "05_DESeq2", "QC", paste0(i,"_correlation.png")), units="in", width=7, height=5, res=300)
    pheatmap::pheatmap(dds_cor, cluster_rows = F, cluster_cols = F, display_numbers = T, main=sprintf("dds %s Correlation",i))
    dev.off()
}

# Variance stabilization for visualization ---------------------
library("vsn")
meanSdPlot(log(counts(dds)+1), ranks = FALSE)

rld <- rlog(dds, blind = FALSE)
saveRDS(rld, file.path(projdir, "05_DESeq2", "pressure_rld.RDS"))


# Custom PCA visualization ---------------------
plotPCACustom = function(rld, color, shape){
  rv = rowVars(assay(rld))
  select = order(rv, decreasing=T)[seq_len(min(500, length(rv)))]
  rld_PCA = prcomp(t(assay(rld)[select,]), center=T, scale=F) # clustering variables not genes, don't scale as rlog already scaled
  # https://support.bioconductor.org/p/100069/
  summary(rld_PCA)
  percent_var = (rld_PCA$sdev^2 / sum(rld_PCA$sdev^2)) * 100
  rld_PCA_df = data.frame(rld_PCA$x, sampleTable)
  dds_pca = ggscatter(rld_PCA_df, x="PC1", y="PC2",
                      color=color, shape=shape,
                      palette = "lancet",
                      size=3,
                      legend = "right") +
    coord_fixed() +
    labs(x = paste0("PC1: ", round(percent_var[1],1), "% variance"), y = paste0("PC2: ", round(percent_var[2],1), "% variance")) +
    theme(text = element_text(family="Arial"), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12)) +
    theme_bw() + labs_pubr()
  return(dds_pca)
}

#plotPCA(rld, intgroup = c("pressure", "BNP"))
plotPCACustom(rld, "pressure", "BNP")
ggsave(file.path(projdir, "05_DESeq2", "QC", "rld_pca.png"), dpi=300, width=5, height=5)

# Heatmaps of count matrix ---------------------
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("pressure","BNP")])
count_hmp = assay(rld)[select,]
rownames(count_hmp) = genes(edb, columns="symbol", filter=GeneIdFilter(rownames(count_hmp)), return.type="data.frame")$symbol
#colnames(count_hmp) = paste(df$condition, df$day, c("R1","R2"), sep="-")
png(file.path(projdir, "05_DESeq2", "QC", "dds_rld_heatmap_first100.png"), units="in", width=7, height=15, res=300)
pheatmap::pheatmap(count_hmp, cluster_rows=TRUE, show_rownames=TRUE, scale="row", cellwidth = 15,
         cluster_cols=F, annotation_col=df, colorRampPalette(rev(brewer.pal(11, "RdBu")))(255))
dev.off()

# Group contrasts ---------------------
resultsNames(dds)
# [1] "Intercept"              "pressure_600Pa_vs_0Pa"  "pressure_1200Pa_vs_0Pa"

BNP_600Pa_vs_0Pa <- results(dds, name = "pressure_600Pa_vs_0Pa", pAdjustMethod = "fdr")
table(BNP_600Pa_vs_0Pa$padj<0.01)
BNP_600Pa_vs_0Pa.df <- BNP_600Pa_vs_0Pa %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(tx2gene %>% dplyr::select(gene_id, gene_name) %>% distinct(), by="gene_id")
write.table(BNP_600Pa_vs_0Pa.df, file.path(projdir, "05_DESeq2", "pressure_600Pa_vs_0Pa.csv"), sep=",", row.names=F, quote=F)

BNP_1200Pa_vs_0Pa <- results(dds, name = "pressure_1200Pa_vs_0Pa", pAdjustMethod = "fdr")
table(BNP_1200Pa_vs_0Pa$padj<0.01)
BNP_1200Pa_vs_0Pa.df <- BNP_1200Pa_vs_0Pa %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(tx2gene %>% dplyr::select(gene_id, gene_name) %>% distinct(), by="gene_id")
write.table(BNP_1200Pa_vs_0Pa.df, file.path(projdir, "05_DESeq2", "pressure_1200Pa_vs_0Pa.csv"), sep=",", row.names=F, quote=F)

make_rld_heatmap <- function(rld, res, tx2gene, genenames) {
    select <- rownames(res)[which(res$padj<0.01)]
    hmp_df <- assay(rld)[select,]
    rownames(hmp_df) <- tx2gene$gene_name[match(rownames(hmp_df), tx2gene$gene_id)]
    pheatmap::pheatmap(hmp_df, cluster_rows=TRUE, show_rownames=genenames, scale="row", cellwidth = 15,
         cluster_cols=F, annotation_col=df, colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255))
}

png(file.path(projdir, "05_DESeq2", "Plots", "dds_rld_heatmap_600Pa_vs_0Pa.png"), units="in", width=7, height=10, res=300)
make_rld_heatmap(rld, BNP_600Pa_vs_0Pa, tx2gene, F)
dev.off()

png(file.path(projdir, "05_DESeq2", "Plots", "dds_rld_heatmap_1200Pa_vs_0Pa.png"), units="in", width=7, height=10, res=300)
make_rld_heatmap(rld, BNP_1200Pa_vs_0Pa, tx2gene, T)
dev.off()