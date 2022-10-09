library(sleuth)
library(dplyr)
library(tidyr)
library(ggplot2)
library(biomaRt)

projdir <- "/nas/homes/benyang/LiuLabCollab"

sampleTable <- data.frame(path = list.files(file.path(projdir, "03_Kallisto"), full.names=T))
sampleTable <- sampleTable %>% 
    mutate(sample = unlist(lapply(basename(path), function(x) unlist(strsplit(x, split="_"))[1])),
            pressure = factor(rep(c("0Pa","600Pa","1200Pa"), each=2), levels=c("0Pa","600Pa","1200Pa")),
            BNP = factor(rep(c("minus","plus"), times=3), levels=c("minus","plus")),
            group = factor(paste(BNP, pressure, sep="_")))
saveRDS(sampleTable, file.path(projdir, "04_sleuth", "sampleTable.RDS"))

mart <- useMart(biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = 'https://useast.ensembl.org')
t2g <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "description"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
saveRDS(t2g, file.path(projdir, "04_sleuth", "t2g.RDS"))

so <- sleuth_prep(
    sample_to_covariates = sampleTable, 
    full_model = ~ group,
    target_mapping = t2g, 
    aggregation_column = "ens_gene",
    num_cores = 30,
    read_bootstrap_tpm = TRUE,
    extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, fit_name = "full")
models(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
saveRDS(so, file.path(projdir, "04_sleuth", "sleuth_obj_reduced_full.RDS"))

sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
table(sleuth_table_gene$qval <= 0.05)
write.table(sleuth_table_gene, file.path(projdir, "04_sleuth", "sleuth_table_gene_reduced_full_lrt.csv"), sep=',', row.names=F, quote=F)

so <- sleuth_fit(so, ~ BNP + pressure, fit_name = "full_flipped")
so <- sleuth_fit(so, ~ pressure, fit_name = "reduced_flipped")
so <- sleuth_lrt(so, 'reduced_flipped', 'full_flipped')
sleuth_table_gene_flipped <- sleuth_results(so, 'reduced_flipped:full_flipped', 'lrt', show_all = FALSE)
table(sleuth_table_gene_flipped$qval <= 0.05)

so$sample_to_covariates <- so$sample_to_covariates %>% mutate(group = paste(BNP, pressure, sep="_"))
so <- sleuth_fit(so, ~ group, fit_name = "full_group")
so <- sleuth_fit(so, ~ 1, fit_name = "reduced_intercept")

#oe <- sleuth_wt(so, which_beta = 'pressure600Pa')

plot_pca(so, text_labels = T, color_by = 'BNP')
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
  trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
  "sample"), offset = 1)
