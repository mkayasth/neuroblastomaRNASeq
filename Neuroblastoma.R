# Loading libraries and Data.
library(edgeR)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggrepel)

TARGET_NBL_htseq_counts <- read.delim("TARGET-NBL.htseq_counts.tsv", header=TRUE)

# removing Ensembl ID data after .
TARGET_NBL_htseq_counts$Ensembl_ID <- gsub("\\..*", "", TARGET_NBL_htseq_counts$Ensembl_ID)

mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version=112)

# Retrieving the Ensembl gene IDs and gene biotype (protein coding)
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
                              filters="biotype",
                              values = "protein_coding",
                              mart = mart)

# Filter TARGET_NBL_htseq_counts to only include protein-coding genes
TARGET_NBL_htseq_counts <- TARGET_NBL_htseq_counts %>%
  filter(Ensembl_ID %in% protein_coding_genes$ensembl_gene_id)


Neuroblastoma_Metadata <- read.table(file = 'Neuroblastoma_Metadata.txt', header = TRUE, sep = '\t')

# Data filtering: only including ALT and Telomerase phenotype. 
metadata_filtered_TMM <- Neuroblastoma_Metadata %>%
  filter(TMM %in% c("ALT", "Telomerase"))

# Data filtering: only include sample IDs present in the other dataset.
metadata_filtered_TMM <- metadata_filtered_TMM %>%
  filter(SampleID %in% colnames(TARGET_NBL_htseq_counts))

rownames(TARGET_NBL_htseq_counts) <- TARGET_NBL_htseq_counts$Ensembl_ID
TARGET_NBL_htseq_counts$Ensembl_ID <- NULL

# Arranging Sample_ID in the same order in both datasets:
# first arranging ALT then Telomerase in the metadata.
metadata_filtered_TMM <- metadata_filtered_TMM %>%
  arrange(factor(TMM, levels = c("ALT", "Telomerase")))

# counts_TMM is a subset of TARGET_NBL_htseq_counts with only people from metadata.
counts_TMM <- TARGET_NBL_htseq_counts[,metadata_filtered_TMM$SampleID]
counts_TMM_backup <- cbind(Ensembl_ID = rownames(counts_TMM), counts_TMM)



# building model matrix.

# First, determining the factors of TMM.
group1 <- as.factor(metadata_filtered_TMM$TMM)

# model matrix ~ without an intercept term.
design <- model.matrix(~group1+0)
design


# creating differential gene expression object.
dge_TMM <- DGEList(counts=counts_TMM,group=group1)

# TMM normalization. Transformation via counts per million function.
dge_TMM <- calcNormFactors(dge_TMM)


# TMM normalization after removing lowly expressed genes with cpm < 1.
keep <- rowSums(cpm(dge_TMM) > 1) >= ceiling(0.05*dim(d)[2])
dge_TMM2 <- dge_TMM[keep, ]

# new DGEList with new count.
dge_TMM3 <- DGEList(counts = dge_TMM2, group = group1)

# recalculating normalization factors.
dge_TMM3 <- calcNormFactors(dge_TMM3)
dge_TMM3_backup <- cbind(Ensembl_ID = rownames(dge_TMM3), dge_TMM3$counts)

# Calculating dispersion and fitting the model.
d <- estimateDisp(dge_TMM3, design, verbose=TRUE)
fit <- glmQLFit(d, design)

# contrast parameter (ALT-Telomerase).
contrast <- makeContrasts(group1ALT-group1Telomerase, levels=design)

# differential expression test.
fit2 <- glmQLFTest(fit, contrast = contrast)

# Adjusted p-value (False discovery rate correction.)
fit2$table$fdr <- p.adjust(fit2$table$PValue, method ="BH")


#Annotating Ensembl -> GeneSymbol.
mart <- useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl', version = 112)

annotation <- getBM(filters = 'ensembl_gene_id',
                    attributes= c("ensembl_gene_id",
                    "hgnc_symbol"),
                    values = rownames(fit2$table),
                mart = mart)


# final dataframe.
fit2$table$ensembl_gene_id <- rownames(fit2$table)

final_version <- fit2$table %>%
  left_join(annotation, by ='ensembl_gene_id')




# filtering for candidate genes.
candidate_genes_ALT <- final_version[c(final_version$fdr <= 0.05 & final_version$logFC >= 0.5),]
candidate_genes_ALT <- candidate_genes_ALT %>%
  filter(hgnc_symbol != "" & !is.na(hgnc_symbol))

candidate_genes2_ALT <- final_version[c(final_version$fdr <= 0.01 & final_version$logFC >= 1.5),]



candidate_genes_Telomerase <- final_version[c(final_version$fdr <= 0.1 & final_version$logFC <= -0.5),]
candidate_genes_Telomerase <- candidate_genes_Telomerase %>%
  filter(hgnc_symbol != "" & !is.na(hgnc_symbol))

candidate_genes2_Telomerase <- final_version[c(final_version$fdr <= 0.01 & final_version$logFC <= -1.5),]

merged_candidates <- rbind(candidate_genes_ALT, candidate_genes_Telomerase)

# Top 10 Genes in terms of fold change.
ALT_significant <- merged_candidates %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 8)

Telomerase_significant <- merged_candidates %>%
  arrange(logFC) %>%
  slice_head(n=8)


final_version <- final_version %>%
  mutate(gene_status = case_when(
    hgnc_symbol %in% ALT_significant$hgnc_symbol ~ "ALT Significant",
    hgnc_symbol %in% Telomerase_significant$hgnc_symbol ~ "Telomerase Significant"))

merged_candidates <- merged_candidates %>%
  mutate(gene_status = case_when(
    hgnc_symbol %in% ALT_significant$hgnc_symbol ~ "ALT Significant",
    hgnc_symbol %in% Telomerase_significant$hgnc_symbol ~ "Telomerase Significant"))

#Volcano plot ~ all samples.
ggplot(data = final_version, aes(x = logFC, y = -log10(fdr), color = gene_status)) +
  geom_point() + 
  scale_color_manual(values = c("ALT Significant" = "red", "Telomerase Significant" = "blue")) +
  theme_minimal() + geom_text_repel(data = bind_rows(ALT_significant, Telomerase_significant), 
                  aes(label = hgnc_symbol), 
                  vjust = 0.5, hjust = 0.5, size = 3,
                  color = "black", box.padding = 0.5,
                  point.padding = 0.5, max.overlaps = Inf
                  )


# Volcano plot ~ candidate samples.
ggplot(data = merged_candidates, aes(x = logFC, y = -log10(fdr), color = gene_status)) +
  geom_point() + 
  scale_color_manual(values = c("ALT Significant" = "red", "Telomerase Significant" = "blue")) +
  theme_minimal() + geom_text_repel(data = bind_rows(ALT_significant, Telomerase_significant), 
                                    aes(label = hgnc_symbol), 
                                    vjust = 0.5, hjust = 0.5, size = 3,
                                    color = "black", , box.padding = 0.5,
                                    point.padding = 0.5, max.overlaps = Inf
  )

# Heatmap.