library(edgeR)
library(tidyverse)
library(biomaRt)
library(ggrepel)
library(ComplexHeatmap)
library(grid)
library(GSVA)
library(ggpubr)
library(survival)
library(survminer)

TARGET_NBL_htseq_counts <- read.delim("TARGET-NBL.htseq_counts.tsv", header=TRUE)

# removing Ensembl ID data after . in the name
TARGET_NBL_htseq_counts$Ensembl_ID <- gsub("\\..*", "", TARGET_NBL_htseq_counts$Ensembl_ID)

mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version=112)

# Retrieving the Ensembl gene IDs and gene biotype (protein coding)
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
                              filters="biotype",
                              values = "protein_coding",
                              mart = mart)

# Filtering TARGET_NBL_htseq_counts to only include protein-coding genes
TARGET_NBL_htseq_counts <- TARGET_NBL_htseq_counts %>%
  filter(Ensembl_ID %in% protein_coding_genes$ensembl_gene_id)


Neuroblastoma_Metadata <- read.table(file = 'Neuroblastoma_Metadata.txt', header = TRUE, sep = '\t')

# Data filtering: only including ALT and Telomerase phenotype. 
metadata_filtered_TMM <- Neuroblastoma_Metadata %>%
  filter(TMM %in% c("ALT", "Telomerase"))

# Data filtering: only include sample IDs present in the other dataset.
metadata_filtered_TMM <- metadata_filtered_TMM %>%
  filter(SampleID %in% colnames(TARGET_NBL_htseq_counts))

#setting Ensembl ID as rowname.
rownames(TARGET_NBL_htseq_counts) <- TARGET_NBL_htseq_counts$Ensembl_ID
TARGET_NBL_htseq_counts$Ensembl_ID <- NULL

# Arranging Sample_ID in the same order in both datasets:
# first arranging ALT then Telomerase in the metadata.
metadata_filtered_TMM <- metadata_filtered_TMM %>%
  arrange(factor(TMM, levels = c("ALT", "Telomerase")))

# counts_TMM is a subset of TARGET_NBL_htseq_counts with patients in order of metadata_filtered_TMM.
counts_TMM <- TARGET_NBL_htseq_counts[,metadata_filtered_TMM$SampleID]

#creating backup with Ensembl ID as column so it is easier to inspect.
counts_TMM_backup <- cbind(Ensembl_ID = rownames(counts_TMM), counts_TMM)



# building model matrix.

# First, determining the factors of TMM.
group1 <- as.factor(metadata_filtered_TMM$TMM)

# model matrix ~ without an intercept term.
design <- model.matrix(~group1+0)



# creating differential gene expression object.
dge_TMM <- DGEList(counts=counts_TMM,group=group1)

# TMM normalization.
dge_TMM <- calcNormFactors(dge_TMM)


# TMM normalization after removing lowly expressed genes with cpm < 1 in 5% of the samples.
keep <- rowSums(cpm(dge_TMM) > 1) >= ceiling(0.05*dim(dge_TMM)[2])
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
candidate_genes_ALT <- final_version[c(final_version$fdr <= 0.02 & final_version$logFC >= 0.9),]
candidate_genes_ALT <- candidate_genes_ALT %>%
  filter(hgnc_symbol != "" & !is.na(hgnc_symbol))

candidate_genes2_ALT <- final_version[c(final_version$fdr <= 0.01 & final_version$logFC >= 1.5),]



candidate_genes_Telomerase <- final_version[c(final_version$fdr <= 0.05 & final_version$logFC <= -0.5),]
candidate_genes_Telomerase <- candidate_genes_Telomerase %>%
  filter(hgnc_symbol != "" & !is.na(hgnc_symbol))

candidate_genes2_Telomerase <- final_version[c(final_version$fdr <= 0.01 & final_version$logFC <= -1.5),]

merged_candidates <- rbind(candidate_genes_ALT, candidate_genes_Telomerase)

# Top 10 Genes in terms of fold change (to label in the volcano plot.
ALT_significant <- merged_candidates %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 10)

Telomerase_significant <- merged_candidates %>%
  arrange(logFC) %>%
  slice_head(n=10)

# final version contains all genes while merged_candidates only contain candidate (ALT and Telomerase genes).
final_version1 <- final_version %>%
  mutate(gene_status = case_when(
    hgnc_symbol %in% candidate_genes_ALT$hgnc_symbol ~ "ALT",
    hgnc_symbol %in% candidate_genes_Telomerase$hgnc_symbol ~ "Telomerase"))


merged_candidates1 <- merged_candidates %>%
  mutate(gene_status = case_when(
    hgnc_symbol %in% candidate_genes_ALT$hgnc_symbol ~ "ALT",
    hgnc_symbol %in% candidate_genes_Telomerase$hgnc_symbol ~ "Telomerase"))

#Volcano plot ~ all samples.
ggplot(data = final_version1, aes(x = logFC, y = -log10(fdr), color = gene_status)) +
  geom_point() + 
  scale_color_manual(values = c("ALT" = "red", "Telomerase" = "blue")) +
  theme_minimal() + geom_text_repel(data = bind_rows(ALT_significant, Telomerase_significant), 
                  aes(label = hgnc_symbol), 
                  vjust = 0.5, hjust = 0.5, size = 3,
                  color = "black", box.padding = 0.5,
                  point.padding = 0.5, max.overlaps = Inf
                  )



# Heatmap ~ candidate genes.

# recalculating normalization factors using scale function.
dge_TMM3_scaled <- t(scale(t(dge_TMM3$counts)))

# only including genes in the merged_candidate dataset.
dge_TMM3_scaled <- dge_TMM3_scaled[rownames(dge_TMM3_scaled) %in% merged_candidates$ensembl_gene_id,]

phenotype_colors <- c("ALT" = "green", "Telomerase" = "black")
row_label_colors <- ifelse(merged_candidates$hgnc_symbol %in% candidate_genes_ALT$hgnc_symbol, 
                           "green", 
                           ifelse(merged_candidates$hgnc_symbol %in% candidate_genes_Telomerase$hgnc_symbol, 
                                  "black", 
                                  "gray"))


# matching the order.
row_order <- match(merged_candidates1$ensembl_gene_id, rownames(dge_TMM3_scaled))
dge_TMM3_scaled_ordered <- dge_TMM3_scaled[row_order, ]

col_order <- match(metadata_filtered_TMM$SampleID, colnames(dge_TMM3_scaled_ordered))
dge_TMM3_scaled_ordered <- dge_TMM3_scaled_ordered[, col_order]

Heatmap(dge_TMM3_scaled_ordered,
        row_labels = merged_candidates1$hgnc_symbol,
        row_names_gp = gpar(fontsize = 4, col = row_label_colors),
        column_names_gp = gpar(fontsize = 0.5),
        top_annotation = HeatmapAnnotation(Condition = anno_simple(metadata_filtered_TMM$TMM,
                                           col = phenotype_colors), which = "column"),
                                           show_column_dend = FALSE,
                                           show_row_dend = FALSE,
                                           column_split = metadata_filtered_TMM$TMM,
                                           row_split = merged_candidates1$gene_status,
        row_order = rownames(dge_TMM3_scaled_ordered),
        column_order = colnames(dge_TMM3_scaled_ordered)

)


# t-testing.

# Initializing an empty data frame for t-test results
t_test_results <- data.frame(
  hgnc_symbol = character(),
  p_value_t_test = numeric(),
  fdr_t_test = numeric(),
  stringsAsFactors = FALSE
)

# storage: all p-values for FDR adjustment later.
all_p_values <- numeric()

# Looping through each gene in merged_candidates
for (i in 1:nrow(merged_candidates)) {
  
  gene_id <- merged_candidates$ensembl_gene_id[i]
  gene_symbol <- merged_candidates$hgnc_symbol[i]
  
  # Extracting the expression values for this gene, keeping it as a matrix
  gene_expression <- dge_TMM3$counts[rownames(dge_TMM3) == gene_id, , drop = FALSE]
  
  # Creating ALT and Telomerase group
  alt_group <- gene_expression[, metadata_filtered_TMM$TMM == "ALT", drop = FALSE]
  telomerase_group <- gene_expression[, metadata_filtered_TMM$TMM == "Telomerase", drop = FALSE]
  
  
  t_test <- t.test(alt_group, telomerase_group)
  p_value_t_test <- t_test$p.value
  all_p_values <- c(all_p_values, p_value_t_test)
  
  
  # Storing the results
  t_test_results <- rbind(t_test_results, data.frame(
    hgnc_symbol = gene_symbol,
    p_value_t_test = p_value_t_test
  ))
}

# Calculating the FDR-adjusted p-values
t_test_results$fdr_t_test <- p.adjust(all_p_values, method = "BH")



# Joining the merged candidate and t-test result table.
merged_candidate_t <- full_join(merged_candidates, t_test_results, by = "hgnc_symbol")



# Running GSVA & ssGSEA.

gene_set_list <- list(candidate_genes_ALT$ensembl_gene_id)

gsvapar <- gsvaParam(dge_TMM3_scaled_ordered, gene_set_list, maxDiff= TRUE)
gsva_result <- gsva(gsvapar)

gsvapar_ssgsea <- ssgseaParam(dge_TMM3_scaled_ordered, 
                            gene_set_list)
ssGSEA_result <- gsva(gsvapar_ssgsea)

gsva_df <- as.data.frame(gsva_result)
ssGSEA_df <- as.data.frame(ssGSEA_result)

gsva_long <- pivot_longer(gsva_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
ssGSEA_long <- pivot_longer(ssGSEA_df, cols = everything(), names_to = "SampleID", values_to = "GSEA_Score")

# to match the order as per the data, need to convert SampleID into factor.
gsva_long$SampleID <- factor(gsva_long$SampleID, levels = colnames(dge_TMM3_scaled_ordered))
ssGSEA_long$SampleID <- factor(ssGSEA_long$SampleID, levels = colnames(dge_TMM3_scaled_ordered))

gsva_long$Color <- ifelse(gsva_long$GSVA_Score > 0, "Positive", "Negative")
ssGSEA_long$Color <- ifelse(ssGSEA_long$GSEA_Score > 0, "Positive", "Negative")

# Creating bar plot for ssGSEA.
ssgsea_annotation <- anno_barplot(
  ssGSEA_long$GSEA_Score, 
  gp = gpar(fill = ifelse(ssGSEA_long$Color == "Positive", "blue", "red")),
  border = FALSE, 
  axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
)

# Creating bar plot for GSVA.
gsva_annotation <- anno_barplot(
  gsva_long$GSVA_Score, 
  gp = gpar(fill = ifelse(gsva_long$Color == "Positive", "blue", "red")),
  border = FALSE,
  axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
)

# Combining the annotations.
complex_annotation <- HeatmapAnnotation(
  ssGSEA = ssgsea_annotation,
  GSVA = gsva_annotation,
  Condition = anno_simple(metadata_filtered_TMM$TMM, col = phenotype_colors),
  which = "column"
)

# Final heatmap + barPlots.
h1 <- Heatmap(
  dge_TMM3_scaled_ordered,
  row_labels = merged_candidates1$hgnc_symbol,
  row_names_gp = gpar(fontsize = 4, col = row_label_colors),
  column_names_gp = gpar(fontsize = 0.5),
  top_annotation = complex_annotation,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_split = metadata_filtered_TMM$TMM,
  row_split = merged_candidates1$gene_status,
  row_order = rownames(dge_TMM3_scaled_ordered),
  column_order = colnames(dge_TMM3_scaled_ordered)
)

draw(h1)




# boxplot - GSVA.
ALT_gsva_df <- gsva_df[, colnames(gsva_df) %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$TMM == "ALT"]]
ALT_gsva_df <- pivot_longer(ALT_gsva_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
ALT_gsva_df <- ALT_gsva_df %>%
  mutate(Phenotype = "ALT")

Tel_gsva_df <- gsva_df[, colnames(gsva_df) %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$TMM == "Telomerase"]]
Tel_gsva_df <- pivot_longer(Tel_gsva_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
Tel_gsva_df <- Tel_gsva_df %>%
  mutate(Phenotype = "Telomerase")

combined_df <- rbind(ALT_gsva_df, Tel_gsva_df)

ggboxplot(combined_df, x = "Phenotype", y = "GSVA_Score", 
          fill = "Phenotype") +
  stat_compare_means(method = "t.test") +
  labs(title = "GSVA Score Comparison")

# boxplot - ssGSEA.
ALT_gsea_df <- ssGSEA_df[, colnames(ssGSEA_df) %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$TMM == "ALT"]]
ALT_gsea_df <- pivot_longer(ALT_gsea_df, cols = everything(), names_to = "SampleID", values_to = "GSEA_Score")
ALT_gsea_df <- ALT_gsea_df %>%
  mutate(Phenotype = "ALT")

Tel_gsea_df <- ssGSEA_df[, colnames(ssGSEA_df) %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$TMM == "Telomerase"]]
Tel_gsea_df <- pivot_longer(Tel_gsea_df, cols = everything(), names_to = "SampleID", values_to = "GSEA_Score")
Tel_gsea_df <- Tel_gsea_df %>%
  mutate(Phenotype = "Telomerase")

combined_df2 <- rbind(ALT_gsea_df, Tel_gsea_df)

ggboxplot(combined_df2, x = "Phenotype", y = "GSEA_Score", 
          fill = "Phenotype") +
  stat_compare_means(method = "t.test") +
  labs(title = "ssGSEA Score Comparison")


###########################

# survival plot.

# First, deleting the Telomerase samples with GSVA score > 0.
combined_df <- combined_df[!(combined_df$GSVA_Score > 0 & combined_df$Phenotype == "Telomerase") & 
           !(combined_df$GSVA_Score < 0 & combined_df$Phenotype == "ALT"), ]

survival_metadata <- metadata_filtered_TMM[metadata_filtered_TMM$SampleID %in% combined_df$SampleID,]
survival_metadata$Vital.Status <- ifelse(survival_metadata$Vital.Status == "Dead", 1, 0) # recoding event to 1 or 0.


# computing survival curve -with overall survival time.
fit <- survfit(Surv(Overall.Survival.Time.in.Days, Vital.Status) ~ TMM, data = survival_metadata)

# plotting the graph.
ggsurvplot(fit,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ncensor.plot = TRUE,
           ggtheme = theme_bw(),
           palette = c("#E7B800", "#2E9FDF"))


## No_TMM is supposed to have high survival rate with significant p value. Checking...
no_tmm_rows <- Neuroblastoma_Metadata[Neuroblastoma_Metadata$TMM_Case == "NO_TMM", ]
survival_metadata2 <- rbind(no_tmm_rows, survival_metadata)
survival_metadata2$Vital.Status <- ifelse(survival_metadata2$Vital.Status %in% c(1, 0), 
                                         survival_metadata2$Vital.Status, 
                                         ifelse(survival_metadata2$Vital.Status == "Dead", 1, 0))

# computing survival curve -with overall survival time. (No_TMM vs TMM and ALT)
survival_metadata2$Vital.Status <- as.numeric(survival_metadata2$Vital.Status)
fit2 <- survfit(Surv(Overall.Survival.Time.in.Days, Vital.Status) ~ TMM_Case, data = survival_metadata2)

# plotting the graph.
ggsurvplot(fit2,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ncensor.plot = TRUE,
           ggtheme = theme_bw(),
           palette = c("#E7B800", "#2E9FDF"))

# comparison of survival probability between different risk group.
fit3 <- survfit(Surv(Overall.Survival.Time.in.Days, Vital.Status) ~ COG.Risk.Group, data = survival_metadata2)      

# plotting the graph.
ggsurvplot(fit3,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ncensor.plot = TRUE,
           ggtheme = theme_bw())

# comparing ALT signatures for different gender.
gsea_male <- ALT_gsea_df[ALT_gsea_df$SampleID %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$Gender == "Male"], ]
gsea_male <- gsea_male %>%
  mutate(Gender = "Male")


gsea_female <- ALT_gsea_df[ALT_gsea_df$SampleID %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$Gender == "Female"], ]
gsea_female <- gsea_female %>%
  mutate(Gender = "Female")

gsea_gender <- rbind(gsea_male, gsea_female)

ggboxplot(gsea_gender, x = "Gender", y = "GSEA_Score", 
          fill = "Gender") +
  stat_compare_means(method = "t.test") +
  labs(title = "GSEA Score Comparison for different Genders in ALT samples")


# comparing ALT signatures for different ploidy.
gsea_hyperploid <- ALT_gsea_df[ALT_gsea_df$SampleID %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$Ploidy == "Hyperdiploid (DI>1)"], ]
gsea_hyperploid <- gsea_hyperploid %>%
  mutate(Ploidy = "Hyperploid")

gsea_diploid <- ALT_gsea_df[ALT_gsea_df$SampleID %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$Ploidy == "Diploid (DI=1)"], ]
gsea_diploid <- gsea_diploid %>%
  mutate(Ploidy = "Diploid")

gsea_ploidy <- rbind(gsea_hyperploid, gsea_diploid)

ggboxplot(gsea_ploidy, x = "Ploidy", y = "GSEA_Score", 
          fill = "Ploidy") +
  stat_compare_means(method = "t.test") +
  labs(title = "GSEA Score Comparison for different Ploidy States in ALT samples")


# comparing ALT signatures for racial characteristics.
gsea_white <- ALT_gsea_df[ALT_gsea_df$SampleID %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$Race == "White"], ]
gsea_white <- gsea_white %>%
  mutate(Race = "White")

gsea_black <- ALT_gsea_df[ALT_gsea_df$SampleID %in% metadata_filtered_TMM$SampleID[metadata_filtered_TMM$Race == "Black or African American"], ]
gsea_black <- gsea_black %>%
  mutate(Race = "Black")

gsea_race <- rbind(gsea_white, gsea_black)

ggboxplot(gsea_race, x = "Race", y = "GSEA_Score", 
          fill = "Race") +
  stat_compare_means(method = "t.test") +
  labs(title = "GSEA Score Comparison for different Races in ALT samples")




