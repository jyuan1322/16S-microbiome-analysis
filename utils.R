# utils.R

library(ShortRead)
library(scales)
library(dplyr)
library(ComplexHeatmap)
library(circlize)  # for colorRamp2

summarize_quality_by_position <- function(fastq_files, output_csv = "combined_quality_summary.csv") {
  quality_list <- list()
  
  for (fq in fastq_files) {
    cat("Processing:", fq, "\n")
    
    fq_reads <- readFastq(fq)

    qs_matrix <- as(quality(fq_reads), "matrix")
    mean_q <- colMeans(qs_matrix)

    # Clean base name
    fq_base <- tools::file_path_sans_ext(basename(fq))
    fq_base <- sub("\\.fastq$", "", fq_base)
    fq_base <- sub("\\.fq$", "", fq_base)
    fq_base <- sub("\\.gz$", "", fq_base)

    quality_list[[fq_base]] <- round(mean_q, 2)
  }

  # Combine into data frame
  quality_df <- as.data.frame(quality_list)
  quality_df$Cycle <- seq_len(nrow(quality_df))
  quality_df <- quality_df[, c("Cycle", setdiff(names(quality_df), "Cycle"))]

  # Mean and standard deviation across samples
  quality_matrix <- as.matrix(quality_df[, -1])
  quality_df$Mean_All <- round(rowMeans(quality_matrix), 2)
  quality_df$SD_All <- round(apply(quality_matrix, 1, sd), 2)

  # Write to CSV
  write.csv(quality_df, output_csv, row.names = FALSE)
  cat("Saved combined quality summary to:", output_csv, "\n")
  
  return(quality_df)
}



# print percent of reads assigned to each taxa level in each file
classification_rate_table <- function(seqtab, taxa) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # Ensure taxa and seqtab match by ASV (column names)
  common_seqs <- intersect(colnames(seqtab), rownames(taxa))
  seqtab <- seqtab[, common_seqs]
  taxa <- taxa[common_seqs, ranks, drop = FALSE]
  
  n_samples <- nrow(seqtab)
  result <- matrix(NA, nrow = n_samples, ncol = length(ranks) + 1)
  rownames(result) <- rownames(seqtab)
  colnames(result) <- c("Total_Reads", ranks)

  for (i in seq_len(n_samples)) {
    sample_counts <- seqtab[i, ]
    total_reads <- sum(sample_counts)
    result[i, "Total_Reads"] <- total_reads
    
    for (rank in ranks) {
      assigned <- !is.na(taxa[, rank])
      assigned_reads <- sum(sample_counts[assigned])
      result[i, rank] <- round(100 * assigned_reads / total_reads, 2)
    }
  }

  return(as.data.frame(result))
}

# Function to compute percent composition at a taxa level
taxa_percent_table <- function(seqtab, taxa, rank = "Genus") {
  # Check the rank exists
  if (!rank %in% colnames(taxa)) {
    stop(paste("Rank", rank, "not found in taxonomy table."))
  }
  
  # Align ASVs between seqtab and taxa
  common_seqs <- intersect(colnames(seqtab), rownames(taxa))
  seqtab <- seqtab[, common_seqs]
  taxa <- taxa[common_seqs, , drop = FALSE]
  
  # Replace NA with "Unassigned" in the selected rank
  tax_level <- taxa[, rank]
  tax_level[is.na(tax_level)] <- "Unassigned"
  
  # Sum reads by taxon at the chosen rank
  tax_counts <- t(rowsum(t(seqtab), group = tax_level))
  
  # Convert to percent per sample
  tax_pct <- tax_counts / rowSums(tax_counts) * 100
  tax_pct <- round(tax_pct, 2)
  
  return(as.data.frame(tax_pct))
}

filter_taxa_by_prevalence <- function(pct_table, min_abundance = 1, min_prevalence = 0.2) {
  # Remove Sample column if present
  sample_col <- "Sample" %in% colnames(pct_table)
  if (sample_col) {
    samples <- pct_table$Sample
    pct_table <- pct_table[, !colnames(pct_table) %in% "Sample"]
  }
  
  # Compute in how many samples each taxon exceeds min_abundance
  prevalence <- colMeans(pct_table > min_abundance)

  # Keep taxa that are above threshold in at least min_prevalence fraction of samples
  keep_taxa <- names(prevalence[prevalence >= min_prevalence])

  # Reconstruct filtered table
  filtered <- pct_table[, keep_taxa, drop = FALSE]

  # Restore Sample column if it was present
  if (sample_col) {
    filtered$Sample <- samples
  }

  return(filtered)
}

plot_taxa_heatmap <- function(pct_table,
                              min_abundance = 1,
                              min_prevalence = 0.2,
                              scale_rows = FALSE,
                              tax_rank = "Genus",
                              cluster_rows = FALSE,
                              cluster_cols = TRUE) {

  # Extract and drop sample column if present
  has_sample <- "Sample" %in% colnames(pct_table)
  if (has_sample) {
    samples <- pct_table$Sample
    pct_table <- pct_table[, !colnames(pct_table) %in% "Sample"]
    rownames(pct_table) <- samples
  }

  # Filter taxa based on abundance and prevalence
  prevalence <- colMeans(pct_table > min_abundance)
  keep_taxa <- names(prevalence[prevalence >= min_prevalence])
  # Always include "Unassigned" if present
  if ("Unassigned" %in% colnames(pct_table)) {
    keep_taxa <- union(keep_taxa, "Unassigned")
  }
  filtered <- pct_table[, keep_taxa, drop = FALSE]

  # Optionally scale rows (z-score normalization)
  if (scale_rows) {
    filtered <- t(scale(t(filtered)))
  }

  # Melt for ggplot
  # df_long <- melt(filtered, varnames = c("Sample", "Taxon"), value.name = "Abundance")
  filtered$Sample <- rownames(filtered)
  filtered$Sample <- sub("^CP101", "", filtered$Sample)
  filtered$Sample <- sub("^(.{2})(.)", "\\1 \\2", filtered$Sample)

  df_long <- melt(filtered, id.vars = "Sample", variable.name = "Taxon", value.name = "Abundance")
  colnames(df_long)[2] <- tax_rank  # Rename for prettier axis label

  # Plot heatmap
  gg <- ggplot(df_long, aes(x = Sample, y = .data[[tax_rank]], fill = Abundance)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = if (scale_rows) 0 else median(df_long$Abundance, na.rm = TRUE),
                         name = if (scale_rows) "Z-score" else "Percent") +
    scale_x_discrete(position = "top") +
    labs(x = "Sample", y = tax_rank,
         title = paste("Heatmap of", tax_rank, "Abundance")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0),
          axis.text.y = element_text(size = 8))

  # Optional clustering (sort rows/columns by hierarchical clustering)
  if (cluster_rows) {
    filtered_clustering <- filtered %>% select(-Sample, Unassigned)
    taxa_order <- hclust(dist(t(filtered_clustering)))$order
    taxa_names_ordered <- colnames(filtered_clustering)[taxa_order]
    taxa_names_ordered <- c(taxa_names_ordered, "Unassigned")
    gg <- gg + scale_y_discrete(limits = taxa_names_ordered)
  }
  if (cluster_cols) {
    sample_order <- hclust(dist(filtered))$order
    gg <- gg + scale_x_discrete(limits = rownames(filtered)[sample_order])
  }

  return(gg)
}


# Function to plot ASV heatmap with optional filtering
# seqtab_filtered: ASV table with samples as rows and ASVs as columns
# meta_data: metadata for samples (partial mayo, treatment)
# min_abundance: minimum read count for an ASV to be included
# min_prevalence: minimum fraction of samples an ASV must be present in
# scale_rows: whether to z-score normalize rows (ASVs)
# cluster_rows: whether to cluster ASVs by hierarchical clustering
# cluster_cols: whether to cluster samples by hierarchical clustering
# Returns a ggplot object
# Note: Requires ggplot2, reshape2, ComplexHeatmap
# If ComplexHeatmap is not available, it falls back to ggplot2 heatmap
plot_asv_heatmap <- function(seqtab_filtered,
                             meta_data,
                             min_abundance = 1,
                             min_prevalence = 0.2,
                             rel_abund = TRUE,
                             cluster_rows = TRUE,
                             cluster_cols = FALSE) {
    library(reshape2)
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)

    # Ensure ASV IDs are columns and sample names are rownames
    asv_table <- as.data.frame(seqtab_filtered)
    asv_table$Sample <- rownames(asv_table)

    # Filter ASVs based on prevalence
    prevalence <- colMeans(asv_table[, !colnames(asv_table) %in% "Sample"] > min_abundance)
    keep_asvs <- names(prevalence[prevalence >= min_prevalence])
    filtered <- asv_table[, c(keep_asvs, "Sample"), drop = FALSE]

    # Join with metadata
    meta <- meta_data[match(filtered$Sample, meta_data$sample_id), ]
    # Reorder by treatment group if clustering is off
    if (!cluster_cols && "treatment" %in% colnames(meta)) {
        treatment_order <- order(meta$treatment)
        filtered <- filtered[treatment_order, ]
        meta <- meta[treatment_order, ]
    }

    # Optionally scale rows (Z-score by ASV)
    mat <- as.matrix(filtered[, keep_asvs])
    # if (scale_rows) {
        # mat <- t(scale(t(mat)))
    if (rel_abund) {
        # mat <- sweep(mat, 1, rowSums(mat), FUN = "/")  # convert to relative abundance
        # Avoid division by zero
        row_totals <- rowSums(mat)
        safe_row_totals <- ifelse(row_totals == 0, 1, row_totals)
        mat <- sweep(mat, 1, safe_row_totals, FUN = "/")
        mat[row_totals == 0, ] <- 0
    }
    rownames(mat) <- filtered$Sample

    # If clustering is enabled, use ComplexHeatmap
    if (cluster_rows || cluster_cols) {
        col_fun <- circlize::colorRamp2(
            breaks = c(min(mat), max(mat)),
            colors = c("white", "red")
        )

        # Define treatment color
        treatment_colors <- c("placebo" = "#1f77b4", "active" = "#ff7f0e")

        # Define annotation
        bottom_annot <- HeatmapAnnotation(
            Treatment = meta$treatment,
            Mayo = meta$total_partial_mayo,
            col = list(
                Treatment = treatment_colors,
                Mayo = circlize::colorRamp2(
                    range(meta$total_partial_mayo, na.rm = TRUE),
                    c("white", "darkgreen")
                )
            ),
            which = "column",
            annotation_name_side = "left",
            annotation_legend_param = list(
                Treatment = list(title = "Treatment"),
                Mayo = list(title = "Mayo Score")
            )
        )

        ht <- Heatmap(
            t(mat),
            name = if (rel_abund) "Rel Abundance" else "Abundance",
            col = col_fun,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_cols,
            show_row_dend = cluster_rows,
            show_column_dend = cluster_cols,
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 8),
            column_names_rot = 90,
            heatmap_legend_param = list(direction = "horizontal"),
            top_annotation = NULL,
            bottom_annotation = bottom_annot
        )
        return(ht)
    }

    # Otherwise, fall back to ggplot2 heatmap
    filtered$Sample <- sub("^CP101", "", filtered$Sample)
    filtered$Sample <- sub("^(.{2})(.)", "\\1 \\2", filtered$Sample)

    df_long <- melt(filtered, id.vars = "Sample",
                    variable.name = "ASV", value.name = "Abundance")

    gg <- ggplot(df_long, aes(x = Sample, y = ASV, fill = Abundance)) +
        geom_tile() +
        scale_fill_gradient(
            low = "white", high = "red",
            name = if (rel_abund) "Rel Abund" else "Abundance"
        ) +
        scale_x_discrete(position = "top") +
        labs(x = "Sample", y = "ASV", title = "Heatmap of ASV Abundance") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0),
            axis.text.y = element_text(size = 6)
        )

    return(gg)
}



plot_asv_heatmap2 <- function(seqtab_filtered,
                             meta_data,
                             min_abundance = 1,
                             min_prevalence = 0.2,
                             rel_abund = TRUE,
                             cluster_rows = TRUE,
                             cluster_cols = FALSE) {
    library(reshape2)
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)
    library(stringr)

    # Ensure ASV IDs are columns and sample names are rownames
    asv_table <- as.data.frame(seqtab_filtered)
    asv_table$Sample <- rownames(asv_table)

    # Filter ASVs based on prevalence
    prevalence <- colMeans(asv_table[, !colnames(asv_table) %in% "Sample"] > min_abundance)
    keep_asvs <- names(prevalence[prevalence >= min_prevalence])
    filtered <- asv_table[, c(keep_asvs, "Sample"), drop = FALSE]
    filtered$SampleSort <- gsub("Day1", "Wk0", gsub("CP101", "", filtered$Sample))

    # Join with metadata
    meta <- meta_data[match(filtered$Sample, meta_data$sample_id), ]

    # Reorder samples: treatment first, then numeric order from sample name
    if (!cluster_cols && "treatment" %in% colnames(meta)) {
        # Extract numeric part from sample name and time point
        all_nums <- str_extract_all(filtered$SampleSort, "\\d+")
        sample_num <- sapply(all_nums, function(v) if(length(v) >= 1) as.numeric(v[1]) else NA)
        time_num <- sapply(all_nums, function(v) if(length(v) >= 2) as.numeric(v[2]) else NA)

        treatment_order <- order(meta$treatment, sample_num, time_num)
        filtered <- filtered[treatment_order, ]
        meta <- meta[treatment_order, ]
    }

    # Optionally scale to relative abundance
    mat <- as.matrix(filtered[, keep_asvs])
    if (rel_abund) {
        row_totals <- rowSums(mat)
        safe_row_totals <- ifelse(row_totals == 0, 1, row_totals)
        mat <- sweep(mat, 1, safe_row_totals, FUN = "/")
        mat[row_totals == 0, ] <- 0
    }
    rownames(mat) <- filtered$Sample

    # If clustering is enabled, use ComplexHeatmap
    if (cluster_rows || cluster_cols) {
        col_fun <- circlize::colorRamp2(
            breaks = c(min(mat), max(mat)),
            colors = c("white", "red")
        )

        treatment_colors <- c("placebo" = "#1f77b4", "active" = "#ff7f0e")

        bottom_annot <- HeatmapAnnotation(
            Treatment = meta$treatment,
            Mayo = meta$total_partial_mayo,
            col = list(
                Treatment = treatment_colors,
                Mayo = circlize::colorRamp2(
                    range(meta$total_partial_mayo, na.rm = TRUE),
                    c("white", "darkgreen")
                )
            ),
            which = "column",
            annotation_name_side = "left",
            annotation_legend_param = list(
                Treatment = list(title = "Treatment"),
                Mayo = list(title = "Mayo Score")
            )
        )

        ht <- Heatmap(
            t(mat),
            name = if (rel_abund) "Rel Abundance" else "Abundance",
            col = col_fun,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_cols,
            show_row_dend = cluster_rows,
            show_column_dend = cluster_cols,
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 8),
            column_names_rot = 90,
            heatmap_legend_param = list(direction = "horizontal"),
            top_annotation = NULL,
            bottom_annotation = bottom_annot
        )
        return(ht)
    }

    # Fallback: ggplot2 heatmap
    df_long <- melt(filtered, id.vars = "Sample",
                    variable.name = "ASV", value.name = "Abundance")

    gg <- ggplot(df_long, aes(x = factor(Sample, levels = unique(filtered$Sample)),
                              y = ASV, fill = Abundance)) +
        geom_tile() +
        scale_fill_gradient(
            low = "white", high = "red",
            name = if (rel_abund) "Rel Abund" else "Abundance"
        ) +
        scale_x_discrete(position = "top") +
        labs(x = "Sample", y = "ASV", title = "Heatmap of ASV Abundance") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0),
            axis.text.y = element_text(size = 6)
        )

    return(gg)
}