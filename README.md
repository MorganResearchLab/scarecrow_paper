# scarecrow

This repo contains R code used in the analysis presented in the [scarecrow](https://github.com/MorganResearchLab/scarecrow) paper. Several of the QC steps are borrowed from the Pachter lab kallisto [tutorial](https://pachterlab.github.io/kallistobustools/tutorials/kb_building_atlas/R/kb_analysis_0_R/).

## Libraries and functions

```R
library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggdist)
library(ggpointdensity)
library(viridis)
library(scales)
library(ggpmisc)

# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "CsparseMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

# Returns a tibble with total UMI counts for each barcode, and
# rank of the total counts, with number 1 for the barcode with the most counts.
get_knee_df <- function(mat) {
  total <- rank <- NULL
  tibble(total = Matrix::colSums(mat),
         rank = row_number(desc(total))) %>%
    distinct() %>%
    dplyr::filter(total > 0) %>%
    arrange(rank)
}

# Minimum total UMI counts for barcode for it to be considered when calculating the inflection point;
# this helps to avoid the noisy part of the curve for barcodes with very few counts.
get_inflection <- function(df, lower = 100) {
  log_total <- log_rank <- total <-  NULL
  df_fit <- df %>%
    dplyr::filter(total >= lower) %>%
    transmute(log_total = log10(total),
              log_rank = log10(rank))
  d1n <- diff(df_fit$log_total)/diff(df_fit$log_rank)
  right.edge <- which.min(d1n)
  10^(df_fit$log_total[right.edge])
}

# Plot a transposed knee plot, showing the inflection point and
# the number of remaining cells after inflection point filtering. It's
# transposed since it's more generalizable to multi-modal data. Taken from the
# BUSpaRse package.
knee_plot <- function(df, inflection, fontsize=10) {
  total <- rank_cutoff <- NULL
  annot <- tibble(inflection = inflection, rank_cutoff = max(df$rank[df$total > inflection]))
  ggplot(df, aes(total, rank)) + geom_path() +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2, color = "gray40") +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2, color = "gray40") +
    geom_text(aes(inflection, rank_cutoff, label = paste(rank_cutoff, "'cells'")), data = annot, vjust = 1, size=3) +
    scale_x_log10(labels = label_comma()) + scale_y_log10(labels = label_comma()) +
    labs(y = "Rank", x = "Total UMIs") +
    annotation_logticks() & theme_bw() &
    theme(text = element_text(size=fontsize))
}

# Plot PCT genes
plot_pct_genes <- function(mat, tr2g, top_n = 20, symbol = "ensembl", fontsize=10) {
  pct_tx <- rowSums(mat)
  gs <- rownames(mat)[order(-pct_tx)]
  df <- as.data.frame(t(mat[gs[1:top_n],]))
  df <- df %>%
    mutate_all(function(x) x/colSums(mat)) %>%
    pivot_longer(everything(), names_to = "gene")
  if (symbol == "ensembl") {
    df <- left_join(df, tr2g, by = "gene")
  } else {
    df <- rename(df, gene_name = gene)
  }
  df %>%
    mutate(gene = fct_reorder(gene_name, value, .fun = median)) %>%
    ggplot(aes(gene, value)) +
    geom_boxplot() +
    labs(x = "", y = "Proportion of total counts") +
    coord_flip() +
    theme_bw() + theme(text = element_text(size=fontsize))
}
```

## Import transcript to gene table

This is taken from the `kallisto` genome indices.

```R
# Load gene to transcript data
tr2g <- read_tsv("~/Documents/scarecrow_test/results/transcripts_to_genes.txt", col_names = c("transcript", "gene", "gene_name"))
tr2g <- distinct(tr2g[, c("gene", "gene_name")])
tr2g[is.na(tr2g$gene_name),]$gene_name <- tr2g[is.na(tr2g$gene_name),]$gene
# uniquify the duplicates (some ensembl gene IDs map to the same gene symbol)
gene_name_dups <- tr2g$gene_name[which(duplicated(tr2g$gene_name))]
for(n in gene_name_dups) {
  tr2g[which(tr2g$gene_name == n),]$gene_name <- paste(tr2g[which(tr2g$gene_name == n),]$gene_name, tr2g[which(tr2g$gene_name == n),]$gene, sep="-")
}
```

## Import kallisto counts matrices

After generating count matrices for each dataset with `kallisto` under different `scarecrow` jitter settings, the `counts_unfiltered/cells*` files from each run are read into a list and named. In the below example the Parse data is imported for count matrices generated at jitter 0, 1, and 2.

```R
setwd("~/Documents/scarecrow/results")
data <- list.files(path = "~/Documents/scarecrow/results", recursive = T, pattern = "cells_x_genes.mtx")

# Parse WTv2 data
raw <- list(read_count_output(dirname(data[grep("Parse-WTv2/J0", data)]), name = "cells_x_genes"),
            read_count_output(dirname(data[grep("Parse-WTv2/J1", data)]), name = "cells_x_genes"),
            read_count_output(dirname(data[grep("Parse-WTv2/J2", data)]), name = "cells_x_genes"))
names(raw) <- c("Parse-WTv2-J0", "Parse-WTv2-J1", "Parse-WTv2-J2")
```

## Processing count matrices

After importing the data, we apply a series of steps to the list of count matrices. These steps are based on the `kallisto` tutorial linked at the start of the page. Briefly, we (1) test of library saturation; (2) generate a knee plot to identify the inflection point in barcode count; (3) filter each count matrix based on its inflection point and return a `Seurat` object; (4) plot mtDNA fraction and RNA counts;  (5) filter data to retain cells whose mtDNA content is below the 99th percentile. The filtering applied here is very rudimentary.

```R
# Test for library saturation
libsat <- lapply(raw, function(res_mat) {
  tot_counts <- colSums(res_mat)
  lib_sat <- tibble(nCount = tot_counts, nGene = colSums(res_mat > 0))

  p1 <- ggplot(lib_sat, aes(nCount, nGene)) +
    geom_pointdensity(size=1, adjust=2) + scale_color_viridis(name="density\n(neighbours)", labels = label_comma()) +
    scale_x_log10(labels = label_comma()) + scale_y_log10(labels = label_comma()) + annotation_logticks() +
    xlab("Number of molecules") + ylab("Number of features") +
    theme_bw() + theme(text = element_text(size=fontsize), legend.position = "right")
  return(list(p1 = p1))
})

# Knee plot for inflection point in barcode count
knees <- lapply(raw, function(res_mat) {
  knee_df <- get_knee_df(res_mat)
  inflection <- get_inflection(knee_df, lower=100)
  return(list(inflection = inflection, knee_plot = knee_plot(knee_df, inflection, fontsize=fontsize)))
})

# Filter and return Seurat Object
mats <- lapply(1:length(raw), function(i) {
  res_mat <- raw[[i]][, colSums(raw[[i]]) > knees[[i]]$inflection]
  res_mat <- res_mat[Matrix::rowSums(res_mat) > 0,]
  tmp <- rownames(res_mat)
  rownames(res_mat) <- tr2g$gene_name[match(rownames(res_mat), tr2g$gene)]
  rownames(res_mat)[which(is.na(rownames(res_mat)))] <- tmp[which(is.na(rownames(res_mat)))]
  CreateSeuratObject(counts = res_mat, project = names(raw)[i], min.cells = 1, min.features = 1)
})

# Plot mt percent and nCount_RNAs
mt_pct <- lapply(mats, function(dat) {
  dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
  p1a <- ggplot(data.frame(X=dat$nFeature_RNA), aes(y=X, slab_color = after_stat(y))) +
    stat_slabinterval(side = "right", lty=0) +
    stat_dotsinterval(slab_shape = 19, quantiles =100, side="left") +
    scale_color_distiller(aesthetics = "slab_color", guide = "colorbar2") +
    ylab("Features per cell") +
    scale_y_continuous(labels = label_comma()) +
    theme_bw() + theme(text = element_text(size=fontsize),
                       axis.text.y = element_blank(), axis.title.y = element_blank(),
                       legend.position = "None") + coord_flip() +
    theme(plot.margin = margin(l = 0.75, r = 0.5, unit = "cm"))

  p1b <- ggplot(data.frame(X=dat$nCount_RNA), aes(y=X, slab_color = after_stat(y))) +
    stat_slabinterval(side = "right", lty=0) +
    stat_dotsinterval(slab_shape = 19, quantiles =100, side="left") +
    scale_color_distiller(aesthetics = "slab_color", guide = "colorbar2") +
    ylab("Molecules per cell") +
    scale_y_continuous(labels = label_comma()) +
    theme_bw() + theme(text = element_text(size=fontsize),
                       axis.text.y = element_blank(), axis.title.y = element_blank(),
                       legend.position = "None") + coord_flip() +
    theme(plot.margin = margin(l = 0.75, r = 0.5, unit = "cm"))
  p1c <- ggplot(data.frame(X=dat$percent.mt), aes(y=X, slab_color = after_stat(y))) +
    stat_slabinterval(side = "right", lty=0) +
    stat_dotsinterval(slab_shape = 19, quantiles =100, side="left") +
    scale_color_distiller(aesthetics = "slab_color", guide = "colorbar2") +
    ylab("% mitochondrial genes") +
    theme_bw() + theme(text = element_text(size=fontsize),
                       axis.text.y = element_blank(), axis.title.y = element_blank(),
                       legend.position = "None") + coord_flip() +
    theme(plot.margin = margin(l = 0.75, r = 0.5, unit = "cm"))
  p1 <- ggarrange(p1a,p1b,p1c,nrow=3)

  p2 <- ggplot(data.frame(X=dat$nCount_RNA, Y=dat$percent.mt), aes(x=X, y=Y)) +
    geom_pointdensity(size=1, adjust=5) + scale_color_viridis(name="density\n(neighbours)", labels = label_comma()) +
    xlab("Molecules per cell") + ylab("% mitochondrial genes") +
    scale_x_continuous(labels = label_comma()) +
    geom_hline(yintercept=quantile(dat$percent.mt,.99), lty=2) +
    annotate("text", x=max(dat$nCount_RNA)/2, y=quantile(dat$percent.mt,.99), hjust=0, vjust=-1, size=3, label = "99th percentile") +
    theme_bw() + theme(text = element_text(size=fontsize), legend.position = "right")

  p3 <- ggplot(data.frame(X=dat$nCount_RNA, Y=dat$nFeature_RNA), aes(x=X, y=Y)) +
    geom_pointdensity(size=1, adjust=5) + scale_color_viridis(name="density\n(neighbours)", labels = label_comma()) +
    xlab("Molecules per cell") + ylab("Features per cell") +
    scale_x_continuous(labels = label_comma()) + scale_y_continuous(labels = label_comma()) +
    geom_hline(yintercept=quantile(dat$nFeature_RNA,.99), lty=2) +
    geom_vline(xintercept=quantile(dat$nCount_RNA,.99), lty=2) +
    annotate("text", x=max(dat$nCount_RNA)/2, y=quantile(dat$nFeature_RNA,.99), hjust=0, vjust=-1, size=3, label = "99th percentile") +
    theme_bw() + theme(text = element_text(size=fontsize), legend.position = "None")

  return(list(p1 = p1, p2 = p2, p3 = p3))
})

# Apply filter and normalize
mats_filtered <- lapply(mats, function(dat) {
  dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
  dat <- subset(dat, subset = percent.mt < quantile(dat$percent.mt,.99))
  NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)
})
```

## Plot

We generate multi-panel plots to illustrate some of the results from the above steps as follows.

```R
for(i in 1:length(raw)){
  ggarrange(ggarrange(libsat[[i]]$p1, knees[[i]]$knee_plot, widths=c(4,4), labels=c("A", "B"), nrow=1),
            ggarrange(mt_pct[[i]]$p1, mt_pct[[i]]$p2, widths=c(3,5), labels=c("C", "D"), nrow=1),
            nrow=2) %>%
    ggsave(width = 8, height = 6, units="in", dpi=300, filename = paste0("plots/", names(raw)[i], "_kallisto.png"))
}
```

Below is an example, showing the plot generated from the Parse data at jitter 2.

<img src="./img/Parse-WTv2-J2_kallisto.png" alt="Parse jitter 2 processing results"/>

Next we compare the mean and median counts per gene between two sets of results. In the below code we compare the first list element to the last, which for the Parse data results is comparing jitter 0 to jitter 2.

```R
# Function to calculate row-wise stats of dgCMatrix
rowstats_dgCMatrix <- function(dgCMatrix) {
  # Convert to triplet form (i,j,x)
  mat_triplet <- summary(dgCMatrix)

  # Convert character '.' to NA if they exist
  if (is.character(mat_triplet$x)) {
    mat_triplet$x[mat_triplet$x == "."] <- NA
    mat_triplet$x <- as.numeric(mat_triplet$x)
  }

  # Filter for non-zero, non-NA values
  valid_values <- mat_triplet[mat_triplet$x != 0 & !is.na(mat_triplet$x), ]

  # Initialize result vector with NAs
  medians <- rep(NA_real_, nrow(dgCMatrix))
  means <- rep(NA_real_, nrow(dgCMatrix))
  min_count <- rep(NA_real_, nrow(dgCMatrix))
  max_count <- rep(NA_real_, nrow(dgCMatrix))

  # If we have any valid values...
  if (nrow(valid_values) > 0) {
    # Split values by row index
    values_by_row <- split(valid_values$x, valid_values$i)
    # Calculate median for each row with valid values
    row_medians <- vapply(values_by_row, median, numeric(1))
    # Calculate mean for each row with valid values
    row_means <- vapply(values_by_row, mean, numeric(1))
    row_min <- vapply(values_by_row, min, numeric(1))
    row_max <- vapply(values_by_row, max, numeric(1))
    # Insert results at the correct row positions
    medians[as.integer(names(row_medians))] <- row_medians
    means[as.integer(names(row_means))] <- row_means
    min_count[as.integer(names(row_min))] <- row_min
    max_count[as.integer(names(row_max))] <- row_max
  }

  return(data.frame(gene=rownames(dgCMatrix), mean=means, median=medians, min=min_count, max=max_count))
}
a <- rowstats_dgCMatrix(mats_filtered[[1]]@assays$RNA$counts)
b <- rowstats_dgCMatrix(mats_filtered[[length(mats_filtered)]]@assays$RNA$counts)
tmp <- merge(a, b, by="gene")

# kallisto
p1 <- ggplot(tmp, aes(x=log(mean.x), y=log(mean.y))) + geom_pointdensity(size=1, adjust=0.5) +
  scale_color_viridis(name="Density\n(neighbours)", labels = label_comma()) +
  xlab("Jitter 0 mean count per gene (log10)") + ylab(paste0("Jitter ", length(raw)-1, " mean count per gene (log10)")) +
  geom_abline(colour="red", lty=2) +
  annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, fill = "white", alpha = 0.7, size = 3,
           label = paste("Above abline:", sum(tmp$mean.y > tmp$mean.x, na.rm=T),
                         "\nBelow abline:", sum(tmp$mean.y < tmp$mean.x, na.rm=T))) +
  theme_bw() + theme(text = element_text(size=fontsize), legend.position = "right")
p2 <- ggplot(tmp, aes(x=log(max.x), y=log(max.y))) + geom_pointdensity(size=1, adjust=0.5) +
  scale_color_viridis(name="Density\n(neighbours)", labels = label_comma()) +
  xlab("Jitter 0 max count per gene (log10)") + ylab(paste0("Jitter ", length(raw)-1, " max count per gene (log10)")) +
  geom_abline(colour="red", lty=2) +
  annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, fill = "white", alpha = 0.7, size = 3,
           label = paste("Above abline:", sum(tmp$max.y > tmp$max.x, na.rm=T),
                         "\nBelow abline:", sum(tmp$max.y < tmp$max.x, na.rm=T))) +
  theme_bw() + theme(text = element_text(size=fontsize), legend.position = "right")

ggarrange(p1,p2, labels=c("A", "B"), nrow=2, common.legend = T, legend = "right") %>%
  ggsave(width = 4, height = 6, units="in", dpi=300, bg = "white",
         filename = paste0("plots/", names(raw)[1], "-", names(raw)[length(raw)], "_counts.png"))
```

This returns the following image:

<img src="./img/Parse-WTv2-J0-Parse-WTv2-J2_counts.png" alt="Parse jitter 0 versus jitter 2 gene counts"/>
