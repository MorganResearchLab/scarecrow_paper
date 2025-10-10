# scarecrow

This repo outlines code and analysis presented in the [scarecrow](https://github.com/MorganResearchLab/scarecrow) paper. BASH scripts mentioned in this doc are provided in the repo at `scarecrow/src/HPC`.

## Data processing
##### BASH

We analysed Parse Evercode (WTv2) (SRA accession: SRR28867557) and Scale Biosciences QuantumnScale (SRA accession: SRR28867558) data. Once downloaded, the data was processed using the default set-based method of barcode matching with `scarecrow`. Details are provided in the `scarecrow` docs for processing the [Evercode](https://github.com/MorganResearchLab/scarecrow/blob/main/docs/example_evercode.md) and [QuantumnScale](https://github.com/MorganResearchLab/scarecrow/blob/main/docs/example_scale.md) data from downloading through to count matrix generation. The Evercode data was processed to allow upto 2 mismatches, for jitter values of 0, 1, and 2. The QuantumnScale data was processed to allow upto 2 mismatches, for jitter values of 0 and 1.


## Process Kallisto count matrices
##### R

The below code block imports the required libraries, adds a couple of functions, imports the gene-transcript table, sets constant parameters and imports the kallisto count matrices for a given dataset into a list. The list is then processed to generate raw metrics plots, outputting matrix dimensions (features x cells) and total UMIs (for Table S2) in the process. The resulting plots (e.g. Figure S2) are then written to file.

```R
# Data wrangling functions
library(dplyr)
library(tidyverse)
library(SingleCellExperiment)
library(mgcv)                    # this is required for the big additive model (BAM) (e.g. middle panel in Figure S8)

# Functions involved in plotting
library(ggplot2)
library(ggpointdensity)
library(ggdensity)
library(ggdist)
library(ggpubr)
library(ggsci)
library(viridis)
library(scales)

# Function to read Kallisto count matrices
read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- Matrix::readMM(paste0(dir, "/", name, ".mtx"))
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
# https://pachterlab.github.io/kallistobustools/tutorials/kb_building_atlas/R/kb_analysis_0_R
get_knee_df <- function(mat) {
  total <- rank <- NULL
  tibble(
    total = Matrix::colSums(mat),
    rank = dplyr::row_number(dplyr::desc(total))
  ) %>%
    distinct() %>%
    dplyr::filter(total > 0) %>%
    arrange(rank)
}

# Minimum total UMI counts for barcode for it to be considered when calculating the inflection point;
# this helps to avoid the noisy part of the curve for barcodes with very few counts.
# https://pachterlab.github.io/kallistobustools/tutorials/kb_building_atlas/R/kb_analysis_0_R
get_inflection <- function(df, lower = 100) {
  log_total <- log_rank <- total <- NULL
  df_fit <- df %>%
    dplyr::filter(total >= lower) %>%
    transmute(
      log_total = log10(total),
      log_rank = log10(rank)
    )
  d1n <- diff(df_fit$log_total) / diff(df_fit$log_rank)
  right.edge <- which.min(d1n)
  10^(df_fit$log_total[right.edge])
}

# Constants
setwd("~/Documents/scarecrow/examples") # folder that contains data
fontsize <- 10                          # font size on plots
nmads <- 3                              # number of median absolute deviations for mtDNA filter
minFeatures <- 1                        # min number of features (genes) per cell
minUMIs <- 100                          # min number of UMIs per cell

# Import gene to transcript table
tr2g <- readr::read_tsv("transcripts_to_genes.txt", col_names = c("transcript", "gene", "gene_name"))
tr2g <- distinct(tr2g[, c("gene", "gene_name")])
tr2g[is.na(tr2g$gene_name),]$gene_name <- tr2g[is.na(tr2g$gene_name),]$gene
# uniquify the duplicates (some ensembl gene IDs map to the same gene symbol)
gene_name_dups <- tr2g$gene_name[which(duplicated(tr2g$gene_name))]
for(n in gene_name_dups) {
  tr2g[which(tr2g$gene_name == n),]$gene_name <- paste(tr2g[which(tr2g$gene_name == n),]$gene_name, tr2g[which(tr2g$gene_name == n),]$gene, sep="-")
}

# Note below, only one dataset is imported for processing at a time, either Parse or Scale. The associated files are imported into a the list 'raw'.

# Parse Evercode Kallisto count matrices
raw <- list(read_count_output("./Parse/kallisto/J0M2/SRR28867558_1_trimmed/counts_unfiltered", name = "cells_x_genes"),
            read_count_output("./Parse/kallisto/J1M2/SRR28867558_1_trimmed/counts_unfiltered", name = "cells_x_genes"),
            read_count_output("./Parse/kallisto/J2M2/SRR28867558_1_trimmed/counts_unfiltered", name = "cells_x_genes"))
names(raw) <- c("Parse-J0", "Parse-J1", "Parse-J2")

# Scale QuantumnScale Kallisto count matrices
raw <- list(read_count_output("./Scale/kallisto/J0M2/SRR28867557_1_trimmed/counts_unfiltered", name = "cells_x_genes"),
            read_count_output("./Scale/kallisto/J1M2/SRR28867557_1_trimmed/counts_unfiltered", name = "cells_x_genes"))
names(raw) <- c("Scale-J0", "Scale-J1")

# Process count matrix
mats <- lapply(raw, function(res_mat) {

  cat(".... pre-filter dimensions\n")
  print(dim(res_mat))
  cat(".... UMI count\n")
  print(sum(res_mat))

  # Generate plot to check for library saturation
  tot_counts <- Matrix::colSums(res_mat)
  lib_sat <- tidyr::tibble(nCount = tot_counts, nGene = Matrix::colSums(res_mat > 0))
  Fig1A <- ggplot(lib_sat, aes(nCount, nGene)) +
    geom_pointdensity(size=1, adjust=2) + scale_color_viridis(name="density\n(neighbours)", labels = label_comma()) +
    scale_x_log10(labels = label_comma()) + scale_y_log10(labels = label_comma()) + annotation_logticks() +
    xlab("Number of UMIs") + ylab("Number of features") +
    theme_bw() + theme(text = element_text(size=fontsize), legend.position = "right")

  # Generate plot to check for inflection point in barcode/UMI count
  knee_df <- get_knee_df(res_mat)
  inflection <- get_inflection(knee_df, lower=100)

  # Knee plot
  total <- rank_cutoff <- NULL
  annot <- tibble(inflection = inflection, rank_cutoff = max(knee_df$rank[knee_df$total > inflection]))
  Fig1B <- ggplot(knee_df, aes(total, rank)) +
    geom_path() +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2, color = "gray40") +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2, color = "gray40") +
    geom_text(aes(inflection, rank_cutoff, label = paste(rank_cutoff, "'cells'")), data = annot, vjust = 1, size = 3) +
    scale_x_log10(labels = label_comma()) +
    scale_y_log10(labels = label_comma()) +
    labs(y = "Rank", x = "Total UMIs") +
    annotation_logticks() + theme_bw() +
    theme(text = element_text(size = fontsize))

  # Filter and return SingleCellExperiment object
  cat("Filtering and outputting SingleCellExperiment (SCE) object\n")
  res <- SingleCellExperiment(assays = list(counts = res_mat[, Matrix::colSums(res_mat) > inflection]))

  # Annotate mitochondrial genes
  mt_genes <- tr2g[grep("^MT-", tr2g$gene_name),]$gene
  res <- scuttle::addPerCellQC(res, subsets=list(Mito=mt_genes[which(mt_genes %in% rownames(res))]))
  qc.percent.mt <- NULL
  if(max(res$subsets_Mito_percent) > 0) {
    cat(".... %MT > 0, running scater::isOutlier for:", nmads,"median absolute deviations\n")
    qc.percent.mt <- scater::isOutlier(res$subsets_Mito_percent, nmads = nmads, type="higher")
  }

  # Plot: Features per cell
  Fig1Ca <- ggplot(data.frame(X=res$detected), aes(y=X, slab_color = after_stat(y))) +
    stat_slabinterval(side = "right", lty=0) +
    stat_dotsinterval(slab_shape = 19, quantiles =100, side="left") +
    scale_color_distiller(aesthetics = "slab_color", guide = "colorbar2") +
    ylab("Features per cell") +
    scale_y_continuous(labels = label_comma()) +
    theme_bw() + theme(text = element_text(size=fontsize),
                       axis.text.y = element_blank(), axis.title.y = element_blank(),
                       legend.position = "None") + coord_flip() +
    theme(plot.margin = margin(l = 0.75, r = 0.5, unit = "cm"))

  # Plot: Molecules per cell
  Fig1Cb <- ggplot(data.frame(X=res$sum), aes(y=X, slab_color = after_stat(y))) +
    stat_slabinterval(side = "right", lty=0) +
    stat_dotsinterval(slab_shape = 19, quantiles =100, side="left") +
    scale_color_distiller(aesthetics = "slab_color", guide = "colorbar2") +
    ylab("Molecules per cell") +
    scale_y_continuous(labels = label_comma()) +
    theme_bw() + theme(text = element_text(size=fontsize),
                       axis.text.y = element_blank(), axis.title.y = element_blank(),
                       legend.position = "None") + coord_flip() +
    theme(plot.margin = margin(l = 0.75, r = 0.5, unit = "cm"))

  # % mitochondrial genes per cell
  Fig1Cc <- ggplot(data.frame(X=res$subsets_Mito_percent), aes(y=X, slab_color = after_stat(y))) +
    stat_slabinterval(side = "right", lty=0) +
    stat_dotsinterval(slab_shape = 19, quantiles =100, side="left") +
    scale_color_distiller(aesthetics = "slab_color", guide = "colorbar2") +
    ylab("% mitochondrial genes") +
    theme_bw() + theme(text = element_text(size=fontsize),
                       axis.text.y = element_blank(), axis.title.y = element_blank(),
                       legend.position = "None") + coord_flip() +
    theme(plot.margin = margin(l = 0.75, r = 0.5, unit = "cm"))
  Fig1C <- ggarrange(Fig1Ca, Fig1Cb, Fig1Cc, nrow=3)
  rm(Fig1Ca, Fig1Cb, Fig1Cc)

  # % mitochondrial genes per cell x molecules per cell
  Fig1D <- ggplot(data.frame(X=res$sum, Y=res$subsets_Mito_percent), aes(x=X, y=Y)) +
    geom_pointdensity(size=1, adjust=5) + scale_color_viridis(name="density\n(neighbours)", labels = label_comma()) +
    xlab("Molecules per cell") + ylab("% mitochondrial genes") +
    scale_x_continuous(labels = label_comma()) +
    geom_hline(yintercept=min(res$subsets_Mito_percent[which(qc.percent.mt==T)]), lty=2) +
    annotate("text", x=max(res$sum)/2, y=min(res$subsets_Mito_percent[which(qc.percent.mt==T)]),
             hjust=0, vjust=-1, size=3, label = paste0("MADS > ", nmads)) +
    theme_bw() + theme(text = element_text(size=fontsize), legend.position = "right")

  # Figure 1
  Fig1 <- ggarrange(ggarrange(Fig1A, Fig1B, widths=c(4,4), labels=c("A", "B"), nrow=1),
                    ggarrange(Fig1C, Fig1D, widths=c(3,5), labels=c("C", "D"), nrow=1),nrow=2)
  rm(Fig1A, Fig1B, Fig1C, Fig1D)

  # Filter data
  cat(".... pre-MT filter dimensions\n")
  print(dim(res))

  if(max(res$subsets_Mito_percent) > 0) {
    cat(".... %MT > 0, running scater::isOutlier for:", nmads,"median absolute deviations\n")
    res <- res[, which(!scater::isOutlier(res$subsets_Mito_percent, nmads = nmads, type="higher"))]
  }
  cat("\t..identifying cells with sufficient UMI count:", minUMIs, "\n")
  sufficient_UMIs <- Matrix::colSums(counts(res)) >= minUMIs
  cat("\t..identifying cells with sufficient feature count:", minFeatures, "\n")
  sufficient_Features <- Matrix::rowSums(counts(res)) >= minFeatures
  res <- res[which(sufficient_Features), which(sufficient_UMIs)]

  cat(".... post filter dimensions\n")
  print(dim(res))

  return(list(mat = res, plot = Fig1))

})
names(mats) <- names(raw)

# Save mats as these will be used for analysis later
save(mats, file=paste0(gsub("-.*", "", names(raw)[1]), ".RData"))

# Output plots to file (e.g. Figure S2)
cat("Outputting metrics plots\n")
lapply(1:length(mats), function(i) {
  ggsave(mats[[i]]$plot, width = 8, height = 6, units="in", dpi=300, filename = paste0("raw_metrics_", names(mats)[i], ".png"))
})
```

The below code block processes the list of matrices `mats` to return the median genes/cell, median UMIs/cell, and max UMIs/gene (for Table S3). The row- and column-wise values are then used to generate plots comparing mean and max gene counts with and without jitter (e.g. Figure S6).

```R
# Median genes/cell (Table S3)
genes_per_cell.median <- lapply(1:length(mats), function(i) median(Matrix::colSums(counts(mats[[i]]$mat)>0)))
# Median umis/cell (Table S3)
umis_per_cell.median <- lapply(1:length(mats), function(i) median(Matrix::colSums(counts(mats[[i]]$mat))))

# Calculate row-wise and col-wise stats of dgCMatrix
sparse_stats_dgCMatrix <- function(dgCMatrix) {
  stopifnot(inherits(dgCMatrix, "dgCMatrix"))

  # Convert to triplet form (i = row, j = col, x = value)
  mat_triplet <- summary(dgCMatrix)

  # Clean up x values
  if (is.character(mat_triplet$x)) {
    mat_triplet$x[mat_triplet$x == "."] <- NA
    mat_triplet$x <- as.numeric(mat_triplet$x)
  }

  # Filter non-zero and non-NA entries
  valid_values <- mat_triplet[mat_triplet$x != 0 & !is.na(mat_triplet$x), ]

  ## Row-wise stats (gene-level)
  row_medians <- row_means <- row_min <- row_max <- rep(NA_real_, nrow(dgCMatrix))
  if (nrow(valid_values) > 0) {
    row_split <- split(valid_values$x, valid_values$i)
    row_medians_partial <- vapply(row_split, median, numeric(1))
    row_means_partial   <- vapply(row_split, mean, numeric(1))
    row_min_partial     <- vapply(row_split, min, numeric(1))
    row_max_partial     <- vapply(row_split, max, numeric(1))
    row_indices <- as.integer(names(row_medians_partial))
    row_medians[row_indices] <- row_medians_partial
    row_means[row_indices]   <- row_means_partial
    row_min[row_indices]     <- row_min_partial
    row_max[row_indices]     <- row_max_partial
  }

  row_stats <- data.frame(
    gene = rownames(dgCMatrix),
    mean = row_means,
    median = row_medians,
    min = row_min,
    max = row_max,
    stringsAsFactors = FALSE
  )

  ## Column-wise stats (cell-level)
  col_medians <- col_means <- col_min <- col_max <- col_total <- rep(NA_real_, ncol(dgCMatrix))
  if (nrow(valid_values) > 0) {
    col_split <- split(valid_values$x, valid_values$j)
    col_medians_partial <- vapply(col_split, median, numeric(1))
    col_means_partial   <- vapply(col_split, mean, numeric(1))
    col_min_partial     <- vapply(col_split, min, numeric(1))
    col_max_partial     <- vapply(col_split, max, numeric(1))
    col_total_partial     <- vapply(col_split, sum, numeric(1))
    col_indices <- as.integer(names(col_medians_partial))
    col_medians[col_indices] <- col_medians_partial
    col_means[col_indices]   <- col_means_partial
    col_min[col_indices]     <- col_min_partial
    col_max[col_indices]     <- col_max_partial
    col_total[col_indices]     <- col_total_partial
  }

  col_stats <- data.frame(
    cell = colnames(dgCMatrix),
    mean = col_means,
    median = col_medians,
    min = col_min,
    max = col_max,
    total = col_total,
    stringsAsFactors = FALSE
  )

  return(list(row_stats = row_stats, col_stats = col_stats))
}
mat_stats <-  lapply(1:length(mats), function(i) sparse_stats_dgCMatrix(counts(mats[[i]]$mat)))
# Max UMIs/gene (Table S3)
umis_per_cell.max <- lapply(1:length(mats), function(i) max(mat_stats[[i]]$row_stats$max))

# Generate plots of mean and max counts per gene with and without jitter (e.g. Figure S6)
tmp <- merge(mat_stats[[1]]$row_stats, mat_stats[[length(mat_stats)]]$row_stats, by="gene")
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
         filename = paste0(names(mats)[1], "-", names(mats)[length(mats)], "_counts.png"))
```

## Evaluating jitter-dependent read-barcode re-assignment
##### R

First we compare the difference in UMI count per unique barcode combination reported by `scarecrow stats` for data processed with and without jitter. Here we process the Scale and Parse datasets independently, comparing jitter 0 to 1 and 0 to 2, respectively. The counts from each jitter are merged on unique barcode and plotted (e.g. Figure S8).

```R
# Parse datasets for jitter = 0 and jitter = 2
file_a <- read.table("./Parse/extracted/J0M2/SRR28867558_1_trimmed.fastq.combined_CB_counts.txt")
file_b <- read.table("./Parse/extracted/J2M2/SRR28867558_1_trimmed.fastq.combined_CB_counts.txt")

# Scale datasets for jitter = 0 and jitter = 1
file_a <- read.table("./Scale/extracted/J0M2/SRR28867557_1_trimmed.fastq.combined_CB_counts.txt")
file_b <- read.table("./Scale/extracted/J1M2/SRR28867557_1_trimmed.fastq.combined_CB_counts.txt")

dataset <- "Parse" # or "Scale"
jitter <- c(0,2)   # Set labels for jitter values of data, e.g. Scale is c(0,1), Parse is c(0,2)

dat <- merge(file_a, file_b, by="V1", all=T) # Merge data on barcode combination
dat[is.na(dat)] <- 0
dat$diff <- dat$V2.y-dat$V2.x
dat$set <- "Intersect"
dat[which(dat$V2.y == 0),]$set <- paste0("J", jitter[1])
dat[which(dat$V2.x == 0),]$set <- paste0("J", jitter[2])
dat$set <- factor(dat$set, levels=c(paste0("J", jitter[1]), "Intersect", paste0("J", jitter[2])))

p1 <- ggplot(dat[which(dat$set==paste0("J", jitter[1])),], aes(x=log10(V2.x), y=diff)) +
  geom_line(colour="red", lty=2, alpha=0.8) + geom_point(size=1, alpha=0.5) +
  xlab(paste("Jitter", jitter[1], "UMIs (log10)")) + ylab("Δ UMI") +
  scale_x_continuous(labels = label_comma()) + scale_y_continuous(labels = label_comma()) +
  theme_bw() + theme(text = element_text(size=fontsize), legend.position = "bottom")

# bam (big additive model)
fit <- bam(diff ~ V2.x + s(V2.x, bs = "tp"), data = dat[which(dat$set=="Intersect"),])
p2 <- ggplot(dat[which(dat$set=="Intersect"),], aes(x=log10(V2.x), y=diff)) +
  geom_line(colour="red", lty=2, alpha=0.8) + geom_point(size=1, alpha=0.5) +
  geom_smooth(method = "bam", formula = y ~ x + s(x, bs="tp"), color = "blue") +
  xlab(paste("Jitter", jitter[1], "UMIs (log10)")) + ylab("Δ UMI") +
  annotate("text", x = -Inf, y = max(dat[which(dat$set=="Intersect"),]$diff),
           label = sprintf("\tSlope: %.3f\n\tIntercept: %.3f", coef(fit)["V2.x"], coef(fit)["(Intercept)"]),
           hjust = 0, vjust = 1, size = 3) +
  scale_x_continuous(labels = label_comma()) + scale_y_continuous(labels = label_comma()) +
  theme_bw() + theme(text = element_text(size=fontsize), legend.position = "bottom")

p3 <- ggplot(dat[which(dat$set==paste0("J", jitter[2])),], aes(x=log10(V2.y), y=diff)) +
  geom_line(colour="red", lty=2, alpha=0.8) + geom_point(size=1, alpha=0.5) +
  xlab(paste("Jitter", jitter[2], "UMIs (log10)")) + ylab("Δ UMI") +
  scale_x_continuous(labels = label_comma()) + scale_y_continuous(labels = label_comma()) +
  theme_bw() + theme(text = element_text(size=fontsize), legend.position = "bottom")

ggarrange(p1,p2,p3, nrow=1, common.legend=F, widths=c(1,3,1)) %>%
  ggsave(width = 12, height = 6, units="in", dpi=300, filename = paste0(dataset, "_", jitter[1], "-", jitter[2], "_UMIs.png"))
```

## Evaluating jitter-dependent read-barcode re-assignment
##### BASH

First, we extract read name and SAM tags (CB, XP, XM) from BAM the files to compare.

```bash
mamba activate samtools_env

PROJECT=./scarecrow/examples/Scale
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 8 --mem 4G --time=12:00:00 \
    ./scarecrow/scripts/samtags.sh --sam ${PROJECT}/extracted/J0M2/SRR28867557_1_trimmed.sam
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 8 --mem 4G --time=12:00:00 \
    ./scarecrow/scripts/samtags.sh --sam ${PROJECT}/extracted/J1M2/SRR28867557_1_trimmed.sam

PROJECT=./scarecrow/examples/Parse
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 8 --mem 4G --time=12:00:00 \
    ./scarecrow/scripts/samtags.sh --sam ${PROJECT}/extracted/J0M2/SRR28867558_1_trimmed.sam
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 8 --mem 4G --time=12:00:00 \
    ./scarecrow/scripts/samtags.sh --sam ${PROJECT}/extracted/J2M2/SRR28867558_1_trimmed.sam
```

Next, we sample n barcodes from the BAM[1] read tags, retrieve read names, and subset these reads from the BAM[2] read tags. The result is that we now have for the same subset of reads, the read tags recorded from each BAM file.

```bash
PROJECT=./scarecrow/examples/Scale
sbatch -p uoa-compute ./scarecrow/scripts/read_reassignment.sh \
    --a ${PROJECT}/extracted/J0M2/SRR28867557_1_trimmed.sam.tags \
    --b ${PROJECT}/extracted/J1M2/SRR28867557_1_trimmed.sam.tags \
    --n 10000

PROJECT=./scarecrow/examples/Parse
sbatch -p uoa-compute ./scarecrow/scripts/read_reassignment.sh \
    --a ${PROJECT}/extracted/J0M2/SRR28867558_1_trimmed.sam.tags \
    --b ${PROJECT}/extracted/J2M2/SRR28867558_1_trimmed.sam.tags \
    --n 10000
```

The resulting files should have a format as follows:

```bash
SRR28867557.20120	CB:Z:TCCGGCTTAT_TCAGCGGTT_TTATCCGGAT	XP:Z:1_1_24	XM:Z:0_0_0
SRR28867557.20159	CB:Z:TCCGGCTTAT_GACTGACGT_AACCTGCGTA	XP:Z:1_1_24	XM:Z:0_0_0
SRR28867557.20165	CB:Z:TCCGGCTTAT_TTCGCGGAT_ACTTGCTAGA	XP:Z:1_1_24	XM:Z:0_0_0
```


## Evaluating jitter-dependent read-barcode re-assignment
##### R

After generating the above SAM tag TSV files for a subset of reads, each dataset is analysed in R independently. The values from these analyses populate Table S4.

```R
# Scale datasets for jitter = 0 and jitter = 1
file_a <- read.table("./Scale/extracted/J0M2/SRR28867557_1_trimmed.sam.tags.subset_10000.tsv")
file_b <- read.table("./Scale/extracted/J1M2/SRR28867557_1_trimmed.sam.tags.subset_10000.tsv")

# Parse datasets for jitter = 0 and jitter = 2
file_a <- read.table("./Parse/extracted/J0M2/SRR28867558_1_trimmed.sam.tags.subset_10000.tsv")
file_b <- read.table("./Parse/extracted/J2M2/SRR28867558_1_trimmed.sam.tags.subset_10000.tsv")

dat <- merge(file_a, file_b, all=T, by="V1") # merge reads on read name
# Extract barcode index positions and mismatch counts
a <- do.call("rbind", lapply(dat$V3.x, function(i) as.numeric(unlist(strsplit(unlist(strsplit(i, "XP:Z:"))[2], "_")))))
b <- do.call("rbind", lapply(dat$V4.x, function(i) as.numeric(unlist(strsplit(unlist(strsplit(i, "XM:Z:"))[2], "_")))))
c <- do.call("rbind", lapply(dat$V3.y, function(i) as.numeric(unlist(strsplit(unlist(strsplit(i, "XP:Z:"))[2], "_")))))
d <- do.call("rbind", lapply(dat$V4.y, function(i) as.numeric(unlist(strsplit(unlist(strsplit(i, "XM:Z:"))[2], "_")))))
barcode_Xvals <- cbind(a, b, c, d)
colnames(barcode_Xvals) <- c(paste("XP", seq(1,ncol(barcode_Xvals)/4), ".x", sep=""),
                             paste("XM", seq(1, ncol(barcode_Xvals)/4), ".x", sep=""),
                             paste("XP", seq(1, ncol(barcode_Xvals)/4), ".y", sep=""),
                             paste("XM", seq(1, ncol(barcode_Xvals)/4), ".y", sep=""))
dat <- cbind(dat, barcode_Xvals)
dat$mismatches.x <- rowSums(dat[, grep("XM[0-9].x", colnames(dat))])
dat$mismatches.y <- rowSums(dat[, grep("XM[0-9].y", colnames(dat))])

# Barcode counts
length(unique(file_a$V2)) # Unique barcodes at jitter = 0 (should be the sampled amount, i.e. 10k)
length(unique(file_b$V2)) # Unique barcodes at jitter > 0

# Barcode re-assignments (FALSE = Different, TRUE = Identical, N/A = N/A)
table(dat$V2.x == dat$V2.y, useNA = "ifany") # counts
(table(dat$V2.x == dat$V2.y, useNA = "ifany") / sum(table(dat$V2.x == dat$V2.y, useNA = "ifany") ))*100 # as percentage

# Count of perfect barcode matches at j > 0 where barcode is different at j0
reassigned <- data.frame(perfect = nrow(dat[which(dat$V2.x != dat$V2.y & dat$mismatches.y == 0),]),
                         imperfect = nrow(dat[which(dat$V2.x != dat$V2.y & dat$mismatches.y > 0),]))
reassigned / sum(reassigned) # counts as percentage

# Mean mismatches
mean(dat$XM1.x) # Mean mismatches for BC1 at jitter = 0
mean(dat$XM2.x) # Mean mismatches for BC2 at jitter = 0
mean(dat$XM3.x) # Mean mismatches for BC3 at jitter = 0
mean(dat$XM1.y, na.rm=T) # Mean mismatches for BC1 at jitter != 0
mean(dat$XM2.y, na.rm=T) # Mean mismatches for BC2 at jitter != 0
mean(dat$XM3.y, na.rm=T) # Mean mismatches for BC3 at jitter != 0
```


## Checking the impact of jitter on kNN clustering
##### BASH

The RData saved for each matrix above is passed to an R script to undertake kNN analysis for a range of k values. The script filters the first and last matrix to retain cells and genes with > 99 counts, and then reduces each matrix to their common cells and genes. The reduced matrices are then processed using a fast nearest-neighbour searching algorithm to record nearest neighbours for a range of k values (2:10, 20, 30, 50). The results are saved as an RData file for downstream analysis.

```bash
sbatch --partition uoa-compute ./scarecrow/scripts/kNN.sh --in ./scarecrow/examples/Parse/Parse.RData # 2720536
sbatch --partition uoa-compute ./scarecrow/scripts/kNN.sh --in ./scarecrow/examples/Scale/Scale.RData # 2720530
```

##### R

The resulting RData processed in R as follows to generate plots (e.g. Figure S10).

```R
# Load kNN RData created in previous setp
load("./Scale_knn.RData")

# Calculate kNN shared neighbour proportions
knn_props <- lapply(1:length(knn_results), function(k) {
  knn_shared <- do.call("rbind", lapply(1:nrow(knn_results[[k]][[1]]$nn.index), function(i) {
    sum(knn_results[[k]][[1]]$nn.index[i,] %in% knn_results[[k]][[2]]$nn.index[i,])
  } ))
  knn_identical <- do.call("rbind", lapply(1:nrow(knn_results[[k]][[1]]$nn.index), function(i) {
    sum(knn_results[[k]][[1]]$nn.index[i,] == knn_results[[k]][[2]]$nn.index[i,])
  } ))
  return(data.frame(k = dim(knn_results[[k]][[1]]$nn.index)[2],
                    shared = knn_shared/dim(knn_results[[k]][[1]]$nn.index)[2],
                    identical = knn_identical/dim(knn_results[[k]][[1]]$nn.index)[2]))
})
knn_props <- do.call("rbind", knn_props)
knn_props_long <- rbind(
  data.frame(k = knn_props$k, value = knn_props$shared, group = "shared"),
  data.frame(k = knn_props$k, value = knn_props$identical, group = "identical")
)

# HDR plot
set.seed(123)  # For reproducibility
knn_props_long$x <- knn_props_long$k + runif(nrow(knn_props_long), -0.00000001, 0.00000001)
knn_props_long$facet <- 1
knn_props_long[which(knn_props_long$k>10),]$facet <- knn_props_long[which(knn_props_long$k>10),]$k
plots <- lapply(unique(knn_props_long$facet), function(f) {
  ggplot(knn_props_long[which(knn_props_long$group=="shared" & knn_props_long$facet == f),], aes(x = x, y = value)) +
    geom_hdr(aes(fill=after_stat(probs)), alpha=1) +
    scale_x_continuous(breaks = seq(1:50)) + ylim(0,1) +
    ylab("Proportion shared") + xlab("k-Nearest neighbours") +
    theme_bw() + theme(text = element_text(size = fontsize), legend.position = "top",
                       plot.margin = margin(l = 0.75, r = 0, unit = "cm"),
                       axis.text.x = element_text(size=6),
                       axis.text.y = element_text(size=6))
})
# Remove y-axis elements for all but the first
plots[-1] <- lapply(plots[-1], function(p) {
  p + xlab("") + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       plot.margin = margin(l = 0, r = 0, unit = "cm"))
})

# Quantile bar plot
qdat <- knn_props %>% group_by(k) %>%
  summarise(
    q25 = quantile(shared, 0.25),
    q50 = quantile(shared, 0.5),  # median
    q75 = quantile(shared, 0.75),
    q100 = quantile(shared, 1.00)
  ) %>%
  pivot_longer(
    cols = starts_with("q"),         # or just: -group
    names_to = "quantile",
    values_to = "value"
  )
qdat$quantile <- as.numeric(gsub("q", "", qdat$quantile))
p2 <- ggplot(qdat, aes(x=quantile, y=value, fill="a")) + geom_col(position="dodge") + facet_grid(.~k) +
  scale_fill_jco() + scale_x_continuous(breaks=c(0,25,50,75,100)) +
  ylab("Proportion of cells") + xlab("Quantile of proportion shared") +
  theme_bw() + theme(text = element_text(size = fontsize),
                     axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6),
                     axis.text.y = element_text(size=6),
                     legend.position="none")

# Combine plots
ggarrange(ggarrange(plotlist=plots, nrow=1, widths=c(1.5,0.1,0.1,0.1), common.legend=T),
          p2, nrow=2, heights = c(2,1))

# Mean proportion shared neighbours by k value
knn_props_long %>% filter(group=="shared") %>% group_by(k) %>% summarize(mean(value))
```


## Exploring the impact of jitter on a typical scRNA-seq workflow
##### R

To explore any possible impact on the application of jitter on barcode matching by processing the data in a typical workflow. We start by loading the RData saved after processing the Kallisto count matrices. This data is further processed by log-normalising the counts, idenentifying the top 10% of highly variable genes, running PCA and clustering followed by doublet-removal.

```R
dataset <- "Parse"
load(paste0(dataset, ".RData"))


# Remove doublets
mats_clean <- lapply(mats, function(dat) {
  set.seed(seed)

  # UMI filter
  sufficient_UMIs <- Matrix::colSums(counts(dat$mat)) >= 500
  dat$mat <- dat$mat[, which(sufficient_UMIs)]

  # Log normalisation
  cat("\t..log-normalising counts data\n")
  dat <- scuttle::logNormCounts(dat$mat)

  # Feature selection
  cat("\t..identifying HVGs\n")
  dec <- scran::modelGeneVar(dat)
  hvg <- scran::getTopHVGs(dec, prop=0.1)

  # PCA
  cat("\t..running PCA\n")
  dat <- scater::runPCA(dat)

  # Clustering
  cat("\t..clustering cells based on PCA using Louvain function\n")
  colLabels(dat) <- scran::clusterCells(dat, use.dimred='PCA', BLUSPARAM=bluster::NNGraphParam(cluster.fun="louvain"))

  # Doublets: Simulation-based (doublet densities)
  cat("\t..computing doublet densiies (scDblFinder)\n")
  set.seed(seed)
  dbl.dens <- scDblFinder::computeDoubletDensity(dat, subset.row=hvg, d=ncol(reducedDim(dat)))
  dat$DoubletScore <- dbl.dens
  dbl.calls <- scDblFinder::doubletThresholding(data.frame(score=dbl.dens), method="griffiths", returnType="call", perSample = FALSE)
  summary(dbl.calls)

  # Doublets: Simulation-based (classification)
  cat("\t..combining doublet density with iterative classification scheme (scDblFinder)\n")
  set.seed(seed)
  dat.dbl <- scDblFinder::scDblFinder(dat, clusters=colLabels(dat), nfeatures = hvg)
  summary(dat.dbl$scDblFinder.class)

  # Union of doublets
  dbl.union <- unique(c(which(dbl.calls=="doublet"), which(dat.dbl$scDblFinder.class=="doublet")))
  cat("\t..UNION of doublet finding methods accounts for", signif(length(dbl.union) / dim(dat)[2], 3)*100, "% of data\n")

  # Return doublet-filtered data
  return(dat[,-dbl.union])
})

# Identify common genes and barcodes between singlet count matrices
mat.a <- mats_clean[[1]]
mat.b <- mats_clean[[length(mats)]] # length(mats) = 2 for Scale; length(mats) = 3 for Parse
barcodes <- colnames(mat.a)[which(colnames(mat.a) %in% colnames(mat.b))]
genes <- rownames(mat.a)[which(rownames(mat.a) %in% rownames(mat.b))]
singlets <- list(mat.a[genes, barcodes], mat.b[genes, barcodes])

# Normalise and cluster
mats_proc <- lapply(singlets, function(dat) {
  set.seed(seed)

  # Log normalisation
  cat("\t..log-normalising counts data\n")
  dat <- scuttle::logNormCounts(dat)

  # PCA
  cat("\t..running PCA\n")
  dat <- scater::runPCA(dat)

  # Clustering
  cat("\t..clustering cells based on PCA using Louvain function\n")
  colLabels(dat) <- scran::clusterCells(dat, use.dimred='PCA', BLUSPARAM=bluster::NNGraphParam(cluster.fun="louvain"))

  # Build neighbor graph
  cat("\t..building nearest-neighbour graph\n")
  g <- scran::buildSNNGraph(dat, use.dimred = "PCA", k = 50)

  # Clustering (Louvain; resolution controlled via cluster parameter)
  cat("\t..clustering\n")
  clusters <- igraph::cluster_louvain(g)$membership
  dat$cluster <- factor(clusters)

  # UMAP
  cat("\t..running UMAP\n")
  dat <- scater::runUMAP(dat, dimred = "PCA")

  return(dat)
})

# Plot UMAPs
p1 <- scater::plotReducedDim(mats_proc[[1]], "UMAP", colour_by = "cluster") + scale_color_jco(name="cluster")
p2 <- scater::plotReducedDim(mats_proc[[2]], "UMAP", colour_by = "cluster") + scale_color_jco(name="cluster")
ggarrange(p1, p2, labels=c("A", "B")) %>%
  ggsave(width = 12, height = 6, units="in", dpi=300, bg = "white",
    filename = paste0(dataset, "_umap.png"))

# Compute adjusted Rand index between the two clustering results
adjRand <- mclust::adjustedRandIndex(factor(mats_proc[[2]]$cluster), factor(mats_proc[[1]]$cluster))
```


## Main paper figure
##### R

The plots on the main figure in the paper were generated based on the Scale data. The UMAP (Fig 1C) generation and adjusted Rand index values are reported from the previous code block. The data underpinning the bar plots (Fig 1B) is available from the supplementary tables (Table S1, Table S2, Table S3). The heatmap of cluster memberships was generated as follows:

```R
# Generate cluster table for plotting heatmap
clusters <- as.data.frame(table(Jitter0 = colData(mats_proc[[1]])$cluster,
                                Jitter1 = colData(mats_proc[[2]])$cluster))
# Normalize
clusters <- clusters %>% group_by(Jitter0) %>% mutate(Prop = Freq / sum(Freq)) %>% ungroup()

# Plot heatmap
ggplot(clusters, aes(x = Jitter1, y = Jitter0, fill = Prop)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#DDDDFF", high = "steelblue")+
  geom_text(aes(label = Freq), size = 2) +  # show actual counts
  labs(x = "Jitter 1 cluster assignments", y = "Jitter 0 cluster assignments", fill = "Proportion") +
  labs(subtitle=paste0("Adjusted Rand index: ", signif(adjRand, 2))) +
  theme_bw() + theme(text = element_text(size=8), plot.subtitle = element_text(size=7))
```
