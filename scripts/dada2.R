#!/usr/bin/env Rscript

library(ggplot2)
library(dada2)

# Debug ---------------------------------------------------------------------------------------------------------------
# fnFs <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/trimmed/DNA-M1-1_R1.fastq"
# fnRs <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/trimmed/DNA-M1-1_R2.fastq"
# filtFs <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/denoising/DNA-M1-1_R1_filtered.fastq"
# filtRs <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/denoising/DNA-M1-1_R2_filtered.fastq"
# asv_table <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/denoising/DNA-M1-1_ASV.fasta"
# errplotF <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/denoising/DNA-M1-1_R1_R1_error_plot.tiff"
# errplotR <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/denoising/DNA-M1-1_R1_R2_error_plot.tiff"
# denoising_R1 <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/denoising/DNA-M1-1_R1_R1_denoise.txt"
# denoising_R2 <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/denoising/DNA-M1-1_R1_R2_denoise.txt"
# report <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/reports/DNA-M1-1_denoising.tsv"
# merged <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/denoising/DNA-M1-1_merging.txt"
# threads <- 6
# sample.names <- "DNA-M1-1"
# max_EE <- 1
# minlength <- 100
# maxlength <- 120

# Get parameters from snakemake ---------------------------------------------------------------------------------------------

# Input
fnFs <- snakemake@input[["r1"]]
fnRs <- snakemake@input[["r2"]]

# Output
filtFs <- snakemake@output[["r1_filt"]]
filtRs <- snakemake@output[["r2_filt"]]
errplotF <- snakemake@output[["errplotF"]]
errplotR <- snakemake@output[["errplotR"]]
denoising_R1 <- snakemake@outpur[["denoiseR1"]]
denoising_R2 <- snakemake@outpur[["denoiseR2"]]
merged <- snakemake@outpur[["merged"]]
asv_table <- snakemake@output[["asv"]]
report <- snakemake@output[["report"]]

# Multi-threading
threads <- snakemake@threads[[1]]

# Parameters
sample.names <- snakemake@params[["sample"]]
max_EE <- snakemake@params[["maxEE"]
minsize <- snakemake@params[["minsize"]]
minlength <- snakemake@params[["minlength"]]
minlength <- snakemake@params[["maxlength"]]


# DADA Workflow -------------------------------------------------------------------------------------------------------------

# Filter reads
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
              maxN=0, maxEE=c(max_EE,max_EE), rm.phix=TRUE,
              compress=TRUE, multithread=threads)

# Learn Error rate
errF <- learnErrors(filtFs, multithread=threads)
errR <- learnErrors(filtRs, multithread=threads)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=threads)
dadaRs <- dada(filtRs, err=errR, multithread=threads)
 
# Merge Reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# ASV table
seqtab <- makeSequenceTable(mergers)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% minlength:maxlength]

# Remove chimeras
seqtab2.nochim <- removeBimeraDenovo(seqtab, method="per-sample", multithread=threads, verbose=TRUE)

# Export sequences and QC reports -------------------------------------------------------------------------------------------

# ASV table
asv <- t(data.frame(seqtab2.nochim))
colnames(asv) <- "count"
asv <- cbind(asv, name = sprintf("ASV_%s", seq(1:dim(asv)[1])))
asv <- cbind(asv, sequence = rownames(asv))
rownames(asv) <- seq(1:dim(asv)[1])
fn <- function(x) paste0(">", x[2], ";size=", x[1], "\n", x[3])
asfasta <- apply(asv, MARGIN=1, fn)
writeLines(asfasta, asv_table)

# Error plots
errFplot <- plotErrors(errF, nominalQ=TRUE) + theme(text=element_text(family="Arial"))
ggsave(errplotF, width = 10, height = 10)
errRplot <- plotErrors(errF, nominalQ=TRUE) + theme(text=element_text(family="Arial"))
ggsave(errplotR, width = 10, height = 10)

# Denoising reports
write.table(dadaFs$clustering, denoising_R1, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(dadaRs$clustering, denoising_R2, quote = FALSE, sep = "\t", row.names = FALSE)

# Merged reads
write.table(mergers, merged, quote = FALSE, sep = "\t", row.names = FALSE)

# Report
getN <- function(x) sum(getUniques(x))
track <- cbind(sample.names,
				out,
				round((out[1]-out[2])/out[1]*100, 2),
				getN(dadaFs),
				getN(dadaRs),
				getN(mergers),
				round((1-getN(mergers)/out[2])*100, 2),
				length(seqtab),
				length(seqtab2),
				round((length(seqtab)-length(seqtab2))/length(seqtab)*100, 2),
				round((sum(seqtab)-sum(seqtab2))/sum(seqtab)*100, 2),
				length(seqtab2.nochim),
				round((length(seqtab2)-length(seqtab2.nochim))/length(seqtab2)*100, 2),
				round((sum(seqtab2)-sum(seqtab2.nochim))/sum(seqtab2)*100, 2),
				sum(seqtab2),
				round(sum(seqtab2)/out[1]*100, 2))
colnames(track) <- c("Sample",
				"Total reads", "Filtered reads",
				"Discarded reads [%]",
				"Denoised R1",
				"Denoised R2",
				"Merged",
				"Merging failures [%]",
				"ASV",
				"Size-filtered ASV",
				"Discarded ASVs [% ASV]",
				"Discarded ASVs [% of reads]",
				"Non-chimeric ASV",
				"Chimeras [% ASV]",
				"Chimeras [% of reads]",
				"Reads in ASVs",
				"Reads in ASVs [%]")
write.table(track, report, quote = FALSE, sep = "\t", row.names = FALSE)