#!/usr/bin/env Rscript

library(ggplot2, quiet=T)
library(dada2, quiet=T)

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
# chimeras_fasta <- "/home/debian/NGS/spezies_indev/tests_results/test/DNA-M1-1/denoising/DNA-M1-1_chiomeras.fasta"
# threads <- 6
# sample.names <- "DNA-M1-1"
# max_EE <- 1
# minlength <- 100
# maxlength <- 120
# chimera <-TRUE

# Get parameters from snakemake ---------------------------------------------------------------------------------------------

# Input
fnFs <- snakemake@input[["r1"]]
fnRs <- snakemake@input[["r2"]]

# Output
filtFs <- snakemake@output[["r1_filt"]]
filtRs <- snakemake@output[["r2_filt"]]
errplotF <- snakemake@output[["errplotF"]]
errplotR <- snakemake@output[["errplotR"]]
denoising_R1 <- snakemake@output[["denoiseR1"]]
denoising_R2 <- snakemake@output[["denoiseR2"]]
merged <- snakemake@output[["merged"]]
asv_table <- snakemake@output[["asv"]]
report <- snakemake@output[["report"]]
chimeras_fasta <- snakemake@output[["chimeras"]]

# Multi-threading
threads <- snakemake@threads[[1]]

# Parameters
sample.names <- snakemake@params[["sample"]]
max_EE <- snakemake@params[["max_EE"]]
minlength <- snakemake@params[["min_length"]]
maxlength <- snakemake@params[["max_length"]]
chimera <- snakemake@params[["chimera"]]

# logging
sink(snakemake@log[[1]], append=FALSE, split=FALSE)


# DADA Workflow -------------------------------------------------------------------------------------------------------------

# Filter reads
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
              maxN=0, maxEE=max_EE, rm.phix=TRUE,
              compress=TRUE, multithread=threads, verbose=TRUE)

# Learn Error rate
errF <- learnErrors(filtFs, multithread=threads, verbose=TRUE)
errR <- learnErrors(filtRs, multithread=threads, verbose=TRUE)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=threads, verbose=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=threads, verbose=TRUE)
 
# Merge Reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# ASV table
seqtab <- makeSequenceTable(mergers)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% minlength:maxlength]

# Remove chimeras
if (chimera == TRUE) {
	seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="per-sample", multithread=threads, verbose=TRUE)
	bimeras <- isBimeraDenovo(seqtab2, multithread=threads, verbose = FALSE)
} else {
	seqtab2.nochim <- seqtab2
}

# Export sequences and QC reports -------------------------------------------------------------------------------------------

# ASV table
asv <- data.frame(seqtab2.nochim)
colnames(asv) <- c("count")
asv <- cbind(asv, name = sprintf("ASV_%s", seq(1:dim(asv)[1])))
asv <- cbind(asv, sequence = rownames(asv))
rownames(asv) <- seq(1:dim(asv)[1])
fn <- function(x) paste0(">", x[2], ";size=", trimws(x[1]), "\n", x[3])
asfasta <- apply(asv, MARGIN=1, fn)
writeLines(asfasta, asv_table)

# Error plots
errFplot <- plotErrors(errF, nominalQ=TRUE)
ggsave(errplotF, width = 10, height = 10)
errRplot <- plotErrors(errF, nominalQ=TRUE)
ggsave(errplotR, width = 10, height = 10)

# Denoising reports
write.table(dadaFs$clustering, denoising_R1, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(dadaRs$clustering, denoising_R2, quote = FALSE, sep = "\t", row.names = FALSE)

# Merged reads
write.table(mergers, merged, quote = FALSE, sep = "\t", row.names = FALSE)

# Chimeras
if (chimera == TRUE) {
	bimeras <- bimeras[bimeras==TRUE]
	asfasta = vector()
	for (i in seq(1:length(bimeras))) {
		asfasta <- c(asfasta, paste0(">bimera_", i, "\n", names(bimeras)[i]))
	}
} else {
	asfasta <- c("Skipped chimera detection")
}
writeLines(asfasta, chimeras_fasta)

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