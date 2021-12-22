#!/usr/bin/env Rscript

library(ggplot2, quiet=T)
library(dada2, quiet=T)

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
sample.names <- snakemake@params[["sample_name"]]
max_EE <- snakemake@params[["max_EE"]]
minlength <- snakemake@params[["min_length"]]
maxlength <- snakemake@params[["max_length"]]
max_mismatch <- snakemake@params[["max_mismatch"]]
chimera <- snakemake@params[["chimera"]]

# logging
log = file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# DADA Workflow -------------------------------------------------------------------------------------------------------------

# Is there a more elegant solution here to handle no results cases?
# can't find a way to properly handle exceptions in R1

tryCatch({
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
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, 
                            maxMismatch=max_mismatch, 
                            returnRejects=TRUE,
                            verbose=TRUE)

    filt_mergers <- mergers[mergers$accept==TRUE]

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
                    getN(filt_mergers),
                    round((1-getN(filt_mergers)/out[2])*100, 2),
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
    
}, error = function(c) {
    
    # Write empty files 
    if (!exists(filtRs)){
        writeLines("", filtRs)
    }
    if (!exists(filtFs)){
        writeLines("", filtFs)
    }
    if (!exists(errplotF)){
        ggplot() +
        theme_void() +
        geom_text(aes(0,0,label='N/A')) +
        xlab(NULL)
        ggsave(errplotF)
    }
    if (!exists(errplotR)){
        ggplot() +
        theme_void() +
        geom_text(aes(0,0,label='N/A')) +
        xlab(NULL)
        ggsave(errplotR)
    }
    if (!exists(denoising_R1)){
        writeLines("sequence	abundance	n0	n1	nunq	pval	birth_from	birth_pval	birth_fold	birth_ham	birth_qave",
                    denoising_R1)
    }
    if (!exists(denoising_R2)){
        writeLines("sequence	abundance	n0	n1	nunq	pval	birth_from	birth_pval	birth_fold	birth_ham	birth_qave",
                    denoising_R2)
    }
    if (!exists(merged)){
        writeLines("sequence	abundance	forward	reverse	nmatch	nmismatch	nindel	prefer	accept",
                    merged)
    }
    if (!exists(asv_table)){
        writeLines("", asv_table)
    }
    if (!exists(report)){
        writeLines("Sample	Total reads	Filtered reads	Discarded reads [%]	Denoised R1	Denoised R2	Merged	Merging failures [%]	ASV	Size-filtered ASV	Discarded ASVs [% ASV]	Discarded ASVs [% of reads]	Non-chimeric ASV	Chimeras [% ASV]	Chimeras [% of reads]	Reads in ASVs	Reads in ASVs [%]", 
                    report)
    }
    if (!exists(chimeras_fasta)){
        writeLines("", chimeras_fasta)
    }
})
