import pandas as pd

shell.executable("bash")

# Rules ------------------------------------------------------------------------

rule unpack_fastq:
    input: 
        r1 = expand("{sample}/trimmed/{sample}_R1.fastq.gz", sample = samples.index),
        r2 = expand("{sample}/trimmed/{sample}_R2.fastq.gz", sample = samples.index),
    output: 
        r1 = temp(expand("{sample}/trimmed/{sample}_R1.fastq", sample = samples.index)),
        r2 = temp(expand("{sample}/trimmed/{sample}_R2.fastq", sample = samples.index)),
    message:
        "Unpacking fastq files"
    shell:
        """
        gzip -kd {input.r1}
        gzip -kd {input.r2}
        """

rule denoise:
    input: 
        r1 = "{sample}/trimmed/{sample}_R1.fastq",
        r2 = "{sample}/trimmed/{sample}_R2.fastq",
    output:
        r1_filt = temp("{sample}/denoising/{sample}_R1_filtered.fasta"),
        r2_filt = temp("{sample}/denoising/{sample}_R2_filtered.fasta"),
        errplotF = "{sample}/denoising/{sample}_R1_errorRates.pdf",
        errplotR = "{sample}/denoising/{sample}_R2_errorRates.pdf",
        denoiseR1 = "{sample}/denoising/{sample}_R1_denoising.txt",
        denoiseR2 = "{sample}/denoising/{sample}_R2_denoising.txt",
        merged = "{sample}/denoising/{sample}_merging.txt",
        asv = "{sample}/denoising/{sample}_ASVs.fasta",
        report = "{sample}/reports/{sample}_denoising.tsv",
        chimeras = "{sample}/denoising/{sample}_chimeras.fasta",
    message:
        "Denoising and filtering {wildcards.sample}"
    threads:
        config["threads_sample"]
    conda:
        "../envs/dada2.yaml"
    params:
        sample = "{sample}",
        max_EE = config["read_filter"]["max_expected_errors"],
        min_length = config["read_filter"]["min_length"],
        max_length = config["read_filter"]["max_length"],
        chimera = config["chimera"],
    log:
        "logs/{sample}_denoising.log"
    script:
        "../scripts/dada.R"

rule collect_denoising_stats:
    input:
        expand("{sample}/reports/{sample}_denoising.tsv", sample = samples.index),
    output:
        "reports/denoising.tsv",
    message:
        "collecting denoising stats"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """
