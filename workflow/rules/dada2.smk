shell.executable("bash")


# Rules -----------------------------------------------------------------------


rule unpack_fastq:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq.gz",
        r2="{sample}/trimmed/{sample}_R2.fastq.gz",
    output:
        r1=temp("{sample}/trimmed/{sample}_R1.fastq"),
        r2=temp("{sample}/trimmed/{sample}_R2.fastq"),
    message:
        "[{wildcards.sample}][denoising] unpacking fastq files"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{sample}/fastq_unpack.log",
    shell:
        """
        exec 2> {log}
        gzip -kd {input.r1}
        gzip -kd {input.r2}
        """


rule denoise:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
    output:
        r1_filt=temp("{sample}/denoising/{sample}_R1_filtered.fasta.gz"),
        r2_filt=temp("{sample}/denoising/{sample}_R2_filtered.fasta.gz"),
        errplotF="{sample}/denoising/{sample}_R1_errorRates.pdf",
        errplotR="{sample}/denoising/{sample}_R2_errorRates.pdf",
        denoiseR1="{sample}/denoising/{sample}_R1_denoising.txt",
        denoiseR2="{sample}/denoising/{sample}_R2_denoising.txt",
        merged="{sample}/denoising/{sample}_merging.txt",
        asv="{sample}/denoising/{sample}_ASVs.fasta",
        report="{sample}/reports/{sample}_denoising.tsv",
        chimeras="{sample}/denoising/{sample}_chimeras.fasta",
    message:
        "[{wildcards.sample}][denoising] denoising and filtering"
    threads: config["threads_sample"]
    conda:
        "../envs/dada2.yaml"
    params:
        sample_name=lambda w, input: w.sample,
        max_EE=config["max_expected_errors"],
        min_length=config["amplicon_min_length"],
        max_length=config["amplicon_max_length"],
        chimera=config["remove_chimera"],
        max_mismatch=config["merging_max_mismatch"],
    log:
        "logs/{sample}/denoising.log",
    script:
        "../scripts/dada.R"


rule collect_denoising_stats:
    input:
        report=expand("{sample}/reports/{sample}_denoising.tsv", sample=samples.index),
    output:
        agg=report(
            "reports/denoising.tsv",
            caption="../report/denoising_stats.rst",
            category="Quality controls",
        ),
    message:
        "[All][denoising] aggregating denoising stats"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/denoising_stats.log",
    shell:
        """
        exec 2> {log}
        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """
