import pandas as pd
import os, json, csv, subprocess

shell.executable("bash")

# Sample getter ----------------------------------------------------------------

def _get_fastq(wildcards,read_pair='fq1'):
    return samples.loc[(wildcards.sample), [read_pair]].dropna()[0]

# Rules primers trimming

rule get_primer_revcomp:
    output:
        primers_rc = temp("common/primer_revcomp.fa"),
    params:
        primers = config["trimming"]["primers_fasta"],
    conda:
        "../envs/seqtk.yaml"
    shell:
        """
        seqtk seq -r {params.primers} > {output.primers_rc}
        """

rule cutadapt:
    input:
        r1 = lambda wildcards: _get_fastq(wildcards, 'fq1'),
        r2 = lambda wildcards: _get_fastq(wildcards, 'fq2'),
        primers_rc = "common/primer_revcomp.fa",
    output:
        r1 = temp("{sample}/trimmed/{sample}_primertrimmed_R1.fastq.gz"),
        r2 = temp("{sample}/trimmed/{sample}_primertrimmed_R2.fastq.gz"),
        trash_R1_5p = temp("{sample}/trimmed/{sample}_noprimer5p_R1.fastq.gz"),
        trash_R2_5p =temp("{sample}/trimmed/{sample}_noprimer5p_R2.fastq.gz"),
        trash_R1_3p = temp("{sample}/trimmed/{sample}_noprimer3p_R1.fastq.gz"),
        trash_R2_3p =temp("{sample}/trimmed/{sample}_noprimer3p_R2.fastq.gz"),
    params:
        error_rate = config["trimming"]["primer_error_rate"],
        primer_3p = config["trimming"]["primers_3end"],
        primers = config["trimming"]["primers_fasta"],
    message:
        "Trimming primers on {wildcards.sample}"
    conda:
        "../envs/cutadapt.yaml"
    log:
        "logs/{sample}_cutadapt.log"
    shell:
        """
        # Simple case only 5p trimming
        if [[ {params.primer_3p} == false ]]
        then
            cutadapt {input.r1} \
                {input.r2} \
                -o {output.r1} \
                -p {output.r2} \
                -g file:{params.primers} \
                -G file:{params.primers} \
                --untrimmed-output {output.trash_R1_5p} \
                --untrimmed-paired-output {output.trash_R2_5p} \
                --error-rate {params.error_rate} \
                2>&1 > {log}
            touch {output.trash_R1_3p}
            touch {output.trash_R2_3p}
            
        # in case trimming of 3p is also nescessary
        else
            cutadapt --interleaved \
                {input.r1} \
                {input.r2} \
                -g file:{params.primers} \
                -G file:{params.primers} \
                --untrimmed-output {output.trash_R1_5p} \
                --untrimmed-paired-output {output.trash_R2_5p} \
                --error-rate {params.error_rate} \
            2>> {log} \
            | cutadapt --interleaved \
                -o {output.r1} \
                -p {output.r2} \
                -a file:{input.primers_rc} \
                -A file:{input.primers_rc} \
                --untrimmed-output {output.trash_R1_3p} \
                --untrimmed-paired-output {output.trash_R2_3p} \
                --error-rate {params.error_rate} \
                - \
            2>&1 >> {log}
        fi
        """
        
rule primer_trimming_stats:
    input:
        before_r1 = lambda wildcards: _get_fastq(wildcards, 'fq1'),
        after_r1 = "{sample}/trimmed/{sample}_primertrimmed_R1.fastq.gz",
        before_r2 = lambda wildcards: _get_fastq(wildcards, 'fq2'),
        after_r2 = "{sample}/trimmed/{sample}_primertrimmed_R2.fastq.gz",
    output:
        report = temp("{sample}/trimmed/{sample}_primer_trimming.tsv"),
    message:
        "Collecting primer trimming statisctics for {wildcards.sample}"
    shell:
        """
        before_r1=$(zcat {input.before_r1} | echo $((`wc -l`/4)))
        after_r1=$(zcat {input.after_r1} | echo $((`wc -l`/4)))
        before_r2=$(zcat {input.before_r2} | echo $((`wc -l`/4)))
        after_r2=$(zcat {input.after_r2} | echo $((`wc -l`/4)))
        
        before=$(( before_r1 + before_r2 ))
        after=$(( after_r1 + after_r2 ))
        
        if [ $after -ne 0 ] 
        then
            perc_discarded=$( python -c "print(f'{{round(100*(1-${{after}}/${{before}}),2)}}')" )
        else
            perc_discarded=0.00
        fi
        
        echo "Sample\tTotal raw reads\tTotal reads after primer trim\tNo primer found [%]" > {output.report}
        echo "{wildcards.sample}\t$before\t$after\t$perc_discarded" >> {output.report}
        """

# Rules quality trimming -------------------------------------------------------

rule run_fastp:
    input:
        r1 = "{sample}/trimmed/{sample}_primertrimmed_R1.fastq.gz",
        r2 = "{sample}/trimmed/{sample}_primertrimmed_R2.fastq.gz",
    output:
        r1 = temp("{sample}/trimmed/{sample}_R1.fastq.gz"),
        r2 = temp("{sample}/trimmed/{sample}_R2.fastq.gz"),
        json = "{sample}/trimmed/{sample}.json",
        html = "{sample}/trimmed/{sample}.html",
    params:
        length_required = config["trimming"]["length_required"],
        qualified_quality_phred = config["trimming"]["qualified_quality_phred"],
        window_size = config["trimming"]["window_size"],
        mean_qual = config["trimming"]["mean_quality"],
    threads:
        config["threads_sample"]
    message:
        "Running fastp on {wildcards.sample}"
    conda:
        "../envs/fastp.yaml"
    log: 
        "logs/{sample}_fastp.log"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -h {output.html} -j {output.json}\
            --length_required {params.length_required} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --cut_by_quality3 \
            --cut_window_size {params.window_size} \
            --cut_mean_quality {params.mean_qual} \
            --disable_adapter_trimming \
            --thread {threads} \
            --report_title 'Sample {wildcards.sample}' \
        > {log} 2>&1
        """

rule parse_fastp:
    input:
        json = "{sample}/trimmed/{sample}.json",
        html = "{sample}/trimmed/{sample}.html",
    output:
        tsv = temp("{sample}/trimmed/{sample}_fastp.tsv"),
    message:
        "Parsing fastp json report for {wildcards.sample}"
    run:
        with open(input.json,'r') as handle:
            data = json.load(handle)
          
        link_path = os.path.join("..", input.html)
        header = "Total bases before quality trim\tTotal reads after quality trim\tTotal bases after quality trim\tQ20 rate after\tQ30 rate after\tDuplication rate\tInsert size peak\tlink_to_report"
        datalist = [data["summary"]["before_filtering"]["total_bases"],
                    data["summary"]["after_filtering"]["total_reads"],
                    data["summary"]["after_filtering"]["total_bases"],
                    data["summary"]["after_filtering"]["q20_rate"],
                    data["summary"]["after_filtering"]["q30_rate"],
                    data["duplication"]["rate"],data["insert_size"]["peak"], 
                    link_path]
                    
        with open (output.tsv,"w") as outfile:
            outfile.write(header+"\n")
            writer=csv.writer(outfile, delimiter='\t')
            writer.writerow(datalist) 

# Reporting rules --------------------------------------------------------------

rule trimming_stats:
    input:
        cutadapt = "{sample}/trimmed/{sample}_primer_trimming.tsv",
        fastp = "{sample}/trimmed/{sample}_fastp.tsv",
    output:
        report = "{sample}/reports/{sample}_trimmed.tsv"
    message:
        "Merging trimming stats for {wildcards.sample}"
    shell:
        """
        paste {input.cutadapt} {input.fastp} > {output.report}
        """

rule collect_trimming_stats:
    input:
        expand("{sample}/reports/{sample}_trimmed.tsv", sample=samples.index),
    output:
        "reports/fastp_stats.tsv",
    message:
        "Collecting fastp stats"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """
