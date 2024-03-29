import pandas as pd

shell.executable("bash")

# Reads merging and quality filtering rules-------------------------------------


rule merge_reads:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq.gz",
        r2="{sample}/trimmed/{sample}_R2.fastq.gz",
    output:
        merged=temp("{sample}/pseudo_reads/{sample}_merged.fastq"),
        notmerged_fwd=temp("{sample}/pseudo_reads/{sample}_notmerged_fwd.fasta"),
        notmerged_rev=temp("{sample}/pseudo_reads/{sample}_notmerged_rev.fasta"),
    threads: config["threads_sample"]
    message:
        "[{wildcards.sample}] merging reads"
    conda:
        "../envs/vsearch.yaml"
    log:
        "logs/{sample}/read_merge.log",
    shell:
        """
        exec 2> {log}

        vsearch --fastq_mergepairs {input.r1} --reverse {input.r2} \
            --threads {threads} --fastqout {output.merged} \
            --fastq_eeout --fastaout_notmerged_fwd {output.notmerged_fwd} \
            --fastaout_notmerged_rev {output.notmerged_rev} \
            --fastq_allowmergestagger > {log} 2>&1
        """


rule qual_stat:
    input:
        merged="{sample}/pseudo_reads/{sample}_merged.fastq",
    output:
        stat="{sample}/pseudo_reads/{sample}_merged_qual.txt",
    message:
        "[{wildcards.sample}] collecting quality statistics"
    conda:
        "../envs/vsearch.yaml"
    log:
        "logs/{sample}/fastq_quality_check.log",
    shell:
        """
        exec 2> {log}

        vsearch --fastq_eestats {input.merged} --output {output.stat}
        """


rule quality_filter:
    input:
        merged="{sample}/pseudo_reads/{sample}_merged.fastq",
    output:
        filtered=temp("{sample}/pseudo_reads/{sample}_filtered.fasta"),
        discarded=temp("{sample}/pseudo_reads/{sample}_discarded.fasta"),
    params:
        minlen=config["amplicon_min_length"],
        maxlen=config["amplicon_max_length"],
        maxee=config["max_expected_errors"],
        maxns=config["max_ns"],
    message:
        "[{wildcards.sample}] quality filtering"
    conda:
        "../envs/vsearch.yaml"
    log:
        "logs/{sample}/fastq_quality_filter.log",
    shell:
        """
        exec 2> {log}

        vsearch --fastq_filter {input.merged} \
        --fastq_maxee {params.maxee} --fastq_minlen {params.minlen} \
        --fastq_maxlen {params.maxlen} --fastq_maxns {params.maxns} \
        --fastaout {output.filtered} \
        --fasta_width 0 \
        --fastaout_discarded {output.discarded} --log {log}
        """


rule dereplicate:
    input:
        filtered="{sample}/pseudo_reads/{sample}_filtered.fasta",
    output:
        derep="{sample}/pseudo_reads/{sample}_derep.fasta",
    message:
        "[{wildcards.sample}] dereplicating"
    conda:
        "../envs/vsearch.yaml"
    log:
        "logs/{sample}/dereplication.log",
    shell:
        """
        exec 2> {log}

        vsearch --derep_fulllength {input.filtered} --strand plus \
        --output {output.derep} --sizeout \
        --relabel {wildcards.sample}_seq --fasta_width 0 --log {log}
        """


rule qc_stats:
    input:
        merged="{sample}/pseudo_reads/{sample}_merged.fastq",
        filtered="{sample}/pseudo_reads/{sample}_filtered.fasta",
        notmerged_fwd="{sample}/pseudo_reads/{sample}_notmerged_fwd.fasta",
        notmerged_rev="{sample}/pseudo_reads/{sample}_notmerged_rev.fasta",
        discarded="{sample}/pseudo_reads/{sample}_discarded.fasta",
        dereplicated="{sample}/pseudo_reads/{sample}_derep.fasta",
    output:
        merging="{sample}/reports/{sample}_merging.tsv",
    message:
        "[{wildcards.sample}] collecting quality filtering summary "
    log:
        "logs/{sample}/qc_stats.log",
    shell:
        """
        exec 2> {log}

        # Parsing fasta/fastq files
        merged=$(grep -c "^@" {input.merged} || true)
        notmerged=$(grep -c "^>" {input.notmerged_fwd} || true)
        total_reads=$(($merged + $notmerged)) 
        filtered=$(grep -c "^>" {input.filtered} || true)
        discarded=$(grep -c "^>" {input.discarded} || true)
        dereplicated=$(grep -c "^>" {input.dereplicated} || true)

        # Calculating fractions
        if [[ $total_reads -eq 0 ]]
        then
            notmerged_perc=0
        else
            notmerged_perc=$(printf %.2f "$((10**3 * (100* $notmerged / $total_reads)))e-3")
        fi

        if [[ $merged -eq 0 ]]
        then
            discarded_perc=0
        else
            discarded_perc=$(printf %.2f "$((10**3 * (100* $discarded / $merged)))e-3")
        fi

        if [[ $filtered -eq 0 ]]
        then
            derep_perc=0
        else
            derep_perc=$(printf %.2f "$((10**3 * (100* $dereplicated / $filtered)))e-3")
        fi

        # Writing report
        echo "Sample\tTotal reads\tPseudo-reads\tMerging failures [%]\tPseudo-reads PF\tDiscarded reads [%]\tUnique sequences\tUnique sequences [%]" > {output.merging}
        echo "{wildcards.sample}\t$total_reads\t$merged\t$notmerged_perc\t$filtered\t$discarded_perc\t$dereplicated\t$derep_perc" >> {output.merging}
        """


rule collect_qc_stats:
    input:
        report=expand("{sample}/reports/{sample}_merging.tsv", sample=samples.index),
    output:
        agg=report(
            "reports/merging_stats.tsv",
            caption="../report/vsearch_qc_stats.rst",
            category="Quality controls",
        ),
    message:
        "[All] collecting quality filtering stats"
    log:
        "logs/all/qc_stats.log",
    shell:
        """
        exec 2> {log}

        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """


# Clustering rules----------------------------


rule cluster:
    input:
        fasta="{sample}/pseudo_reads/{sample}_derep.fasta",
    output:
        centroids=temp("{sample}/clustering/{sample}_clusters.fasta"),
        stat="{sample}/clustering/{sample}_clustertab.txt",
    params:
        clusterID=config["cluster_identity"],
    conda:
        "../envs/vsearch.yaml"
    threads: config["threads"]
    message:
        "[{wildcards.sample}] distance greedy clustering sequences"
    log:
        "logs/{sample}/clustering.log",
    shell:
        """
        exec 2> {log}

        vsearch --cluster_size {input.fasta} --threads {threads} \
        --id {params.clusterID} --strand plus --sizein --sizeout \
        --fasta_width 0 --centroids {output.centroids} --uc {output.stat} \
        --log {log}
        """


rule sort_otu:
    input:
        fasta="{sample}/clustering/{sample}_clusters.fasta",
    output:
        sorted="{sample}/clustering/{sample}_clusters_sorted.fasta",
    params:
        min_size=config["cluster_minsize"],
    conda:
        "../envs/vsearch.yaml"
    threads: config["threads"]
    message:
        "[{wildcards.sample}] sorting centroids and size-filtering"
    log:
        "logs/{sample}/sort_otu.log",
    shell:
        """
        exec 2> {log}

        vsearch --sortbysize {input.fasta} --threads {threads} \
        --sizein --sizeout \
        --fasta_width 0 --minsize {params.min_size} --output {output.sorted} \
        --log {log}
        """


# Chimera detection-------------


rule chimera_denovo:
    # should be skipped if config["remove_chimera"] is False
    input:
        sorted="{sample}/clustering/{sample}_clusters_sorted.fasta",
    output:
        nonchim="{sample}/clustering/{sample}_clusters_nonchimera.fasta",
        chimera="{sample}/clustering/{sample}_clusters_chimeras.txt",
    conda:
        "../envs/vsearch.yaml"
    message:
        "[{wildcards.sample}] de novo chimera detection"
    log:
        "logs/{sample}/denovo_chimera.log",
    shell:
        """
        exec 2> {log}

        vsearch --uchime_denovo {input.sorted} \
        --sizein --sizeout --fasta_width 0 \
        --qmask none --nonchimeras {output.nonchim} \
        --uchimeout {output.chimera} | tee {log} 2>&1
        """


rule relabel_otu:
    input:
        fasta="{sample}/clustering/{sample}_clusters_nonchimera.fasta"
        if config["remove_chimera"]
        else "{sample}/clustering/{sample}_clusters_sorted.fasta",
    output:
        renamed="{sample}/clustering/{sample}_OTUs.fasta",
    conda:
        "../envs/vsearch.yaml"
    message:
        "[{wildcards.sample}] relabelling OTUs"
    log:
        "logs/{sample}/relabel_otus.log",
    shell:
        """
        exec 2> {log}

        vsearch --fastx_filter {input.fasta} \
        --sizein --sizeout --fasta_width 0 \
        --relabel OTU_ --fastaout {output.renamed}
        """


rule create_otu_tab:
    input:
        fasta="{sample}/clustering/{sample}_OTUs.fasta",
    output:
        tab="{sample}/clustering/{sample}_OTUs.txt",
    message:
        "[{wildcards.sample}] export OTU table"
    log:
        "logs/{sample}/out_tab.log",
    shell:
        """
        exec 2> {log}

        grep "^>" {input.fasta} \
            | sed -e 's/;size=/\t/' \
            | tr -d '>' > {output.tab}
        """


rule clustering_stats:
    input:
        derep="{sample}/pseudo_reads/{sample}_derep.fasta",
        clusters="{sample}/clustering/{sample}_clusters.fasta",
        sizefilt="{sample}/clustering/{sample}_clusters_sorted.fasta",
        non_chimera="{sample}/clustering/{sample}_OTUs.fasta",
    output:
        report="{sample}/reports/{sample}_clustering.tsv",
    message:
        "[{wildcards.sample}] collecting clustering stats"
    log:
        "logs/{sample}/cluastering_stats.log",
    shell:
        """
        exec 2> {log}

        # Collecting counts
        uniques_seq=$(grep -c "^>" {input.derep} || true)
        uniques_reads=$(grep "^>" {input.derep} | awk -F '=' '{{s+=$2}}END{{print s}}' || true)

        clusters_seq=$(grep -c "^>" {input.clusters} || true) 
        # number of reads in clusters is same as number of reads in derep

        size_filt_seq=$(grep -c "^>" {input.sizefilt} || true)
        size_filt_reads=$(grep "^>" {input.sizefilt} | awk -F '=' '{{s+=$2}}END{{print s}}' || true)

        non_chimera_seq=$(grep -c "^>" {input.non_chimera} || true)
        non_chimera_reads=$(grep "^>" {input.non_chimera} | awk -F '=' '{{s+=$2}}END{{print s}}' || true)

        # Get percents
        discarded_seq=$(($clusters_seq - $size_filt_seq))
        if [[ $clusters_seq -eq 0 ]]
        then
            discarded_perc_clust=0
        else
            discarded_perc_clust=$(printf %.2f "$((10**3 * (100* $discarded_seq / $clusters_seq)))e-3")
        fi

        discarded_reads=$(($uniques_reads - $size_filt_reads))
        if [[ $uniques_reads -eq 0 ]]
        then
            discarded_perc_reads=0
            clustered_perc=0
        else
            discarded_perc_reads=$(printf %.2f "$((10**3 * (100* $discarded_reads / $uniques_reads)))e-3")
            clustered_perc=$(printf %.2f "$((10**3 * (100* $non_chimera_reads / $uniques_reads)))e-3")
        fi

        chim_seq=$(($size_filt_seq - $non_chimera_seq))
        if [[ $discarded_perc_reads -eq 0 ]]
        then
            chim_seq_perc=0
        else
            chim_seq_perc=$(printf %.2f "$((10**3 * (100* $chim_seq / $size_filt_seq)))e-3")
        fi

        chim_reads=$(($size_filt_reads - $non_chimera_reads))
        if [[ $size_filt_reads -eq 0 ]]
        then
            chim_reads_perc=0
        else
            chim_reads_perc=$(printf %.2f "$((10**3 * (100* $chim_reads / $size_filt_reads)))e-3")
        fi

        # Writting report
        echo "Sample\tUnique sequences\tClusters\tClusters above size filter\tDiscarded clusters[% of clusters]\tDiscarded clusters[% of reads]\tNon-chimeric clusters (OTU)\tChimeras [% of clusters]\tChimeras [% of reads]\tPseudo-reads clustered\tPseudo-reads clustered [%]" > {output.report}
        echo "{wildcards.sample}\t$uniques_seq\t$clusters_seq\t$size_filt_seq\t$discarded_perc_clust\t$discarded_perc_reads\t$non_chimera_seq\t$chim_seq_perc\t$chim_reads_perc\t$non_chimera_reads\t$clustered_perc" >> {output.report}
        """


rule collect_clustering_stats:
    input:
        report=expand("{sample}/reports/{sample}_clustering.tsv", sample=samples.index),
    output:
        agg=report(
            "reports/clustering_stats.tsv",
            caption="../report/clustering_stats.rst",
            category="Quality controls",
        ),
    message:
        "[All] collecting clustering stats"
    log:
        "logs/all/clustering_stats.log",
    shell:
        """
        exec 2> {log}

        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """
