import pandas as pd
import os, json, csv, subprocess

shell.executable("bash")
    
# Settings ---------------------------
 
workdir: config["workdir"]

try:
    __version__ = subprocess.check_output(["git", "describe"], cwd= workflow.basedir).strip().decode("utf-8")
except subprocess.CalledProcessError:
    __version__ = "Unknown"

samples = pd.read_csv(config["samples"], index_col="sample", sep = "\t", engine="python")
samples.index = samples.index.astype('str', copy=False) # in case samples are integers, need to convert them to str

# Functions ------------------------------------------

def _get_fastq(wildcards,read_pair='fq1'):
    return samples.loc[(wildcards.sample), [read_pair]].dropna()[0]
    
# Rules ------------------------------
 
rule all:
    input: 
        # VSEARCH
        expand("{sample}/{sample}.derep.fasta", sample = samples.index),
        "VSEARCH/otus.fasta",
        expand("{sample}/{sample}_otutab.tsv", sample = samples.index),
        # Taxonomy
        "taxonomy/assignment_report.tsv",
        expand("{sample}/{sample}_composition.tsv", sample = samples.index),
        # Krona
        expand("{sample}/{sample}_krona_table.txt", sample = samples.index),
        "reports/Krona_chart.html",
        # Sample reports
        expand("trimmed/reports/{sample}.tsv", sample = samples.index),
        expand("{sample}/{sample}_qc_filtering_report.tsv", sample = samples.index),
        expand("{sample}/sequence_quality.stats", sample = samples.index),
        expand("{sample}/{sample}_mapping_report.tsv", sample = samples.index),
        expand("{sample}/{sample}_taxonomy_stats.tsv", sample = samples.index),
        expand("{sample}/{sample}_composition.tsv", sample = samples.index),
        expand("{sample}/{sample}_result_summary.tsv", sample = samples.index),
        expand("{sample}/{sample}_summary.tsv", sample = samples.index),
        # Global reports
        "reports/fastp_stats.tsv",
        "reports/qc_filtering_stats.tsv",
        "reports/clustering_stats.tsv",
        "reports/mapping_stats.tsv",
        "reports/assignment_stats.tsv",
        "reports/taxonomy_stats.tsv",
        "reports/summary.tsv",
        "reports/result_summary.tsv",
        "reports/software_versions.tsv",
        "reports/db_versions.tsv",
        "reports/centroid_size.tsv",
        # Markdown
        "reports/report.html",

# Fastp rules----------------------------
 
rule run_fastp:
    input:
        r1 = lambda wildcards: _get_fastq(wildcards, 'fq1'),
        r2 = lambda wildcards: _get_fastq(wildcards, 'fq2')
    output:
        r1 = temp("trimmed/{sample}_R1.fastq.gz"),
        r2 = temp("trimmed/{sample}_R2.fastq.gz"),
        json = "trimmed/reports/{sample}.json",
        html = "trimmed/reports/{sample}.html"
    params:
        length_required = config["fastp"]["length_required"],
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"],
        window_size = config["fastp"]["window_size"],
        mean_qual = config["fastp"]["mean_quality"]
    threads: config["threads_sample"]
    message: "Running fastp on {wildcards.sample}"
    conda: "envs/fastp.yaml"
    log: 
        "logs/{sample}_fastp.log"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json}\
        --length_required {params.length_required} --qualified_quality_phred {params.qualified_quality_phred} -3 -W {params.window_size} -M {params.mean_qual} --detect_adapter_for_pe --thread {threads} --report_title 'Sample {wildcards.sample}' |\
        tee {log} 2>&1"

rule parse_fastp:
    input:
        json = "trimmed/reports/{sample}.json",
        html = "trimmed/reports/{sample}.html"
    output:
        tsv = "trimmed/reports/{sample}.tsv"
    message: "Parsing fastp json report"
    run:
        with open(input.json,'r') as handle:
            data = json.load(handle)
          
        link_path = os.path.join("..", input.html)
        header = "Sample\tTotal reads before\tTotal bases before\tTotal reads after\tTotal bases after\tQ20 rate after\tQ30 rate after\tDuplication rate\tInsert size peak\tlink_to_report"
        datalist = [wildcards.sample, data["summary"]["before_filtering"]["total_reads"],data["summary"]["before_filtering"]["total_bases"],data["summary"]["after_filtering"]["total_reads"],data["summary"]["after_filtering"]["total_bases"],data["summary"]["after_filtering"]["q20_rate"],data["summary"]["after_filtering"]["q30_rate"],data["duplication"]["rate"],data["insert_size"]["peak"], link_path]
        with open (output.tsv,"w") as outfile:
            outfile.write(header+"\n")
            writer=csv.writer(outfile, delimiter='\t')
            writer.writerow(datalist) 

rule collect_fastp_stats:
    input:
        expand('trimmed/reports/{sample}.tsv', sample=samples.index)
    output:
        "reports/fastp_stats.tsv"
    message: "Collecting fastp stats"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """
 
# Reads merging and quality filtering rules----------------------------

rule merge_reads:
    input:
        r1 = "trimmed/{sample}_R1.fastq.gz",
        r2 = "trimmed/{sample}_R2.fastq.gz"
    output:
        merged = temp("{sample}/{sample}.merged.fastq"),
        notmerged_fwd = temp("{sample}/{sample}.notmerged.fwd.fasta"),
        notmerged_rev = temp("{sample}/{sample}.notmerged.rev.fasta")
    threads: config["threads_sample"]
    message: "Merging reads on {wildcards.sample}"
    conda: "envs/vsearch.yaml"
    log:
        "logs/{sample}_merge.log"
    shell:
        "vsearch --fastq_mergepairs {input.r1} --reverse {input.r2} --threads {threads} --fastqout {output.merged} \
        --fastq_eeout --fastaout_notmerged_fwd {output.notmerged_fwd} --fastaout_notmerged_rev {output.notmerged_rev} --log {log}"    

rule qual_stat:
    input: 
        merged = "{sample}/{sample}.merged.fastq"
    output:
        stat = "{sample}/sequence_quality.stats"
    message: "Collecting quality statistics for {wildcards.sample}"
    conda: "envs/vsearch.yaml"
    shell:
        "vsearch --fastq_eestats {input.merged} --output {output.stat}"
        
rule quality_filter: 
    input: 
        merged = "{sample}/{sample}.merged.fastq"
    output:
        filtered = temp("{sample}/{sample}_filtered.fasta"),
        discarded = temp("{sample}/{sample}_discarded.fasta")
    params:
        minlen= config["read_filter"]["min_length"],
        maxlen = config["read_filter"]["max_length"],
        maxee = config["read_filter"]["max_expected_errors"]
    message: "Quality filtering {wildcards.sample}"
    conda: "envs/vsearch.yaml"
    log:
        "logs/{sample}_filter.log"
    shell:
        "vsearch --fastq_filter {input.merged} --fastq_maxee {params.maxee} --fastq_minlen {params.minlen} --fastq_maxlen {params.maxlen}\
        --fastq_maxns 0 --fastaout {output.filtered} --fasta_width 0 --fastaout_discarded {output.discarded} --log {log}"


rule dereplicate:
    input: 
        filtered = "{sample}/{sample}_filtered.fasta"
    output:
        derep = "{sample}/{sample}.derep.fasta"
    message: "Dereplicating {wildcards.sample}"
    conda: "envs/vsearch.yaml"
    log:
        "logs/{sample}_derep.log"
    shell:
        "vsearch --derep_fulllength {input.filtered} --strand plus --output {output.derep} --sizeout --relabel {wildcards.sample}_seq --fasta_width 0 --log {log}"

rule qc_stats:
    input:
        merged = "{sample}/{sample}.merged.fastq",
        filtered = "{sample}/{sample}_filtered.fasta",
        notmerged_fwd = "{sample}/{sample}.notmerged.fwd.fasta",
        notmerged_rev = "{sample}/{sample}.notmerged.rev.fasta",
        discarded = "{sample}/{sample}_discarded.fasta",
        dereplicated = "{sample}/{sample}.derep.fasta"
    output:
        "{sample}/{sample}_qc_filtering_report.tsv"
    message: "Collecting quality filtering summary for {wildcards.sample}"
    shell:
        """
        # Parsing fasta/fastq files
        merged=$(grep -c "^@" {input.merged})
        notmerged=$(grep -c "^>" {input.notmerged_fwd})
        filtered=$(grep -c "^>" {input.filtered})
        discarded=$(grep -c "^>" {input.discarded})
        dereplicated=$(grep -c "^>" {input.dereplicated})
        # Calculating fractions
        reads_total=$(($merged + $notmerged))
        notmerged_perc=$(echo "scale=2;(100* $notmerged / $reads_total)" | bc)
        discarded_perc=$(echo "scale=2;(100* $discarded / $merged)" | bc)
        kept=$(echo "scale=2;(100* $filtered / $reads_total)" | bc)
        # Writing report
        echo "Sample\tTotal reads\tMerged reads\tMerging failures\tMerging failures [%]\tQuality filtered reads\tDiscarded reads\tDiscarded reads [%]\tNumber of unique sequences\tReads kept [%]" > {output}
        echo "{wildcards.sample}\t$reads_total\t$merged\t$notmerged\t$notmerged_perc\t$filtered\t$discarded\t$discarded_perc\t$dereplicated\t$kept" >> {output}
        """
        
rule collect_qc_stats:
    input:
        expand("{sample}/{sample}_qc_filtering_report.tsv", sample = samples.index)
    output:
        "reports/qc_filtering_stats.tsv"
    message: "Collecting quality filtering stats"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """
        
# Clustering rules----------------------------

rule merge_samples:
    input:
        expand("{sample}/{sample}.derep.fasta", sample = samples.index)
    output:
        temp("VSEARCH/all.fasta")
    message: "Merging samples"
    shell:
        """
        cat {input} > {output}
        """

rule derep_all:
    input:
        "VSEARCH/all.fasta"
    output:
        temp("VSEARCH/all.derep.fasta")
    conda: "envs/vsearch.yaml"
    message: "Dereplicating"
    log: 
        "logs/derep_all.log"
    shell:
        "vsearch --derep_fulllength {input} --sizein --sizeout --fasta_width 0 --output {output} --log {log}"
        
rule cluster:
    input: 
        "VSEARCH/all.derep.fasta"
    output:
        "VSEARCH/centroids.fasta"
    params:
        clusterID = config["cluster"]["cluster_identity"]
    conda: "envs/vsearch.yaml"
    threads: config["threads"]
    message: "Clustering sequences"
    log:
        "logs/clustering.log"
    shell:
        "vsearch --cluster_size {input} --threads {threads} --id {params.clusterID} --strand plus --sizein --sizeout --fasta_width 0 --centroids {output} --log {log}"
        
rule sort_all:
    input: 
        "VSEARCH/centroids.fasta"
    output:
        "VSEARCH/sorted.fasta"
    params:
        min_size= config["cluster"]["cluster_minsize"]
    conda: "envs/vsearch.yaml"
    threads: config["threads"]
    message: "Sorting centroids and removing singleton"
    log:
        "logs/sort_all.log"
    shell:
        "vsearch --sortbysize {input} --threads {threads} --sizein --sizeout --fasta_width 0 --minsize {params.min_size} --output {output} --log {log}"
  
rule centroid_histogram:
    input:
        "VSEARCH/centroids.fasta"
    output:
        "reports/centroid_size.tsv"
    message:
        "Collecting centroid size distribution"
    shell:
        """
        echo "sample\tseq_id\tsize" > {output} 
        grep "^>" {input} | sed -e "s/_seq/\t/" -e "s/;size=/\t/" | tr -d ">"  >> {output}   
        """  
        
# Chimera detection-------------

if config["chimera"]["method"] == "denovo_reference":
    include: "rules/chimera_denovo_and_ref.rule" 
elif config["chimera"]["method"] == "denovo":
    include: "rules/chimera_denovo_only.rule" 
elif config["chimera"]["method"] == "reference":
    include: "rules/chimera_ref_only.rule"
elif config["chimera"]["method"] == "None":
    include: "rules/no_chimera.rule"
    
# Reads mapping rules----------------------------

rule map_sample:
    input:
        fasta = "{sample}/{sample}.derep.fasta",
        db = "VSEARCH/otus.fasta"
    output:
        otu = "{sample}/{sample}_otutab.tsv",
    params:
        clusterID = config["cluster"]["cluster_identity"]
    threads: config["threads_sample"]
    conda: "envs/vsearch.yaml"
    message: "Mapping {wildcards.sample} to OTUs"
    log:
        "logs/{sample}_map_reads.log"
    shell:
        """
        vsearch --usearch_global {input.fasta} --threads {threads} --db {input.db} --id {params.clusterID}\
        --strand plus --sizein --sizeout --fasta_width 0 --qmask none --dbmask none --otutabout {output.otu} --log {log}    

        tail -n +2 {output.otu} | sort -k 2,2nr -o {output.otu} 
        """

rule mapping_stats:
    input: 
        fasta = "{sample}/{sample}.derep.fasta",
        otu = "{sample}/{sample}_otutab.tsv"
    output:
        "{sample}/{sample}_mapping_report.tsv"
    message: "Collecting mapping summary for {wildcards.sample}"
    shell:
        """
        # Collecting counts
        nreads=$(grep "^>" {input.fasta} | awk -F '=' '{{s+=$2}}END{{print s}}')
        nmapped=$(awk '{{s+=$2}}END{{print s}}' {input.otu})
        notu=$(grep -c "OTU_" {input.otu})
        max=$(head -n 1 {input.otu} | cut -f 2)
        min=$(tail -n 1 {input.otu} | cut -f 2)
        
        # Calculating fractions
        map_perc=$(echo "scale=2;(100* $nmapped / $nreads)" | bc)
        
        #Writting to file
        echo "Sample\tRead number\tReads mapped\tReads mapped [%]\tOTU number\tMax count\tMin count" > {output}
        echo "{wildcards.sample}\t$nreads\t$nmapped\t$map_perc\t$notu\t$max\t$min" >> {output}
        """

rule collect_mapping_stats:
    input:
        expand("{sample}/{sample}_mapping_report.tsv", sample = samples.index)
    output:
        "reports/mapping_stats.tsv"
    message: "Collecting mapping stats"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """
        
# Taxonomic assignement rules----------------------------

if config["workflow"]["taxonomy"] == "blast":
    include: "rules/blast.rule"
elif config["workflow"]["taxonomy"] == "sintax":
    include: "rules/sintax.rule"
        
rule tax_stats:
    input:
        "{sample}/{sample}_composition.tsv" 
    output:
        "{sample}/{sample}_taxonomy_stats.tsv"
    message: "Collecting taxonomy stats for {wildcards.sample}"
    shell:
        """
        echo "Sample\tNo Blast hit\tSpecy consensus\tGenus consensus\tFamily consensus\tHigher rank consensus" > {output}
        
        nohits=$(grep -c "-" {input} || true)
        spec=$(grep -c "species" {input} || true)
        gen=$(grep -c "genus" {input} || true)
        fam=$(grep -c "family" {input} || true)
        other=$(( $(grep -c "OTU_" {input} || true) - $nohits - $spec - $gen - $fam ))
        
        echo "{wildcards.sample}\t$nohits\t$spec\t$gen\t$fam\t$other" >> {output}
        """
        
rule collect_tax_stats:
    input:
        samples = expand("{sample}/{sample}_taxonomy_stats.tsv", sample = samples.index),
        all = "reports/assignment_stats.tsv"
    output:
        "reports/taxonomy_stats.tsv"
    message: "Collecting blast statistics"
    shell:
        """              
        # Summary
        nohits=$(grep -c "-" <(tail -n +2 {input.all}) || true)
        spec=$(grep -c "species" <(tail -n +2 {input.all}) || true)
        gen=$(grep -c "genus" <(tail -n +2 {input.all}) || true)
        fam=$(grep -c "family" <(tail -n +2 {input.all}) || true)
        other=$(( $(grep -c "OTU_" <(tail -n +2 {input.all}) || true) - $nohits - $spec - $gen - $fam ))
        
        echo "All\t$nohits\t$spec\t$gen\t$fam\t$other" >> {output}
        
        # Per sample
        for i in {input.samples}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        
        # Insert Header 
        sed -i "1 i\Sample\tNo Blast hit\tSpecy consensus\tGenus consensus\tFamily consensus\tHigher rank consensus" {output}
        """

rule summarize_results:
    input:
        compo = "{sample}/{sample}_composition.tsv"
    output:
        report = "{sample}/{sample}_result_summary.tsv"
    message:
        "Summarizing results for {wildcards.sample}"
    run:
        df = pd.read_csv(input.compo, sep="\t", header=0)
        groups = df.groupby(['Consensus', 'Rank', 'Taxid'])['Count'].sum().sort_values(ascending=False).to_frame().reset_index()
        groups['perc']= round(groups['Count']/groups['Count'].sum() *100, 2)
        groups.insert(0, 'Sample', wildcards.sample)
        groups.rename(columns={"perc":"Percent of total"}, index={"-": "No match"}, inplace = True)
        groups.to_csv(output.report, sep="\t", index = False)
        
rule collect_results:
    input:
        expand("{sample}/{sample}_result_summary.tsv", sample = samples.index)
    output:
        "reports/result_summary.tsv"
    message: "Collecting results"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """

# Krona rules----------------------------

rule krona_table:
    input:
        "{sample}/{sample}_result_summary.tsv"
    output:
        "{sample}/{sample}_krona_table.txt"
    params: 
        names = config["blast"]["names_dmp"],
        nodes = config["blast"]["nodes_dmp"]
    message:
        "Exporting {wildcards.sample} in Krona input format"
    script:
        "scripts/krona_table.py"
        
rule krona:
    input:
        expand("{sample}/{sample}_krona_table.txt", sample = samples.index)
    output:
        "reports/Krona_chart.html"
    params:
        samples.index
    message:
        "Producing graphical summary result"
    conda:
        "envs/krona.yaml"
    shell:
        """
        i=0
        for file in {input}
        do
            file_list[$i]="${{file}},$(echo ${{file}} | cut -d"/" -f1)"
            ((i+=1))
        done
        
        ktImportText -o {output} ${{file_list[@]}}
        """
    
# Report rules----------------------------

rule summary_report:
    input:
        fastp = "trimmed/reports/{sample}.tsv",
        filter = "{sample}/{sample}_qc_filtering_report.tsv",
        map = "{sample}/{sample}_mapping_report.tsv",
        tax = "{sample}/{sample}_taxonomy_stats.tsv"
    output:
        "{sample}/{sample}_summary.tsv"
    message: "Summarizing statistics for {wildcards.sample}"
    shell:
        """
        echo "Sample\tQ30 rate\tFiltered reads\tFiltered reads [%]\tMapped reads [%]\tOTU number\tSpecy consensus\tGenus consensus\tHigher rank\tNo consensus" > {output}
        
        Q30=$(tail -n +2 {input.fastp} | cut -d'\t' -f7)
        fil_reads=$(tail -n +2 {input.filter} | cut -d'\t' -f6)
        fil_perc=$(tail -n +2 {input.filter} | cut -d'\t' -f10) 
        mapped=$(tail -n +2 {input.map} | cut -d'\t' -f4)
        otu=$(tail -n +2 {input.map} | cut -d'\t' -f5)
        spec=$(tail -n +2 {input.tax} | cut -d'\t' -f3)
        gen=$(tail -n +2 {input.tax} | cut -d'\t' -f4)
        high=$(tail -n +2 {input.tax} | cut -d'\t' -f5)
        noc=$(tail -n +2 {input.tax} | cut -d'\t' -f2)
        
        echo "{wildcards.sample}\t$Q30\t$fil_reads\t$fil_perc\t$mapped\t$otu\t$spec\t$gen\t$high\t$noc" >> {output}
        """
        
rule collect_summaries:        
    input:
        expand("{sample}/{sample}_summary.tsv", sample= samples.index)
    output:
        "reports/summary.tsv"
    message: "Collecting summary reports"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """

rule report_all:
    input:
        summary = "reports/summary.tsv",
        fastp = "reports/fastp_stats.tsv",
        qc_filtering = "reports/qc_filtering_stats.tsv",
        clustering = "reports/clustering_stats.tsv",
        centroids = "reports/centroid_size.tsv",
        mapping = "reports/mapping_stats.tsv",
        blast = "reports/assignment_stats.tsv",
        blast_rep = "taxonomy/assignment_report.tsv",
        taxonomy = "reports/taxonomy_stats.tsv",
        result = "reports/result_summary.tsv",
        db = "reports/db_versions.tsv",
        soft = "reports/software_versions.tsv",
        merged_qc = expand("{sample}/sequence_quality.stats", sample = samples.index),
        krona = "reports/Krona_chart.html"
    params:
        workdir = config["workdir"],
        max_ee = config["read_filter"]["max_expected_errors"],
        max_len = config["read_filter"]["max_length"],
        min_len = config["read_filter"]["min_length"],
        min_cluster_size = config["cluster"]["cluster_minsize"],
        taxonomy = config["workflow"]["taxonomy"],
        version = __version__
    output:
        "reports/report.html"
    conda:
        "envs/rmarkdown.yaml"
    log:
        "logs/rmarkdown.log"
    message: "Generating html report"
    script:
        "scripts/write_report.Rmd"
          
rule software_versions:
    output:
        "reports/software_versions.tsv"
    message: "Collecting software versions"
    shell:
        """
        echo "Software\tVersion" > {output}
        paste <(echo "fastp") <(grep fastp= {workflow.basedir}/envs/fastp.yaml | cut -d "=" -f2) >> {output}
        paste <(echo "blast") <(grep blast= {workflow.basedir}/envs/blast.yaml | cut -d "=" -f2) >> {output}
        paste <(echo "vsearch") <(grep vsearch= {workflow.basedir}/envs/vsearch.yaml | cut -d "=" -f2) >> {output}
        """

rule database_version:
    output:
        "reports/db_versions.tsv"
    message: "Collecting databases versions"
    params:
        chimera = config["chimera"]["chimera_DB"],
        blast = config["blast"]["blast_DB"],
        taxdb = config["blast"]["taxdb"],
        taxdump = config["blast"]["nodes_dmp"]
    shell:
        """
        echo "Database\tLast modified\tFull path" > {output}      
        paste <(echo "Chimera") <(date +%F -r {params.chimera}) <(echo {params.chimera}) >> {output}
        paste <(echo "BLAST") <(date +%F -r {params.blast}) <(echo {params.blast}) >> {output}
        paste <(echo "taxdb") <(date +%F -r {params.taxdb}/taxdb.bti) <(echo {params.chimera}/taxdb[.bti/.btd]) >> {output}
        paste <(echo "taxdump") <(date +%F -r {params.taxdump}) <(echo $(dirname {params.taxdump})/[names.dmp/nodes.dmp]) >> {output}
        """
        
# Workflow----------------------------

onstart:
    print("\nYou are using FooDMe version: {}".format(__version__))

onsuccess:
    print("\nWorkflow finished, no error")
    
onerror:
    print("\nAn error occured, please consult the latest log file in {}".format(os.path.join(config["workdir"], ".snakemake", "log")))