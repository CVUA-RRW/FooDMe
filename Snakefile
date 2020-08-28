import pandas as pd
import os, subprocess

shell.executable("bash")
    
# Settings ------------------------------------------------------------------------------------------------------------------
 
workdir: config["workdir"]

samples = pd.read_csv(config["samples"], index_col="sample", sep = "\t", engine="python")
samples.index = samples.index.astype('str', copy=False) # in case samples are integers, need to convert them to str

# Functions -----------------------------------------------------------------------------------------------------------------

def git_version():
    try:
        __version__ = subprocess.check_output(["git", "describe", "--always"], cwd= workflow.basedir).strip().decode("utf-8")
    except subprocess.CalledProcessError:
        __version__ = "Unknown"
    finally:
        return(__version__)
    
# Input rule ----------------------------------------------------------------------------------------------------------------
 
rule all:
    input: 
        # Summaries
        "reports/summary.tsv",
        "reports/db_versions.tsv",
        "reports/software_versions.tsv",
        
        # krona
        expand("{sample}/reports/{sample}_krona_chart.html", sample = samples.index),
        "reports/Krona_chart.html"
        
        # Markdown reports
        # expand("{sample}/reports/{sample}_report.html", sample = samples.index)
        # "reports/report.html"

# Includes ------------------------------------------------------------------------------------------------------------------
 
include: "rules/fastp.rule"
include: "rules/vsearch.rule" if config["cluster"]["method"] == "otu" else "rules/dada2.rule" 
include: "rules/blast.rule"

# Krona rules ---------------------------------------------------------------------------------------------------------------

rule krona_table:
    input:
        "{sample}/reports/{sample}_composition.tsv"
    output:
        "{sample}/krona/{sample}_krona_table.txt"
    params: 
        lineage = config["taxonomy"]["rankedlineage_dmp"],
        nodes = config["taxonomy"]["nodes_dmp"]
    message:
        "Exporting {wildcards.sample} in Krona input format"
    script:
        "scripts/krona_table.py"
    
rule krona:
    input:
        "{sample}/krona/{sample}_krona_table.txt"
    output:
        "{sample}/reports/{sample}_krona_chart.html"
    message: "producing graphical summary for {wildcards.sample}"
    conda:
        "envs/krona.yaml"
    shell:
        "ktImportText -o {output} {input}"
        
rule krona_all:
    input:
        expand("{sample}/krona/{sample}_krona_table.txt", sample = samples.index)
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
    
# Summary report rules ------------------------------------------------------------------------------------------------------

rule summary_report:
    input:
        fastp = "{sample}/reports/{sample}_trimmed.tsv",
        merging = "{sample}/reports/{sample}_merging.tsv" if config["cluster"]["method"] == "otu" else "{sample}/reports/{sample}_denoising.tsv",
        clustering = "{sample}/reports/{sample}_clustering.tsv" if config["cluster"]["method"] == "otu" else "{sample}/reports/{sample}_denoising.tsv",
        tax = "{sample}/reports/{sample}_taxonomy_assignement_stats.tsv"
    output:
        "{sample}/reports/{sample}_summary.tsv"
    message: "Summarizing statistics for {wildcards.sample}"
    params:
        method = config["cluster"]["method"]
    shell:
        """
        if [[ {params.method} == "otu" ]] 
        then
            echo "Sample\tQ30 rate\tInsert size peak\tRead number\tPseudo-reads\tReads in OTU\tOTU number\tSpecies consensus\tGenus consensus\tHigher rank consensus\tNo match" > {output}
            
            Q30=$(tail -n +2 {input.fastp} | cut -d'\t' -f7)
            size=$(tail -n +2 {input.fastp} | cut -d'\t' -f9)
            reads=$(tail -n +2 {input.merging} | cut -d'\t' -f2)
            pseudo=$(tail -n +2 {input.merging} | cut -d'\t' -f5)
            clustered=$(tail -n +2 {input.clustering} | cut -d'\t' -f10)
            otu=$(tail -n +2 {input.tax} | cut -d'\t' -f2)
            spec=$(tail -n +2 {input.tax} | cut -d'\t' -f5)
            gen=$(tail -n +2 {input.tax} | cut -d'\t' -f7)
            high=$(($(tail -n +2 {input.tax} | cut -d'\t' -f9) + $(tail -n +2 {input.tax} | cut -d'\t' -f11)))
            none=$(tail -n +2 {input.tax} | cut -d'\t' -f3)

            echo "{wildcards.sample}\t$Q30\t$size\t$reads\t$pseudo\t$clustered\t$otu\t$spec\t$gen\t$high\t$none" >> {output}
        else
            echo "Sample\tQ30 rate\tInsert size peak\tRead number\tPseudo-reads\tReads in ASV\tASV number\tSpecies consensus\tGenus consensus\tHigher rank consensus\tNo match" > {output}
        
            Q30=$(tail -n +2 {input.fastp} | cut -d'\t' -f7)
            size=$(tail -n +2 {input.fastp} | cut -d'\t' -f9)
            reads=$(tail -n +2 {input.clustering} | cut -d'\t' -f2)
            pseudo=$(tail -n +2 {input.clustering} | cut -d'\t' -f6)
            clustered=$(tail -n +2 {input.clustering} | cut -d'\t' -f16)
            otu=$(tail -n +2 {input.tax} | cut -d'\t' -f2)
            spec=$(tail -n +2 {input.tax} | cut -d'\t' -f5)
            gen=$(tail -n +2 {input.tax} | cut -d'\t' -f7)
            high=$(($(tail -n +2 {input.tax} | cut -d'\t' -f9) + $(tail -n +2 {input.tax} | cut -d'\t' -f11)))
            none=$(tail -n +2 {input.tax} | cut -d'\t' -f3)
            
            echo "{wildcards.sample}\t$Q30\t$size\t$reads\t$pseudo\t$clustered\t$otu\t$spec\t$gen\t$high\t$none" >> {output}
        fi
        """

rule collect_summaries:
    input:
        expand("{sample}/reports/{sample}_summary.tsv", sample= samples.index)
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
        
# Rmarkdown reports rules ---------------------------------------------------------------------------------------------------
######################################################
# FIXME: report all + generate single sample reports #
######################################################

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
        taxonomy = config["taxonomy"]["method"],
        version = git_version()
    output:
        "reports/report.html"
    conda:
        "envs/rmarkdown.yaml"
    log:
        "logs/rmarkdown.log"
    message: "Generating html report"
    script:
        "scripts/write_report.Rmd"

# Version reports -----------------------------------------------------------------------------------------------------------

rule software_versions:
    output:
        "reports/software_versions.tsv"
    message: "Collecting software versions"
    params:
        method = config["cluster"]["method"]
    shell:
        """
        echo "Software\tVersion" > {output}
        
        paste <(echo "fastp") <(grep fastp= {workflow.basedir}/envs/fastp.yaml | cut -d "=" -f2) >> {output}
        
        if [[ {params.method} == "otu" ]] 
        then
            paste <(echo "vsearch") <(grep vsearch= {workflow.basedir}/envs/vsearch.yaml | cut -d "=" -f2) >> {output}
        else
            paste <(echo "dada2") <(grep bioconductor-dada2= {workflow.basedir}/envs/dada2.yaml | cut -d "=" -f2) >> {output}
        fi
        
        paste <(echo "blast") <(grep blast= {workflow.basedir}/envs/blast.yaml | cut -d "=" -f2) >> {output}
        """

rule database_version:
    output:
        "reports/db_versions.tsv"
    message: "Collecting databases versions"
    params:
        blast = config["blast"]["blast_DB"],
        taxdb = config["blast"]["taxdb"],
        taxdump_nodes = config["taxonomy"]["nodes_dmp"],
        taxdump_lin = config["taxonomy"]["rankedlineage_dmp"]
    shell:
        """
        echo "Database\tLast modified\tFull path" > {output}      
        
        paste <(echo "BLAST") <(date +%F -r {params.blast}) <(echo {params.blast}) >> {output}
        
        paste <(echo "taxdb.bti") <(date +%F -r {params.taxdb}/taxdb.bti) <(echo {params.taxdb}/taxdb.bti) >> {output}
        paste <(echo "taxdb.btd") <(date +%F -r {params.taxdb}/taxdb.btd) <(echo {params.taxdb}/taxdb.btd) >> {output}
        
        paste <(echo "taxdump lineages") <(date +%F -r {params.taxdump_lin}) <(echo {params.taxdump_lin}) >> {output}
        paste <(echo "taxdump nodes") <(date +%F -r {params.taxdump_nodes}) <(echo {params.taxdump_nodes}) >> {output}
        """
        
# Workflow ------------------------------------------------------------------------------------------------------------------

onstart:
    print("\nYou are using FooDMe version: {}".format(git_version()))

onsuccess:
    print("\nWorkflow finished, no error")
    
onerror:
    print("\nAn error occured, please consult the latest log file in {}".format(os.path.join(config["workdir"], ".snakemake", "log")))