import pandas as pd

shell.executable("bash")

# Krona rules ---------------------------------------------------------------------------------------------------------------

rule krona_table:
    input:
        "{sample}/reports/{sample}_composition.tsv",
    output:
        "{sample}/krona/{sample}_krona_table.txt",
    params:
        lineage = config["taxonomy"]["rankedlineage_dmp"],
        nodes = config["taxonomy"]["nodes_dmp"],
    message:
        "Exporting {wildcards.sample} in Krona input format"
    script:
        "../scripts/krona_table.py"

rule krona:
    input:
        "{sample}/krona/{sample}_krona_table.txt",
    output:
        "{sample}/reports/{sample}_krona_chart.html",
    message:
        "producing graphical summary for {wildcards.sample}"
    conda:
        "../envs/krona.yaml"
    shell:
        "ktImportText -o {output} {input}"

rule krona_all:
    input:
        expand("{sample}/krona/{sample}_krona_table.txt", sample = samples.index),
    output:
        "reports/krona_chart.html",
    params:
        samples.index,
    message:
        "Producing graphical summary result"
    conda:
        "../envs/krona.yaml"
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
        tax = "{sample}/reports/{sample}_taxonomy_assignement_stats.tsv",
    output:
        "{sample}/reports/{sample}_summary.tsv",
    message:
        "Summarizing statistics for {wildcards.sample}"
    params:
        method = config["cluster"]["method"],
    shell:
        """
        if [[ {params.method} == "otu" ]] 
        then
            echo "Sample\tQ30 rate\tInsert size peak\tRead number\tPseudo-reads\tReads in OTU\tOTU number\t(Sub-)Species consensus\tGenus consensus\tHigher rank consensus\tNo match" > {output}
            
            Q30=$(tail -n +2 {input.fastp} | cut -d'\t' -f9)
            size=$(tail -n +2 {input.fastp} | cut -d'\t' -f11)
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
            echo "Sample\tQ30 rate\tInsert size peak\tRead number\tPseudo-reads\tReads in ASV\tASV number\t(Sub-)Species consensus\tGenus consensus\tHigher rank consensus\tNo match" > {output}
        
            Q30=$(tail -n +2 {input.fastp} | cut -d'\t' -f9)
            size=$(tail -n +2 {input.fastp} | cut -d'\t' -f11)
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
        expand("{sample}/reports/{sample}_summary.tsv", sample= samples.index),
    output:
        "reports/summary.tsv",
    message:
        "Collecting summary reports"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """
        
# Rmarkdown reports rules ---------------------------------------------------------------------------------------------------

rule report_sample: 
    input:
        summary = "{sample}/reports/{sample}_summary.tsv",
        fastp = "{sample}/reports/{sample}_trimmed.tsv",
        qc_filtering = "{sample}/reports/{sample}_merging.tsv" if config["cluster"]["method"] == "otu" else "{sample}/reports/{sample}_denoising.tsv",
        clustering = "{sample}/reports/{sample}_clustering.tsv" if config["cluster"]["method"] == "otu" else "{sample}/reports/{sample}_denoising.tsv",
        blast_rep = "{sample}/reports/{sample}_blast_stats.tsv",
        taxonomy = "{sample}/reports/{sample}_taxonomy_assignement_stats.tsv",
        result = "{sample}/reports/{sample}_composition.tsv",
        db = "reports/db_versions.tsv",
        soft = "reports/software_versions.tsv",
        krona = "{sample}/reports/{sample}_krona_chart.html",
    params:
        method = config["cluster"]["method"],
        workdir = config["workdir"],
        version = git_version(),
        sample = "{sample}",
    output:
        "{sample}/reports/{sample}_report.html",
    conda:
        "../envs/rmarkdown.yaml"
    message:
        "Generating html report for {wildcards.sample}"
    script:
        "../scripts/write_report.Rmd"

rule report_all: 
    input:
        summary = "reports/summary.tsv",
        fastp = "reports/fastp_stats.tsv",
        qc_filtering = "reports/merging_stats.tsv" if config["cluster"]["method"] == "otu" else "reports/denoising.tsv",
        clustering = "reports/clustering_stats.tsv" if config["cluster"]["method"] == "otu" else "reports/denoising.tsv",
        blast_rep = "reports/blast_stats.tsv",
        taxonomy = "reports/taxonomy_assignement_stats.tsv",
        result = "reports/composition_summary.tsv",
        db = "reports/db_versions.tsv",
        soft = "reports/software_versions.tsv",
        krona = "reports/krona_chart.html",
    params:
        method = config["cluster"]["method"],
        workdir = config["workdir"],
        version = git_version(),
        sample = "all",
    output:
        "reports/report.html",
    conda:
        "../envs/rmarkdown.yaml"
    message:
        "Generating global html report"
    script:
        "../scripts/write_report.Rmd"

# Version reports -----------------------------------------------------------------------------------------------------------

rule software_versions:
    output:
        "reports/software_versions.tsv",
    message:
        "Collecting software versions"
    params:
        method = config["cluster"]["method"],
    shell:
        """
        echo "Software\tVersion" > {output}
        
        paste <(echo "fastp") <(grep fastp= {workflow.basedir}/envs/fastp.yaml | cut -d "=" -f2) >> {output}
        
        paste <(echo "cutadapt") <(grep cutadapt= {workflow.basedir}/envs/cutadapt.yaml | cut -d "=" -f2) >> {output}
        
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
        "reports/db_versions.tsv",
    message:
        "Collecting databases versions"
    params:
        blast = config["blast"]["blast_DB"],
        taxdb = config["blast"]["taxdb"],
        taxdump_nodes = config["taxonomy"]["nodes_dmp"],
        taxdump_lin = config["taxonomy"]["rankedlineage_dmp"],
    shell:
        """
        echo "Database\tLast modified\tFull path" > {output}      
        
        paste <(echo "BLAST") <(date +%F -r {params.blast}.nto) <(echo {params.blast}) >> {output}
        
        paste <(echo "taxdb.bti") <(date +%F -r {params.taxdb}/taxdb.bti) <(echo {params.taxdb}/taxdb.bti) >> {output}
        paste <(echo "taxdb.btd") <(date +%F -r {params.taxdb}/taxdb.btd) <(echo {params.taxdb}/taxdb.btd) >> {output}
        
        paste <(echo "taxdump lineages") <(date +%F -r {params.taxdump_lin}) <(echo {params.taxdump_lin}) >> {output}
        paste <(echo "taxdump nodes") <(date +%F -r {params.taxdump_nodes}) <(echo {params.taxdump_nodes}) >> {output}
        """
