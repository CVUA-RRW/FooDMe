import pandas as pd

shell.executable("bash")

# Krona rules -----------------------------------------------------------------


rule krona_table:
    input:
        compo="{sample}/reports/{sample}_composition.tsv",
        tax="common/taxonomy.json",
    output:
        krt="{sample}/krona/{sample}_krona_table.txt",
    message:
        "Exporting {wildcards.sample} in Krona input format"
    conda:
        "../envs/taxidtools.yaml"
    log:
        "logs/{sample}/krona_table.log",
    script:
        "../scripts/krona_table.py"


rule krona:
    input:
        table="{sample}/krona/{sample}_krona_table.txt",
    output:
        graph=report(
            "{sample}/reports/{sample}_krona_chart.html",
            caption="../report/krona.rst",
            category="Results",
            subcategory="{wildcards.sample}",
        ),
    message:
        "Producing graphical summary for {wildcards.sample}"
    conda:
        "../envs/krona.yaml"
    log:
        "logs/{sample}/krona.log",
    shell:
        "ktImportText -o {output.graph} {input.table} 2> {log}"


rule krona_all:
    input:
        report=expand("{sample}/krona/{sample}_krona_table.txt", sample=samples.index),
    output:
        agg=report(
            "reports/krona_chart.html",
            caption="../report/krona_glob.rst",
            category="Results",
            subcategory="Global",
        ),
    params:
        samples.index,
    message:
        "Producing graphical summary result"
    conda:
        "../envs/krona.yaml"
    log:
        "logs/all/krona.log",
    shell:
        """
        exec 2> {log}
        i=0
        for file in {input.report}
        do
            file_list[$i]="${{file}},$(echo ${{file}} | cut -d'/' -f1)"
            ((i+=1))
        done

        ktImportText -o {output.agg} ${{file_list[@]}}
        """


# Summary report rules --------------------------------------------------------


rule summary_report:
    input:
        fastp="{sample}/reports/{sample}_trimmed.tsv",
        merging="{sample}/reports/{sample}_merging.tsv"
        if config["cluster"]["method"] == "otu"
        else "{sample}/reports/{sample}_denoising.tsv",
        clustering="{sample}/reports/{sample}_clustering.tsv"
        if config["cluster"]["method"] == "otu"
        else "{sample}/reports/{sample}_denoising.tsv",
        tax="{sample}/reports/{sample}_taxonomy_assignement_stats.tsv",
        compo="{sample}/reports/{sample}_composition.tsv",
    output:
        report="{sample}/reports/{sample}_summary.tsv",
    message:
        "Summarizing statistics for {wildcards.sample}"
    params:
        method=config["cluster"]["method"],
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{sample}/summary_report.log",
    shell:
        """
        exec 2> {log}

        if [[ {params.method} == "otu" ]] 
        then
            echo "Sample\tQ30 rate\tInsert size peak\tRead number\tPseudo-reads\tReads in OTU\tOTU number\tAssigned reads\t(Sub-)Species consensus\tGenus consensus\tHigher rank consensus\tNo match" > {output.report}

            Q30=$(tail -n +2 {input.fastp} | cut -d'\t' -f9)
            size=$(tail -n +2 {input.fastp} | cut -d'\t' -f11)
            reads=$(tail -n +2 {input.merging} | cut -d'\t' -f2)
            pseudo=$(tail -n +2 {input.merging} | cut -d'\t' -f5)
            clustered=$(tail -n +2 {input.clustering} | cut -d'\t' -f10)
            otu=$(tail -n +2 {input.tax} | cut -d'\t' -f2)
            assigned=$(tail -n +2 {input.compo} | awk '$2 != "No match"' | cut -d'\t' -f5 | awk '{{s+=$1}}END{{print s}}')
            spec=$(tail -n +2 {input.tax} | cut -d'\t' -f5)
            gen=$(tail -n +2 {input.tax} | cut -d'\t' -f7)
            high=$(($(tail -n +2 {input.tax} | cut -d'\t' -f9) + $(tail -n +2 {input.tax} | cut -d'\t' -f11)))
            none=$(tail -n +2 {input.tax} | cut -d'\t' -f3)

            echo "{wildcards.sample}\t$Q30\t$size\t$reads\t$pseudo\t$clustered\t$otu\t$assigned\t$spec\t$gen\t$high\t$none" >> {output.report}
        else
            echo "Sample\tQ30 rate\tInsert size peak\tRead number\tPseudo-reads\tReads in ASV\tASV number\tAssigned reads\t(Sub-)Species consensus\tGenus consensus\tHigher rank consensus\tNo match" > {output.report}

            Q30=$(tail -n +2 {input.fastp} | cut -d'\t' -f9)
            size=$(tail -n +2 {input.fastp} | cut -d'\t' -f11)
            reads=$(tail -n +2 {input.clustering} | cut -d'\t' -f2)
            pseudo=$(tail -n +2 {input.clustering} | cut -d'\t' -f6)
            clustered=$(tail -n +2 {input.clustering} | cut -d'\t' -f16)
            otu=$(tail -n +2 {input.tax} | cut -d'\t' -f2)
            assigned=$(tail -n +2 {input.compo} | awk '$2 != "No match"' | cut -d'\t' -f5 | awk '{{s+=$1}}END{{print s}}')
            spec=$(tail -n +2 {input.tax} | cut -d'\t' -f5)
            gen=$(tail -n +2 {input.tax} | cut -d'\t' -f7)
            high=$(($(tail -n +2 {input.tax} | cut -d'\t' -f9) + $(tail -n +2 {input.tax} | cut -d'\t' -f11)))
            none=$(tail -n +2 {input.tax} | cut -d'\t' -f3)

            echo "{wildcards.sample}\t$Q30\t$size\t$reads\t$pseudo\t$clustered\t$otu\t$assigned\t$spec\t$gen\t$high\t$none" >> {output.report}
        fi
        """


rule collect_summaries:
    input:
        report=expand("{sample}/reports/{sample}_summary.tsv", sample=samples.index),
    output:
        agg=report(
            "reports/summary.tsv",
            caption="../report/summary.rst",
            category="Quality controls",
        ),
    message:
        "Aggregating summary reports"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/collect_summaries.log",
    shell:
        """
        exec 2> {log}
        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """


# Rmarkdown reports rules -----------------------------------------------------


rule report_sample:
    input:
        summary="{sample}/reports/{sample}_summary.tsv",
        fastp="{sample}/reports/{sample}_trimmed.tsv",
        qc_filtering="{sample}/reports/{sample}_merging.tsv"
        if config["cluster"]["method"] == "otu"
        else "{sample}/reports/{sample}_denoising.tsv",
        clustering="{sample}/reports/{sample}_clustering.tsv"
        if config["cluster"]["method"] == "otu"
        else "{sample}/reports/{sample}_denoising.tsv",
        blast_rep="{sample}/reports/{sample}_blast_stats.tsv",
        taxonomy="{sample}/reports/{sample}_taxonomy_assignement_stats.tsv",
        result="{sample}/reports/{sample}_composition.tsv",
        db="reports/db_versions.tsv",
        soft="reports/software_versions.tsv",
        krona="{sample}/reports/{sample}_krona_chart.html",
    params:
        method=config["cluster"]["method"],
        workdir=config["workdir"],
        version=version,
        sample=lambda w, input: w.sample,
    output:
        report=report(
            "{sample}/reports/{sample}_report.html",
            caption="../report/markdown_sample.rst",
            category="Results",
            subcategory="{wildcards.sample}",
        ),
    conda:
        "../envs/rmarkdown.yaml"
    log:
        "logs/{sample}/report.log",
    message:
        "Generating html report for {wildcards.sample}"
    script:
        "../scripts/write_report.Rmd"


rule report_all:
    input:
        summary="reports/summary.tsv",
        fastp="reports/fastp_stats.tsv",
        qc_filtering="reports/merging_stats.tsv"
        if config["cluster"]["method"] == "otu"
        else "reports/denoising.tsv",
        clustering="reports/clustering_stats.tsv"
        if config["cluster"]["method"] == "otu"
        else "reports/denoising.tsv",
        blast_rep="reports/blast_stats.tsv",
        taxonomy="reports/taxonomy_assignement_stats.tsv",
        result="reports/composition_summary.tsv",
        db="reports/db_versions.tsv",
        soft="reports/software_versions.tsv",
        krona="reports/krona_chart.html",
    params:
        method=config["cluster"]["method"],
        workdir=config["workdir"],
        version=version,
        sample="all",
    output:
        report=report(
            "reports/report.html",
            caption="../report/markdown_glob.rst",
            category="Results",
            subcategory="Global",
        ),
    conda:
        "../envs/rmarkdown.yaml"
    log:
        "logs/all/report.log",
    message:
        "Generating global html report"
    script:
        "../scripts/write_report.Rmd"


# Version reports -------------------------------------------------------------


rule software_versions:
    output:
        report="reports/software_versions.tsv",
    message:
        "Collecting software versions"
    params:
        dir=f"{workflow.basedir}/workflow/envs/",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/common/software_version.log",
    shell:
        """
        echo "Package\tVersion" > {output.report}
        for env in $(ls {params.dir}*.yaml)
          do
          cat $env \
            | tr "\n" "@" \
            | sed -E 's/(.*)dependencies:(.*)/\2/' \
            | sed -E 's/\s{2}-\s{1}/\n/' \
            | tr -d " " \
            | tr "@" "\n" \
            | awk 'NF' \
            | tr "=" "\t" \
            | sed -E 's/^-//' \
            >> {output.report}
        done
        """


rule database_version:
    output:
        report="reports/db_versions.tsv",
    message:
        "Collecting databases versions"
    params:
        blast=config["blast"]["blast_DB"],
        taxdb=config["blast"]["taxdb"],
        taxdump_nodes=config["taxonomy"]["nodes_dmp"],
        taxdump_lin=config["taxonomy"]["rankedlineage_dmp"],
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/common/db_versions.log",
    shell:
        """
        exec 2> {log}

        echo "Database\tLast modified\tFull path" \
            > {output.report}

        paste \
            <(echo "BLAST") \
            <(date +%F -r {params.blast}.nto) \
            <(echo {params.blast}) \
            >> {output.report}

        paste \
            <(echo "taxdb.bti") \
            <(date +%F -r {params.taxdb}/taxdb.bti) \
            <(echo {params.taxdb}/taxdb.bti) \
            >> {output.report}

        paste \
            <(echo "taxdb.btd") \
            <(date +%F -r {params.taxdb}/taxdb.btd) \
            <(echo {params.taxdb}/taxdb.btd) \
            >> {output.report}

        paste \
            <(echo "taxdump lineages") \
            <(date +%F -r {params.taxdump_lin}) \
            <(echo {params.taxdump_lin}) \
            >> {output.report}

        paste \
            <(echo "taxdump nodes") \
            <(date +%F -r {params.taxdump_nodes}) \
            <(echo {params.taxdump_nodes}) \
            >> {output.report}
        """
