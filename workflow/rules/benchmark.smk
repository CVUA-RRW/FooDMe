shell.executable("bash")

# Use {bchmk_sample} with index 'benchmark_index' instead of {sample} to make explicit that we are
# working with 'benchmark.index' here and not with 'samples.index'
# See common.smk for the index specifications


rule confusion_matrix:
    input:
        compo="{bchmk_sample}/reports/{bchmk_sample}_composition.tsv",
        truth=config["benchmark_reference"],
        tax="common/taxonomy.json",
    output:
        confmat="{bchmk_sample}/benchmarking/{bchmk_sample}_confusion_matrix.tsv",
    params:
        threshold=config["benchmark_threshold"],
        target_rank=config["benchmark_rank"],
        sample=lambda w: w.bchmk_sample,
    message:
        "[{wildcards.bchmk_sample}][benchamrking] calculating confusion table"
    conda:
        "../envs/taxidtools.yaml"
    log:
        "logs/{bchmk_sample}/confusion_matrix.log",
    script:
        "../scripts/confusion_matrix.py"


rule collect_confusion_matrices:
    input:
        report=expand(
            "{bchmk_sample}/benchmarking/{bchmk_sample}_confusion_matrix.tsv",
            bchmk_sample=benchmark_index,
        ),
    output:
        agg=report(
            "benchmarking/confusion_matrix.tsv",
            caption="../report/confusion_matrix.rst",
            category="Benchmarking",
        ),
    message:
        "[All][benchamrking] collecting confusion tables"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/confusion_matrix.log",
    shell:
        """
        exec 2> {log}
        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """


rule yields:
    input:
        summary="{bchmk_sample}/reports/{bchmk_sample}_summary.tsv",
    output:
        yields="{bchmk_sample}/benchmarking/{bchmk_sample}_yield.tsv",
    message:
        "[{wildcards.bchmk_sample}][benchamrking] calculating yield"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{bchmk_sample}/yield.log",
    shell:
        """
        exec 2> {log}

        # Read data in
        IFS='\t'
        read -r sample q30 size read_in merged clustered n_cluster assigned sp gn hr nm <<< $(cat {input.summary} | tail -n +2)

        # Calculate yields
        if [[ $read_in -eq 0 ]]
        then
            merged_perc=0
            clustered_perc=0
            assigned_perc=0
        else
            merged_perc=$(printf %.2f "$((10**3 * (100* $merged / $read_in)))e-3")
            clustered_perc=$(printf %.2f "$((10**3 * (100* $clustered / $read_in)))e-3")
            assigned_perc=$(printf %.2f "$((10**3 * (100* $assigned / $read_in)))e-3")
        fi

        # Write report
        echo "Sample\tTotal reads [%]\tMerged reads [%]\tClustered reads [%]\tAssigned reads [%]" > {output.yields}
        echo "$sample\t100.00\t$merged_perc\t$clustered_perc\t$assigned_perc" >> {output.yields}
        """


rule collect_yield:
    input:
        report=expand(
            "{bchmk_sample}/benchmarking/{bchmk_sample}_yield.tsv",
            bchmk_sample=benchmark_index,
        ),
    output:
        agg=report(
            "benchmarking/yield.tsv",
            caption="../report/yields.rst",
            category="Benchmarking",
        ),
    message:
        "[All][benchamrking] collecting yield reports"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/yield.log",
    shell:
        """
        exec 2> {log}
        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """


rule metrics_sample:
    input:
        confmat="{bchmk_sample}/benchmarking/{bchmk_sample}_confusion_matrix.tsv",
    output:
        metrics="{bchmk_sample}/benchmarking/{bchmk_sample}_metrics.tsv",
    params:
        sample=lambda w: w.bchmk_sample,
    message:
        "[{wildcards.bchmk_sample}][benchamrking] calculating metrics for {wildcards.bchmk_sample}"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{bchmk_sample}/metrics.log",
    script:
        "../scripts/benchmark_metrics.py"


rule metrics_global:
    input:
        confmat="benchmarking/confusion_matrix.tsv",
    output:
        metrics="aggregated_samples/benchmarking/aggregated_metrics.tsv",
    params:
        sample="all",
    message:
        "[All][benchamrking] calculating metrics for aggregated samples"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/agg_metrics.log",
    script:
        "../scripts/benchmark_metrics.py"


rule collect_benchmarking_metrics:
    input:
        agg_metrics="aggregated_samples/benchmarking/aggregated_metrics.tsv",
        metrics=expand(
            "{bchmk_sample}/benchmarking/{bchmk_sample}_metrics.tsv",
            bchmk_sample=benchmark_index,
        ),
    output:
        agg=report(
            "benchmarking/metrics.tsv",
            caption="../report/metrics.rst",
            category="Benchmarking",
        ),
    message:
        "[All][benchamrking] collecting metrics reports"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/metrics.log",
    shell:
        """
        exec 2> {log}
        cat {input.agg_metrics} > {output.agg}
        for i in {input.metrics}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """


rule prcurve:
    input:
        confmat="benchmarking/confusion_matrix.tsv",
    output:
        pr_curve="benchmarking/pr_curve.tsv",
    message:
        "[All][benchamrking] calculating precision-recall curve"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/pr_curve.log",
    script:
        "../scripts/pr-curve.py"


rule benchmarking_report:
    input:
        confmat="benchmarking/confusion_matrix.tsv",
        metrics="benchmarking/metrics.tsv",
        pr_curve="benchmarking/pr_curve.tsv",
        yields="benchmarking/yield.tsv",
        db="reports/db_versions.tsv",
        soft="reports/software_versions.tsv",
    output:
        report=report(
            "benchmarking/benchmarking_report.html",
            caption="../report/bchm_markdown.rst",
            category="Benchmarking",
        ),
    params:
        workdir=config["workdir"],
        version=version,
        logo=f"{workflow.basedir}/../docs/logo.png",
    message:
        "[All][benchamrking] creating benchmarking report"
    conda:
        "../envs/rmarkdown.yaml"
    log:
        "logs/all/benchmark_report.log",
    script:
        "../scripts/benchmark_report.Rmd"
