shell.executable("bash")

# Use {bchmk_sample} instead of {sample} to make explicit that we are 
# working with 'benchmark.index' here and not with 'samples.index'
# See common.smk for the index specifications


rule confusion_matrix:
    input:
        compo="{bchmk_sample}/reports/{bchmk_sample}_composition.tsv",
        truth=lambda wildcards: get_sample_reference(wildcards),
        tax = "common/taxonomy.json",
    output:
        confmat="{bchmk_sample}/benchmarking/{bchmk_sample}_confusion_matrix.tsv",
    params:
        threshold=config["benchmark"]["threshold"],
        target_rank=config["benchmark"]["target_rank"],
        sample={wildcards.bchmk_sample},
    message:
        "Finding out the truth for {wildcards.bchmk_sample}"
    conda:
        "../envs/taxidtools.yaml"
    log:
        "logs/{bchmk_sample}/confusion_matrix.log"
    script:
        "../scripts/confusion_matrix.py"


rule collect_confusion_matrices:
    input:
        report=expand("{bchmk_sample}/benchmarking/{bchmk_sample}_confusion_matrix.tsv", sample=benchmark.index),
    output:
        agg=report(
            "benchmarking/confusion_matrix.tsv",
            caption="../report/confusion_matrix.rst",
            category="Benchmarking",
        ),
    message:
        "Rassembling the truth"
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


# rule runtime:


# rule yield:


# rule PRcurve:


# rule PRcurve_all:


# rule metrics:
    # precision
    # recall
    # PR-AUC
    # L2 distance
    # Mean absolute error
    # yield
