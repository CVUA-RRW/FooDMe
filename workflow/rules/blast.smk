shell.executable("bash")


#  Rules DB Masking -----------------------------------------------------------


rule prep_taxonomy:
    output:
        tax="common/taxonomy.json",
    params:
        nodes=config["nodes_dmp"],
        rankedlineage=config["rankedlineage_dmp"],
        taxid=config["taxid_filter"],
    message:
        "[Common][assignement] preparing taxonomy definitions"
    conda:
        "../envs/taxidtools.yaml"
    log:
        "logs/common/taxonomy_prep.log",
    script:
        "../scripts/filter_taxonomy.py"


rule get_taxid_from_db:
    output:
        taxlist="common/taxid_list.txt",
    params:
        blast_DB=config["blast_DB"],
        taxdb=config["taxdb"],
    message:
        "[Common][assignement] collecting BLAST database entries"
    conda:
        "../envs/blast.yaml"
    log:
        "logs/common/taxid_from_db.log",
    shell:
        """
        exec 2> {log}

        export BLASTDB={params.taxdb}

        blastdbcmd -db {params.blast_DB} -tax_info -outfmt %T \
        > {output.taxlist}
        """


rule create_blast_mask:
    input:
        taxlist="common/taxid_list.txt",
        tax="common/taxonomy.json",
    output:
        mask="common/taxid_mask.txt",
    message:
        "[Common][assignement] masking BLAST Database"
    params:
        taxid=config["taxid_filter"],
    message:
        "Preparing list of searchable taxids"
    conda:
        "../envs/taxidtools.yaml"
    log:
        "logs/common/blast_mask.log",
    script:
        "../scripts/make_blast_mask.py"


rule apply_blocklist:
    input:
        taxids=get_mask(),
        blocklist=get_blocklist(),
    output:
        mask="common/blast_mask.txt",
    message:
        "[Common][assignement] applying taxid blocklist"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/common/blocklist.log",
    script:
        "../scripts/apply_blocklist.py"


rule no_masking:
    output:
        mask=temp("common/nomask"),
    message:
        "[Common][assignement] Skipping BLAST database masking"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/common/blocklist.log",
    shell:
        """
        touch {output.mask} 2> {log}
        """


rule no_blocklist:
    output:
        block=temp("common/noblock"),
    message:
        "[Common][assignement] skipping taxid blocklist cration"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/common/blocklist.log",
    shell:
        """
        touch {output.block} > {log}
        """


# Rules Blast -----------------------------------------------------------------


rule blast_otus:
    input:
        query="{sample}/clustering/{sample}_OTUs.fasta"
        if config["cluster_method"] == "otu"
        else "{sample}/denoising/{sample}_ASVs.fasta",
        mask="common/blast_mask.txt",
    output:
        report="{sample}/taxonomy/{sample}_blast_report.tsv",
    params:
        blast_DB=config["blast_DB"],
        taxdb=config["taxdb"],
        e_value=config["blast_evalue"],
        perc_identity=config["blast_identity"],
        qcov=config["blast_qcov"],
    threads: config["threads_sample"]
    message:
        "[{wildcards.sample}][assignement] BLASTing clusters against local database"
    conda:
        "../envs/blast.yaml"
    log:
        "logs/{sample}/blast.log",
    shell:
        """
        export BLASTDB={params.taxdb}

        if [ {input.mask} = "common/blast_mask.txt" ]
        then
            masking="-taxidlist common/blast_mask.txt"
        else
            masking=""
        fi

        blastn -db {params.blast_DB} \
            -query {input.query} \
            -out {output.report} \
            -task 'megablast' \
            -evalue {params.e_value} \
            -perc_identity {params.perc_identity} \
            -qcov_hsp_perc {params.qcov} $masking \
            -outfmt '6 qseqid sseqid evalue pident bitscore sacc staxid length mismatch gaps stitle' \
            -num_threads {threads} \
        2> {log} 

        sed -i '1 i\query\tsubject\tevalue\tidentity\tbitscore\tsubject_acc\tsubject_taxid\talignment_length\tmismatch\tgaps\tsubject_name' {output.report}
        """


rule filter_blast_acc:
    input:
        report="{sample}/taxonomy/{sample}_blast_report.tsv",
    output:
        report="{sample}/taxonomy/{sample}_blast_report_prefiltered.tsv",
    params:
        acc_list=config["seq_blocklist"],
    message:
        "[{wildcards.sample}][assignement] remove negative sequences from BLAST results"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{sample}/filter_blast_acc.log",
    shell:
        """
        exec 2> {log}
        if [ -s {input.report} ]; then
          grep -v -f {params.acc_list} {input.report} > {output.report}
        else
          touch {output.report}
        fi
        """


rule filter_blast_bitscores:
    input:
        report=lambda wildcards: get_acc_blocklist(wildcards),
    output:
        filtered="{sample}/taxonomy/{sample}_blast_report_filtered.tsv",
    params:
        bit_diff=config["bit_score_diff"],
    message:
        "[{wildcards.sample}][assignement] filtering BLAST results"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{sample}/filter_blast_bitscore.log",
    script:
        "../scripts/filter_blast.py"


rule find_consensus:
    input:
        blast="{sample}/taxonomy/{sample}_blast_report_filtered.tsv",
        tax="common/taxonomy.json",
    output:
        consensus="{sample}/taxonomy/{sample}_consensus_table.tsv",
    params:
        min_consensus=config["min_consensus"],
    message:
        "[{wildcards.sample}][assignement] consensus taxonomy determination"
    log:
        "logs/{sample}/find_consensus.log",
    conda:
        "../envs/taxidtools.yaml"
    script:
        "../scripts/min_consensus_filter.py"


# Rules reports ----------------------------------------------------------------


rule blast_stats:
    input:
        otus="{sample}/clustering/{sample}_OTUs.fasta"
        if config["cluster_method"] == "otu"
        else "{sample}/denoising/{sample}_ASVs.fasta",
        blast=lambda wildcards: get_acc_blocklist(wildcards),
        filtered="{sample}/taxonomy/{sample}_blast_report_filtered.tsv",
        lca="{sample}/taxonomy/{sample}_consensus_table.tsv",
    output:
        "{sample}/reports/{sample}_blast_stats.tsv",
    params:
        bit_diff=config["bit_score_diff"],
    message:
        "[{wildcards.sample}][assignement] collecting BLAST stats"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{sample}/blast_stats.log",
    shell:
        """
        exec 2> {log}

        if [ -s {input.blast} ]
        then
            # Get list of all OTUs
            OTUs=$(grep "^>" {input.otus} | cut -d";" -f1 | tr -d '>' | sort -u)

            for otu in $OTUs
            do
                size=$(grep -E "^>${{otu}}\>" {input.otus}  | cut -d"=" -f2)
                bhits=$(grep -c -E "^${{otu}};" {input.blast} || true)
                if [ $bhits -eq 0 ]
                then
                    # When there is no blast hit
                    echo "{wildcards.sample}\t$otu\t$size\t0\t0\t0\t0\t0\t-\t-\t-\t- (1.0)\t../{input.blast}\t../{input.filtered}" >> {output}
                else
                    # Otherwise collect and print stats to file
                    bit_best=$(grep -E "^${{otu}};" {input.blast} | cut -f5 | cut -d. -f1 | sort -rn | head -n1)
                    bit_low=$(grep -E "^${{otu}};" {input.blast} | cut -f5 | cut -d. -f1 | sort -n | head -n1)
                    bit_thr=$(($bit_best - {params.bit_diff}))
                    shits=$(grep -c -E "^${{otu}}\>" {input.filtered})
                    cons=$(grep -E "^${{otu}}\>" {input.lca} | cut -d'\t' -f2-5)

                    echo "{wildcards.sample}\t$otu\t$size\t$bhits\t$bit_best\t$bit_low\t$bit_thr\t$shits\t$cons\t../{input.blast}\t../{input.filtered}" >> {output}
                fi
            done
            # Sort by size and add header (just to get hits on top)
            sort -k3,3nr -o {output} {output}
            sed -i '1 i\Sample\tQuery\tCount\tBlast hits\tBest bit-score\tLowest bit-score\tBit-score threshold\tSaved Blast hits\tConsensus\tRank\tTaxid\tDisambiguation\tlink_report\tlink_filtered' {output}

        else
            echo "{wildcards.sample}\t-\t-\t0\t0\t0\t0\t0\t-\t-\t-\t-\t../{input.blast}\t../{input.filtered}" > {output}
            sed -i '1 i\Sample\tQuery\tCount\tBlast hits\tBest bit-score\tLowest bit-score\tBit-score threshold\tSaved Blast hits\tConsensus\tRank\tTaxid\tDisambiguation\tlink_report\tlink_filtered' {output}
        fi
        """


rule collect_blast_stats:
    input:
        report=expand("{sample}/reports/{sample}_blast_stats.tsv", sample=samples.index),
    output:
        agg=report(
            "reports/blast_stats.tsv",
            caption="../report/blast_stats.rst",
            category="Quality controls",
        ),
    message:
        "[All][assignement] aggregating BLAST stats"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/blast_stats.log",
    shell:
        """
        head -n 1 {input.report[0]} > {output.agg}
        for i in {input.report}; do 
          cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """


rule tax_stats:
    input:
        "{sample}/reports/{sample}_blast_stats.tsv",
    output:
        "{sample}/reports/{sample}_taxonomy_assignement_stats.tsv",
    message:
        "[{wildcards.sample}][assignement] collecting taxonomy assignement stats"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{sample}/taxonomy_stats.log",
    shell:
        """
        exec 2> {log}

        echo "Sample\tQuery\tUnknown sequences\tUnknown sequences [%]\t(Sub-)Species consensus\t(Sub-)Species consensus [%]\tGenus consensus\tGenus consensus [%]\tFamily consensus\tFamily consensus [%]\tHigher rank consensus\tHigher rank consensus [%]" > {output}

        all=$(grep -c -E "OTU_|ASV_" <(tail -n +2 {input}) || true)
        nohits=$(grep -c "[[:blank:]]-[[:blank:]]" {input} || true)
        spec=$(grep -c "species" {input} || true)
        gen=$(grep -c "genus" {input} || true)
        fam=$(grep -c "family" {input} || true)
        other=$(( $all - $nohits - $spec - $gen - $fam ))

        if [ $all -ne 0 ]
        then
            nohits_perc=$(printf %.2f "$((10**3 * (100* $nohits / $all)))e-3")
            spec_perc=$(printf %.2f "$((10**3 * (100* $spec / $all)))e-3")
            gen_perc=$(printf %.2f "$((10**3 * (100* $gen / $all)))e-3")
            fam_perc=$(printf %.2f "$((10**3 * (100* $fam / $all)))e-3")
            other_perc=$(printf %.2f "$((10**3 * (100* $other / $all)))e-3")

            echo "{wildcards.sample}\t$all\t$nohits\t$nohits_perc\t$spec\t$spec_perc\t$gen\t$gen_perc\t$fam\t$fam_perc\t$other\t$other_perc" >> {output}

        else
            echo "{wildcards.sample}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" >> {output}
        fi
        """


rule collect_tax_stats:
    input:
        report=expand(
            "{sample}/reports/{sample}_taxonomy_assignement_stats.tsv",
            sample=samples.index,
        ),
    output:
        agg=report(
            "reports/taxonomy_assignement_stats.tsv",
            caption="../report/taxonomic_ass_stats.rst",
            category="Quality controls",
        ),
    message:
        "[All][assignement] collecting taxonomy assignement stats"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/taxonomy_stats.log",
    shell:
        """
        exec 2> {log}
        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """


rule summarize_results:
    input:
        compo="{sample}/reports/{sample}_blast_stats.tsv",
    output:
        report=report(
            "{sample}/reports/{sample}_composition.tsv",
            caption="../report/compo_sample.rst",
            category="Results",
            subcategory="{wildcards.sample}",
        ),
    params:
        sample_name=lambda w, input: w.sample,
    message:
        "[{wildcards.sample}][assignement] summarizing results"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{sample}/summarize_results.log",
    script:
        "../scripts/summarize_results.py"


rule collect_results:
    input:
        report=expand("{sample}/reports/{sample}_composition.tsv", sample=samples.index),
    output:
        agg=report(
            "reports/composition_summary.tsv",
            caption="../report/compo_glob.rst",
            category="Results",
            subcategory="Global",
        ),
    message:
        "[All][assignement] aggregating compositions"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/all/collect_results.log",
    shell:
        """
        exec 2> {log}

        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """
