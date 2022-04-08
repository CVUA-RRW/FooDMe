import pandas as pd
from os import stat

shell.executable("bash")


#  Rules DB Masking -----------------------------------------------------------


rule prep_taxonomy:
    output:
        tax="common/taxonomy.json",
    params:
        nodes=config["taxonomy"]["nodes_dmp"],
        rankedlineage=config["taxonomy"]["rankedlineage_dmp"],
        taxid=config["blast"]["taxid_filter"],
    message:
        "Preparing taxonomy definitions"
    conda:
        "../envs/taxidtools.yaml"
    script:
        "../scripts/filter_taxonomy.py"


rule get_taxid_from_db:
    output:
        taxlist="common/taxid_list.txt",
    params:
        blast_DB=config["blast"]["blast_DB"],
        taxdb=config["blast"]["taxdb"],
    message:
        "Collecting BLAST database entries"
    conda:
        "../envs/blast.yaml"
    shell:
        """
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
        "Masking BLAST Database"
    params:
        taxid=config["blast"]["taxid_filter"],
    message:
        "Preparing list of searchable taxids"
    conda:
        "../envs/taxidtools.yaml"
    script:
        "../scripts/make_blast_mask.py"


rule apply_blocklist:
    input:
        taxids=get_mask(),
        blocklist=get_blocklist(),
    output:
        mask="common/blast_mask.txt",
    message:
        "Applying taxid blocklist"
    script:
        "../scripts/apply_blocklist.py"


rule no_masking:
    output:
        mask=temp("common/nomask"),
    message:
        "Skipping BLAST database masking"
    shell:
        """
        touch {output.mask}
        """


rule no_blocklist:
    output:
        block=temp("common/noblock"),
    message:
        "Skipping taxid blocklist cration"
    shell:
        """
        touch {output.block}
        """


# Rules Blast -----------------------------------------------------------------


rule blast_otus:
    input:
        query="{sample}/clustering/{sample}_OTUs.fasta"
        if config["cluster"]["method"] == "otu"
        else "{sample}/denoising/{sample}_ASVs.fasta",
        mask="common/blast_mask.txt",
    output:
        report="{sample}/taxonomy/{sample}_blast_report.tsv",
    params:
        blast_DB=config["blast"]["blast_DB"],
        taxdb=config["blast"]["taxdb"],
        e_value=config["blast"]["e_value"],
        perc_identity=config["blast"]["perc_identity"],
        qcov=config["blast"]["qcov"],
    threads: config["threads"]
    message:
        "BLASTing {wildcards.sample} clusters against local database"
    conda:
        "../envs/blast.yaml"
    log:
        "logs/{sample}_blast.log",
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
        > logs/blast.log 2>&1

        sed -i '1 i\query\tsubject\tevalue\tidentity\tbitscore\tsubject_acc\tsubject_taxid\talignment_length\tmismatch\tgaps\tsubject_name' {output.report}
        """


rule filter_blast:
    input:
        report="{sample}/taxonomy/{sample}_blast_report.tsv",
    output:
        filtered="{sample}/taxonomy/{sample}_blast_report_filtered.tsv",
    params:
        bit_diff=config["blast"]["bit_score_diff"],
    message:
        "Filtering BLAST results for {wildcards.sample}"
    run:
        if stat(input["report"]).st_size == 0:
            with open(output["filtered"], "w") as fout:
                fout.write(
                    "query\tsubject\tevalue\tidentity\tbitscore\tsubject_acc\tsubject_taxid\talignment_length\tmismatch\tgaps\tsubject_name"
                )
        else:
            df = pd.read_csv(input["report"], sep="\t", header=0)

            if df.empty:
                df.to_csv(output["filtered"], sep="\t", header=True, index=False)

            else:
                sd = dict(tuple(df.groupby("query")))

                dfout = pd.DataFrame()

                for key, val in sd.items():
                    dfout = pd.concat([dfout, 
                                       val[val["bitscore"] >= max(val["bitscore"]) - params.bit_diff]]
                    )

                dfout["query"] = dfout["query"].str.split(";").str[0]

                dfout.to_csv(output["filtered"], sep="\t", header=True, index=False)


rule find_consensus:
    input:
        blast="{sample}/taxonomy/{sample}_blast_report_filtered.tsv",
        tax="common/taxonomy.json",
    output:
        consensus="{sample}/taxonomy/{sample}_consensus_table.tsv",
    params:
        min_consensus=config["taxonomy"]["min_consensus"],
    message:
        "Consensus taxonomy determination"
    conda:
        "../envs/taxidtools.yaml"
    script:
        "../scripts/min_consensus_filter.py"


# Rules reports ----------------------------------------------------------------


rule blast_stats:
    input:
        otus="{sample}/clustering/{sample}_OTUs.fasta"
        if config["cluster"]["method"] == "otu"
        else "{sample}/denoising/{sample}_ASVs.fasta",
        blast="{sample}/taxonomy/{sample}_blast_report.tsv",
        filtered="{sample}/taxonomy/{sample}_blast_report_filtered.tsv",
        lca="{sample}/taxonomy/{sample}_consensus_table.tsv",
    output:
        "{sample}/reports/{sample}_blast_stats.tsv",
    params:
        bit_diff=config["blast"]["bit_score_diff"],
    message:
        "Collecting BLAST stats for {wildcards.sample}"
    shell:
        """
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
        agg=report("reports/blast_stats.tsv",
                   caption="../report/blast_stats.rst",
                   category="Quality controls"),
    message:
        "Aggregating BLAST stats"
    shell:
        """
        cat {input.report[0]} | head -n 1 > {output.agg}
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
        "Collecting taxonomy assignement stats for {wildcards.sample}"
    shell:
        """
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
        agg=report("reports/taxonomy_assignement_stats.tsv",
                   caption="../report/taxonomic_ass_stats.rst",
                   category="Quality controls")
    message:
        "Collecting taxonomy assignement stats"
    shell:
        """
        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """


rule summarize_results:
    input:
        compo="{sample}/reports/{sample}_blast_stats.tsv",
    output:
        report=report("{sample}/reports/{sample}_composition.tsv",
                      caption="../report/compo_sample.rst",
                      category="Results",
                      subcategory="{wildcards.sample}"),
    message:
        "Summarizing results for {wildcards.sample}"
    run:
        df = pd.read_csv(input.compo, sep="\t", header=0).fillna(0)

        # Empty input case
        if len(df["Query"]) == 1 and df["Query"].head(1).item() == "-":
            with open(output.report, "w") as fout:
                fout.write(
                    "Sample\tConsensus\tRank\tTaxid\tCount\tDisambiguation\tPercent of total"
                )
        else:
            groups = df.groupby(["Consensus", "Rank", "Taxid"]).agg(
                {"Count": "sum", "Disambiguation": concatenate_uniq}
            )
            groups = groups.sort_values("Count", ascending=False).reset_index()
            assigned, notassigned = (
                groups[groups["Consensus"] != "-"],
                groups[groups["Consensus"] == "-"],
            )
            assigned["perc"] = round(groups["Count"] / groups["Count"].sum() * 100, 2)
            notassigned["perc"] = "-"
            groups = pd.concat([assigned, notassigned])
            groups.insert(0, "Sample", wildcards.sample)
            groups.rename(columns={"perc": "Percent of total"}, inplace=True)
            groups["Consensus"].replace({"-": "No match"}, inplace=True)
            groups["Taxid"].replace({0: "-"}, inplace=True)
            groups.to_csv(output.report, sep="\t", index=False)


rule collect_results:
    input:
        report=expand("{sample}/reports/{sample}_composition.tsv", sample=samples.index),
    output:
        agg=report("reports/composition_summary.tsv",
                   caption="../report/compo_glob.rst",
                   category="Results",
                   subcategory="Global")
    message:
        "Aggregating compositions"
    shell:
        """
        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """
