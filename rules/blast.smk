import pandas as pd
from os import stat

shell.executable("bash")

# Functions --------------------------------------------------------------------

def get_mask():
    if config["blast"]["taxid_filter"] == "None":
        return "common/nomask"
    else:
        return "common/taxid_mask.txt"

def concatenate_uniq(entries):
    s = str(entries.to_list()[0])
    uniq = list(set(s.split("; ")))
    uniq_no_zero = ["-" if e=="0" else e for e in uniq]
    return "; ".join(uniq_no_zero)

#  Rules DB Masking ------------------------------------------------------------

rule get_taxid_from_db:
    output:
        "common/taxid_list.txt",
    params:
        blast_DB = config["blast"]["blast_DB"],
        taxdb = config["blast"]["taxdb"],
    conda:
        "../envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        blastdbcmd -db {params.blast_DB} -tax_info -outfmt %T > {output}
        """

rule create_blast_mask:
    input:
        "common/taxid_list.txt",
    output:
        "common/taxid_mask.txt",
    params:
        taxid = config["blast"]["taxid_filter"],
        lineage = config["taxonomy"]["rankedlineage_dmp"],
        nodes = config["taxonomy"]["nodes_dmp"],
    message:
        "Preparing list of searchable taxids"
    script:
        "../scripts/make_blast_mask.py"

rule no_masking:
    output:
        temp("common/nomask"),
    shell:
        """
        touch {output}
        """

# Rules Blast ------------------------------------------------------------------

rule blast_otus:
    input: 
        query = "{sample}/clustering/{sample}_OTUs.fasta" if config["cluster"]["method"] == "otu" else "{sample}/denoising/{sample}_ASVs.fasta",
        mask = get_mask(),
    output:
        "{sample}/taxonomy/{sample}_blast_report.tsv",
    params:
        blast_DB = config["blast"]["blast_DB"],
        taxdb = config["blast"]["taxdb"],
        e_value = config["blast"]["e_value"],
        perc_identity = config["blast"]["perc_identity"],
        qcov = config["blast"]["qcov"],
    threads:
        config["threads"]
    message:
        "BLASTing {wildcards.sample} clusters against local database"
    conda:
        "../envs/blast.yaml"
    log:
        "logs/{sample}_blast.log"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        if [ {input.mask} = "common/taxid_mask.txt" ]
        then
            masking="-taxidlist common/taxid_mask.txt"
        else
            masking=""
        fi
        
        blastn -db {params.blast_DB} \
            -query {input.query} \
            -out {output} \
            -task 'megablast' \
            -evalue {params.e_value} \
            -perc_identity {params.perc_identity} \
            -qcov_hsp_perc {params.qcov} $masking \
            -outfmt '6 qseqid sseqid evalue pident bitscore sacc staxid length mismatch gaps stitle' \
            -num_threads {threads} \
        > logs/blast.log 2>&1
        
        sed -i '1 i\query\tsubject\tevalue\tidentity\tbitscore\tsubject_acc\tsubject_taxid\talignment_length\tmismatch\tgaps\tsubject_name' {output}
        """

rule filter_blast:
    input:
        "{sample}/taxonomy/{sample}_blast_report.tsv",
    output:
        "{sample}/taxonomy/{sample}_blast_report_filtered.tsv",
    params:
        bit_diff= config["blast"]["bit_score_diff"],
    message:
        "Filtering BLAST results for {wildcards.sample}"
    run:
        if stat(input[0]).st_size == 0:
            with open(output[0], 'w') as fout:
                fout.write("query\tsubject\tevalue\tidentity\tbitscore\tsubject_acc\tsubject_taxid\talignment_length\tmismatch\tgaps\tsubject_name")
        else:
            df = pd.read_csv(input[0], sep = '\t', header= 0)
            
            if df.empty:
                df.to_csv(output[0], sep='\t', header=True, index=False)
            
            else:
                sd = dict(tuple(df.groupby('query')))

                dfout = pd.DataFrame()

                for key, val in sd.items():
                    dfout = dfout.append( val[val['bitscore'] >= max(val['bitscore']) - params.bit_diff] )
                
                dfout['query'] = dfout['query'].str.split(';').str[0]

                dfout.to_csv(output[0], sep='\t', header=True, index=False)

rule blast2lca:
    input:
        "{sample}/taxonomy/{sample}_blast_report_filtered.tsv",
    output:
        "{sample}/taxonomy/{sample}_consensus_table.tsv" ,
    params:
        lineage = config["taxonomy"]["rankedlineage_dmp"],
        nodes = config["taxonomy"]["nodes_dmp"],
    message:
        "Lowest common ancestor determination"
    script:
        "../scripts/blast_to_lca.py"

# Rules reports ----------------------------------------------------------------

rule blast_stats:
    input:
        otus = "{sample}/clustering/{sample}_OTUs.fasta" if config["cluster"]["method"] == "otu" else "{sample}/denoising/{sample}_ASVs.fasta",
        blast = "{sample}/taxonomy/{sample}_blast_report.tsv",
        filtered = "{sample}/taxonomy/{sample}_blast_report_filtered.tsv",
        lca = "{sample}/taxonomy/{sample}_consensus_table.tsv",
    output:
        "{sample}/reports/{sample}_blast_stats.tsv",
    params:
        bit_diff= config["blast"]["bit_score_diff"],
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
                    echo "{wildcards.sample}\t$otu\t$size\t0\t0\t0\t0\t0\t-\t-\t-\t-" >> {output}
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
            echo "{wildcards.sample}\t-\t-\t0\t0\t0\t0\t0\t-\t-\t-\t-" > {output}
            sed -i '1 i\Sample\tQuery\tCount\tBlast hits\tBest bit-score\tLowest bit-score\tBit-score threshold\tSaved Blast hits\tConsensus\tRank\tTaxid\tDisambiguation\tlink_report\tlink_filtered' {output}
        fi
        """

rule collect_blast_stats:
    input:
        expand("{sample}/reports/{sample}_blast_stats.tsv", sample = samples.index),
    output:
        "reports/blast_stats.tsv",
    message:
        "Collecting BLAST stats"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
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
       expand("{sample}/reports/{sample}_taxonomy_assignement_stats.tsv", sample = samples.index),
    output:
        "reports/taxonomy_assignement_stats.tsv",
    message:
        "Collecting taxonomy assignement stats"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """

rule summarize_results:
    input:
        compo = "{sample}/reports/{sample}_blast_stats.tsv",
    output:
        report = "{sample}/reports/{sample}_composition.tsv",
    message:
        "Summarizing results for {wildcards.sample}"
    run:
        df = pd.read_csv(input.compo, sep="\t", header=0).fillna(0)

        # Empty input case
        if len(df["Query"]) == 1 and  df["Query"].head(1).item() == "-":
            with open(output.report, 'w') as fout:
                fout.write("Sample\tConsensus\tRank\tTaxid\tCount\tDisambiguation\tPercent of total")
        else:
            groups = df.groupby(['Consensus', 'Rank', 'Taxid']).agg({'Count': 'sum', 'Disambiguation': concatenate_uniq})
            groups = groups.sort_values("Count", ascending = False).reset_index()
            groups['perc'] = round(groups['Count']/groups['Count'].sum() *100, 2)
            groups.insert(0, 'Sample', wildcards.sample)
            groups.rename(columns = {"perc":"Percent of total"},inplace = True)
            groups["Consensus"].replace({"-": "No match"}, inplace = True)
            groups["Taxid"].replace({0: "-"}, inplace = True)
            groups.to_csv(output.report, sep = "\t", index = False)
        
rule collect_results:
    input:
        expand("{sample}/reports/{sample}_composition.tsv", sample = samples.index),
    output:
        "reports/composition_summary.tsv",
    message:
        "Collecting results"
    shell:
        """
        cat {input[0]} | head -n 1 > {output}
        for i in {input}; do 
            cat ${{i}} | tail -n +2 >> {output}
        done
        """
