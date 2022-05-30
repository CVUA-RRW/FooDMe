#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import numpy as np
import pandas as pd
import taxidTools as txd


sys.stderr = open(snakemake.log[0], "w")


def closest_node(taxid, ref_taxids, tax):
    """
    Find closest taxid in a list of refs
    
    returns closest reference
    """
    lcas = [tax.lca([taxid, tx]) for tx in ref_taxids]
    dists = [tax.distance(taxid, tx.taxid) for tx in lcas]
    return ref_taxids[np.argmin(dists)]


def format_rank(taxid, tax, ranks): # ranks can pass custom classification
    l=txd.Lineage(taxid)
    l.filter(ranks)
    # Missing ranks ranks are replaced by DummyNodes
    return [node for node in l 
            if not isinstance(node, txd.DummyNode)][0]


def main(
        compo,
        truth,
        output,
        sample,
        threshold,
        target_rank,
        taxonomy
        ):
    ### Part 1: Data prep --------------------------------------------

    # Ensure target_rank is linean
    # Waiting to think of a smartest solution
    # best I found so far is to use a normalized list of ranks 
    # to ensure that input rank threhsold is valid and consistent with the reuslts.
    # This limits the choice to 'species' and 'genus' in practice (or higher)
    # but enforces consistency with NCBI taxonomy.
    # --> NCBI only Linnean is conserved between all organisms
    # and the other ranks or 'randomly' present.
    # Solution to think about:
    # hard-code an extendended rank list to use to filter the ranks?
    try:
        assert(target_rank in txd.linne())
        ranks =  txd.linne()
    except AssertionError:
        raise ValueError(f"Parameter 'target_rank' must be in {txd.linne()}")

    tax = txd.load(taxonomy)

    # Load data
    compo_tbl = pd.read_csv(compo, sep="\t")
    truth_tbl = pd.read_csv(truth, sep="\t")

    # Filter truth table
    truth_tbl = truth_tbl.loc[truth_tbl['sample'] == sample]

    # Filter and normalize compo to range [0, 1] (input in %)
    compo_tbl = compo_tbl.loc[
        (compo_tbl['Taxid'] != "-") & (compo_tbl['Percent of assigned'] >= threshold)
        ]
    compo_tbl['Percent of assigned'] = compo_tbl['Percent of assigned']/100

    # Preformating tables
    pred = pd.DataFrame(compo_tbl, 
        columns=['Sample', 'Taxid', 'Count', 'Percent of assigned'], 
        copy = True
        ).rename({
            'Sample': 'sample', 
            'Taxid': 'taxid', 
            'Count': 'read_count', 
            'Percent of assigned': 'pred_ratio'
            }, axis = 1)

    exp = pd.DataFrame(truth_tbl, 
        columns=['sample', 'taxid', 'proportion'], 
        copy = True
        ).rename({'proportion': 'exp_ratio'}, axis = 1)

    ### Part2: Matching expectations to reality ---------------------

    # Find matches 
    matches = [
        closest_node(taxid, exp.loc[:,'taxid'].to_list(), tax) 
        for taxid in pred.loc[:,'taxid'].to_list()
        ]
    pred['ref_match'] = matches

    # Get Lca of exp and pred matches to know what ranks correspond
    norms = [
        format_rank(tax.get(str(tx)), ranks) if tax.isAncestorOf(m,tx)  # speed up if same lineage
        else format_rank(tax.lca([tx,m]), ranks) 
        for tx,m in zip(pred.loc[:,'taxid'], pred.loc[:,'ref_match']) 
        ]
    pred['match_rank'] = [node.rank for node in norms]

    # Post format finished tables
    pred = pred.set_index(
        'taxid'
        ).filter(
            items=['pred_ratio', "ref_match", "match_rank"]
         )
    pred['predicted'] = 1

    exp = exp.set_index('taxid'
        ).filter(items=['exp_ratio']
        )
    exp['expected'] = 1

    ### Part 3: building confusion table ------------------------

    # Index of max accepted rank
    max_index = ranks.index(target_rank)

    # checking correctness of results
    pred['expected'] = [
        0 if ranks.index(mrk) > max_index
        else 1
        for mrk in pred['match_rank']
        ]

    # reindex predicted using ref_match, keep taxid as info and sum ratios below the rank limit
    # This is to avoid merging data on assignement with too high rank for duplicated taxa
    pred['ref_match'] = [
        match_val if ranks.index(match_rank) <= max_index
        else taxid_val
        for taxid_val, match_val, match_rank in zip(pred.index, pred['ref_match'], pred['match_rank'])
        ]

    conftable = pred.reset_index(drop=True).groupby(
                            ["ref_match", "predicted", "expected"]).sum()
    conftable = conftable.reset_index().rename(columns={'ref_match': 'Taxid'})
    conftable = conftable.astype('float64')
    exp = exp.astype('float64')

    conftable = pd.merge(conftable, 
        exp,
        how='outer')

    # Missing values are missagnignements!
    conftable = conftable.fillna(0)

    contable.to_csv(output, sep="\t", header=True, index=False)


if __name__ == '__main__':
    main(
        compo = snakemake.input['compo'],
        truth = snakemake.input['truth'],
        output = snakemake.output['confmat'],
        sample = snakemake.params['params'],
        threshold = snakemake.params['threshold'],
        target_rank = snakemake.params['target_rank'],
        taxonomy = snakemake.input['tax']
        )
