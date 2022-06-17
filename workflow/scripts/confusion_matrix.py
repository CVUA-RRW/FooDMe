#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


sys.stderr = open(snakemake.log[0], "w")


import numpy as np
import pandas as pd
import taxidTools as txd


def closest_node(taxid, ref_taxids, tax):
    """
    Find closest taxid in a list of refs
    
    returns closest reference
    """
    lcas = [tax.lca([taxid, tx]) for tx in ref_taxids]
    dists = [tax.distance(taxid, tx.taxid) for tx in lcas]
    return ref_taxids[np.argmin(dists)]


def format_rank(taxid, tax, ranks): # ranks can pass custom classification
    """
    Return next corresponding rank for given taxid and rank list
    """
    l=txd.Lineage(taxid)
    l.filter(ranks)
    # Missing ranks ranks are replaced by DummyNodes
    try:
        next_node = [node for node in l if not isinstance(node, txd.DummyNode)][0]
        
        return next_node.rank
    except IndexError:
        return "root"

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
        ranks =  txd.linne()+['root']
    except AssertionError:
        raise ValueError(f"Parameter 'target_rank' must be in {txd.linne()}")

    tax = txd.load(taxonomy)

    threshold=threshold/100 # Everything in the interval [0,1]

    # Load data
    compo_tbl = pd.read_csv(compo, sep="\t")
    truth_tbl = pd.read_csv(truth, sep="\t")

    # Filter truth table
    truth_tbl = truth_tbl.loc[truth_tbl['sample'] == sample]

    # Filter and normalize compo to range [0, 1] (input in %)
    compo_tbl = compo_tbl.loc[(compo_tbl['Taxid'] != "-")]
    compo_tbl['Percent of assigned'] = compo_tbl['Percent of assigned'].astype(float)/100

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

    ### Part2: Checking the acceptability of predictions ---------------------

    # Find matches 
    pred['ref_match'] = pred.apply(
        lambda x: closest_node(int(x['taxid']), exp.loc[:,'taxid'].astype(int).to_list(), tax), 
        axis=1,
    )

    # Get Lca of exp and pred matches to know what ranks correspond
    pred['match_rank'] = pred.apply(
        lambda x: format_rank(tax.lca([int(x['taxid']),int(x['ref_match'])]), tax, ranks),
        axis=1,
    )

    # Drop unnescessary columns
    pred = pred.set_index(
        'taxid'
        ).filter(
            items=['pred_ratio', "ref_match", "match_rank"]
         )

    # Counts as predicted only if quantif above threshold and rank under limit
    # Index of max accepted rank
    max_index = ranks.index(target_rank)
    
    pred['predicted'] = pred.apply(
        lambda x: 1 if (x['pred_ratio']>=threshold) and (ranks.index(x['match_rank'])<= max_index) else 0,
        axis=1
    )

    exp = exp.set_index('taxid'
        ).filter(items=['exp_ratio']
        )
    exp['expected'] = 1

    ### Part 3: building confusion table ------------------------

    # reindex predicted using ref_match, keep taxid as info and sum ratios below the rank limit
    # This is to avoid merging data on assignement with too high rank for duplicated taxa
    pred['ref_match'] = [
        match_val if ranks.index(match_rank) <= max_index
        else taxid_val
        for taxid_val, match_val, match_rank in zip(pred.index, pred['ref_match'], pred['match_rank'])
        ]

    conftable = pred.reset_index(
                    drop=True
                ).rename(
                    columns={'ref_match': 'Taxid'}
                ).groupby(
                    ["Taxid", "predicted", "match_rank"]
                ).sum()
    conftable = conftable.reset_index(["match_rank", "predicted"])
    conftable = conftable.astype({ "match_rank": str, 'predicted': int, 'pred_ratio': float})
    exp = exp.astype({'expected': int, 'exp_ratio': float})

    conftable = pd.merge(
        conftable, 
        exp,
        how='outer',
        left_index=True,
        right_index=True)

    # Missing values are missagnignements!
    conftable = conftable.fillna(0)
    conftable = conftable.reset_index().rename(columns={'index': 'Taxid'})
    conftable["Sample"] = sample
    conftable = conftable[['Sample', 'Taxid',  "match_rank", 'predicted', 'expected', 'pred_ratio', 'exp_ratio']]

    conftable.to_csv(output, sep="\t", header=True, index=False)


if __name__ == '__main__':
    main(
        compo = snakemake.input['compo'],
        truth = snakemake.input['truth'],
        output = snakemake.output['confmat'],
        sample = snakemake.params['sample'],
        threshold = snakemake.params['threshold'],
        target_rank = snakemake.params['target_rank'],
        taxonomy = snakemake.input['tax']
        )
