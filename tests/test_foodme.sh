BASEDIR=$(dirname "$0")
cd $BASEDIR
cd ..

python foodme.py -l tests/samples.tsv \
    -d tests/output \
    --rankedlineage_dmp tests/minitaxdump/minirankedlineage.dmp \
    --nodes_dmp tests/minitaxdump/mininodes.dmp \
    --taxdb tests/miniblast \
    --blastdb tests/miniblast/miniblast \
    --primers_fasta tests/primers/16S.fa \
    --denoise \
    --trim_3end \
    -T 1 \
    -c ~/anaconda3/envs