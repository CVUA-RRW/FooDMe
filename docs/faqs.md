# Frequently Asked Questions

**Q: Do I need to install mamba?**

**A:** Installing mamba is optional. Simply using conda should be enough, although considerably slower, to 
create environments. If using only conda you may need to add the argument `--conda-frontend conda`
when executing snakemake.

**Q: Can I use FooDMe for single end reads?**

**A:** FooDMe is designed for Illumina paired-end sequencing and single read sequencing is not supported.

**Q: Are IonTorrent data supported?**

**A:** Not yet. But it should be possible with some modifications. [Get in touch](https://github.com/CVUA-RRW/FooDMe/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=) if this is something you need.

**Q: Can I use FooDMe for bacteria metabarcoding?**

**A:** You could but it is not recommended. There are pipeline more suitables for micro-organisms out there.

**Q: Can I use nucleotide databases other than RefSeq or BLAST NT (e.g. BOLD)?**

**A:** Any nucleotide database can be used, however they need to be formatted by the BLAST+ tools suit to be usable in FooDMe.
Consult the [BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK569841/) for instructions.

**Q: In the reports, some decimals are point-separated and others are comma-separated. Why is that and how do I fix it?**

**A:** The comma-spearated values are generated in some cases depending on your system's langugae (e.g. French, German,...).
To enforce point-separated values, change your system locale by modifying the `~/.bashrc` file with `LANG=C`. 
