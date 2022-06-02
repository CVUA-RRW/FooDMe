# Frequently Asked Questions

**Q: Do I need to install mamba?**

**A:** No you dont´t. Simply using conda should be enough, although considerably slower, to 
create environments. If using only conda you may need to add the argument `--conda-frontend conda`
when executing snakemake.

**Q: Can I use FooDMe for single end reads?**

**A:** No. If that is something you need, [get in touch](https://github.com/CVUA-RRW/FooDMe/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=), it may be possible to implement.

**Q: Are IonTorrent data supported?**

**A:** Not yet. But it should be possible with some modifications. [Get in touch](https://github.com/CVUA-RRW/FooDMe/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=) if this is something you need.

**Q: Can I use FooDMe for bacteria metabarcoding?**

**A:** You could but it is not recommended. There are pipeline more suitables for micro-organisms out there.

**Q: In the reports, some decimals are point-separated and others are comma-separated. Why is that and how do I fix it?**

**A:** The comma-spearated values are generated in some cases depending on your system's langugae (e.g. French, German,...).
To enforce point-separated values, change your system locale by modifying the `~/.bashrc` file with `LANG=C`. 
