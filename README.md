# 5utr ['suter']
Annotation resource/tool for 5'UTR variants

Contact: Sander Pajusalu sander.pajusalu@yale.edu

***This tool is in active development***

5utr is an [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) plugin aiming to provide different annotations relevant to 5' UTR (untranslated region) variants.

*Current limitations*
- Works only on single nucleotide substitutions (SNVs)
- Some annotations are available only for the nearest 100bp to CDS.

## Usage

`./vep -i variants.vcf.gz  --cache --offline --assembly GRCh37 --dir_plugin /path/to/Sutr_dir --no_stats --plugin Sutr,/path/to/Sutr_resource_table_b37.tsv.gz,ALL -o out.file.txt`

`./vep -i variants.vcf.gz  --cache --offline --assembly GRCh38 --dir_plugin /path/to/Sutr_dir --no_stats --plugin Sutr,/path/to/Sutr_resource_table_hg38.tsv.gz,ALL -o out.file.txt`

Instaed of 'ALL' you can specify columns from Sutr.tsv.gz to be included in the annotation. See VEP manual for general usage instructions (e.g. using and downloading cache, additional annotations like gnomAD and ClinVar etc.)

The table is in nested format, meaning that if the same position is covered by multiple transcripts, the respective fields are comma-separated lists, where the values correspond to the respective transcripts marked in "transcript" field in the same order.

## Annotations

The annotations files can be accessed [here](https://owncloud.ut.ee/owncloud/index.php/s/2o4FWdRWtxAsjMM), be sure to download right genome assembly matching with your vcf.

The following annotations are available:

**chr, pos, REF, ALT** - locus and allele

**transcript**

**distance** - distance from main start codon of CDS. -1 marks the nearest base in 5'UTR, +1 is the first base of CDS (start codon).

**AUG_event** - whether the variant creates a novel AUG (*AUG_gain*), removes the AUG in reference sequence (*AUG_loss*) or changes the position of AUG at the same locus (*AUG_change*)

**var_frame, AUGgain_stop, var_kozak** - these fields are relevant for AUG_gain and AUG_change (for the created AUG) - **var_frame** specifies if it is *in-frame*  or *out-of-frame* with canonical start codon (CDS); **AUGgain_stop** indicates whether there is a stop before canonical start codon; **var_kozak** specifies the strength of kozak surrounding the novel AUG (*weak, moderate, strong*)
**ref_frame, AUGloss_stop, ref_kozak** - the same as previous, but now these define properties for AUG in the reference sequence which was removed by the variant.

**STOP_event, STOPgain_start, STOPloss_start** - whether the variant creates a novel STOP codon (*STOP_gain*), removes the STOP in reference sequence (*STOP_loss*) or changes the position of STOP at the same locus (*STOP_change*). **STOP_gain_start** and **STOPloss_start** indicate whether there is a start codon somewhere upstream of the created/removed STOP, but located still in the same 5'UTR. Only the ones *with_start* should be biologically relevant.

**kozak_strength, kozak_cat** - **kozak_strength** specifiec whether the variant makes kozak sequence stronger or weaker or remains the same. **kozak_cat** higlights whether the kozak change affects main AUG (canonical start codon) or upstream AUG, and whether the affected nucleotide is at relative position -3 (m3) or +4 (p4).

**TE_log2fold** - change in predicted translational efficiency caused by the variant in log2fold scale. 0 indicates no change, -1 indicates 2-fold reduction of translation and +1 2-fold increase in translation. The predictions are calculated by the method and model published by [Sample et al. in Nature Biotech 2019](https://www.nature.com/articles/s41587-019-0164-5)

**Delta_dsRNA, Delta_G4** - change in minimum free energy (secondary structure) affecting double-stranded RNA or G4 quadruplex structures respectively. Negative values show that the RNA secondary structure has become more stable due to the variant, which is associated with reduced translation. The calculations were done for the nearest 100bps in 5'UTRs and follow the methods published by [Murat et al, in Genome Biology 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1602-2)
