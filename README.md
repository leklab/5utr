# 5utr ['suter']
Annotation resource/tool for 5'UTR variants

***This tool is in active development***

5utr is an [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) plugin aiming to provide different annotations relevant to 5' UTR (untranslated region) variants.

*Current limitations*
- Works only on single nucleotide substitutions (SNVs)
- Is limited to canonical transcripts and 50 bp range from the main ATG (start codon)
- Only GRCh37 reference build is implemented

## Usage

`./vep -i variants.vcf.gz  --cache --offline --assembly GRCh37 --dir_plugin /path/to/Sutr_dir --no_stats --plugin Sutr,deltas_vep.tsv.gz -o out.file.txt`

## Annotations

### delta_TE

Change in predicted translational efficiency from the [Sample et al 2019 publication](https://doi.org/10.1101/310375). The values were calculated using the [published model](https://github.com/pjsample/human_5utr_modeling) for all possible SNVs in 50 bp range as well as the respective reference sequences without any scaling. Then the delta score was calculated, by subtracting the variant value from reference value. The lower scores are indicative for lower translational efficiency for variant if compared to the reference.

### delta_dsRNA & delta_rG4

Change in minimum free energy (MFE) for variant if compared to the reference. Calculated as described in [Murat et al. Genome Biology 2018 19:229](https://doi.org/10.1186/s13059-018-1602-2) using ViennaRNA package using the following commands:

dsRNA: `RNAfold -g whole_utrs_allsnvs50bp.txt > whole_utrs_allsnvs50bp_G4.txt`

G4: `RNAfold whole_utrs_allsnvs50bp.txt > whole_utrs_allsnvs50bp_dsRNA.txt`

rG4 scores were calculated then by subtracting G4 results from dsRNA (rG4 = G4 - dsRNA).

Then the delta scores were calculated by subtracting the variant values from the respective reference values. The lower delta scores (negative) would mean that MFE is decreased by variant, and the secondary structure is more stable thus indicative for lower translational activity, as demonstrated in Murat et al. publication. 
