# 5utr ['suter']
Annotation resource/tool for 5'UTR variants

***This tool is in active development***

5utr is an [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) plugin aiming to provide different annotations relevant to 5' UTR (untranslated region) variants.

*Current limitations*
- Works only on single nucleotide substitutions (SNVs)
- Only GRCh37 reference build is implemented

## Usage

`./vep -i variants.vcf.gz  --cache --offline --assembly GRCh37 --dir_plugin /path/to/Sutr_dir --no_stats --plugin Sutr,Sutr.tsv.gz,ALL -o out.file.txt`

Instaed of 'ALL' you can specify columns from Sutr.tsv.gz to be included in the annotation.

## Annotations

To be updated. If interested please contact me at sander.pajusalu@yale.edu
