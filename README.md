# 5utr ['suter']
Annotation resource/tool for 5'UTR variants

***This tool is in active development***

5utr is an [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) plugin aiming to provide different annotations relevant to 5' UTR (untranslated region) variants.

*Current limitations*
- Works only on single nucleotide substitutions (SNVs)
- Some annotations are available only for the nearest 100bp to CDS.

## Usage

`./vep -i variants.vcf.gz  --cache --offline --assembly GRCh37 --dir_plugin /path/to/Sutr_dir --no_stats --plugin Sutr,/path/to/Sutr_resource_table_b37.tsv.gz,ALL -o out.file.txt`

`./vep -i variants.vcf.gz  --cache --offline --assembly GRCh38 --dir_plugin /path/to/Sutr_dir --no_stats --plugin Sutr,/path/to/Sutr_resource_table_hg38.tsv.gz,ALL -o out.file.txt`

Instaed of 'ALL' you can specify columns from Sutr.tsv.gz to be included in the annotation. See VEP manual for general usage instructions (e.g. using and downloading cache, additional annotations like gnomAD and ClinVar etc.)

## Annotations

The annotations files can be accessed [here](https://owncloud.ut.ee/owncloud/index.php/s/2o4FWdRWtxAsjMM), be sure to download right genome assembly matching with your vcf.
