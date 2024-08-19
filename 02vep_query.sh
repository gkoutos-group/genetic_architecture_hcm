#Exectuted using VEP version 110
#Requires cached version of human genome VEP annotations 110 GRCh38
#Sources some annotations from VEP
mkdir vep_output

inp=$1

vep \
    --cache --dir_cache ~/.cache --offline \
    --dir_plugins /rds/projects/r/russdr-bb-data/scripts_progs/tools/plugins \
    --format vcf -i $inp \
    --show_ref_allele --sift p --polyphen p --nearest symbol --distance bp --clin_sig_allele 1 --hgvs \
    --fasta reference_files/hg38.fa \
    --fields Uploaded_variation,Location,Allele,REF_ALLELE,Consequence,NEAREST,HGVSc,HGVSp \
    --tab --force_overwrite --fork 20 -o stdout | grep -v '##' | gzip -c > $(echo $inp | cut -d'.' -f1)_vep.tsv.gz


