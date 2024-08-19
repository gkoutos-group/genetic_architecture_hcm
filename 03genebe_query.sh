#script to run genebe using python script wrapper and join together with vep output

source genebe/bin/activate

inp=$1
ivout=$2

python 03genebe.py $inp

Rscript 03vepiv.R ${ivout}_annotations_genebe.tsv.gz ${ivout}_vep.tsv.gz

