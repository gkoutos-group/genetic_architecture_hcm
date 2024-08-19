cohort_name=$1
input_ids_file=$2

#carry out some qc, ensure alignment and convert vcf into minor allele count matrix
#qc - 10.1186/1471-2105-15-125
#select only genes of interest including only those with 90% having depth greater than 10 and not missing genomtype
#set any that fit the criteria of having Depth <= 8 or Genotype Quality <= 10 as missing
#normalise data and check against reference
#take only alleles that appear a minimum of 1 time
#update chromosome names to remove 'chr' and set ID column
#convert to a table with position columns and GT for all individuals
#a series of sed commands to translate diploid alt allele counts into total alt allele counts e.g. 1/0 = 1
bcftools view -R reference_files/gene_list_both.txt -S ${input_ids_file} -c 1 ${cohort_name}.vcf.gz |\
    bcftools +setGT - -- -t q -n . -i 'FORMAT/DP <= 8 | FORMAT/GQ <= 10' |\
    bcftools norm -m -any --fasta-ref reference_files/hg38.fa - |\
    bcftools view --min-ac 1 - |\
    bcftools annotate --rename-chrs reference_files/chr_names.txt --set-id '%CHROM\_%POS\_%REF\_%ALT' - |\
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' -H - |\
    sed 's/\.\/\./NA/g' | \
    sed 's/0\/0/0/g' | \
    sed 's/\.\/0/0/g' | \
    sed 's/0\/\./0/g' | \
    sed 's/1\/0/1/g' | \
    sed 's/0\/1/1/g' | \
    sed 's/1\/\./1/g' | \
    sed 's/\.\/1/1/g' | \
    sed 's/1\/1/2/g' | \
    sed 's/\.|\./NA/g' | \
    sed 's/0|0/0/g' | \
    sed 's/\.|0/0/g' | \
    sed 's/0|\./0/g' | \
    sed 's/1|0/1/g' | \
    sed 's/0|1/1/g' | \
    sed 's/1|\./1/g' | \
    sed 's/\.|1/1/g' | \
    sed 's/1|1/2/g' | \
    gzip -c > ${cohort_name}.tsv.gz


