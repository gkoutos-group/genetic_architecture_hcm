#python wrapper script to run genebe for all positions
#input user details in variables below to use 

import genebe as gnb
import pandas as pd
import argparse

user_name=<INSERT USER NAME>
api_key_user=<INSERT USER API KEY>


parser = argparse.ArgumentParser()
parser.add_argument("file")
args = parser.parse_args()

file_out = file_out.replace('.tsv.gz','_annotations_genebe.tsv.gz')

dat = pd.read_csv(args.file, sep = '\t')

vars = dat['[3]ID'].tolist()

input_variants =  list()

for str in vars:
   new_string =  str.replace('_','-')
   input_variants.append(new_string)

df=gnb.annotate_variants_list_to_dataframe(input_variants,flatten_consequences=True, genome='hg38', username = user_name, api_key = api_key_user,use_netrc=False, endpoint_url='https://api.genebe.net/cloud/api-public/v1/variants', use_ensembl=True, use_refseq=True)

df = df[df['pos'].notna()]

df['pos'] = df['pos'].astype('int')

cols = ['chr','pos','ref','alt','gene_symbol','acmg_classification','acmg_criteria']

df = df[cols]

df.rename({'chr':'CHROM','pos':'POS','ref':'REF','alt':'ALT','gene_symbol':'Ref.Gene','acmg_classification':'ACMG','acmg_criteria':'ACMG_criteria'}, axis=1, inplace=True)

df.to_csv(file_out, sep = '\t', index = False)
