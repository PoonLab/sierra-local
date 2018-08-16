import pandas as pd
import argparse
from pathlib import Path
import os

'''
HOW TO USE THIS DATA PREP SCRIPT:
1. Create /hivdb/hivd-data/ folder located in the highest-level sierralocal directory.
2. From HIVdb, get the PR, RT, IN genotype-treatment correlation datasets (as text files) and save as .tsv files as PR.tsv etc.
3. Run this script.
'''

hivdb_data_path = Path(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'hivdb', 'hivdb-data')))

#specify a gene
parser = argparse.ArgumentParser(description='Create test FASTA files.')
parser.add_argument('-gene', type=str, nargs='+', help='List one or more genes: IN, PR, RT')
parser.add_argument('-size', nargs='?', default=1000, const=1000, type=int, help='Batch size.')
args = parser.parse_args()

# Maximum number of sequences per FASTA file. sierrapy MAY be limiting this to 2500 per query
block_size = args.size

for g in args.gene:
	# read the gene's TSV file from HIVdb
	df = pd.read_csv(hivdb_data_path/'{}.tsv'.format(g), sep='\t')

	# iterate over the file, in BLOCKS
	sum_count = 0
	for k,j in enumerate(range(0, len(df), block_size)):
		count = 0
		#write each block to a unique FASTA file
		with open(hivdb_data_path/'{}-{}.fa'.format(g, k), 'w') as output:
			for i in range(min(block_size, len(df)-(k*block_size))):
				if df.iloc[k*block_size + i].loc['NASeq'] == df.iloc[k*block_size + i].loc['NASeq']: # skips NA sequences
					header = str(df.iloc[k*block_size + i].loc['AccessionID'])+'.'+str(df.iloc[k*block_size + i].loc['IsolateName'])+'.'+str(df.iloc[k*block_size + i].loc['Subtype']+'.'+str(count))
					output.write('>'+header+'\n'+str(df.iloc[k*block_size + i].loc['NASeq'])+'\n')
					count += 1
		sum_count += count
		print("{} sequences written to file {}-{}.fa".format(count, g, k))
	print("{} sequences written in total for gene {}.".format(sum_count, g))