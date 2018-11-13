import subprocess
from sierralocal import main
import os
from pathlib import Path
import time
import random

# set True to run timing analysis
run_sierra_local = True
run_sierra_py = False

# correspond to the number of files produced by data_prep.py
gene_ranges = {
	'IN':range(0, 13),
	'PR':range(0, 106),
	'RT':range(0, 113)
}

sample_size = 10

#sample set of files
#can set a seed
random.seed(1)
samples = {}
for g in gene_ranges.keys():
	samples.update({g:random.sample(gene_ranges[g], sample_size)})

print(samples)

# sierralocal
if run_sierra_local:
	for g in samples.keys():
		for i in samples[g]:
			main.score(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'hivdb', 'hivdb-data', '{}-{}.fa'.format(g,i))))


#sierrapy
if run_sierra_py:
	# PR_samples = [42, 61, 97, 11, 12, 96, 105, 18, 2, 89]
	# RT_samples = [88, 38, 78, 10, 37, 76, 101, 60, 99, 17]
	# IN_samples = [0, 8, 3, 4, 1, 10, 6, 2, 7, 11] #duplicates in 8
	for g in samples.keys():
		for i in samples[g]:
			start = time.time()
			subprocess.run(["sierrapy", "fasta",
				os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'hivdb', 'hivdb-data', '{}-{}.fa'.format(g, i))),
				"-o",
				os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'hivdb', 'hivdb-data', '{}-{}.fa'.format(g, i)))])
			end = time.time()
			print(end - start)