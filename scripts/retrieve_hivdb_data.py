"""
Use urllib2 module to stream Stanford HIVdb published data set and
filter a limited number of records to a local text file.
Use iteration so we don't necessarily have to download the entire file.

This accomplishes the same task as scripts/data_prep.py but I just wanted
to do it my way :-)
"""

from urllib import request
from csv import DictReader
import codecs
import argparse

__author__ = '@ArtPoon'

parser = argparse.ArgumentParser()
parser.add_argument('gene', choices=['PR', 'RT', 'INT'],
                  help='Target gene (PR, RT or INT).')
parser.add_argument('outfile', type=argparse.FileType('w'),
                    help='Destination file for sequences.')
parser.add_argument('-size', type=int, default=100,
                    help='Number of sequences to retrieve.')
parser.add_argument('-step', type=int, default=10,
                    help='Interval size between output sequences.')

args = parser.parse_args()


# Stanford University HIV Drug Resistance Database, Genotype-Rx Datasets
base_url = 'https://hivdb.stanford.edu/modules/lookUpFiles/geno-rx-datasets'

handle = request.urlopen('{}/{}.txt'.format(base_url, args.gene))
rows = DictReader(codecs.iterdecode(handle, 'utf-8'), delimiter='\t')

noutput = 0
for ln, row in enumerate(rows):
    accno = row['AccessionID']
    if accno == '':
        # skip records without a GenBank accession number (no sequence)
        continue

    if ln % args.step == 0:
        args.outfile.write('>{}.{}.{}.{}\n{}\n'.format(
            accno,
            row['IsolateName'],
            row['Subtype'],
            noutput,
            row['NASeq']
        ))
        noutput += 1

    if noutput == args.size:
        # exit condition
        break

