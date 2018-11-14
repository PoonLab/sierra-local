"""
A script for processing FASTA files containing sequences from Stanford's
HIVdb database using multiple cores.  This was used for comparing outputs
to the SierraPy web service client.
"""

from sierralocal.main import sierralocal
from glob import glob
import os
from multiprocessing import Pool
import sys



# settings for integrase (IN)
#xml = 'sierralocal/data/HIVDB_8.6.1.ca26d72f.xml'
#version = 'v8.6.1'

def process(file):
    outfile = file.replace('.fa', '-rs3-'+version+'.json')
    if os.path.exists(outfile):
        sys.stdout.write('skipping '+file+'\n')
        return
    sierralocal(fasta=file, outfile=outfile, xml=xml)

if __name__ == '__main__':
    # assumes we're calling this script from the repository root
    files = glob('hivdb/hivdb-data/PR*.fa')
    xml = 'sierralocal/data/HIVDB_8.5.d926dfff.xml'
    version = 'v8.5'
    with Pool(processes=11) as pool:
        pool.map(process, files)

    files = glob('hivdb/hivdb-data/RT*.fa')
    with Pool(processes=11) as pool:
        pool.map(process, files)

    files = glob('hivdb/hivdb-data/IN*.fa')
    xml = 'sierralocal/data/HIVDB_8.6.1.ca26d72f.xml'
    version = 'v8.6.1'
    with Pool(processes=11) as pool:
        pool.map(process, files)
