from gotoh2.aligner import Aligner
import unittest
import seqUtils
import score_alg
from hivdb import HIVdb
import os
import argparse
import sys

g2 = Aligner()

#argv parsing
parser = argparse.ArgumentParser(description='Score a sequence.')
parser.add_argument('-i', nargs='+', dest='inputfiles', type=str, help='List of input files.')
#parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), help='List of output files. Defaults to same names as input files.', required=False)
args = parser.parse_args()

def align(reference, query):
    """
    'Procrustean' alignment of query nucleotide sequence to reference sequence

    :param query: Nucleotide sequence as a string
    :param reference: Reference sequence as a string
    :return: Aligned sequence fitted wrt relative insertions to reference
    """

    # Use gotoh2 to align the query to reference sequence
    query = query.strip()
    assert type(query) is str, 'Query must be a string.'
    aligned_query, aligned_reference, aligned_score = g2.align(
        query, reference)
    fitted_query = ''
    ''' 
    TODO:   flag indels not in sets of three
                exceptions for indels in proximity summing to 3n
    '''
    # Remove insertions in the sequence relative to reference sequence
    for idx, char in enumerate(aligned_reference):
        # Iff the query does not have insertion relative to reference, add
        # the nucleotide/relative deletion to the output string
        if char != '-':
            fitted_query += aligned_query[idx]

    return fitted_query


def translate(seq, offset=0, resolve=False, return_list=False):
    """
    Translate nucleotide sequence into amino acid sequence.
    Synonymous nucleotide mixtures are resolved to the corresponding residue.
    Nonsynonymous nucleotide mixtures are encoded with '?' 

    :param seq: nucleotide sequence to translate
    :param offset: offset by X shifts sequence to the right by X bases
    :param resolve: Boolean to resolve ambiguous AA for purposes of aligning AA
    :param return_list: Boolean to return as a list

    :return: string (default) or list of the translated amino acid sequence
    """
    aa_seq = seqUtils.translate_nuc(seq, offset, resolve, return_list)
    return aa_seq

def import_fasta(filename):
    '''
    Returns a dictionary from FASTA file.

    :param filename: string of the filename
    :return: dictionary of header:sequence pairs
    '''
    with open(filename, 'r' ) as file:
        fasta = seqUtils.parse_fasta(file)
    return fasta

def main():
    #Some of Tammy's code
    path = os.getcwd() + '/HIVDB.xml'
    algorithm = HIVdb(path)
    definitions = algorithm.parse_definitions(algorithm.root)
    # print(definitions['gene'])
    # print(definitions['level'])
    # print(definitions['drugclass'])
    # print(definitions['globalrange'])
    # print(definitions['comment'])
    database = algorithm.parse_drugs(algorithm.root)
    # print(database)
    comments = algorithm.parse_comments(algorithm.root)
    # print(comments)

    # Parsing and processing the nucleotide sequence file
    # Scoring each sequence
    # Output to a JSON
    for file in args.inputfiles:
        fasta = import_fasta(file)
        nuc_seq = fasta.values()[0]
        reference = nuc_seq #placeholder reference
        aligned_nuc_seq = align(reference, nuc_seq)
        aa_seq = translate(aligned_nuc_seq)
        scores = score_alg.score_drugs(database, aa_seq)
        print(scores)

        with open(file.split('.fasta')[0]+'.out','w') as outfile:
            outfile.write(str(scores))


class TestAlignments(unittest.TestCase):
    def test(self):
        self.assertEqual(align('TAGACTTACCAC', 'ACTTAGCAT'), '-A--CTTAGCAT')
        self.assertEqual(align('ACTTAGCAT', 'TAGACTTACCAC'), 'ACTTACCAC')

    # Check out this error later
    # def test_empty(self):
    #    self.assertEqual(align('-----------','ACGTACGTACGT'),'')

    def test_deletions(self):
        self.assertEqual(
            align('AGTACGCTCGTAGCAT', 'AGTACTAGCAT'), 'AGTAC-----TAGCAT')

    def test_insertions(self):
        self.assertEqual(
            align('AGTACTAGCAT', 'ACTACGCTCGTAGCAT'), 'ACTACTAGCAT')


if __name__ == '__main__':
    # unittest.main()
    main()
