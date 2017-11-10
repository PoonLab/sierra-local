from gotoh2.aligner import Aligner
import unittest
import seqUtils
import score_alg
from hivdb import HIVdb
import os
import argparse
import sys
import reference
import hxb2

g2 = Aligner()

def align(reference, query):
    """
    'Procrustean' alignment of query nucleotide sequence to HXB2 reference sequence

    :param query: Nucleotide sequence as a string
    :param reference: Reference sequence as a string
    :return: Aligned sequence fitted wrt relative insertions to reference
    """

    # Use gotoh2 to align the query to reference sequence
    query = query.strip()
    assert type(query) is str, 'Query must be a string.'

    #Gets the best alignment based on the alignments to three HXB2 genes
    maxscore = 0
    for gene in [hxb2.integrase,hxb2.pro,hxb2.rt]:
        temp = g2.align(query, gene)
        if temp[2] > maxscore:
            aligned_query, aligned_reference, aligned_score = temp
            maxscore = temp[2]

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
    # CLI arg parsing
    parser = argparse.ArgumentParser(
        description='Use the HIVdb algorithm for mutation-based resistance scoring of sequences.'
    )
    parser.add_argument('-i', '--i', nargs='+', dest='inputfiles', type=str, help='List of input files.')

    args = parser.parse_args(['-i','testsequences.fasta'])
    #Some of Tammy's code
    path = os.getcwd() + '/HIVDB.xml'
    algorithm = HIVdb(path)
    definitions = algorithm.parse_definitions(algorithm.root)
    #print(definitions['gene'])
    #print(definitions['level'])
    #print(definitions['drugclass'])
    #print(definitions['globalrange'])
    # print(definitions['comment'])
    database = algorithm.parse_drugs(algorithm.root)
    #print(database)
    comments = algorithm.parse_comments(algorithm.root)
    #print(comments)

    # Parsing and processing the nucleotide sequence file
    # Scoring each sequence
    for file in args.inputfiles:
        fastadict = import_fasta(file)
        for header,fasta in fastadict.items():
            nuc_seq = fasta
            aligned_nuc_seq = align(hxb2.integrase, nuc_seq)
            #print aligned_nuc_seq
            aa_seq = translate(aligned_nuc_seq)
            print header,aa_seq
            scores = score_alg.score_drugs(database, aa_seq)
            print(scores)

            #TODO: output to a formatted JSON consistent with the sierra-client JSON output

            #with open(file.split('.fasta')[0]+'.out','w') as outfile:
            #    outfile.write(str(scores))

if __name__ == '__main__':
    main()
