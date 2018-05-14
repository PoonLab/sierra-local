import score_alg
from hivdb import HIVdb
import os
import argparse
import sys
from nucaminohook import NucAminoAligner
import jsonwriter
import re
import string

cwd = os.getcwd()

def parse_args():
    '''
    CLI argument parser. Current options include input FASTA files only
    :return: args object
    '''
    parser = argparse.ArgumentParser(
        description='Local execution of Stanford HIVdb algorithm for mutation-based resistance scoring of sequences.'
    )
    parser.add_argument('fasta', nargs='+', type=str, help='List of input files.')
    parser.add_argument('-o', dest='outfile', default=None, type=str, help='Output filename.')
    parser.add_argument('-xml', default='./data/HIVDB.xml', 
                        help='Path to HIVDB algorithm XML file, which can be downloaded using the provided script update_HIVDB.py')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    path = args.xml  #'./data/HIVDB.xml'
    algorithm = HIVdb(path)
    definitions = algorithm.parse_definitions(algorithm.root)
    database = algorithm.parse_drugs(algorithm.root)
    comments = algorithm.parse_comments(algorithm.root)
    print('writing output to {}'.format(args.outfile))

    for file in args.fasta:
        names, scores, ordered_mutation_list, genes, sequence_lengths = scorefile(file, database)
        if args.outfile == None:
            outputname = os.path.splitext(file)[0] + '-local.json'
        else:
            outputname = args.outfile
        jsonwriter.write_to_json(outputname, names, scores, genes, ordered_mutation_list, sequence_lengths)

def parse_fasta(handle):
    """
    Parse open file as FASTA, return dictionary of
    headers and sequences as key-value pairs.
    """
    names = []
    sequences = []
    sequence = ''
    for i in handle:
        if i[0] == '$':  # skip h info
            continue
        elif i[0] == '>' or i[0] == '#':
            if len(sequence) > 0:
                names.append(h)
                sequences.append(sequence)
                sequence = ''
                h = i.strip('\n')[1:]
            else:
                h = i.strip('\n')[1:]
        else:
            sequence += i.strip('\n').upper()
    names.append(h)
    sequences.append(sequence)
    return names, sequences


def scorefile(file, database):
    '''
    Returns a set of corresponding names, scores, and ordered mutations for a given FASTA file containing pol sequences
    :param file: the FASTA file name containing arbitrary number of sequences and headers
    :param database: the HIVdb drug scores and notations
    :return: list of names, list of scores, list of ordered mutations
    '''
    aligner = NucAminoAligner(file)
    aligner.align_file()
    names, genes, muts = aligner.get_mutations(aligner.gene_map())
    ordered_mutation_list = []
    scores = []
    sequence_lengths = []

    with open(file, 'r') as fastafile:
        headers, sequence_list = parse_fasta(fastafile)
        sequence_lengths = [len(s.replace('N', '')) / 3 for s in sequence_list]
    with open('sequence_out.txt','w') as out:
        out.write('\n'.join([str(sequence_lengths[i]) + " " + names[i] + " " + str(genes[i]) + sequence_list[i][1:30] +'...'+sequence_list[i][-15:-1] for i in range(1, len(sequence_lengths))]))

    if len(names) != len(muts):
        print('help')

    for index, query in enumerate(names):
        #print("scoring",query)
        ordered_mutation_list.append(sorted(zip(muts[index].keys(), [x[1] for x in muts[index].values()], [x[0] for x in muts[index].values()])))
        scores.append(score_alg.score_drugs(database, muts[index]))
    return names, scores, ordered_mutation_list, genes, sequence_lengths


if __name__ == '__main__':
    main()
