from sierralocal import score_alg
from sierralocal.hivdb import HIVdb
import os
import argparse
from sierralocal.nucaminohook import NucAminoAligner
from sierralocal.jsonwriter import JSONWriter
from pathlib import Path
import time

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
    parser.add_argument('-xml', default=str(Path('.')/'data'/'HIVDB.xml'), 
                        help='Path to HIVDB algorithm XML file, which can be downloaded using the provided script update_HIVDB.py')
    parser.add_argument('-skipalign', action='store_true', help='Skip NucAmino alignment if TSV file already present.')
    parser.add_argument('-cleanup', action='store_true', help='Deletes NucAmino alignment file after processing.')
    args = parser.parse_args()
    return args

def score(filename):
    """
    Functionality as a Python module. Can import this function from sierralocal
    """
    assert os.path.isfile(filename), "Provided file not found."
    # path = Path('./sierralocal/data/HIVDB_8.4.d926dfff.xml')
    algorithm = HIVdb()
    definitions = algorithm.parse_definitions(algorithm.root)
    database = algorithm.parse_drugs(algorithm.root)
    comments = algorithm.parse_comments(algorithm.root)
    time_start = time.time()
    count = 0
    sequence_headers, sequence_scores, ordered_mutation_list, file_genes, sequence_lengths, file_firstlastNA, file_trims = scorefile(filename, database, False)
    count += len(sequence_headers)
    print("{} sequences found in file {}.".format(len(sequence_headers), filename))
    output_file = os.path.splitext(filename)[0] + '-local.json'
    writer = JSONWriter()
    writer.write_to_json(output_file, sequence_headers, sequence_scores, file_genes, ordered_mutation_list, sequence_lengths, file_firstlastNA, file_trims)
    time_end = time.time()
    print("Time elapsed: {:{prec}} seconds ({:{prec}} it/s)".format(time_end - time_start, count/(time_end - time_start), prec='.5'))
    os.remove(os.path.splitext(filename)[0] + '.tsv')


def main():
    args = parse_args()
    path = args.xml  #'./data/HIVDB.xml'
    # algorithm = HIVdb(path)
    algorithm = HIVdb()
    definitions = algorithm.parse_definitions(algorithm.root)
    database = algorithm.parse_drugs(algorithm.root)
    comments = algorithm.parse_comments(algorithm.root)
    time_start = time.time()
    count = 0
    for input_file in args.fasta:
        sequence_headers, sequence_scores, ordered_mutation_list, file_genes, sequence_lengths, file_firstlastNA, file_trims = scorefile(input_file, database, args.skipalign)
        count += len(sequence_headers)
        print("{} sequences found in file {}.".format(len(sequence_headers), input_file))
        if args.outfile == None:
            output_file = os.path.splitext(input_file)[0] + '-local.json'
        else:
            output_file = args.outfile
        writer = JSONWriter()
        writer.write_to_json(output_file, sequence_headers, sequence_scores, file_genes, ordered_mutation_list, sequence_lengths, file_firstlastNA, file_trims)
    time_end = time.time()
    print("Time elapsed: {:{prec}} seconds ({:{prec}} it/s)".format(time_end - time_start, count/(time_end - time_start), prec='.5'))
    if args.cleanup:
        os.remove(os.path.splitext(input_file)[0] + '.tsv')


def scorefile(input_file, database, skipalign):
    '''
    Returns a set of corresponding names, scores, and ordered mutations for a given FASTA file containing pol sequences
    :param file: the FASTA file name containing arbitrary number of sequences and headers
    :param database: the HIVdb drug scores and notations
    :return: list of names, list of scores, list of ordered mutations
    '''
    aligner = NucAminoAligner()
    if not skipalign:
        aligner.align_file(input_file)
    sequence_headers, file_genes, file_mutations, file_firstlastNA, file_trims = aligner.get_mutations(aligner.gene_map(), input_file)
    ordered_mutation_list = []
    sequence_scores = []
    sequence_lengths = []

    with open(input_file, 'r') as fastafile:
        sequence_list = get_input_sequences(fastafile)
        sequence_lengths = [len(s.replace('N', '')) / 3 for s in sequence_list]

    for index, query in enumerate(sequence_headers):
        ordered_mutation_list.append(sorted(zip(file_mutations[index].keys(), [x[1] for x in file_mutations[index].values()], [x[0] for x in file_mutations[index].values()])))
        sequence_scores.append(score_alg.score_drugs(database, file_mutations[index]))
    return sequence_headers, sequence_scores, ordered_mutation_list, file_genes, sequence_lengths, file_firstlastNA, file_trims

def get_input_sequences(handle):
    """
    Parse open file as FASTA, return dictionary of
    headers and sequences as key-value pairs.
    """
    sequences = []
    sequence = ''
    for i in handle:
        if i[0] == '>' or i[0] == '#':
            if len(sequence) > 0:
                sequences.append(sequence)
                sequence = ''
        else:
            sequence += i.strip('\n').upper()
    sequences.append(sequence)
    return sequences

if __name__ == '__main__':
    main()
