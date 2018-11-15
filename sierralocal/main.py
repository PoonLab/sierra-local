from sierralocal import score_alg
from sierralocal.hivdb import HIVdb
import os
import argparse
from sierralocal.nucaminohook import NucAminoAligner
from sierralocal.jsonwriter import JSONWriter
from sierralocal.utils import get_input_sequences
from pathlib import Path
import time
import sys


def score(filename, xml_path=None, forceupdate=False):
    """
    Functionality as a Python module. Can import this function from sierralocal
    """

    algorithm = HIVdb(xml_path, forceupdate)
    definitions = algorithm.parse_definitions(algorithm.root)
    database = algorithm.parse_drugs(algorithm.root)
    comments = algorithm.parse_comments(algorithm.root)

    time_start = time.time()

    sequence_headers, sequence_scores, ordered_mutation_list, file_genes, sequence_lengths, \
    file_trims, subtypes = scorefile(filename, algorithm)

    count = len(sequence_headers)

    print("{} sequences found in file {}.".format(len(sequence_headers), filename))
    output_file = os.path.splitext(filename)[0] + '-local.json'
    writer = JSONWriter(algorithm)
    writer.write_to_json(output_file, sequence_headers, sequence_scores, file_genes,
                         ordered_mutation_list, sequence_lengths, file_trims, subtypes)
    time_end = time.time()
    print("Time elapsed: {:{prec}} seconds ({:{prec}} it/s)".format(time_end - time_start, count/(time_end - time_start), prec='.5'))
    # cleanup is default action
    #os.remove(os.path.splitext(filename)[0] + '.tsv')


def scorefile(input_file, algorithm):
    '''
    Returns a set of corresponding names, scores, and ordered mutations for a given FASTA file containing pol sequences
    :param file: the FASTA file name containing arbitrary number of sequences and headers
    :param database: the HIVdb drug scores and notations
    :return: list of names, list of scores, list of ordered mutations
    '''
    aligner = NucAminoAligner()
    result = aligner.align_file(input_file)

    print('Aligned '+input_file)
    sequence_headers, file_genes, file_mutations, file_trims, subtypes = \
        aligner.get_mutations(result)

    ordered_mutation_list = []
    sequence_scores = []
    sequence_lengths = []

    # iteration over records in file
    for index, query in enumerate(sequence_headers):
        genes = file_genes[index]
        mutations = file_mutations[index]

        scores = []
        mutation_lists = []
        length_lists = []

        # iterate by gene
        for idx, gene_info in enumerate(genes):
            gene, firstAA, lastAA, firstNA, lastNA = gene_info

            length_lists.append(lastNA - firstNA + 1)

            # convert format
            mutation_lists.append(
                sorted(zip(
                    mutations[idx].keys(),  # position
                    [x[1] for x in mutations[idx].values()],  # aa
                    [x[0] for x in mutations[idx].values()]   # wt
                ))
            )
            scores.append(score_alg.score_drugs(algorithm, gene, mutations[idx]))

        ordered_mutation_list.append(mutation_lists)
        sequence_scores.append(scores)
        sequence_lengths.append(length_lists)

    return sequence_headers, sequence_scores, ordered_mutation_list, file_genes, \
           sequence_lengths, file_trims, subtypes


def sierralocal(fasta, outfile, xml=None, cleanup=False, forceupdate=False):
    """
    Contains all initializing and processing calls.

    :param fasta:  relative or absolute paths to FASTA file to process; multiple files may be
                   passed as a list object
    :param outfile:  file path to write JSON results
    :param xml:  <optional> path to local copy of HIVdb algorithm XML file
    :param skipalign:  <optional> to save time, skip NucAmino alignment step (reuse TSV output)
    :param forceupdate:  <optional> forces sierralocal to update its local copy of the HIVdb algorithm

    :return:  a tuple of (number of records processed, time elapsed initializing algorithm)
    """

    # initialize algorithm and jsonwriter
    time0 = time.time()
    algorithm = HIVdb(path=xml, forceupdate=forceupdate)
    definitions = algorithm.parse_definitions(algorithm.root)
    database = algorithm.parse_drugs(algorithm.root)
    comments = algorithm.parse_comments(algorithm.root)
    writer = JSONWriter(algorithm)
    time_elapsed = time.time() - time0

    # accommodate single file path argument
    if type(fasta) is str:
        fasta = [fasta]

    # begin processing
    count = 0
    for input_file in fasta:
        prefix = os.path.splitext(input_file)[0]

        # process and score file
        sequence_headers, sequence_scores, ordered_mutation_list, file_genes, sequence_lengths, \
        file_trims, subtypes = scorefile(input_file, algorithm)

        count += len(sequence_headers)
        print("{} sequences found in file {}.".format(len(sequence_headers), input_file))

        # output results for the file
        if outfile == None:
            output_file = prefix+'_results.json'
        else:
            output_file = outfile

        writer.write_to_json(output_file, sequence_headers, sequence_scores, file_genes, ordered_mutation_list,
                             sequence_lengths, file_trims, subtypes)

        if cleanup:
            # delete alignment file
            os.remove(prefix+'.tsv')

    return count, time_elapsed


def parse_args():
    """
    CLI argument parser. Current options include input FASTA files only
    :return: args object
    """
    parser = argparse.ArgumentParser(
        description='Local execution of Stanford HIVdb algorithm for mutation-based resistance scoring of sequences.'
    )
    parser.add_argument('fasta', nargs='+', type=str, help='List of input files.')
    parser.add_argument('-o', dest='outfile', default=None, type=str, help='Output filename.')
    parser.add_argument('-xml', default=None,
                        help='Path to HIVDB algorithm XML file, which can be downloaded using the provided script updater.py')
    parser.add_argument('--cleanup', action='store_true',
                        help='Deletes NucAmino alignment file after processing.')
    parser.add_argument('--forceupdate', action='store_true',
                        help='Forces update of HIVdb algorithm. Requires network connection.')
    args = parser.parse_args()
    return args


def main():
    """
    Main function called from CLI.
    """
    args = parse_args()

    # check that the FASTA file exists
    if not os.path.exists(args.fasta):
        print("Error: there is no file {}".format(args.fasta))
        sys.exit()

    time_start = time.time()
    count, time_elapsed = sierralocal(args.fasta, args.outfile, args.xml,
                                      cleanup=args.cleanup, forceupdate=args.forceupdate)
    time_diff = time.time() - time_start

    print("Time elapsed: {:{prec}} seconds ({:{prec}} it/s)".format(
        time_diff,
        count / (time_diff-time_elapsed),  # adjust for XML processing time
        prec='.5'
    ))


if __name__ == '__main__':
    main()
