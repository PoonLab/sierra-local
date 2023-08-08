import os
import sys
import time
import argparse

from sierralocal import score_alg
from sierralocal.hivdb import HIVdb
from sierralocal.jsonwriter import JSONWriter
from sierralocal.nucaminohook import NucAminoAligner


def score(filename, xml_path=None, tsv_path=None, forceupdate=False, do_subtype=False, program='post'): # pragma: no cover
    """
    Functionality as a Python module. Can import this function from sierralocal.
    @param filename: str, Path to FASTA file containing sequences
    @param xml_path: str <optional>, Path to ASI2 XML file
    @param tsv_path: (optional) str, Path to tab-separated APOBEC DRM file
    @param forceupdate: bool, DEPRECATED. Uses Selenium to retrieve ASI2 and TSV files.
    @param do_subtype: bool, ???
    """
    algorithm = HIVdb(asi2=xml_path, apobec=tsv_path, forceupdate=forceupdate)
    time_start = time.time()

    sequence_headers, sequence_scores, ordered_mutation_list, file_genes, \
    sequence_lengths, file_trims, subtypes, na_sequence = scorefile(filename, algorithm, do_subtype)

    count = len(sequence_headers)

    print("{} sequences found in file {}.".format(len(sequence_headers), filename))
    output_file = os.path.splitext(filename)[0] + '-local.json'
    writer = JSONWriter(algorithm)
    writer.write_to_json(output_file, sequence_headers, sequence_scores, file_genes,
                         ordered_mutation_list, sequence_lengths, file_trims, subtypes, na_sequence)
    time_end = time.time()
    print("Time elapsed: {:{prec}} seconds ({:{prec}} it/s)".format(
        time_end - time_start, count/(time_end - time_start), prec='.5'))
    # cleanup is default action
    #os.remove(os.path.splitext(filename)[0] + '.tsv')


def scorefile(input_file, algorithm, do_subtype=False, program='post'):
    """
    Returns a set of corresponding names, scores, and ordered
    mutations for a given FASTA file containing pol sequences
    @param input_file: str, the FASTA file name containing arbitrary
    number of sequences and headers
    @param algorithm: sierralocal.hivdb.HIVdb, the HIVdb drug scores and notations
    @param do_subtype: bool <optional>, ???
    @return: list of names, list of scores, list of ordered mutations, list of NA sequence
    """
    aligner = NucAminoAligner(algorithm, program=program)
    result = aligner.align_file(input_file, program=program)

    print('Aligned ' + input_file)
    sequence_headers, file_genes, file_mutations, file_trims, subtypes = \
        aligner.get_mutations(result, do_subtype=do_subtype)

    ordered_mutation_list = []
    sequence_scores = []
    sequence_lengths = []
    na_sequence = {}

    for index, value in enumerate(result):
        na_sequence[value['Name']] = value['Sequence']

    # iteration over records in file
    for index, query in enumerate(sequence_headers):
        genes = file_genes[index]
        mutations = file_mutations[index]

        scores = []
        mutation_lists = []
        length_lists = []

        # iterate by gene
        for idx, gene_info in enumerate(genes):
            gene, first_aa, last_aa, first_na, last_na = gene_info

            length_lists.append(last_na - first_na + 1)

            # convert format
            mutation_lists.append(
                sorted(zip(
                    mutations[idx].keys(),  # position
                    [x[1] for x in mutations[idx].values()],  # aa
                    [x[0] for x in mutations[idx].values()],  # wt
                    [x[2] for x in mutations[idx].values()]   # text
                ))
            )
            scores.append(score_alg.score_drugs(algorithm,
                                                gene,
                                                mutations[idx]))

        ordered_mutation_list.append(mutation_lists)
        sequence_scores.append(scores)
        sequence_lengths.append(length_lists)

    return sequence_headers, sequence_scores, ordered_mutation_list, \
           file_genes, sequence_lengths, file_trims, subtypes, na_sequence

def sierralocal(fasta, outfile, xml=None, json=None, cleanup=False, forceupdate=False,
                program='post', do_subtype=False): # pragma: no cover
    """
    Contains all initializing and processing calls.

    @param fasta:  relative or absolute paths to FASTA file to process; multiple files may be
                   passed as a list object
    @param outfile:  file path to write JSON results
    @param xml: <optional> str, path to local copy of HIVdb algorithm XML file
    @param json: <optional> str, path to local copy of HIVdb algorithm APOBEC DRM file
    @param cleanup:  <optional> bool, to delete alignment file
    @param forceupdate: <optional> bool, forces sierralocal to update its local copy of the HIVdb algorithm
    @return: tuple, a tuple of (number of records processed, time elapsed initializing algorithm)
    """

    # initialize algorithm and jsonwriter
    time0 = time.time()
    algorithm = HIVdb(asi2=xml, apobec=json, forceupdate=forceupdate)
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
        sequence_headers, sequence_scores, ordered_mutation_list, file_genes, \
        sequence_lengths, file_trims, subtypes, na_sequence = scorefile(input_file, algorithm,
                                                                        program=program, do_subtype=do_subtype)

        count += len(sequence_headers)
        print("{} sequences found in file {}.".format(len(sequence_headers), input_file, na_sequence))

        # output results for the file
        if outfile == None:
            output_file = prefix + '_results.json'
        else:
            output_file = outfile

        writer.write_to_json(output_file, sequence_headers, sequence_scores,
                             file_genes, ordered_mutation_list, sequence_lengths,
                             file_trims, subtypes, na_sequence)

        if cleanup:
            # delete alignment file
            os.remove(prefix+'.tsv')

    return count, time_elapsed


def parse_args(): # pragma: no cover
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
                        help='<optional> Path to HIVdb ASI2 XML file')
    parser.add_argument('-json', default=None,
                        help='<optional> Path to JSON HIVdb APOBEC DRM file')
    parser.add_argument('--cleanup', action='store_true',
                        help='Deletes NucAmino alignment file after processing.')
    parser.add_argument('--forceupdate', action='store_true',
                        help='Forces update of HIVdb algorithm. Requires network connection.')
    parser.add_argument('-alignment', default='post', choices=['post', 'nuc'],
                        help='Alignment program to use, "post" for post align and "nuc" for nucamino')
    parser.add_argument('-subtype', action='store_true', default=False,
                        help='Subtype while running sierralocal'
                        )
    args = parser.parse_args()
    return args


def main(): # pragma: no cover
    """
    Main function called from CLI.
    """
    args = parse_args()

    # check that FASTA files in list all exist
    for file in args.fasta:
        if not os.path.exists(file):
            print("Error: there is no file {}".format(file))
            sys.exit()

    time_start = time.time()
    count, time_elapsed = sierralocal(args.fasta, args.outfile, xml=args.xml,
                                      json=args.json, cleanup=args.cleanup, forceupdate=args.forceupdate,
                                      program=args.alignment, do_subtype=args.subtype)
    time_diff = time.time() - time_start

    print("Time elapsed: {:{prec}} seconds ({:{prec}} it/s)".format(
        time_diff,
        count / (time_diff-time_elapsed),  # adjust for XML processing time
        prec='.5'
    ))


if __name__ == '__main__':
    main()
