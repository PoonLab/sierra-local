import os
import sys
import time
import argparse
import json
from pathlib import Path
import csv

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
    sequence_lengths, file_trims, subtypes, na_sequence, ambiguous, names = scorefile(filename, algorithm, do_subtype)

    count = len(sequence_headers)

    print("{} sequences found in file {}.".format(len(sequence_headers), filename))
    output_file = os.path.splitext(filename)[0] + '-local.json'
    writer = JSONWriter(algorithm)
    writer.write_to_json(output_file, sequence_headers, sequence_scores, file_genes, ordered_mutation_list,
                        sequence_lengths, file_trims, subtypes, na_sequence, ambiguous, names)
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

    # hold all NNN sequences and void them when scoring drugs
    ambiguous = {}  # {sequence name: {gene : set(positions of NNN)}}
    gene_order = {}


    print('Aligned ' + input_file)
    sequence_headers, file_genes, file_mutations, file_trims, subtypes = \
        aligner.get_mutations(result, do_subtype=do_subtype)


    for ind, gene in enumerate(file_genes):
        seq_n = sequence_headers[ind]
        gene_order[seq_n] = []
        
        for ind2, protein in enumerate(gene):
            gene_order[seq_n].append(protein[0])

    for ind, sequence in enumerate(file_mutations):
        seq_n = sequence_headers[ind]
        ambiguous[seq_n] = {}
        for ind2, muts in enumerate(sequence):
            ambiguous[seq_n][gene_order[seq_n][ind2]] = set()
            for position, AA in muts.items():
                # >4 results in X AA, which sierrapy seems to ignore for SDRMs
                if len(AA[1]) > 4:
                    ambiguous[seq_n][gene_order[seq_n][ind2]].add(position)

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
           file_genes, sequence_lengths, file_trims, subtypes, na_sequence, ambiguous, gene_order

def sierralocal(fasta, outfile, xml=None, json=None, cleanup=False, forceupdate=False,
                apobec_csv=None, unusual_csv=None, sdrms_csv=None, mutation_csv=None,
                updater_outdir=None, program='post', do_subtype=False): # pragma: no cover
    """
    Contains all initializing and processing calls.

    @param fasta:  relative or absolute paths to FASTA file to process; multiple files may be
                   passed as a list object
    @param outfile:  file path to write JSON results
    @param xml: <optional> str, path to local copy of HIVdb algorithm XML file
    @param json: <optional> str, path to local copy of HIVdb algorithm APOBEC DRM file
    @param cleanup:  <optional> bool, to delete alignment file
    @param forceupdate: <optional> bool, forces sierralocal to update its local copy of the HIVdb algorithm
    @param apobec_csv: str <optional>, Path to CSV APOBEC csv file (default: apobecs.csv)
    @param unusual_csv: str <optional>, Path to CSV file to determine if is unusual (default: rx-all_subtype-all.csv)
    @param sdrms_csv: str <optional>, Path to CSV file to determine SDRM mutations (default: sdrms_hiv1.csv)
    @param mutation_csv: str <optional>, Path to CSV file to determine mutation type (default: mutation-type-pairs_hiv1.csv)
    @return: tuple, a tuple of (number of records processed, time elapsed initializing algorithm)
    """

    # initialize algorithm and jsonwriter
    time0 = time.time()
    algorithm = HIVdb(asi2=xml, apobec=json, forceupdate=forceupdate, updater_outdir=updater_outdir)
    writer = JSONWriter(algorithm, apobec_csv, unusual_csv, sdrms_csv, mutation_csv)
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
        sequence_lengths, file_trims, subtypes, na_sequence, ambiguous, names = scorefile(input_file, algorithm,
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
                             file_trims, subtypes, na_sequence, ambiguous, names)

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
                        help='<optional> Path to HIVdb ASI2 XML file (default: HIVDB_9.4.xml)')
    parser.add_argument('-json', default=None,
                        help='<optional> Path to JSON HIVdb APOBEC DRM file')
    parser.add_argument('--cleanup', action='store_true',
                        help='Deletes NucAmino alignment file after processing.')
    parser.add_argument('--forceupdate', action='store_true',
                        help='Forces update of HIVdb algorithm. Requires network connection.')
    parser.add_argument('-alignment', default='post', choices=['post', 'nuc'],
                        help='Alignment program to use, "post" for post align and "nuc" for nucamino')
    parser.add_argument('-apobec_csv', default=None,
                        help='<optional> Path to CSV APOBEC csv file (default: apobecs.csv)')
    parser.add_argument('-unusual_csv', default=None,
                        help='<optional> Path to CSV file to determine if is unusual (default: rx-all_subtype-all.csv)')
    parser.add_argument('-sdrms_csv', default=None,
                        help='<optional> Path to CSV file to determine SDRM mutations (default: sdrms_hiv1.csv)')
    parser.add_argument('-mutation_csv', default=None,
                        help='<optional> Path to CSV file to determine mutation type (default: mutation-type-pairs_hiv1.csv)')
    parser.add_argument('-updater_outdir', default=None,
                        help='<optional> Path to folder to store updated files from updater (default: sierralocal/data folder))')

    args = parser.parse_args()
    return args


def check_input(apobec_path, unusual_path, sdrms_path, mutation_path):
    """
    Check if the input for the files are valid based on the first row of the csv.

    apobec_path: path to apobec_drms.csv
    unusual_path: path to rx-all_subtype-all.csv
    sdrms_path: path to sdrms_hiv1.csv
    mutation_path: path to mutation-type-pairs_hiv1.csv
    """
    exp = {
        "apobec_csv": ["gene", "position", "aa"],
        "unusual_csv": ["gene", "position", "aa", "percent", "count", "total", "reason", "isUnusual"],
        "sdrms_csv": ["drug_class", "gene", "position", "aa"],
        "mutation_csv": ["strain", "gene", "drugClass", "position", "aas", "mutationType", "isUnusual"],
    }

    paths = {
        "apobec_csv": apobec_path,
        "unusual_csv": unusual_path,
        "sdrms_csv": sdrms_path,
        "mutation_csv": mutation_path,
    }
    for key, path in paths.items():
        if path is None:
            continue
        try:
            with open(path, newline="", encoding="utf-8-sig") as f:
                reader = csv.reader(f)
                header = next(reader)
        except Exception as e:
            sys.exit(f"Could not open {key} file '{path}': {e}")

        if header != exp[key]:
            print(
                f"Invalid header in {key} file '{path}'.\n"
                f"Expected: {exp[key]}\nFound:    {header}"
            )
            sys.exit()


def main(): # pragma: no cover
    """
    Main function called from CLI.
    """
    args = parse_args()

    # check for valid file inputs
    check_input(args.apobec_csv, args.unusual_csv, args.sdrms_csv, args.mutation_csv)

    mod_path = Path(os.path.dirname(__file__))

    if args.updater_outdir:
        target_dir = args.updater_outdir
    else:
        target_dir = os.path.join(mod_path, "data")
    
    # Create directory if it doesn't exist
    os.makedirs(target_dir, exist_ok=True)

    # check that FASTA files in list all exist
    for file in args.fasta:
        if not os.path.exists(file):
            print("Error: there is no file {}".format(file))
            sys.exit()

    time_start = time.time()
    count, time_elapsed = sierralocal(args.fasta, args.outfile, xml=args.xml,
                                      json=args.json, cleanup=args.cleanup, forceupdate=args.forceupdate,
                                      apobec_csv=args.apobec_csv, unusual_csv=args.unusual_csv, 
                                      sdrms_csv=args.sdrms_csv, mutation_csv=args.mutation_csv, updater_outdir=target_dir,
                                      program=args.alignment)
    time_diff = time.time() - time_start

    print("Time elapsed: {:{prec}} seconds ({:{prec}} it/s)".format(
        time_diff,
        count / (time_diff-time_elapsed),  # adjust for XML processing time
        prec='.5'
    ))


if __name__ == '__main__':
    main()
