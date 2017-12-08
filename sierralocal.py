import score_alg
from hivdb import HIVdb
import os
import argparse
import sys
from nucaminohook import NucAminoAligner
import jsonwriter
from seqUtils import parse_fasta

cwd = os.getcwd()

def parse_args():
    '''
    CLI argument parser. Current options include input FASTA files only
    :return: args object
    '''
    parser = argparse.ArgumentParser(description='Use the HIVdb algorithm for mutation-based resistance scoring of sequences.')
    parser.add_argument('fasta', nargs='+', type=str, help='List of input files.')
    parser.add_argument('-o', dest='outfile', default=None, type=str, help='Output filename.')
    #args = parser.parse_args(['testsequences.fasta'])
    args = parser.parse_args()
    return args

def gene_map():
    pol_start = 2085
    pol_nuc_map = {
        'PR':(2253,2549),
        'RT':(2550,3869),
        'IN':(4230,5096)
    }
    convert = lambda x:(x-pol_start)/3
    pol_aa_map = {}
    for key, val in pol_nuc_map.items():
        pol_aa_map[key] = (convert(val[0]), convert(val[1]))
    return pol_aa_map

def main():
    args = parse_args()
    path = cwd + '/HIVDB.xml'
    algorithm = HIVdb(path)
    definitions = algorithm.parse_definitions(algorithm.root)
    database = algorithm.parse_drugs(algorithm.root)
    comments = algorithm.parse_comments(algorithm.root)
    print args.outfile

    for file in args.fasta:
        names, scores, ordered_mutation_list, genes = scorefile(file, database)
        if args.outfile == None:
            outputname = file.replace('.fasta','-local.json')
        else:
            outputname = args.outfile
        jsonwriter.write_to_json(outputname, names, scores, genes, ordered_mutation_list)
    

def scorefile(file, database):
    '''
    Returns a set of corresponding names, scores, and ordered mutations for a given FASTA file containing pol sequences
    :param file: the FASTA file name containing arbitrary number of sequences and headers
    :param database: the HIVdb drug scores and notations
    :return: list of names, list of scores, list of ordered mutations
    '''
    aligner = NucAminoAligner(file)
    aligner.align_file()
    names, genes, muts = aligner.get_mutations(gene_map())
    ordered_mutation_list = []
    scores = []
    for index, query in enumerate(names):
        print "scoring",query
        #raw_sequence = parse_fasta(open(file,'r')).items()[index][1]
        ordered_mutation_list.append(sorted(zip(muts[index].keys(), [x[1] for x in muts[index].values()], [x[0] for x in muts[index].values()])))
        scores.append(score_alg.score_drugs(database, muts[index]))
    return names, scores, ordered_mutation_list, genes


if __name__ == '__main__':
    main()
