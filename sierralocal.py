import score_alg
from hivdb import HIVdb
import os
import argparse
import sys
from nucaminohook import NucAminoAligner
import jsonwriter

cwd = os.getcwd()

def parse_args():
    parser = argparse.ArgumentParser(description='Use the HIVdb algorithm for mutation-based resistance scoring of sequences.')
    parser.add_argument('fasta', nargs='+', type=str, help='List of input files.')
    args = parser.parse_args(['KU127836.fasta'])
    return args

def main():
    args = parse_args()
    path = cwd + '/HIVDB.xml'
    algorithm = HIVdb(path)
    definitions = algorithm.parse_definitions(algorithm.root)
    database = algorithm.parse_drugs(algorithm.root)
    comments = algorithm.parse_comments(algorithm.root)
    scores = scorefile(args, database)
    #print scores, mutationlist
    jsonwriter.write_to_json(scores, ['RT'])
    

def scorefile(args, database):
    # Scoring each file
    for file in args.fasta:
        aligner = NucAminoAligner(file)
        aligner.align_file()
        names, genes, muts = aligner.get_mutations()

        for idx, mutationlist in enumerate(muts):
            #print mutationlist
            scores = score_alg.score_drugs(database, mutationlist)
            return scores


if __name__ == '__main__':
    main()
