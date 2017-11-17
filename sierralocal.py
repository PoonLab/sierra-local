import score_alg
from hivdb import HIVdb
import os
import argparse
import sys
from nucaminohook import NucAminoAligner

cwd = os.getcwd()

def parse_args():
    parser = argparse.ArgumentParser(description='Use the HIVdb algorithm for mutation-based resistance scoring of sequences.')
    parser.add_argument('fasta', nargs='+', type=str, help='List of input files.')
    args = parser.parse_args(['testsequences.fasta'])
    return args

def main():
    args = parse_args()
    print args
    path = cwd + '/HIVDB.xml'
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
    scorefile(args, database)

def scorefile(args, database):
    # Scoring each sequence
    for file in args.fasta:
        aligner = NucAminoAligner(file)
        aligner.align_file()
        names, res = aligner.get_mutations()

        for idx, mutationlist in enumerate(res):
            print names[idx]
            scores = score_alg.score_drugs(database, mutationlist)
            print(scores)

            #TODO: output to a formatted JSON consistent with the sierra-client JSON output

            #with open(file.split('.fasta')[0]+'.out','w') as outfile:
            #    outfile.write(str(scores))

if __name__ == '__main__':
    main()
