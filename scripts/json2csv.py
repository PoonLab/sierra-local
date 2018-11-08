import json
import csv
import argparse

def json2csv(infile, outfile, write_header=True):
    results = json.load(infile)
    writer = csv.writer(outfile)

    first_row = True
    for result in results:
        if first_row and write_header:
            # use this row to generate labels
            drugs = []
            for d in result['drugResistance']:
                for d2 in d['drugScores']:
                    drugs.append(d2['drug']['name'])

            writer.writerow(['name', 'subtype'] + drugs)
            first_row = False

        scores = []
        for d in result['drugResistance']:
            for d2 in d['drugScores']:
                scores.append(d2['score'])

        writer.writerow(
            [result['inputSequence']['header'], result['subtypeText']] + scores
        )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert JSON generated by running a Sierra processor "
                    "into a more accessible CSV format, retaining only "
                    "drug-specific resistance scores.")

    parser.add_argument('json', help='<input> Sierra JSON results file',
                        type=argparse.FileType('r'))
    parser.add_argument('csv', help='<output> CSV file',
                        type=argparse.FileType('w'))
    args = parser.parse_args()

    json2csv(infile=args.json, outfile=args.csv)