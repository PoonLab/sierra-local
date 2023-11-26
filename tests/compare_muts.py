import argparse
import json

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process JSON files and generate drug resistance scores.')

    parser.add_argument('--input1', required=True, help='Path to the first input JSON file')
    parser.add_argument('--input2', required=True, help='Path to the second input JSON file')
    parser.add_argument('--output1', required=True, help='Path to the first output JSON file')
    parser.add_argument('--output2', required=True, help='Path to the second output JSON file')

    return parser.parse_args()

def process_json(input_path):
    with open(input_path, 'r') as file:
        data = json.load(file)

    result_dict = {}
    for item in data:
        name = item['inputSequence']['header']
        result_dict[name] = {}

        for genes in item['drugResistance']:
            gene = genes['gene']['name']
            result_dict[name][gene] = {}

            for drugs in genes['drugScores']:
                drug = drugs['drug']['name']
                result_dict[name][gene][drug] = float(drugs['score'])

    return result_dict

def main():
    args = parse_arguments()

    # Process the first input file
    dic_py = process_json(args.input1)

    # Process the second input file
    dic_local = process_json(args.input2)

    # Write results to the first output file
    with open(args.output1, 'w') as file:
        json.dump(dic_py, file, indent=4)

    # Write results to the second output file
    with open(args.output2, 'w') as file:
        json.dump(dic_local, file, indent=4)

if __name__ == "__main__":
    main()