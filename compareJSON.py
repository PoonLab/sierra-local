import json
import argparse

parser = argparse.ArgumentParser(description='Compare sierra JSON outputs')
parser.add_argument('files', nargs=2, type=str, help='sierrapy output first, sierralocal output second')
#args = parser.parse_args(['IN-sierra.json', 'IN-local.json'])
args = parser.parse_args()

def parse_results(j):
    '''
    Parse results from a json
    '''
    output = {}
    for sequence_results in j:
        sequence_dict = dict(sequence_results)
        output[sequence_dict['inputSequence']['header']] = []
        drugResistance = sequence_dict['drugResistance'][0]
        for drugScore in drugResistance['drugScores']:
            text = drugScore['text']
            score = drugScore['score']
            displayAbbr = drugScore['drug']['displayAbbr']
            partialScores = drugScore['partialScores']
            mutation = drugScore
            row = [text, score, displayAbbr, partialScores]
            output[sequence_dict['inputSequence']['header']].append(row)
    return output

def main():
    out = open('comparejson.txt','w')
    # Sierrapy goes first
    with open(args.files[0],'r') as spy:
        j = json.load(spy)
    # sierra-local goes second
    with open(args.files[1],'r') as slocal:
        k = json.load(slocal)

    headers = ['text','score','drug']

    #Parse SierraPy results
    sequence_sierrapy_results = parse_results(j)

    #Parse SierraLocal results
    sequence_sierralocal_results = parse_results(k)

    perfect = 0
    for header in sequence_sierrapy_results:
        # Sort by accession number
        sierrapy = sorted(sequence_sierrapy_results[header], key=lambda x : x[2])
        sierralocal = sorted(sequence_sierralocal_results[header], key=lambda x : x[2])

        py_to_local_count = 0
        py_to_local_total = 0
        local_to_py_count = 0
        local_to_py_total = 0

        #Compare sierrapy to sierralocal
        for index, drugresistance in enumerate(sierrapy):
            if drugresistance not in sierralocal:
                out.write(str(header) +' '+ str(drugresistance[2]) +' '+ str(drugresistance[1]) +' '+'not in sierralocal output'+'\n')
                if drugresistance[1] != sierralocal[index][1]:
                    out.write("Score mismatch\n")
                if drugresistance[3] != sierralocal[index][3]:
                    out.write("partialScore mismatch\n")
            else:
                py_to_local_count += 1
            py_to_local_total += 1

        #compare sierralocal to sierrapy
        for index, drugresistance in enumerate(sierralocal):
            if drugresistance not in sierrapy:
                out.write(str(header) +' '+ str(drugresistance[2]) +' '+ str(drugresistance[1]) +' '+'not in sierrapy output'+'\n')
                if drugresistance[1] != sierrapy[index][1]:
                    out.write("Score mismatch\n")
                if drugresistance[3] != sierrapy[index][3]:
                    out.write("partialScore mismatch\n")

            else:
                local_to_py_count += 1
            local_to_py_total += 1
        #print str((py_to_local_count+local_to_py_count)/float(py_to_local_total+local_to_py_total)*100)
        if py_to_local_count == py_to_local_total and local_to_py_count == local_to_py_total:
            perfect += 1
        '''
        print 'py_to_local_count', py_to_local_count
        print 'py_to_local_total', py_to_local_total
        print 'local_to_py_count', local_to_py_count
        print 'local_to_py_total', local_to_py_total
        '''
    out.write(str(perfect) + '/' + str(len(sequence_sierrapy_results))+'\n')
    out.write(str(float(perfect)*100/len(sequence_sierrapy_results))+'%')
    out.close()
if __name__ == '__main__':
    main()