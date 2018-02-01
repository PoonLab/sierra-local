import json
import argparse

parser = argparse.ArgumentParser(description='Compare sierra JSON outputs')
parser.add_argument('files', nargs=2, type=str, help='sierrapy output first,sierralocal output second')
#args = parser.parse_args(['IN-sierra.json', 'IN-local.json'])
args = parser.parse_args()

'''

This is a testing script for validating sierralocal against sierrapy

Instructions:
1. Generate sierrapy output
2. Generate sierralocal output using same data
3. Run this file using the two output files (in that order) as args

'''

def parse_results(j):
    '''
    Parse results from a json
    '''
    output = {}
    for sequence_results in j:
        sequence_dict = dict(sequence_results)
        output[sequence_dict['inputSequence']['header']] = {}
        drugResistance = sequence_dict['drugResistance'][0]
        for drugScore in drugResistance['drugScores']:
            output[sequence_dict['inputSequence']['header']][drugScore['drug']['name']] = {}
            output[sequence_dict['inputSequence']['header']][drugScore['drug']['name']]['text'] = drugScore['text']
            output[sequence_dict['inputSequence']['header']][drugScore['drug']['name']]['score'] = drugScore['score']
            output[sequence_dict['inputSequence']['header']][drugScore['drug']['name']]['partialScores'] = {}
            for index, pscore in enumerate(drugScore['partialScores']):
                output[sequence_dict['inputSequence']['header']][drugScore['drug']['name']]['partialScores'][str(index)] = {}
                output[sequence_dict['inputSequence']['header']][drugScore['drug']['name']]['partialScores'][str(index)]['score'] = pscore['score']
                output[sequence_dict['inputSequence']['header']][drugScore['drug']['name']]['partialScores'][str(index)]['mutations'] = {}
                for mutation in pscore['mutations']:
                    output[sequence_dict['inputSequence']['header']][drugScore['drug']['name']]['partialScores'][str(index)]['mutations'][mutation['text']] = mutation['comments'][0]['text']
    return output

def main():
    out = open('comparejson.txt','w')
    # Sierrapy goes first
    with open(args.files[0],'r') as spy:
        j = json.load(spy)
    # sierra-local goes second
    with open(args.files[1],'r') as slocal:
        k = json.load(slocal)

    #Parse SierraPy results into a dictionary
    dict1 = parse_results(j)

    #Parse SierraLocal results
    dict2 = parse_results(k)
    #print dict1
    #print dict2

    errors = 0
    totalchecks = 0
    totalmatches = 0

    for header in dict1:
        drug_errors = 0
        perfect = True
        if not dict2.has_key(header):
            errors += 1
            out.write(header+' '+"missing from sierralocal\n")
            perfect = False
            continue
        for drug in dict1[header]:
            totalchecks += 1
            if not dict2[header].has_key(drug):
                errors += 1
                out.write(drug+' '+"missing from sierralocal for"+' '+header+'\n')
                perfect = False
                continue
            if dict1[header][drug]['text'] != dict2[header][drug]['text']:
                perfect = False
                out.write("level mismatch:"+' '+drug+' '+dict1[header][drug]['text']+' '+ dict2[header][drug]['text']+' '+ header+'\n')
            if dict1[header][drug]['score'] != dict2[header][drug]['score']:
                perfect = False
                out.write("score mismatch:"+' '+drug+' '+str(dict1[header][drug]['score'])+' '+ str(dict2[header][drug]['score'])+' '+ header+'\n')
            if perfect:
                totalmatches += 1
            #for pscore in dict1[header][drug]['partialScores']:
            #    print dict1[header][drug]['partialScores'][pscore]
        if not perfect:
            errors += 1



    out.write("Errors: "+str(errors))
    out.write("\nCount: "+str(len(dict1)))
    out.write("\n"+str(float(len(dict1)-errors)/len(dict1)))
    out.write("\n"+str(totalmatches)+'/'+str(totalchecks))
    out.write("\n"+str(float(totalmatches)/totalchecks))
    out.close()
if __name__ == '__main__':
    main()