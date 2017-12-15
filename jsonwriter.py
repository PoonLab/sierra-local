import json
import xml.etree.ElementTree as xml
import os
import re
import hashlib
from hivdb import HIVdb
import csv

names = {}
names['3TC'] = 'LMV'
cwd = os.getcwd()
path = cwd + '/HIVDB.xml'
algorithm = HIVdb(path)
root = xml.parse(path).getroot()
version = root.find('ALGVERSION').text
version_date = root.find('ALGDATE').text
definitions = algorithm.parse_definitions(algorithm.root)
levels = definitions['level']
globalrange = definitions['globalrange']
database = algorithm.parse_drugs(algorithm.root)
comments = algorithm.parse_comments(algorithm.root)
with open('apobec-drms.5b7e1215.tsv','r') as csvfile:
    ApobecDRMs = list(csv.reader(csvfile, delimiter='\t'))
with open('INSTI-comments.csv','r') as INSTI_file:
    INSTI_comments = dict(csv.reader(INSTI_file,delimiter='\t'))
with open('PI-comments.csv','r') as PI_file:
    PI_comments = dict(csv.reader(PI_file,delimiter='\t'))
with open('RT-comments.csv','r') as RT_file:
    RT_comments = dict(csv.reader(RT_file,delimiter='\t'))


def validationresults():
    validationResults = {}
    validationResults['message'] = None
    validationResults['level'] = None
    return [validationResults]

def drugresistance(scores,genes):
    drugResistance = {}
    drugResistance['version'] = {}
    drugResistance['version']['text'] = version
    drugResistance['version']['publishDate'] = version_date
    drugScores = []
    # TODO: find a way to restrict drug classes, e.g. only NRTI/NNRTI for RT gene sequence...
    for gene in genes:
        for drugclass in definitions['gene'][gene]:
            classlist = definitions['drugclass'][drugclass]
            for drug in classlist:
                # Only record data if the score is non-zero
                #if float(scores[drug][0]) == 0.0:
                #    continue
                drugScore = {}

                #Infer resistance level text from the score and globalrange
                resistancelevel = -1
                for key in globalrange:
                    maximum = float(globalrange[key]['max'])
                    minimum = float(globalrange[key]['min'])
                    if minimum <= scores[drug][0] <= maximum:
                        resistancelevel = str(key)
                        break
                resistance_text = levels[resistancelevel]

                drugScore['text'] = resistance_text.keys()[0] #resistance level
                drugScore['drugClass'] = {'name':drugclass} #e.g. NRTI
                drugScore['score'] = float(scores[drug][0]) # score for this paritcular drug
                drugScore['drug'] = {}
                drugScore['drug']['displayAbbr'] = drug
                if names.has_key(drug):
                    drugScore['drug']['name'] = names[drug]
                else:
                    drugScore['drug']['name'] = drug
                drugScore['partialScores'] = []
                # create partial score, for each mutation, datastructure
                for index,pscore in enumerate(scores[drug][1]):
                    pscore = float(pscore)
                    if not pscore == 0.0:
                        pscoredict = {}
                        pscoredict['score'] = pscore
                        pscoredict['mutations'] = []
                        for combination in scores[drug][2][index]:
                            # find the mutation classification "type" based on the gene
                            type_ = drugclass
                            pos = re.findall(u'([0-9]+)',combination)[0]
                            muts = re.findall(u'(?<=[0-9])([A-Za-z])+',combination)[0]
                            if gene == 'IN':
                                for key in INSTI_comments:
                                    if pos in key and muts in key:
                                        type_ = INSTI_comments[key]
                                        break
                            elif gene == 'PR':
                                for key in PI_comments:
                                    if pos in key and muts in key:
                                        type_ = PI_comments[key]
                                        break
                            elif gene == 'RT':
                                for key in PI_comments:
                                    if pos in key and muts in key:
                                        type_ = RT_comments[key]
                                        break
                            # do the other stuff
                            mut = {}
                            mut['text'] = combination
                            mut['comments'] = [{
                                'text' : findComment(gene, combination, comments, definitions['comment']),
                                'type' : type_
                            }]
                            mut['primaryType'] = type_
                            pscoredict['mutations'].append(mut)                    
                        drugScore['partialScores'].append(pscoredict)
                drugScores.append(drugScore)
    drugResistance['drugScores'] = drugScores
    drugResistance['gene'] = {'name':genes[0]}
    return [drugResistance]

def alignedgenesequences(ordered_mutation_list, genes):
    dic = {}
    dic['prettyPairwise'] = {}
    mutationline = []
    prev = 0
    for tuple in ordered_mutation_list:
        mutationline = mutationline + [' - ']*(tuple[0]-1-prev) + [str(tuple[1]).center(3)]
        prev = tuple[0]
    dic['prettyPairwise']['mutationLine'] = mutationline
    dic['prettyPairwise']['alignedNAsLine'] = []
    dic['prettyPairwise']['refAALine'] = []
    dic['prettyPairwise']['positionLine'] = []
    dic['lastAA'] = None
    dic['firstAA'] = None
    dic['mutations'] = []
    for a in ordered_mutation_list:
        mutdict = {}
        mutdict['isInsertion'] = None
        mutdict['isDeletion'] = None
        mutdict['consensus'] = a[2]
        mutdict['AAs'] = a[1]
        mutdict['isApobecDRM'] = isApobecDRM(genes[0], a[2], a[0], a[1])
        mutdict['position'] = a[0]
        dic['mutations'].append(mutdict)
    return [dic]

def inputsequence(name):
    out = {
        'header' : name,
        'SHA512' : ''
    }
    return out

def write_to_json(filename, names, scores, genes, ordered_mutation_list):
    '''
    The main function to write passed result to a JSON file
    :param filename: the filename to write the JSON to
    :param names: list of sequence headers
    :param scores: list of sequence scores
    :param genes: list of genes in pol
    :param ordered_mutation_list: ordered list of mutations in the query sequence relative to reference
    '''
    out = []
    for index, score in enumerate(scores):
        data = {}
        data['subtypeText'] = 'NULL'
        data['validationResults'] = validationresults()
        data['drugResistance'] = drugresistance(score, genes[index])
        data['alignedGeneSequences'] = alignedgenesequences(ordered_mutation_list[index], genes[index])
        data['inputSequence'] = inputsequence(names[index])
        out.append(data)

    with open(filename,'w') as outfile:
        json.dump(out, outfile, indent=2)
        print "Writing JSON to file {}".format(filename)

def findComment(gene, mutation, comments, details):
    trunc_mut = re.findall(r'\d+\D',mutation)[0] #163K
    pos = re.findall(u'([0-9]+)',trunc_mut)[0]
    muts = re.findall(u'(?<=[0-9])([A-Za-z])+',trunc_mut)[0]
    for g, mutationdict in comments.items():
        for item in mutationdict.keys():
            if pos in item and muts in item:
                full_mut = mutationdict[item]
                if details.has_key(full_mut) and g == gene:
                    return details[full_mut]['1']
                else:
                    print full_mut
                    print pos, muts, trunc_mut 

def isApobecDRM(gene, consensus, position, AA):
    for row in ApobecDRMs[1:]:
        if [gene, consensus, position, AA] == row:
            return True
    return False