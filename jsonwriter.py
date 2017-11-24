import json
import xml.etree.ElementTree as xml
import os
import re
import hashlib
from hivdb import HIVdb

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
#print definitions['comment']
database = algorithm.parse_drugs(algorithm.root)
comments = algorithm.parse_comments(algorithm.root)
#print comments

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
                if float(scores[drug][0]) == 0.0:
                    continue
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
                            mut = {}
                            mut['text'] = combination
                            mut['comments'] = [{
                                'text' : findComment(combination, comments, definitions['comment']),
                                'type' : drugclass
                            }]
                            mut['primaryType'] = drugclass
                            pscoredict['mutations'].append(mut)                    
                        drugScore['partialScores'].append(pscoredict)
                drugScores.append(drugScore)
    drugResistance['drugScores'] = drugScores
    return [drugResistance]

def alignedgenesequences():
    return None

def inputsequence():
    out = {
        'header' : '',
        'SHA512' : ''
    }
    return out

def write_to_json(scores, genes):
    data = {}
    data['subtypeText'] = 'NULL'
    data['validationResults'] = validationresults()
    data['drugResistance'] = drugresistance(scores, genes)
    data['alignedGeneSequences'] = alignedgenesequences()
    data['inputSequence'] = inputsequence()

    with open('test.json','w') as outfile:
        json.dump([data], outfile, indent=2)

def findComment(mutation, comments, details):
    # TODO: use the gene information to make this more accurate
    trunc_mut = re.findall(r'\d+\D+',mutation)[0]
    for gene, mutationdict in comments.items():
        for item in mutationdict.keys():
            if trunc_mut in item:
                full_mut = mutationdict[item]
                if details.has_key(full_mut):
                    return details[full_mut]['1']