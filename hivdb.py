import xml.etree.ElementTree as xml
import re
import os

class HIVdb():
    def __init__(self, path):
        self.root = xml.parse(path).getroot()
        self.algname = self.root.find('ALGNAME').text
        self.version = self.root.find('ALGVERSION').text
        self.version_date = self.root.find('ALGDATE').text


    def parse_definitions(self, root):
        self.definitions = {
            'gene': {},  # gene target names and drug classes
            'drugclass': {},  # maps drug class to drugs
            'globalrange': {},  # maps from score to level
            'level': {},  # maps from level to S/I/R symbols
            'comment': {}  # maps comments to id string
        }

        for element in root.getchildren():
            if element.tag == 'GENE_DEFINITION':
                gene = element.find('NAME').text
                drug_classes = element.find('DRUGCLASSLIST').text.split(', ')
                self.definitions['gene'].update({gene: drug_classes})

            elif element.tag == 'DRUGCLASS':
                drugclass = element.find('NAME').text
                druglist = element.find('DRUGLIST').text.split(',')
                self.definitions['drugclass'].update({drugclass: druglist})

            elif element.tag == 'GLOBALRANGE':
                globalrange = element.find('GLOBALRANGE').text.split(',')
                for item in globalrange:
                    order = int(re.split('=>', item)[1].strip('() '))  # str containing order number: '1'
                    range = re.split('=>', item)[0].strip('() ')  # str containing the range: '-INF TO 9'
                    min = re.split('TO', range)[0].strip()  # str containing min val in range: '-INF'
                    max = re.split('TO', range)[1].strip()  # str containing max val in range: '9'
                    # convert_to_num converts a string to integer, and also 'INF' and '-INF' to infinite representations
                    def convert_to_num(s):
                        if s == '-INF':
                            return float('-inf')
                        elif s == 'INF':
                            return float('inf')
                        else:
                            return int(s)
                    min = convert_to_num(min)
                    max = convert_to_num(max)
                    self.definitions['globalrange'].update({order: [min, max]})  # TODO: convert to numerical

            elif element.tag == 'LEVEL_DEFINTIION':
                order = element.find('ORDER').text
                original = element.find('ORIGINAL').text
                sir = element.find('SIR').text
                self.definitions['level'].update({order: {original: sir}})

            elif element.tag == 'COMMENT_DEFINITIONS':
                id = element.find('COMMENT_STRING').text
                comment = element.find('TEXT').text
                # sort_tag = element.find('SORT_TAG').text  # TODO: is this necessary?? (it is always 1)
                self.definitions['comment'].update({id: comment})


    """ parse_drugs iterates through each drug in HIVDB, 
        parses condition for a specific drug, 
        and assigns a library of the drug resistant mutation conditions to the dictionary of drugs
        
        @param root: algorithm root
    """
    def parse_drugs(self, root):
        self.drugs = {}
        #count =0
        # can this be eliminated? May not need it
        for element in root.getchildren():
            if element.tag == 'DRUG':
                drug = element.find('NAME').text                            # name of the drug
                fullname = element.find('FULLNAME').text                    # full name of the drug
                condition = element.find('RULE').find('CONDITION').text     # drug conditions
                cond_dict = self.parse_condition(condition)                 # dictionary of parsed drug conditions
                #count += 1
                #if count == 10:
                    #print(drug, self.drms)
                self.drugs[drug] = self.drugs[fullname] = cond_dict         # is this needed? if self.drugs = {} isn't needed then scrap this


    """parse_condition function takes a given condition (one of four types)
        'MAXAND' condition: MAX ((41L AND 215ACDEILNSV) => 5, (41L AND 215FY) => 15)
        'MAX' condition: MAX ( 219E => 5, 219N => 5, 219Q => 5, 219R => 5 )
        'AND' condition: (67EGN AND 215FY AND 219ENQR) => 5
        'single-drm' condition: 62V => 5
        
    @param condition: given drug condition to parse
    """
    def parse_condition(self, condition):
        # drug resistant mutation (DRM)
        mutation_list = ((condition.split('(',1)[1].rstrip(')')) + ',').split('\n')    # NOTE: \n may not be stable or transferable on other platforms
        self.drms= []                            # library of drms for the drug of the given condition

        for drm in mutation_list:
            if drm.strip().startswith('MAX'):
                max_lib = []
                max_chunks = re.findall('\(?\d+[A-z]+[\s?AND\s?\d+\w]+\)?\s?=>\s?\d+', drm)
                iter = 0
                for chunk in max_chunks:
                    # for both MAX conditions, need to create a mini-library that will keep all of the individual DRMs together
                    self._parse_scores(max_lib, drm, chunk, iter)
                    iter += 1               # iter is created to keep track of which DRM associated with which index in the list of scores in _parse_scores function
                self.drms.append(max_lib)   # finally append this mini-library to the DRMs library

            else:
                iter = 0
                self._parse_scores(self.drms, drm, drm, iter)

        return self.drms


    """ _parse_scores function is a helper function to parse_condition.
        Parses the residues, amino acids, and scores associated with a particular DRM,
        then updates the specified dictionary library 
        
        @param drm_lib: given library to be updated with DRMs
        @param drm: full name of original drug resistant mutation                           
        @param chunk: drm_group of one of the condition types   (kinda vague...will fix)
        @param iter: index to keep track of which DRM is associated with respective index in the extracted list of scores
    """                                                                                                                     #TODO: combine drm and iter so less params
    def _parse_scores(self, drm_lib, drm, chunk, iter):
        scores = re.findall('[0-9]+(?=\W)', drm.strip())   # extract scores in same order as grouped drm tuples; stored in (indexable) list
        rANDr = re.findall('\d+[A-z]+[\s?AND\s?\d+\w]+', chunk)

        for combo_group in rANDr:
            mut_list = []
            residueAA = re.findall('[0-9]+[A-z]+', combo_group.strip())  # TODO: needs testing
            for mutation in residueAA:
                residue = int(re.findall('[\d]+(?!\d)(?=\w)', mutation)[0])
                aa = str(re.findall('[0-9]+([A-z]+)', mutation)[0])
                mut_list.append(tuple((residue, aa)))

            # populate the drms library with the new drm condition
            drm_lib.append({'group': mut_list, 'value': int(scores[iter])})
            # wipe out scores stored in variable scores for next batch
            scores[:] = []


    """
    score_drugs function first checks if the drug is in the HIVdb
    if found, calculates score with a given drug and sequence according to Stanford algorithm
    
    @param drugname: name of the drug you want the score for
    @param sequence: user provided sequence of type str (tolerates whitespace on either side, will strip it out later)
    @return score: calculated drm mutation score
    """
    def score_drugs(self, drugname, sequence):
        FOUND = False
        if drugname not in self.drugs.keys(): print("Drugname: " +  drugname + " not found.")
        else: FOUND = True
        if FOUND:
            score = 0

            for condition in self.drugs[drugname]:
                #print(condition)
                # list of potential scores
                candidates = [0]

                # if every item in condition['group'] is satisfied
                    # append to list of potential scores
                value = condition['value']
                #print(value)
                residueAAtuples = condition['group']
                                                                 # TODO: Currently works only for single drm and combo, not the MAX conditions yet (deal with the mini-libraries)
                #print(len(residueAAtuples))
                count = 0
                for i in residueAAtuples:
                    residue = i[0]
                    aa = i[1]
                    if sequence[residue-1] != aa: continue
                    else: count += 1
                    if count == len(residueAAtuples):
                        candidates.append(value)

                # take the max of what's in the list of potential scores and update total score
                # doesn't matter for the single drm or combo condition because they only have one associated value anyways
                score += max(candidates)

            return score








def main():
    alg = HIVdb("/home/tng92/git/sierra-local/HIVDB.xml")
    alg.parse_definitions(alg.root)
    alg.parse_drugs(alg.root)
    #print(alg.definitions)
    #print(alg.drugs.keys())
    alg.score_drugs("etravirine", 'DAAAAAGAAELGAAAATCTQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')




main()

