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


    def parse_drugs(self, root):
        self.drugs = {}                                                     # can this be eliminated? May not need it
        for element in root.getchildren():
            if element.tag == 'DRUG':
                drug = element.find('NAME').text                            # name of the drug
                fullname = element.find('FULLNAME').text                    # full name of the drug
                condition = element.find('RULE').find('CONDITION').text     # drug conditions
                cond_dict = self.parse_condition(condition)                 # dictionary of parsed drug conditions
                #print(condition, '\n\n')
                self.drugs[drug] = self.drugs[fullname] = cond_dict         # is this needed? if self.drugs = {} isn't needed then scrap this


    def parse_condition(self, condition):
        # drug resistant mutation (DRM)
        mutation_list = condition.split('(',1)[1].rstrip(')').split('\n')       # NOTE: \n may not be stable or transferable on other platforms
        self.drms= []                            # list of DRM groups of 1 or more DRM tuples indicating residue and aa(s), as well as corresponding value
        for drm in mutation_list:
            drm_tuples = []
            drm_group = drm.strip().rstrip(',')

            # will apply to all max conditions and mixed conditions b/c when we score we're going to 'take the max' of every condition regardless
            # applies to 'AND' condition now too! yay
            # applies to a 'MAX'-only condition or a single drm condition
            # in a 'MAX' condition, creates a dictionary of a list of drm_group tuples and a list of score values
            # if a single drm condition, creates a dictionary of one drm_group tuple and one score value
            rANDr = re.findall('([0-9]+[A-z]+(\s?AND\s?[0-9]+[A-z]+)*)', drm_group)     # match one or more residueAA in combination(s) with a single associated score
            for combo_group in rANDr:
                residueAA = re.findall('[0-9]+[A-z]+', combo_group)             # TODO: this section has been refactored but needs testing
                for mutation in residueAA:
                    residue = re.findall('[\d]+(?!\d)(?=\w)', mutation)
                    aa = re.findall('[0-9]+([A-z]+)', mutation)
                    drm_tuples.append((residue, aa))

            scores = re.findall('[0-9]+(?=\W)', drm_group)            # extracts scores in same order as grouped drm tuples; stored in an (indexable) list
            self.parse_indiv_cond(self.drms, drm_tuples, scores)

        return self.drms

    ###
        # parse_indiv_cond function takes the list associated with the algorithm
        # populates the library with a given drug resistant mutation condition
        # @param self.drms; list library we are populating with each sub-condition
        # @param drm_tuples; given set of one or more separated (residue, aa) tuples
        # @params scores; score associated with each tuple/ combination of tuples
    ###
    def parse_indiv_cond(self, drms_lib, drm_tuples, score):                            # TODO: new function is created here, can possibly scrap the rest of the parse functions; needs testing
        iter = 0
        for mutation in drm_tuples:
            drms_lib.append({'group': mutation, 'value': score[iter]})
            iter += 1


    ###
    # score_drugs function first checks if the drug is in the HIVdb
    # if found, calculates score with a given drug and sequence according to Stanford algorithm
    # @param drugname; name of the drug you want the score for
    # @param sequence; user provided sequence of type str (tolerates whitespace on either side, will strip it out later)
    # @return score; calculated drm mutation score
    ###
    def score_drugs(self, drugname, sequence):
        FOUND = False
        if drugname not in self.drugs.keys(): print("Drugname: " +  drugname + " not found.")
        else: FOUND = True
        if FOUND:
            #calculate scores
            for i in self.drms:
                
                #don't forget to take the max of everything every time; doesn't matter for the single drm or combination condition because they only have one associated value anyways







def main():
    alg = HIVdb("/home/tng92/git/sierra-local/HIVDB.xml")
    alg.parse_definitions(alg.root)
    alg.parse_drugs(alg.root)
    #print(alg.definitions)
    #print(alg.drugs.keys())
    print(alg.score_drugs("etravirine", 'DAAAAAGAAELGAAAATCTQAAA'))


main()

