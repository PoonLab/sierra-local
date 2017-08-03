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
        self.drugs = {}
        for element in root.getchildren():
            if element.tag == 'DRUG':
                drug = element.find('NAME').text
                fullname = element.find('FULLNAME').text
                rule = element.find('RULE')
                condition = rule.find('CONDITION').text
                condition = self._partition_scores(condition)
                self.drugs[drug] = self.drugs[fullname] = condition


    def _partition_scores(self, condition):
        # drug resistant mutation (DRM)
        mutation_list = condition.lstrip('SCORE FROM(').rstrip(')').split('\n')
        self.drm_scores = {
            'single_drm_dict': {},
            'max_dict': {},
            'combo_dict': {}
        }
        for drm in mutation_list:
            single_drm = drm.strip().rstrip(',')
            if single_drm.startswith('MAX'):
                possibilities = re.findall(r'[\S]+[\s]*[=][>][\s][\d]+', single_drm)
                # TODO: deal with the case when they incorporate boolean 'AND' into the 'MAX' scorelist of scoreitems
                for aa in possibilities:
                    mutation = aa.split('=>')
                    drm = mutation[0].strip()
                    score = int(mutation[1].strip())
                    self.drm_scores['max_dict'].update({drm: score})
            elif single_drm.find('AND') != -1:
                combinations = re.split('(| AND |) => ', single_drm)
                for aa in combinations[:-1]:
                    drm = aa
                    score = int(combinations[len(combinations) - 1])
                    self.drm_scores['combo_dict'].update({drm: score})
            else:
                score_cond = single_drm.split()
                drm = score_cond[0].strip()
                score = int(score_cond[2].strip())
                self.drm_scores['single_drm_dict'].update({drm: score})
        return self.drm_scores


    def score_drugs(self, drugname):
        FOUND = False
        if drugname not in self.drugs.keys():
            print("Drugname: " +  drugname + " not found.")
        else:
            #calculating score
            FOUND = True



def main():
    alg = HIVdb("/home/tng92/git/sierra-local/HIVDB.xml")
    alg.parse_definitions(alg.root)
    alg.parse_drugs(alg.root)
    #print(alg.definitions)
    #print(alg.drugs.keys())
    #print(alg.drugs[('ETR', 'etravirine')])
    alg.score_drugs("ETR")
    alg.score_drugs("etravirine")


main()

