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
                condition = self.partition_scores(condition)
                print(condition, '\n\n')
                self.drugs[drug] = self.drugs[fullname] = condition


    def partition_scores(self, condition):
        # drug resistant mutation (DRM)
        mutation_list = condition.lstrip('SCORE FROM(').rstrip(')').split('\n')
        drm_scores = {
            'single_drm_dict': {},
            'max_dict': {},
            'combo_dict': {},
            'mixed_dict': {}
        }
        for drm in mutation_list:
            drm_group = drm.strip().rstrip(',')
            if drm_group.startswith('MAX') and drm_group.find('AND') != -1:
                mixed_drms = self.parse_mixed_condition(drm_group)
                for aa in mixed_drms:
                    mixed_group = {}
                    mutation = aa.split('=>')
                    score = int(mutation[1].strip())
                    drm = str(mutation[0].strip())
                    if drm.find('AND') != -1:
                        combinations = self.parse_combo_condition(drm)
                        for choice in combinations:
                            drm = choice
                            mixed_group.update({drm: score})
                    else:
                        drm = mutation[0].strip()
                        mixed_group.update({drm:score})
                    drm_scores['mixed_dict'].update({drm_group: mixed_group})


            elif drm_group.startswith('MAX'):
                max_group = {}
                max_drms = self.parse_max_condition(drm_group)
                for aa in max_drms:
                    mutation = aa.split('=>')
                    drm = mutation[0].strip()
                    score = int(mutation[1].strip())
                    max_group.update({drm: score})
                drm_scores['max_dict'].update({drm_group: max_group})

            elif drm_group.find('AND') != -1:
                mutation = drm_group.split('=>')
                score = int(mutation[1].strip())
                combo_drms = self.parse_combo_condition(drm_group)
                combo_group = {}
                for aa in combo_drms[:-1]:
                    drm = aa
                    combo_group.update({drm: score})
                drm_scores['combo_dict'].update({drm_group: combo_group})

            else:    # parsing for a single drm condition
                score_cond = drm_group.split()
                drm = score_cond[0].strip()
                score = int(score_cond[2].strip())
                drm_scores['single_drm_dict'].update({drm: score})
        return drm_scores


    def parse_max_condition(self, drm_group):
        regex = '[\S]+[\s]*=>[\s]*[\d]+'
        max_drms = re.findall(regex, drm_group)
        return(max_drms)

    def parse_combo_condition(self, drm_group):
        combinations = re.split('[(]|[\s]*AND[\s]*|[)][\s]*', drm_group)
        return combinations

    def parse_mixed_condition(self, drm_group):
        regex = '[\S]+[\s]*=>[\s]*[\d]+ | [(]{1}[\d]+[\S]+[\s]*AND[\s]*[\S]+[\s]*=>[\s]*[\d]+'
        mixed_drms = re.findall(regex, drm_group)
        return(mixed_drms)


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

