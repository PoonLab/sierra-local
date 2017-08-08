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
                #print(condition, '\n\n')
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

            # if this is a mixed condition with both conditional 'AND' as well as an overall 'MAX'
            # storing each drug resistant mutation group as a nested dictionary within the mixed dictionary
            if drm_group.startswith('MAX') and drm_group.find('AND') != -1:
                mixed_drms = self.parse_mixed_condition(drm_group)
                drm_scores['mixed_dict'].update({drm_group: mixed_drms})

            # if this is a 'MAX' condition, store each DRM group as a nested dictionary within max_dict
            elif drm_group.startswith('MAX'):
                max_drms = self.parse_max_condition(drm_group)
                drm_scores['max_dict'].update({drm_group: max_drms})

            # if this is strictly and 'AND' condition, store individual DRMs into nested dictionary within combo_dict
            elif drm_group.find('AND') != -1:
                mutation = drm_group.split('=>')
                score = int(mutation[1].strip())
                combo_drms = self.parse_combo_condition(drm_group, score)
                drm_scores['combo_dict'].update({drm_group: combo_drms})

            # parsing for a single drm condition, no helper function necessary
            else:
                score_cond = drm_group.split()
                drm = score_cond[0].strip()
                score = int(score_cond[2].strip())
                drm_scores['single_drm_dict'].update({drm: score})

        return drm_scores


    # example condition: 'MAX ( 65E => 10, 65N => 30, 65R => 45 )'
    # @return {'65E': 10, '65N': 30, '65R': 45}
    def parse_max_condition(self, drm_group):
        max_group = {}
        regex = '[\S]+[\s]*=>[\s]*[\d]+'
        max_drms = re.findall(regex, drm_group)

        for aa in max_drms:
            mutation = aa.split('=>')
            drm = mutation[0].strip()
            score = int(mutation[1].strip())
            max_group.update({drm: score})
        return(max_group)


    # example condition: '(40F AND 41L AND 210W AND 215FY) => 5'
    # @return {'40F': 5, '41L': 5, '210W': 5, '215FY': 5}
    def parse_combo_condition(self, drm_group, score):
        combo_group = {}
        regex = '[\d]+[A-Za-z]+'
        combinations = re.findall(regex, drm_group)

        for aa in combinations:
            drm = aa
            combo_group.update({drm: score})
        return combo_group


    # example condition: 'MAX ((210W AND 215ACDEILNSV) => 5, (210W AND 215FY) => 10)'
    # @return {'(210W AND 215ACDEILNSV)': {'210W': 5, '215ACDEILNSV': 5}
    def parse_mixed_condition(self, drm_group):
        mixed_group = {}  # nested dictionary
        regex = '[\S]+[\s]*=>[\s]*[\d]+|[(]{1}[\d]+[\S]+[\s]*AND[\s]*[\S]+[\s]*=>[\s]*[\d]+'
        mixed_drms = re.findall(regex, drm_group)

        for aa in mixed_drms:
            mutation = aa.split('=>')
            score = int(mutation[1].strip())
            drm = str(mutation[0].strip())
            # parse each combination condition within the mixed condition and update nested dictionary
            if drm.find('AND') != -1:
                combinations = self.parse_combo_condition(aa, score)
                mixed_group.update({aa: combinations})
            # for the possibility that there may be a single DRM mixed with the 'AND' conditionals
            else:
                mixed_group.update({aa: score})

        return(mixed_group)


    # every dictionary one-level deep within self.drugs will be run through, checking with the user-given sequence
    # stores 4 individuals scores, which are then summed and returned
    def score_drugs(self, drugname, sequence):
        FOUND = False
        if drugname not in self.drugs.keys(): print("Drugname: " +  drugname + " not found.")
        else: FOUND = True
        if FOUND:
            #calculate scores
            single_score = self.score_single(drugname, sequence)
            max_score = self.score_max(drugname, sequence)
            combo_score = self.score_combo(drugname, sequence)
            mixed_score = self.score_mixed(drugname, sequence)

            print(max_score)
            print(mixed_score)
            total = single_score + max_score + combo_score + mixed_score
            return total


    def score_single(self, drugname, sequence):
        single_drm_scores = self.drugs[drugname]['single_drm_dict']
        score = 0
        # preliminary testing
        # single_drm_scores = {'15Q': 10, '7G': 10, '5Y': 15}
        for drm in single_drm_scores.keys():
            residue = int(re.findall('[\d]+', drm)[0])
            aa = re.findall('[0-9]+([A-Za-z])', drm)
            if sequence[residue - 1] in aa:              # TODO: check and throw IndexOutofRange Error (then apply to all score functions)
                score += single_drm_scores[drm]
        return score


    def score_max(self, drugname, sequence):
        max_drm_scores = self.drugs[drugname]['max_dict']
        scores = [0]
        # preliminary testing
        # max_drm_scores = {'MAX ( 1D => 10, 17E => 10, 19F => 15)': {'1D': 10, '17E': 10, '19F': 15}}
        for drm in max_drm_scores.keys():
            key_val_pairs = self.parse_max_condition(drm)
            for key in key_val_pairs.keys():
                residue = int(re.findall('[\d]+', key)[0])
                aa = re.findall('[0-9]+([A-Za-z]+)', key)
                if sequence[residue - 1] in aa:
                    scores.append(key_val_pairs[key])

        return(max(scores))


    def score_combo(self, drugname, sequence):
        combo_drm_scores = self.drugs[drugname]['combo_dict']
        score = 0
        # preliminary testing
        # combo_drm_scores = {'(10E AND 18C AND 20Q) => 15': {'10E': 15, '18C': 15, '20Q': 15}}
        for drm in combo_drm_scores.keys():
            count = 0
            mutation = drm.split('=>')
            score = int(mutation[1].strip())
            key_val_pairs = self.parse_combo_condition(drm, score)
            for key in key_val_pairs.keys():
                residue = int(re.findall('[\d]+', key)[0])
                aa = re.findall('[0-9]+([A-Za-z]+)', key)
                if sequence[residue - 1] not in aa:
                    # print(sequence[residue -1], " is not in ", key_val_pairs.keys())
                    continue
                count += 1
                if count == len(key_val_pairs.keys()):
                    score += (key_val_pairs[key])
        return(score)

    def score_mixed(self, drugname, sequence):
        mixed_drm_scores = self.drugs[drugname]['mixed_dict']
        total_scores = [0]
        # preliminary testing
        # mixed_drm_scores= {'MAX ((17T AND 18C) => 5, (19F AND 18C) => 15)': {'(17T AND 18C) => 5': {'17T': 5, '18C': 5},'(19F AND 18C) => 15': {'19F': 15, '18C': 15}}}
        for drm in  mixed_drm_scores.keys():
            key_val_pairs = self.parse_mixed_condition(drm)

            for key in key_val_pairs:
                count = 0
                mutation = key.split('=>')
                score = int(mutation[1].strip())

                if key.find('AND') != -1:
                    split_keys = re.findall('[0-9]+[A-Za-z]+', key)

                    for i in split_keys:
                        residue = int(re.findall('[\d]+(?!\d)(?=\w)', i)[0])
                        aa = re.findall('[0-9]+([A-Za-z]+)', i)

                        if sequence[residue - 1] not in aa:
                            # print(sequence[residue -1], " is not present in sequence.")
                            continue
                        count += 1

                        if count == len(split_keys):
                            total_scores.append(score)
                else:
                    aa = re.findall('[0-9]+([A-Za-z]+)', key)
                    residue = int(re.findall('[\d]+', key)[0])
                    if sequence[residue - 1] in aa:
                        total_scores.append(score)

        return(max(total_scores))







def main():
    alg = HIVdb("/home/tng92/git/sierra-local/HIVDB.xml")
    alg.parse_definitions(alg.root)
    alg.parse_drugs(alg.root)
    #print(alg.definitions)
    #print(alg.drugs.keys())
    print(alg.score_drugs("etravirine", 'DAAAAAGAAELGAAAATCTQAAA'))


main()

