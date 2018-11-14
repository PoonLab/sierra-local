import xml.etree.ElementTree as xml
import requests
import urllib.request
import glob, os
from pathlib import Path
import re
import sierralocal.updater as updater


class HIVdb():
    """
    Define a class for handling network transactions with the Stanford HIVdb
    webserver, to retrieve the rules-based prediction algorithm as ASI XML,
    and convert this information into Python objects.
    """
    def __init__(self, path=None, forceupdate=False):
        file_found = False
        file_newest = False  # FIXME: this appears to be unused
        self.xml_filename = None
        self.BASE_URL = 'https://hivdb.stanford.edu'

        # args force update flag is TRUE
        if forceupdate:
            self.xml_filename = updater.update_HIVDB()

        if path is None:
            # If user has not specified XML path
            # Iterate over possible HIVdb ASI files matching the glob pattern
            dest = str(Path(os.path.dirname(__file__))/'data'/'HIVDB*.xml')
            print("searching path " + dest)
            files = glob.glob(dest)

            # find the newest XML that can be parsed
            intermed = []
            for file in files:
                version = re.search("HIVDB_([0-9]\.[0-9.]+)\.", file).group(1)
                intermed.append((version, file))
            intermed.sort(reverse=True)

            for version, file in intermed:
                try:
                    xml.parse(file)
                    file_found = True  # if no exception thrown
                    self.xml_filename = file
                    break
                except:
                    print('Failed to parse XML file. Please post an issue at http://github.com/PoonLab/sierra-local/issues.')
                    raise
        else:
            # The user has specified XML path
            if os.path.isfile(path):
                # Ensure is a file
                try:
                    # Ensure is XML parseable
                    xml.parse(path)
                    file_found = True
                    self.xml_filename = path
                except:
                    raise
            else:
                print("HIVDB XML cannot be found at user specified path {}".format(path))
                raise

        # Parseable XML file not found. Update from web
        if not file_found:
            print("Error: could not find local copy of HIVDB XML, attempting download...")
            self.xml_filename = updater.update_HIVDB()

        dest = str(Path(os.path.dirname(__file__))/'data'/'apobec.tsv')
        print(dest)
        if not os.path.isfile(dest):
            print("Error: could not retrieve APOBEC DRM data")
            updater.updateAPOBEC()

        # Set algorithm metadata
        self.root = xml.parse(str(self.xml_filename)).getroot()
        self.algname = self.root.find('ALGNAME').text
        self.version = self.root.find('ALGVERSION').text
        self.version_date = self.root.find('ALGDATE').text
        print("HIVdb version",self.version)


    def parse_definitions(self, root):
        """ parse_definitions function meant to assemble definitions into nested dictionary and return

            @param root: algorithm root
            @return definitions: dictionary of gene, level, druglcass, globalrange, and comment dictionaries
        """
        self.definitions = {
            'gene': {},  # gene target names and drug classes
            'level': {},  # maps from level to S/I/R symbols
            'drugclass': {},  # maps drug class to drugs
            'globalrange': {},  # maps from score to level
            'comment': {}  # maps comments to id string
        }
        # Convert list of elements from class 'xml.etree.ElementTree.Element' to type 'str'
        element_list = list(map(lambda x: xml.tostring(x).strip().decode("utf-8"), root.getchildren()))
        # Find the index of element 'DEFINITIONS' so that it's children may be iterated over to parse definitions
        def_index = [i for i, item in enumerate(element_list) if re.search('<DEFINITIONS>', item)]
        def_ind = def_index[0]  # un-list the index of 'DEFINITIONS' element

        definitions = root.getchildren()[def_ind]
        comment_definitions = definitions.getchildren()[-1]  # TODO: swap out hard-coded index with variable

        globalrange = definitions.find('GLOBALRANGE').text.split(',')
        default_grange = self.parse_globalrange(self.definitions['globalrange'], globalrange)

        for element in definitions.getchildren():
            if element.tag == 'GENE_DEFINITION':
                gene = element.find('NAME').text
                drug_classes = element.find('DRUGCLASSLIST').text.split(',')
                self.definitions['gene'].update({gene: drug_classes})

            elif element.tag == 'LEVEL_DEFINITION':
                order = element.find('ORDER').text
                original = element.find('ORIGINAL').text
                sir = element.find('SIR').text
                self.definitions['level'].update({order: {original: sir}})

            elif element.tag == 'DRUGCLASS':
                name = element.find('NAME').text
                druglist = element.find('DRUGLIST').text.split(',')
                self.definitions['drugclass'].update({name: druglist})

            elif element.tag == 'COMMENT_DEFINITIONS':
                for comment_str in comment_definitions.getchildren():
                    id = comment_str.attrib['id']
                    comment = comment_str.find('TEXT').text
                    sort_tag = comment_str.find('SORT_TAG').text
                    self.definitions['comment'].update({id: {sort_tag: comment}})

        return(self.definitions)



    def parse_globalrange(self, grange_dict, scorerange):
        """ parse_globalrange function meant to assemble a scorerange into a dictionary and return

            @param grange_dict: dictionary for storing the scoring range specifications on
            @param scorerange: text to parse out the scorerange for
            @return grange_dict: dictionary for storing scorerange of the parsed range
        """
        for item in scorerange:
            order = int(re.split('=>', item)[1].strip('() '))  # str containing order number: '1'
            range = re.split('=>', item)[0].strip('() ')  # str containing the range: '-INF TO 9'
            min = re.split('TO', range)[0].strip()  # str containing min val in range: '-INF'
            max = re.split('TO', range)[1].strip()  # str containing max val in range: '9'

            # convert_to_num converts strings to integers, and also 'INF' and '-INF' to their
            # numerical representations
            def convert_to_num(s):
                if s == '-INF':
                    return float('-inf')
                elif s == 'INF':
                    return float('inf')
                else:
                    return int(s)

            min = convert_to_num(min)
            max = convert_to_num(max)
            grange_dict.update({order: {'min': min, 'max': max}})
        return grange_dict


    def parse_drugs(self, root):
        """ parse_drugs iterates through each drug in HIVDB,
            parses condition for a specific drug,
            and assigns a library of the drug resistant mutation conditions to the dictionary of drugs
            Also includes score ranges associated with the drug, (which is most often the default globalrange)

            @param root: algorithm root
            @return self.drugs: populated dictionary of drugs, associated with their (global) score ranges
        """
        self.drugs = {}

        for element in root.getchildren():
            if element.tag == 'DRUG':
                drug = element.find('NAME').text                            # drug name
                fullname = element.find('FULLNAME').text                    # drug full name
                condition = element.find('RULE').find('CONDITION').text     # drug conditions
                cond_dict = self.parse_condition(condition)                 # dictionary of parsed drug conditions

                scorerange = list(element.find('RULE').find('ACTIONS').find('SCORERANGE'))[0]
                if scorerange.find('USE_GLOBALRANGE') == None:
                    self.drugs[drug] = self.drugs[fullname] = (cond_dict, self.definitions['globalrange'])    #default
                else:
                    sep_dict = {}
                    globalrange = self.parse_globalrange(sep_dict, scorerange)
                    self.drugs[drug] = self.drugs[fullname] = (cond_dict, globalrange)

        return self.drugs


    def parse_condition(self, condition):
        """ parse_condition function takes a given condition (one of four types)
            'MAXAND' condition: MAX ((41L AND 215ACDEILNSV) => 5, (41L AND 215FY) => 15)
            'MAX' condition: MAX ( 219E => 5, 219N => 5, 219Q => 5, 219R => 5 )
            'AND' condition: (67EGN AND 215FY AND 219ENQR) => 5
            'single-drm' condition: 62V => 5

            @param condition: given drug condition to parse
            @return self.drms: list library updated with all the DRM conditions associated with given drug condition
        """
        # drug resistant mutation (DRM)
        mutation_list = ((condition.split('(',1)[1].rstrip(')')) + ',').split('\n')    # NOTE: \n may not be stable or transferable on other platforms
        self.drms= []                            # library of drms for the drug of the given condition

        for drm in mutation_list:
            if drm.strip().startswith('MAX'):
                max_lib = []
                max_chunks = re.findall('(\(?\d+[A-z]+[\s?AND\s?\d+\w]+\)?\s?=>\s?\d+|\(?\d+[A-z]+[\s?AND\s?\d+\w]+\)?\s?=>\s?-\d+)', drm)
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

                                                                                                                    #TODO: combine drm and iter so less params
    def _parse_scores(self, drm_lib, drm, chunk, iter):
        """ _parse_scores function is a helper function to parse_condition.
            Parses the residues, amino acids, and scores associated with a particular DRM,
            then updates the specified list library

            @param drm_lib: given library to be updated with DRMs
            @param drm: full name of original drug resistant mutation
            @param chunk: drm_group of one of the condition types   (kinda vague...will fix)
            @param iter: index to keep track of which DRM is associated with respective index in the extracted list of scores
        """
        scores = re.findall('([0-9]+(?=\W)|-[0-9]+(?=\W))', drm.strip())   # extract scores in same order as grouped drm tuples; stored in (indexable) list
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


    def parse_comments(self, root):
        """ parse_comments function retrieves all mutation comments and stores in a dictionary
            Dictionary is further separated into protease (PR), integrase (IN), and reverse transcriptase (RT) dictionaries

            @param root: algorithm root
            @return self.comments: a dictionary with key = mutation and value = mutation comment reference for PR, IN, RT
        """
        self.comments = {}

        for element in root.getchildren():
            if element.tag == 'MUTATION_COMMENTS':

                for gene in element.getchildren():
                    name = gene.find('NAME').text
                    gene_dict = {}

                    for rule in gene.findall('RULE'):
                        condition = rule.find('CONDITION').text
                        reference = rule.find('ACTIONS').find('COMMENT').get('ref')
                        gene_dict.update({condition: reference})

                    self.comments.update({name: gene_dict})

        return(self.comments)


    #def retrieve_comment(self, genetype, condition, reference):
        #if condition.contains('i') or condition.contains('d'):
            # condition can include an insertion or deletion


def main():
    hivdb = HIVdb()

if __name__ == "__main__":
    main()
