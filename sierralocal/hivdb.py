import os
import re
import sys
import glob
from pathlib import Path
import xml.etree.ElementTree as xml


class HIVdb():
    """
    Define a class for handling network transactions with the Stanford HIVdb
    webserver, to retrieve the rules-based prediction algorithm as ASI XML,
    and convert this information into Python objects.
    """
    def __init__(self, asi2=None, apobec=None, forceupdate=False):
        self.xml_filename = None
        self.tsv_filename = None
        self.BASE_URL = 'https://hivdb.stanford.edu'

        if forceupdate:
            # DEPRECATED, requires selenium, chrome and chromedriver
            import sierralocal.updater as updater
            self.xml_filename = updater.update_HIVDB()
            self.tsv_filename = updater.updateAPOBEC()
        else:
            self.set_hivdb_xml(asi2)
            self.set_apobec_tsv(apobec)

        # Set algorithm metadata
        self.root = xml.parse(str(self.xml_filename)).getroot()
        self.algname = self.root.find('ALGNAME').text
        self.version = self.root.find('ALGVERSION').text
        self.version_date = self.root.find('ALGDATE').text
        print("HIVdb version", self.version)

    def set_hivdb_xml(self, path):
        """
        Assigns user specified XML file path to self.xml_filename.
        Finds a newest XML that can be parsed and assigns.
        @param path: str, path to XML file.
        """
        file_found = False

        # If user has not specified XML path
        # Iterate over possible HIVdb ASI files matching the glob pattern
        if path is None:
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
                    print('Failed to parse XML file. Please post an issue at '
                          'http://github.com/PoonLab/sierra-local/issues.')
                    raise

        # The user has specified XML path
        else:
            if os.path.isfile(path):    # Ensure is a file
                try:
                    xml.parse(path)     # Ensure is XML parseable
                    file_found = True
                    self.xml_filename = path
                except:
                    raise
            else:
                raise FileNotFoundError(
                        "HIVDB XML cannot be found at user specified "
                        "path {}".format(path))

        # Parseable XML file not found. Update from web
        if not file_found:
            print("Error: could not find local copy of HIVDB XML.")
            print("Manually download from https://hivdb.stanford.edu/page/"
                  "release-notes/#algorithm.updates")
            sys.exit()
            # self.xml_filename = updater.update_HIVDB()

    def set_apobec_tsv(self, path):
        """
        Attempt to locate a local APOBEC DRM file (tsv format)
        @param path: str, path to TSV file
        """
        if path is None:
            dest = str(Path(os.path.dirname(__file__))/'data'/'apobec*.tsv')
            print("searching path {}".format(dest))
            files = glob.glob(dest)
            for file in files:
                # no version numbering, take first hit
                # TODO: some basic format check on TSV file
                if os.path.isfile(file):
                    self.tsv_filename = file
                    return
        else:
            self.tsv_filename = path
            return

        # if we end up here, no local files found
        print("Error: could not locate local APOBEC DRM data file.")
        print("Manually download from https://hivdb.stanford.edu/page/"
              "release-notes/#data.files")
        sys.exit()
        # self.tsv_filename = updater.updateAPOBEC()

    def parse_definitions(self, root):
        """
        parse_definitions function meant to assemble definitions into
        nested dictionary and return.
        @param root: xml.etree.ElementTree.ElementTree, algorithm root
        @return self.definitions: dict, dictionary with keys = root xml tags
        """
        self.definitions = {
            'gene': {},         # gene target names and drug classes
            'level': {},        # maps from level to S/I/R symbols
            'drugclass': {},    # maps drug class to drugs
            'globalrange': {},  # maps from score to level
            'comment': {}       # maps comments to id string
        }

        # Convert list of items of 'xml.etree.ElementTree.Element' type
        # to list of items of 'str' type
        element_list = \
            list(map(lambda x: xml.tostring(x).strip().decode("utf-8"),
                     root.findall('./')))

        # Find the index of element 'DEFINITIONS' so that it's
        # children may be iterated over to parse definitions
        def_index = [
            i for i, item in enumerate(element_list)
            if re.search('<DEFINITIONS>', item)]
        def_ind = def_index[0]  # un-list the index of 'DEFINITIONS' element

        definitions = root.findall('./')[def_ind]
        comment_definitions = definitions.findall('./')[-1]  # TODO: swap out hard-coded index with variable

        globalrange = definitions.find('GLOBALRANGE').text.split(',')
        default_grange = \
            self.parse_globalrange(self.definitions['globalrange'],
                                   globalrange)

        for element in definitions.findall('./'):
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
                for comment_str in comment_definitions.findall('./'):
                    id = comment_str.attrib['id']
                    comment = comment_str.find('TEXT').text
                    sort_tag = comment_str.find('SORT_TAG').text
                    self.definitions['comment'].update({id: {sort_tag: comment}})

        return self.definitions

    def parse_globalrange(self, grange_dict, scorerange):
        """
        parse_globalrange function meant to assemble a scorerange 
        into a dictionary and return
        @param grange_dict: dict, dictionary for storing the 
        scoring range specifications on
        @param scorerange: list, list of strings text to parse out
        the scorerange for. Eg: ['(-INF TO 9 => 1', '  10 TO 14 => 2',
                                 '  15 TO 29 => 3', '  30 TO 59 => 4',
                                 '  60 TO INF => 5)']
        @return grange_dict: dict, dictionary for storing scorerange 
        of the parsed range
        """
        for item in scorerange:
            order = int(re.split('=>', item)[1].strip('() '))  # str containing order number: '1'
            range = re.split('=>', item)[0].strip('() ')  # str containing the range: '-INF TO 9'
            min = re.split('TO', range)[0].strip()  # str containing min val in range: '-INF'
            max = re.split('TO', range)[1].strip()  # str containing max val in range: '9'

            # convert_to_num converts strings to integers, and also 
            # 'INF' and '-INF' to their numerical representations
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
        """ 
        parse_drugs iterates through each drug in HIVDB, parses condition 
        for a specific drug, and assigns a library of the drug resistant 
        mutation conditions to the dictionary of drugs. Also includes scoreranges
        associated with the drug, (which is most often the default globalrange)
        @param root: xml.etree.ElementTree.ElementTree, algorithm root
        @return self.drugs: dict, populated dictionary of drugs, associated 
        with their (global) score ranges
        """
        self.drugs = {}

        for element in root.findall('./'):
            if element.tag == 'DRUG':
                drug = element.find('NAME').text                 # drug name
                fullname = element.find('FULLNAME').text         # drug full name
                condition = element.find('RULE/CONDITION').text  # drug conditions
                cond_dict = self.parse_condition(condition)      # dictionary of parsed drug conditions

                # scorerange = list(element.find('RULE').find('ACTIONS').find('SCORERANGE'))[0]
                # same as above but PEP 8 compliant
                scorerange = list(element.find('RULE/ACTIONS/SCORERANGE'))[0]
                if scorerange.find('USE_GLOBALRANGE') is None:
                    self.drugs[drug] = self.drugs[fullname] = \
                        (cond_dict, self.definitions['globalrange'])    #default
                else:
                    sep_dict = {}
                    globalrange = self.parse_globalrange(sep_dict, scorerange)
                    self.drugs[drug] = self.drugs[fullname] = \
                        (cond_dict, globalrange)

        return self.drugs

    def parse_condition(self, condition):
        """
        parse_condition function takes a given condition (one of four types)
        'MAXAND' condition: MAX ((41L AND 215ACDEILNSV) => 5, (41L AND 215FY) => 15)
        'MAX' condition: MAX ( 219E => 5, 219N => 5, 219Q => 5, 219R => 5 )
        'AND' condition: (67EGN AND 215FY AND 219ENQR) => 5
        'single-drm' condition: 62V => 5
        @param condition: str, given drug condition to parse
        @return self.drms: list, library updated with all the DRM conditions
        associated with given drug condition
        """
        # drug resistant mutation (DRM)
        mutation_list = (
                (condition.split('(', 1)[1].rstrip(')')) + ','
            ).split('\n')   # NOTE: \n may not be stable or transferable on other platforms
        self.drms= []       # library of drms for the drug of the given condition

        for drm in mutation_list:
            if drm.strip().startswith('MAX'):
                max_lib = []
                max_chunks = re.findall('(\(?\d+[A-z]+[\s?AND\s?\d+\w]+\)?\s?=>\s?\d+|\(?\d+[A-z]+[\s?AND\s?\d+\w]+\)?\s?=>\s?-\d+)', drm)
                iter = 0
                for chunk in max_chunks:
                    # for both MAX conditions, need to create a mini-library
                    # that will keep all of the individual DRMs together
                    self._parse_scores(max_lib, drm, chunk, iter)
                    iter += 1               # iter is created to keep track of which DRM associated with which index in the list of scores in _parse_scores function
                self.drms.append(max_lib)   # finally append this mini-library to the DRMs library

            else:
                iter = 0
                self._parse_scores(self.drms, drm, drm, iter)

        return self.drms    #TODO: combine drm and iter so less params

    def _parse_scores(self, drm_lib, drm, chunk, iter):
        """
        _parse_scores function is a helper function to parse_condition.
        Parses the residues, amino acids, and scores associated with a 
        particular DRM, then updates the specified list library
        @param drm_lib: list, given library to be updated with DRMs
        @param drm: str, full name of original drug resistant mutation
        @param chunk: str, drm group of one of the condition types. Eg: '65R => 45'
        @param iter: int, index to keep track of which DRM is associated 
        with respective index in the extracted list of scores
        """
        scores = re.findall('([-]?[0-9]+(?=\W)|[-]?[0-9]+(?=$))', drm.strip())   # extract scores in same order as grouped drm tuples; stored in (indexable) list
        r_and_r = re.findall('\d+[A-z]+[\s?AND\s?\d+\w]+', chunk)

        for combo_group in r_and_r:
            mut_list = []
            residue_amino_acid = re.findall('[0-9]+[A-z]+', combo_group.strip())  # TODO: needs testing
            for mutation in residue_amino_acid:
                residue = int(re.findall('[\d]+(?!\d)(?=\w)', mutation)[0])
                amino_acid = str(re.findall('[0-9]+([A-z]+)', mutation)[0])
                mut_list.append(tuple((residue, amino_acid)))

            # populate the drms library with the new drm condition
            drm_lib.append({'group': mut_list, 'value': int(scores[iter])})
            # wipe out scores stored in variable scores for next batch
            scores[:] = []

    def parse_comments(self, root):
        """ 
        parse_comments function retrieves all mutation comments and stores
        in a dictionary, which is further separated into protease (PR),
        integrase (IN), and reverse transcriptase (RT) dictionaries.
        @param root: xml.etree.ElementTree.ElementTree, algorithm root
        @return self.comments: dict, a dictionary with key = mutation and
        value = mutation comment reference for PR, IN, RT
        """
        self.comments = {}

        for element in root.findall('./'):
            if element.tag == 'MUTATION_COMMENTS':

                for gene in element.findall('.//GENE'):
                    name = gene.find('NAME').text
                    gene_dict = {}

                    for rule in gene.findall('RULE'):
                        condition = rule.find('CONDITION').text
                        reference = rule.find('ACTIONS/COMMENT').get('ref')
                        gene_dict.update({condition: reference})

                    self.comments.update({name: gene_dict})

        return self.comments


def main():
    # simply attempt to load XML and TSV files
    HIVdb()


if __name__ == "__main__":
    main()
