import xml.etree.ElementTree as xml

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
            'level': {}  # maps from level to S/I/R symbols
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
                globalrange = element.find('NAME').text
                self.definitions['globalrange'].update({globalrange})
            elif element.tag == 'LEVEL_DEFINTIION':
                level = element.find('NAME').text
                order = element.find('ORDER').text.split(',')  # TODO: check level parameter structures
                original = element.find('ORIGINAL').text.split(',')
                sir = element.find('SIR').text.split(',')

    def parse_drugs(self, root):
        self.drugs = {
            'name': {},
            'fullname': {},
            'rule': {
                'condition': {},
                'actions': {
                    'scorerange': False
                }
            }
        }
        for element in root.getchildren():



