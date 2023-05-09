import unittest
from xml.etree.ElementTree import parse
from sierralocal.score_alg import score_drugs, score_single
from sierralocal.hivdb import HIVdb


class TestScoreAlg(unittest.TestCase):
    def setUp(self):
        self.algorithm = HIVdb()
        self.root = parse(self.algorithm.xml_filename)
        self.algorithm.parse_definitions(self.root)
        self.algorithm.parse_drugs(self.root)

    def testScoreSingle(self):
        # Setting params
        drugname = 'ABC'
        seq_mutations = {41: ['FY', 'L'], 65: ['FY', 'E']}

        exp_tuple = (15, [5, 10], [['FY41L'], ['FY65E']])
        res_tuple = score_single(self.algorithm, drugname, seq_mutations)
        
        self.assertEqual(exp_tuple, res_tuple)

        # Setting params
        drugname = 'BIC'
        seq_mutations = {51: ['F', 'Y']}

        exp_tuple = (10, [10], [['F51Y']])
        res_tuple = score_single(self.algorithm, drugname, seq_mutations)

        self.assertEqual(exp_tuple, res_tuple)


    def testScoreDrugs(self):
        # Setting params
        gene = 'IN'
        seq_mutations = {51: ['F', 'Y'], 66: ['N', 'K']}

        exp_dict = {'BIC': (25, [10, 15], [['F51Y'], ['N66K']]),
                    'CAB': (35, [15, 20], [['F51Y'], ['N66K']]),
                    'DTG': (25, [10, 15], [['F51Y'], ['N66K']]),
                    'EVG': (75, [15, 60], [['F51Y'], ['N66K']]),
                    'RAL': (75, [15, 60], [['F51Y'], ['N66K']])}
        
        res_scores = score_drugs(self.algorithm, gene, seq_mutations)

        self.assertEqual(exp_dict, res_scores)


if __name__ == '__main__':
    unittest.main()