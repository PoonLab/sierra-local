import unittest
from xml.etree.ElementTree import parse
from sierralocal.main import *
from sierralocal.hivdb import HIVdb

class TestScore(unittest.TestCase):
    def setUp(self):
        self.algorithm = HIVdb()
        self.root = parse(self.algorithm.xml_filename)
        self.definitions = self.algorithm.parse_definitions(self.root)
        self.drugs = self.algorithm.parse_drugs(self.root)

    def testScorefile(self):
        # Setting params
        input_file = r'tests\hxb2-pr.fa'

        exp_sequence_headers = ['HXB2-PR', 'shift1', 'shift2',
                                'plus1', 'plus_codon', 'del1_after_3codons',
                                'insAAA_after_3codons']
        exp_sequence_scores = [[{'ATV/r': (0, [], []),
                                 'DRV/r': (0, [], []),
                                 'FPV/r': (0, [], []),
                                 'IDV/r': (0, [], []),
                                 'LPV/r': (0, [], []),
                                 'NFV': (0, [], []),
                                 'SQV/r': (0, [], []),
                                 'TPV/r': (0, [], [])}],
                               [{'ATV/r': (0, [], []),
                                 'DRV/r': (0, [], []),
                                 'FPV/r': (0, [], []),
                                 'IDV/r': (0, [], []),
                                 'LPV/r': (0, [], []),
                                 'NFV': (0, [], []),
                                 'SQV/r': (0, [], []),
                                 'TPV/r': (0, [], [])}],
                               [{'ATV/r': (0, [], []),
                                 'DRV/r': (0, [], []),
                                 'FPV/r': (0, [], []),
                                 'IDV/r': (0, [], []),
                                 'LPV/r': (0, [], []),
                                 'NFV': (0, [], []),
                                 'SQV/r': (0, [], []),
                                 'TPV/r': (0, [], [])}],
                               [{'ATV/r': (0, [], []),
                                 'DRV/r': (0, [], []),
                                 'FPV/r': (0, [], []),
                                 'IDV/r': (0, [], []),
                                 'LPV/r': (0, [], []),
                                 'NFV': (0, [], []),
                                 'SQV/r': (0, [], []),
                                 'TPV/r': (0, [], [])}],
                               [{'ATV/r': (0, [], []),
                                 'DRV/r': (0, [], []),
                                 'FPV/r': (0, [], []),
                                 'IDV/r': (0, [], []),
                                 'LPV/r': (0, [], []),
                                 'NFV': (0, [], []),
                                 'SQV/r': (0, [], []),
                                 'TPV/r': (0, [], [])}],
                               [{'ATV/r': (0, [], []),
                                 'DRV/r': (0, [], []),
                                 'FPV/r': (0, [], []),
                                 'IDV/r': (0, [], []),
                                 'LPV/r': (0, [], []),
                                 'NFV': (0, [], []),
                                 'SQV/r': (0, [], []),
                                 'TPV/r': (0, [], [])}],
                               [{'ATV/r': (0, [], []),
                                 'DRV/r': (0, [], []),
                                 'FPV/r': (0, [], []),
                                 'IDV/r': (0, [], []),
                                 'LPV/r': (0, [], []),
                                 'NFV': (0, [], []),
                                 'SQV/r': (0, [], []),
                                 'TPV/r': (0, [], [])}]]
        exp_ordered_mutation_list = [[[(3, 'V', 'I'), (37, 'S', 'N')]],
                                     [[(3, 'V', 'I'), (37, 'S', 'N')]],
                                     [[(3, 'V', 'I'), (37, 'S', 'N')]],
                                     [[(3, 'V', 'I'), (37, 'S', 'N')]],
                                     [[(3, 'V', 'I'), (37, 'S', 'N')]],
                                     [[(3, 'V', 'I'), (4, 'X', 'T'), (37, 'S', 'N')]],
                                     [[(3, 'V', 'I'), (37, 'S', 'N')]]]
        exp_file_genes = [[('PR', 1, 99, 1, 294)],
                          [('PR', 2, 99, 3, 293)],
                          [('PR', 2, 99, 2, 292)],
                          [('PR', 1, 99, 2, 295)],
                          [('PR', 1, 99, 1, 297)],
                          [('PR', 1, 99, 1, 293)],
                          [('PR', 1, 99, 1, 297)]]
        exp_sequence_lengths = [[294], [291], [291], [294], [297], [293], [297]]
        exp_file_trims = [[(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)]]
        exp_subtypes = ['', '', '', '', '', '', '']

        res_sequence_headers, res_sequence_scores, res_ordered_mutation_list, res_file_genes, \
        res_sequence_lengths, res_file_trims, res_subtypes = scorefile(input_file=input_file, algorithm=self.algorithm)

        self.assertEqual(exp_sequence_headers, res_sequence_headers)
        self.assertEqual(exp_sequence_scores, res_sequence_scores)
        self.assertEqual(exp_ordered_mutation_list, res_ordered_mutation_list)
        self.assertEqual(exp_file_genes, res_file_genes)
        self.assertEqual(exp_sequence_lengths, res_sequence_lengths)
        self.assertEqual(exp_file_trims, res_file_trims)
        self.assertEqual(exp_subtypes, res_subtypes)


if __name__ == '__main__':
    unittest.main()