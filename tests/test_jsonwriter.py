import os
import unittest

from sierralocal.hivdb import HIVdb
from sierralocal.jsonwriter import JSONWriter


class TestJsonWriter(unittest.TestCase):
    def setUp(self):
        self.algorithm = HIVdb()
        self.writer = JSONWriter(self.algorithm)

    # These are written for HIVdb version 9.4, the text may change
    def testFormatValidationResults(self):
        # Setting params
        validated = [('WARNING', "The ('RT', 1, 99, 1, 294) sequence contains just 194 codons, which is not sufficient for a comprehensive interpretation.")]

        exp_validation = [{'level': 'WARNING',
                           'message': "The ('RT', 1, 99, 1, 294) sequence contains just 194 codons, which is not sufficient for a comprehensive interpretation."}]
        res_validation = self.writer.format_validation_results(validated)

        self.assertEqual(exp_validation, res_validation)

        # Setting params
        validated = [('CRITICAL', 'Unable to process sequence.')]

        exp_validation = [{'level': 'CRITICAL',
                           'message': 'Unable to process sequence.'}]
        res_validation = self.writer.format_validation_results(validated)

        self.assertEqual(exp_validation, res_validation)

        # Setting params
        validated = [('WARNING', "The ('RT', 1, 99, 1, 294) sequence had 1 amino acid trimmed from its 5\u2032-end due to poor quality."),
                     ('WARNING', "The ('RT', 1, 99, 1, 294) sequence had 1 amino acid trimmed from its 3\u2032-end due to poor quality.")]

        exp_validation = [{'level': 'WARNING',
                           'message': "The ('RT', 1, 99, 1, 294) sequence had 1 amino acid trimmed from its 5\u2032-end due to poor quality."},
                          {'level': 'WARNING',
                           'message': "The ('RT', 1, 99, 1, 294) sequence had 1 amino acid trimmed from its 3\u2032-end due to poor quality."}]
        res_validation = self.writer.format_validation_results(validated)

        self.assertEqual(exp_validation, res_validation)

    def testFormatDrugResistance(self):
        # Setting params

        scores = {'BIC': (25, [10, 15], [['F51Y'], ['N66K']]),
                  'DTG': (25, [10, 15], [['F51Y'], ['N66K']]),
                  'EVG': (75, [15, 60], [['F51Y'], ['N66K']]),
                  'RAL': (75, [15, 60], [['F51Y'], ['N66K']]),
                  'CAB': (75, [15, 60], [['F51Y'], ['N66K']])}
        gene = 'IN'

        exp_drug_resistance = \
                        {'drugScores': [{'drug': {'displayAbbr': 'BIC', 'name': 'BIC'},
                                         'drugClass': {'name': 'INSTI'},
                                         'level': 3,
                                         'partialScores': [{'mutations': [{'comments': [{'text': 'V151I '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'accessory '
                                                                                                 'INSTI '
                                                                                                 'selected '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'occurs '
                                                                                                 'in '
                                                                                                 '1% '
                                                                                                 'to '
                                                                                                 '3% '
                                                                                                 'of '
                                                                                                 'viruses '
                                                                                                 'from '
                                                                                                 'ART-naive '
                                                                                                 'persons '
                                                                                                 'depending '
                                                                                                 'on '
                                                                                                 'subtype. '
                                                                                                 'Alone, '
                                                                                                 'it '
                                                                                                 'appears '
                                                                                                 'to '
                                                                                                 'have '
                                                                                                 'little '
                                                                                                 'or '
                                                                                                 'no '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'INSTI '
                                                                                                 'susceptibility. '
                                                                                                 'V151L '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'reduces '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG '
                                                                                                 'by '
                                                                                                 '15 '
                                                                                                 'to '
                                                                                                 '20-fold '
                                                                                                 'and '
                                                                                                 'to '
                                                                                                 'CAB '
                                                                                                 'and '
                                                                                                 'DTG '
                                                                                                 'by '
                                                                                                 'about '
                                                                                                 '3-fold. '
                                                                                                 'V151A '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'minimally '
                                                                                                 'reduced '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG. '
                                                                                                 '$listMutsIn{151(NOT '
                                                                                                 'IAL)} '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'unusual '
                                                                                                 'mutation '
                                                                                                 'at '
                                                                                                 'this '
                                                                                                 'position.',
                                                                                         'type': 'Accessory'}],
                                                                           'primaryType': 'Accessory',
                                                                           'text': 'F51Y'}],
                                                            'score': 10.0},
                                                           {'mutations': [{'comments': [{'text': 'T66K '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'uncommon '
                                                                                                 'non-polymorphic '
                                                                                                 'INSTI-selected '
                                                                                                 'mutation. '
                                                                                                 'It '
                                                                                                 'is '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'high-level '
                                                                                                 'EVG '
                                                                                                 'resistance, '
                                                                                                 'intermediate/high-level '
                                                                                                 'RAL '
                                                                                                 'resistance, '
                                                                                                 'and '
                                                                                                 'low-level '
                                                                                                 'DTG '
                                                                                                 'and '
                                                                                                 'CAB '
                                                                                                 'resistance. '
                                                                                                 'Its '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'BIC '
                                                                                                 'is '
                                                                                                 'not '
                                                                                                 'known.',
                                                                                         'type': 'Major'}],
                                                                           'primaryType': 'Major',
                                                                           'text': 'N66K'}],
                                                            'score': 15.0}],
                                         'score': 25.0,
                                         'text': 'Low-Level Resistance'},
                                        {'drug': {'displayAbbr': 'CAB', 'name': 'CAB'},
                                         'drugClass': {'name': 'INSTI'},
                                         'level': 5,
                                         'partialScores': [{'mutations': [{'comments': [{'text': 'V151I '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'accessory '
                                                                                                 'INSTI '
                                                                                                 'selected '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'occurs '
                                                                                                 'in '
                                                                                                 '1% '
                                                                                                 'to '
                                                                                                 '3% '
                                                                                                 'of '
                                                                                                 'viruses '
                                                                                                 'from '
                                                                                                 'ART-naive '
                                                                                                 'persons '
                                                                                                 'depending '
                                                                                                 'on '
                                                                                                 'subtype. '
                                                                                                 'Alone, '
                                                                                                 'it '
                                                                                                 'appears '
                                                                                                 'to '
                                                                                                 'have '
                                                                                                 'little '
                                                                                                 'or '
                                                                                                 'no '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'INSTI '
                                                                                                 'susceptibility. '
                                                                                                 'V151L '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'reduces '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG '
                                                                                                 'by '
                                                                                                 '15 '
                                                                                                 'to '
                                                                                                 '20-fold '
                                                                                                 'and '
                                                                                                 'to '
                                                                                                 'CAB '
                                                                                                 'and '
                                                                                                 'DTG '
                                                                                                 'by '
                                                                                                 'about '
                                                                                                 '3-fold. '
                                                                                                 'V151A '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'minimally '
                                                                                                 'reduced '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG. '
                                                                                                 '$listMutsIn{151(NOT '
                                                                                                 'IAL)} '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'unusual '
                                                                                                 'mutation '
                                                                                                 'at '
                                                                                                 'this '
                                                                                                 'position.',
                                                                                         'type': 'Accessory'}],
                                                                           'primaryType': 'Accessory',
                                                                           'text': 'F51Y'}],
                                                            'score': 15.0},
                                                           {'mutations': [{'comments': [{'text': 'T66K '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'uncommon '
                                                                                                 'non-polymorphic '
                                                                                                 'INSTI-selected '
                                                                                                 'mutation. '
                                                                                                 'It '
                                                                                                 'is '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'high-level '
                                                                                                 'EVG '
                                                                                                 'resistance, '
                                                                                                 'intermediate/high-level '
                                                                                                 'RAL '
                                                                                                 'resistance, '
                                                                                                 'and '
                                                                                                 'low-level '
                                                                                                 'DTG '
                                                                                                 'and '
                                                                                                 'CAB '
                                                                                                 'resistance. '
                                                                                                 'Its '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'BIC '
                                                                                                 'is '
                                                                                                 'not '
                                                                                                 'known.',
                                                                                         'type': 'Major'}],
                                                                           'primaryType': 'Major',
                                                                           'text': 'N66K'}],
                                                            'score': 60.0}],
                                         'score': 75.0,
                                         'text': 'High-Level Resistance'},
                                        {'drug': {'displayAbbr': 'DTG', 'name': 'DTG'},
                                         'drugClass': {'name': 'INSTI'},
                                         'level': 3,
                                         'partialScores': [{'mutations': [{'comments': [{'text': 'V151I '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'accessory '
                                                                                                 'INSTI '
                                                                                                 'selected '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'occurs '
                                                                                                 'in '
                                                                                                 '1% '
                                                                                                 'to '
                                                                                                 '3% '
                                                                                                 'of '
                                                                                                 'viruses '
                                                                                                 'from '
                                                                                                 'ART-naive '
                                                                                                 'persons '
                                                                                                 'depending '
                                                                                                 'on '
                                                                                                 'subtype. '
                                                                                                 'Alone, '
                                                                                                 'it '
                                                                                                 'appears '
                                                                                                 'to '
                                                                                                 'have '
                                                                                                 'little '
                                                                                                 'or '
                                                                                                 'no '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'INSTI '
                                                                                                 'susceptibility. '
                                                                                                 'V151L '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'reduces '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG '
                                                                                                 'by '
                                                                                                 '15 '
                                                                                                 'to '
                                                                                                 '20-fold '
                                                                                                 'and '
                                                                                                 'to '
                                                                                                 'CAB '
                                                                                                 'and '
                                                                                                 'DTG '
                                                                                                 'by '
                                                                                                 'about '
                                                                                                 '3-fold. '
                                                                                                 'V151A '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'minimally '
                                                                                                 'reduced '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG. '
                                                                                                 '$listMutsIn{151(NOT '
                                                                                                 'IAL)} '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'unusual '
                                                                                                 'mutation '
                                                                                                 'at '
                                                                                                 'this '
                                                                                                 'position.',
                                                                                         'type': 'Accessory'}],
                                                                           'primaryType': 'Accessory',
                                                                           'text': 'F51Y'}],
                                                            'score': 10.0},
                                                           {'mutations': [{'comments': [{'text': 'T66K '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'uncommon '
                                                                                                 'non-polymorphic '
                                                                                                 'INSTI-selected '
                                                                                                 'mutation. '
                                                                                                 'It '
                                                                                                 'is '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'high-level '
                                                                                                 'EVG '
                                                                                                 'resistance, '
                                                                                                 'intermediate/high-level '
                                                                                                 'RAL '
                                                                                                 'resistance, '
                                                                                                 'and '
                                                                                                 'low-level '
                                                                                                 'DTG '
                                                                                                 'and '
                                                                                                 'CAB '
                                                                                                 'resistance. '
                                                                                                 'Its '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'BIC '
                                                                                                 'is '
                                                                                                 'not '
                                                                                                 'known.',
                                                                                         'type': 'Major'}],
                                                                           'primaryType': 'Major',
                                                                           'text': 'N66K'}],
                                                            'score': 15.0}],
                                         'score': 25.0,
                                         'text': 'Low-Level Resistance'},
                                        {'drug': {'displayAbbr': 'EVG', 'name': 'EVG'},
                                         'drugClass': {'name': 'INSTI'},
                                         'level': 5,
                                         'partialScores': [{'mutations': [{'comments': [{'text': 'V151I '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'accessory '
                                                                                                 'INSTI '
                                                                                                 'selected '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'occurs '
                                                                                                 'in '
                                                                                                 '1% '
                                                                                                 'to '
                                                                                                 '3% '
                                                                                                 'of '
                                                                                                 'viruses '
                                                                                                 'from '
                                                                                                 'ART-naive '
                                                                                                 'persons '
                                                                                                 'depending '
                                                                                                 'on '
                                                                                                 'subtype. '
                                                                                                 'Alone, '
                                                                                                 'it '
                                                                                                 'appears '
                                                                                                 'to '
                                                                                                 'have '
                                                                                                 'little '
                                                                                                 'or '
                                                                                                 'no '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'INSTI '
                                                                                                 'susceptibility. '
                                                                                                 'V151L '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'reduces '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG '
                                                                                                 'by '
                                                                                                 '15 '
                                                                                                 'to '
                                                                                                 '20-fold '
                                                                                                 'and '
                                                                                                 'to '
                                                                                                 'CAB '
                                                                                                 'and '
                                                                                                 'DTG '
                                                                                                 'by '
                                                                                                 'about '
                                                                                                 '3-fold. '
                                                                                                 'V151A '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'minimally '
                                                                                                 'reduced '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG. '
                                                                                                 '$listMutsIn{151(NOT '
                                                                                                 'IAL)} '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'unusual '
                                                                                                 'mutation '
                                                                                                 'at '
                                                                                                 'this '
                                                                                                 'position.',
                                                                                         'type': 'Accessory'}],
                                                                           'primaryType': 'Accessory',
                                                                           'text': 'F51Y'}],
                                                            'score': 15.0},
                                                           {'mutations': [{'comments': [{'text': 'T66K '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'uncommon '
                                                                                                 'non-polymorphic '
                                                                                                 'INSTI-selected '
                                                                                                 'mutation. '
                                                                                                 'It '
                                                                                                 'is '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'high-level '
                                                                                                 'EVG '
                                                                                                 'resistance, '
                                                                                                 'intermediate/high-level '
                                                                                                 'RAL '
                                                                                                 'resistance, '
                                                                                                 'and '
                                                                                                 'low-level '
                                                                                                 'DTG '
                                                                                                 'and '
                                                                                                 'CAB '
                                                                                                 'resistance. '
                                                                                                 'Its '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'BIC '
                                                                                                 'is '
                                                                                                 'not '
                                                                                                 'known.',
                                                                                         'type': 'Major'}],
                                                                           'primaryType': 'Major',
                                                                           'text': 'N66K'}],
                                                            'score': 60.0}],
                                         'score': 75.0,
                                         'text': 'High-Level Resistance'},
                                        {'drug': {'displayAbbr': 'RAL', 'name': 'RAL'},
                                         'drugClass': {'name': 'INSTI'},
                                         'level':5,
                                         'partialScores': [{'mutations': [{'comments': [{'text': 'V151I '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'accessory '
                                                                                                 'INSTI '
                                                                                                 'selected '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'occurs '
                                                                                                 'in '
                                                                                                 '1% '
                                                                                                 'to '
                                                                                                 '3% '
                                                                                                 'of '
                                                                                                 'viruses '
                                                                                                 'from '
                                                                                                 'ART-naive '
                                                                                                 'persons '
                                                                                                 'depending '
                                                                                                 'on '
                                                                                                 'subtype. '
                                                                                                 'Alone, '
                                                                                                 'it '
                                                                                                 'appears '
                                                                                                 'to '
                                                                                                 'have '
                                                                                                 'little '
                                                                                                 'or '
                                                                                                 'no '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'INSTI '
                                                                                                 'susceptibility. '
                                                                                                 'V151L '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'that '
                                                                                                 'reduces '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG '
                                                                                                 'by '
                                                                                                 '15 '
                                                                                                 'to '
                                                                                                 '20-fold '
                                                                                                 'and '
                                                                                                 'to '
                                                                                                 'CAB '
                                                                                                 'and '
                                                                                                 'DTG '
                                                                                                 'by '
                                                                                                 'about '
                                                                                                 '3-fold. '
                                                                                                 'V151A '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'extremely '
                                                                                                 'rare '
                                                                                                 'mutation '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'minimally '
                                                                                                 'reduced '
                                                                                                 'susceptibility '
                                                                                                 'to '
                                                                                                 'RAL '
                                                                                                 'and '
                                                                                                 'EVG. '
                                                                                                 '$listMutsIn{151(NOT '
                                                                                                 'IAL)} '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'unusual '
                                                                                                 'mutation '
                                                                                                 'at '
                                                                                                 'this '
                                                                                                 'position.',
                                                                                         'type': 'Accessory'}],
                                                                           'primaryType': 'Accessory',
                                                                           'text': 'F51Y'}],
                                                            'score': 15.0},
                                                           {'mutations': [{'comments': [{'text': 'T66K '
                                                                                                 'is '
                                                                                                 'an '
                                                                                                 'uncommon '
                                                                                                 'non-polymorphic '
                                                                                                 'INSTI-selected '
                                                                                                 'mutation. '
                                                                                                 'It '
                                                                                                 'is '
                                                                                                 'associated '
                                                                                                 'with '
                                                                                                 'high-level '
                                                                                                 'EVG '
                                                                                                 'resistance, '
                                                                                                 'intermediate/high-level '
                                                                                                 'RAL '
                                                                                                 'resistance, '
                                                                                                 'and '
                                                                                                 'low-level '
                                                                                                 'DTG '
                                                                                                 'and '
                                                                                                 'CAB '
                                                                                                 'resistance. '
                                                                                                 'Its '
                                                                                                 'effect '
                                                                                                 'on '
                                                                                                 'BIC '
                                                                                                 'is '
                                                                                                 'not '
                                                                                                 'known.',
                                                                                         'type': 'Major'}],
                                                                           'primaryType': 'Major',
                                                                           'text': 'N66K'}],
                                                            'score': 60.0}],
                                         'score': 75.0,
                                         'text': 'High-Level Resistance'}],
                         'gene': {'name': 'IN'},
                         'version': {'publishDate': '2022-12-07', 'text': '9.4'}}
                
        res_drug_resistance = self.writer.format_drug_resistance(scores, gene)
        
        self.assertEqual(exp_drug_resistance, res_drug_resistance)

    def testFormatAlignedGeneSequences(self):

        # Setting params
        omlist = [(3, 'V', 'I', 'V'), (37, 'S', 'N', 'S')]
        gene = 'PR'
        nalist = (1, 99)
        
        exp_records = {'firstAA': 1,
                       'lastAA': 99,
                       'gene': {'name': 'PR',
                                'length': None},
                       'mutations': [{"consensus": "I",
                                      "position": 3,
                                      "AAs": "V",
                                      "isInsertion": False,
                                      "isDeletion": False,
                                      "isApobecMutation": False,
                                      "isApobecDRM": False,
                                      "isUnusual": False,
                                      "isSDRM": False,
                                      "hasStop": False,
                                      "primaryType": "Other",
                                      "text": "I3V"},
                                     {"consensus": "N",
                                      "position": 37,
                                      "AAs": "S",
                                      "isInsertion": False,
                                      "isDeletion": False,
                                      "isApobecMutation": False,
                                      "isApobecDRM": False,
                                      "isUnusual": False,
                                      "isSDRM": False,
                                      "hasStop": False,
                                      "primaryType": "Other",
                                      "text": "N37S"}],
                       'SDRMs': []}
        res_records = self.writer.format_aligned_gene_sequences(omlist, gene, nalist)
        self.assertEqual(exp_records, res_records)

        omlist = [(3, 'V', 'I', 'V'), (4, 'X', 'T', 'X'), (37, 'S', 'N', 'S')]
        gene = 'PR'
        nalist = (2, 99)  
        exp_records = {'firstAA': 2,
                       'lastAA': 99,
                       'gene': {'name': 'PR',
                                'length': None},
                       'mutations': [{"consensus": "I",
                                      "position": 3,
                                      "AAs": "V",
                                      "isInsertion": False,
                                      "isDeletion": False,
                                      "isApobecMutation": False,
                                      "isApobecDRM": False,
                                      "isUnusual": False,
                                      "isSDRM": False,
                                      "hasStop": False,
                                      "primaryType": "Other",
                                      "text": "I3V"},
                                     {'consensus': 'T',
                                      'position': 4,
                                      'AAs': 'X',
                                      "isInsertion": False,
                                      "isDeletion": False,
                                      "isApobecMutation": False,
                                      "isApobecDRM": False,
                                      "isUnusual": True,
                                      "isSDRM": False,
                                      "hasStop": False,
                                      "primaryType": "Other",
                                      'text': 'T4X'},
                                     {"consensus": "N",
                                      "position": 37,
                                      "AAs": "S",
                                      "isInsertion": False,
                                      "isDeletion": False,
                                      "isApobecMutation": False,
                                      "isApobecDRM": False,
                                      "isUnusual": False,
                                      "isSDRM": False,
                                      "hasStop": False,
                                      "primaryType": "Other",
                                      "text": "N37S"}],
                       'SDRMs': []}
        res_records = self.writer.format_aligned_gene_sequences(omlist, gene, nalist)
        self.assertEqual(exp_records, res_records)

    def testFormatInputSequence(self):
        # Setting params
        headers = ['HXB2-PR', 'shift1', 'shift2',
                  'plus1', 'plus_codon', 'del1_after_3codons',
                  'insAAA_after_3codons']
        sha512 = 'f90ddd77e400dfe6a3fcf479b0' \
              '0b1ee29e7015c5bb8cd70f5f15' \
              'b4886cc339275ff553fc8a053f' \
              '8ddc7324f45168cffaf81f8c3a' \
              'c93996f6536eef38e5e40768'
        
        for header in headers:
            self.assertEqual({'header': header, 'SHA512': sha512}, self.writer.format_input_sequence(header, ' '))

    def testWriteToJson(self):
        # Setting params
        sequence_headers = ['HXB2-PR', 'shift1', 'shift2',
                            'plus1', 'plus_codon', 'del1_after_3codons',
                            'insAAA_after_3codons']
        sequence_scores = \
                    [[{'ATV/r': (0, [], []),
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
        ordered_mutation_list = [[[(3, 'V', 'I', 'V'), (37, 'S', 'N', 'S')]],
                                 [[(1, 'X', 'P', 'X'), (3, 'V', 'I', 'V'), (37, 'S', 'N', 'S')]],
                                 [[(3, 'V', 'I', 'V'), (37, 'S', 'N', 'S')]],
                                 [[(3, 'V', 'I', 'V'), (37, 'S', 'N', 'S')]],
                                 [[(3, 'V', 'I', 'V'), (37, 'S', 'N', 'S')]],
                                 [[(1, 'X', 'P', 'X'), (2, 'S', 'Q', 'S'), (3, 'G', 'I', 'G'), (4, 'P', 'T', 'P'), (37, 'S', 'N', 'S')]],
                                 [[(1, 'Q', 'P', 'Q'), (2, 'V', 'Q', 'V'), (3, 'K', 'I', 'K'), (37, 'S', 'N', 'S')]]]
        file_genes = [[('PR', 1, 99, 1, 294)],
                      [('PR', 2, 99, 3, 293)],
                      [('PR', 2, 99, 2, 292)],
                      [('PR', 1, 99, 2, 295)],
                      [('PR', 1, 99, 1, 297)],
                      [('PR', 1, 99, 1, 293)],
                      [('PR', 1, 99, 1, 297)]]
        sequence_lengths = [[294], [291], [291], [294], [297], [293], [297]]
        file_trims = [[(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)]]
        subtypes = ['', '', '', '', '', '', '']
        na_sequence = {'HXB2-PR':'',
                       'shift1':'',
                       'shift2':'',
                       'plus1':'',
                       'plus_codon': '',
                       'del1_after_3codons': '',
                       'insAAA_after_3codons': ''}

        file_path = r'tests/hxb2-pr-local.json'
        self.assertFalse(os.path.exists(file_path))
        self.writer.write_to_json(file_path, sequence_headers, sequence_scores, file_genes,
                                  ordered_mutation_list, sequence_lengths, file_trims, subtypes,
                                  na_sequence)
        self.assertTrue(file_path)
        os.remove(file_path)

    def testValidateSequence(self):
        # Setting params
        file_genes = [[('PR', 1, 99, 1, 294)],
                      [('PR', 2, 99, 3, 293)],
                      [('PR', 2, 99, 2, 292)],
                      [('PR', 1, 99, 2, 295)],
                      [('PR', 1, 99, 1, 297)],
                      [('PR', 1, 99, 1, 293)],
                      [('PR', 1, 99, 1, 297)]]
        sequence_lengths = [[294], [291], [291], [294], [297], [293], [297]]
        file_trims = [[(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)]]
        
        for genes, sequence_length, file_trim in zip(file_genes, sequence_lengths, file_trims):
            self.assertEqual([], self.writer.validate_sequence(genes, sequence_length, file_trim))

        # Setting params
        genes = [('RT', 1, 99, 1, 294)]
        sequence_length = [194]
        file_trim = [(0, 0)]

        exp_validation = [('WARNING', "The ('RT', 1, 99, 1, 294) sequence contains just 194 codons, which is not sufficient for a comprehensive interpretation.")]
        res_validation = self.writer.validate_sequence(genes, sequence_length, file_trim)

        self.assertEqual(exp_validation, res_validation)

        # Setting params
        genes = [()]
        sequence_length = [85]
        file_trim = [(0, 0)]

        exp_validation = [('CRITICAL', 'Unable to process sequence.')]
        res_validation = self.writer.validate_sequence(genes, sequence_length, file_trim)

        self.assertEqual(exp_validation, res_validation)

        # Setting params
        genes = [('RT', 1, 99, 1, 294)]
        sequence_length = [295]
        file_trim = [(1, 1)]

        exp_validation = [('WARNING', "The ('RT', 1, 99, 1, 294) sequence had 1 amino acid trimmed from its 5\u2032-end due to poor quality."),
                          ('WARNING',"The ('RT', 1, 99, 1, 294) sequence had 1 amino acid trimmed from its 3\u2032-end due to poor quality.")]
        res_validation = self.writer.validate_sequence(genes, sequence_length, file_trim)

        self.assertEqual(exp_validation, res_validation)

    def testFindComment(self):
        # Setting params
        gene = 'IN'
        mutation = 'F51Y'
        comments = self.algorithm.comments
        details = self.algorithm.definitions['comment']

        exp_text = \
            'V151I is an accessory INSTI selected mutation that occurs in 1% to 3% of ' \
            'viruses from ART-naive persons depending on subtype. Alone, it appears to ' \
            'have little or no effect on INSTI susceptibility. V151L is an extremely rare ' \
            'mutation that reduces susceptibility to RAL and EVG by 15 to 20-fold and to ' \
            'CAB and DTG by about 3-fold. V151A is an extremely rare mutation associated ' \
            'with minimally reduced susceptibility to RAL and EVG. $listMutsIn{151(NOT ' \
            'IAL)} is an unusual mutation at this position.'

        res_text = self.writer.find_comment(gene, mutation, comments, details)

        self.assertEqual(exp_text, res_text)

    def testIsApobecDrm(self):
        # Setting params
        gene = 'IN'
        consensus = 'G'
        position = 163
        amino_acids = 'TRAG'

        result = self.writer.is_apobec_drm(gene, consensus, position, amino_acids)
        self.assertTrue(result)

        # Setting params
        gene = 'RT'
        consensus = 'D'
        position = 67
        amino_acids = 'TRAG'

        result = self.writer.is_apobec_drm(gene, consensus, position, amino_acids)
        self.assertFalse(result)

    def testIsSdrm(self):
        # Setting params
        gene = 'IN'
        consensus = 'G'
        position = 163
        amino_acids = 'TRAG'

        result = self.writer.is_sdrm(gene, position, amino_acids)
        self.assertFalse(result)

        # Setting params
        gene = 'RT'
        consensus = 'D'
        position = 67
        amino_acids = 'TRAG'

        result = self.writer.is_sdrm(gene, position, amino_acids)
        self.assertTrue(result)

    def testIsUnusual(self):
        # Setting params
        gene = 'PR'
        consensus = 'I'
        position = 3
        amino_acids = 'V'
        text = 'V'

        result = self.writer.is_unusual(gene, position, amino_acids, text)
        self.assertFalse(result)

        # Setting params
        gene = 'PR'
        position = 37
        amino_acids = 'S'
        consensus = 'N'
        text = 'S'

        result = self.writer.is_unusual(gene, position, amino_acids, text)
        self.assertFalse(result)

    def testPrimaryType(self):
        # Setting params
        gene = 'PR'
        consensus = 'I'
        position = 3
        amino_acids = 'V'
        text = 'V'

        result = self.writer.primary_type(gene, position, amino_acids)
        self.assertEqual(result, 'Other')

        # Setting params
        gene = 'PR'
        position = 37
        amino_acids = 'S'
        consensus = 'N'
        text = 'S'

        result = self.writer.primary_type(gene, position, amino_acids)
        self.assertTrue(result, 'Other')

    def testApobecMutation(self):
        # Setting params
        gene = 'PR'
        consensus = 'I'
        position = 3
        amino_acids = 'V'
        text = 'V'

        result = self.writer.is_apobec_mutation(gene, position, amino_acids)
        self.assertFalse(result)

        # Setting params
        gene = 'PR'
        position = 37
        amino_acids = 'S'
        consensus = 'N'
        text = 'S'

        result = self.writer.is_apobec_mutation(gene, position, amino_acids)
        self.assertFalse(result)


if __name__ == '__main__':
    unittest.main()