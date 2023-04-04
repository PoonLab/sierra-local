import os
import unittest
from math import inf
from xml.etree.ElementTree import parse, ParseError

from sierralocal.hivdb import HIVdb

class TestHIVDb(unittest.TestCase):
    def setUp(self):
        self.algorithm = HIVdb()
        self.root = parse(self.algorithm.xml_filename)       

    def testSetHivdbXml(self):
        exp_xml_filename = os.path.abspath(r'sierralocal\data\HIVDB_8.8.a126e04c.xml')
        res_xml_filename = os.path.abspath(self.algorithm.xml_filename)
        self.assertEqual(exp_xml_filename, res_xml_filename)

        # Setting params
        xml_path = r'sierralocal\data\hivfacts\data\algorithms\HIVDB_9.3.xml'
        algorithm = HIVdb(asi2=xml_path)

        exp_xml_filename = os.path.abspath(xml_path)
        res_xml_filename = os.path.abspath(algorithm.xml_filename)

        self.assertEqual(exp_xml_filename, res_xml_filename)

        with self.assertRaises(FileNotFoundError) as err:
            invalid_path = r'does\not\exist.txt'
            algorithm = HIVdb(asi2=invalid_path)
        
        with self.assertRaises(ParseError):
            invalid_xml = r'sierralocal\data\RT-comments.csv'
            algorithm = HIVdb(asi2=invalid_xml)

    def testSetApobecJson(self):
        exp_json_filename = os.path.abspath(r'sierralocal\data\apobec_drms.c9583ac2.json')
        res_json_filename = os.path.abspath(self.algorithm.json_filename)
        self.assertEqual(exp_json_filename, res_json_filename)

        # Setting params
        json_path = r'sierralocal\data\apobec_drms.c9583ac2.json'

        algorithm = HIVdb(apobec=json_path)
        exp_json_filename = os.path.abspath(json_path)
        res_json_filename = os.path.abspath(algorithm.json_filename)
        self.assertEqual(exp_json_filename, res_json_filename)

    def testParseDefinitions(self):
        exp_def_dict = \
            {'comment': {'IN118ACDEFHIKLMNPQSTVWY_-': {'1': 'G118R is an extremely rare '
                                                            'non-polymorphic mutation '
                                                            'selected in patients '
                                                            'receiving RAL and DTG and in '
                                                            'vitro by DTG. It causes '
                                                            'intermediate reductions in '
                                                            'RAL and EVG susceptibility '
                                                            'and low-level reductions in '
                                                            'DTG and BIC susceptibility. '
                                                            '$listMutsIn{118(NOT R)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN118R': {'1': 'G118R is a rare non-polymorphic mutation '
                                            'selected in patients receiving RAL and DTG and '
                                            'in vitro by DTG. It causes a variable reduction '
                                            'in susceptibility to each of the INSTIs of about '
                                            '5-fold.'},
                            'IN119R': {'1': 'S119R is a polymorphic mutation that is weakly '
                                            'selected by INSTIs usually in combination with '
                                            'several major INSTI-associated DRMs. Alone, it '
                                            'has little, if any effect, on INSTI '
                                            'susceptibility . '},
                            'IN121ACDEGHIKLMNPQRSTVW_-': {'1': 'F121Y is a non-polymorphic '
                                                            'mutation selected in vitro by '
                                                            'RAL and EVG. It has been '
                                                            'reported rarely in patients '
                                                            'receiving RAL. It causes '
                                                            'intermediate to high-level '
                                                            'reductions in RAL and EVG '
                                                            'susceptibility but does not '
                                                            'appear to reduce DTG or BIC '
                                                            'susceptibility. '
                                                            '$listMutsIn{121(NOT Y)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN121Y': {'1': 'F121Y is a non-polymorphic mutation selected in '
                                            'vitro by RAL and EVG. It has been reported '
                                            'rarely in patients receiving RAL. It causes '
                                            'intermediate to high-level reductions in RAL and '
                                            'EVG susceptibility but does not appear to reduce '
                                            'DTG or BIC susceptibility.'},
                            'IN128T': {'1': 'A128T is a relatively nonpolymorphic possible '
                                            'INSTI-selected mutation, which does not appear '
                                            'to reduce INSTI susceptibility '},
                            'IN138CFGHILMNPQRSVWY_-': {'1': 'E138K/A/T are non-polymorphic '
                                                            'mutations selected in patients '
                                                            'receiving RAL, EVG and rarely '
                                                            'DTG. They usually occur in '
                                                            'combination with Q148 mutations. '
                                                            'In combination with Q148 they '
                                                            'cause high-level resistance to '
                                                            'RAL and EVG and intermediate '
                                                            'reductions in susceptibility to '
                                                            'DTG and BIC. E138D is a '
                                                            'polymorphism that occurs in 1% '
                                                            'to 2% of viruses from '
                                                            'INSTI-naive patients. It does '
                                                            'not appear to be selected by '
                                                            'INSTIs or to reduce INSTI '
                                                            'susceptibility. '
                                                            '$listMutsIn{138(NOT KATD)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN138D': {'1': 'E138D is a polymorphism that occurs in 1% to 2% '
                                            'of viruses from INSTI-naive patients. It does '
                                            'not appear to be selected by INSTIs or to reduce '
                                            'INSTI susceptibility.'},
                            'IN138KAT': {'1': 'E138K/A are non-polymorphic mutations selected '
                                            'in patients receiving RAL, EVG, and DTG. They '
                                            'usually occur in combination with Q148 '
                                            'mutations. Alone they do not reduce INSTI '
                                            'susceptibility. However, when they occur in '
                                            'combination with Q148 mutations, they are '
                                            'associated with high-level resistance to RAL '
                                            'and EVG and intermediate reductions in DTG and '
                                            'BIC susceptibility. E138T is an uncommon '
                                            'nonpolymorphic INSTI-selected mutation that '
                                            'appears to have an effect similar to E138K/A.'},
                            'IN140DEFHIKLMNPQRTVWY_-': {'1': 'G140S/A/C are non-polymorphic '
                                                            'mutations that usually occur '
                                                            'with Q148 mutations. Alone, '
                                                            'they have minimal effects on '
                                                            'INSTI susceptibility. However, '
                                                            'in combination with Q148 '
                                                            'mutations they are associated '
                                                            'with high-level resistance to '
                                                            'RAL and EVG and intermediate '
                                                            'reductions in DTG and BIC '
                                                            'susceptibility. '
                                                            '$listMutsIn{140(NOT SAC)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN140SAC': {'1': 'G140S/A/C are non-polymorphic mutations that '
                                            'usually occur with Q148 mutations. Alone, they '
                                            'have minimal effects on INSTI susceptibility. '
                                            'However, in combination with Q148 mutations '
                                            'they are associated with high-level resistance '
                                            'to RAL and EVG and intermediate reductions in '
                                            'DTG and BIC susceptibility.'},
                            'IN142T': {'1': 'P142T is a rare nonpolymorphic mutation selected '
                                            'in vitro and in vivo by RAL. It usually occurs '
                                            'in combination with other INSTI-resistance '
                                            'mutations. Its effect on susceptibility has not '
                                            'been studied.'},
                            'IN143CR': {'1': 'Y143C/R are non-polymorphic mutations '
                                            'associated with high-level RAL resistance. '
                                            'Alone, they have minimal effects on EVG '
                                            'susceptibility. However, they are associated '
                                            'with intermediate to high-level reductions in '
                                            'EVG susceptibility when they occur in '
                                            'combination with one or more accessory '
                                            'INSTI-resistance mutations. Y143 mutations do '
                                            'not reduce DTG or BIC susceptibility.'},
                            'IN143DEFILMNPQTVW_-': {'1': 'Y143C/R are non-polymorphic '
                                                        'mutations associated with '
                                                        'high-level RAL resistance. Alone, '
                                                        'they have minimal effects on EVG '
                                                        'susceptibility. However, they are '
                                                        'associated with intermediate to '
                                                        'high-level reductions in EVG '
                                                        'susceptibility when they occur in '
                                                        'combination with one or more '
                                                        'accessory INSTI-resistance '
                                                        'mutations. Y143 mutations do not '
                                                        'reduce DTG or BIC susceptibility. '
                                                        '$listMutsIn{143(NOT CRHKGSA)} is an '
                                                        'unusual mutation at this position.'},
                            'IN143H': {'1': 'Y143C/R are non-polymorphic mutations associated '
                                            'with high-level RAL resistance. Alone, they have '
                                            'minimal effects on EVG susceptibility. However, '
                                            'they are associated with intermediate to '
                                            'high-level reductions in EVG susceptibility when '
                                            'they occur in combination with one or more '
                                            'accessory INSTI-resistance mutations. Y143 '
                                            'mutations do not reduce DTG or BIC '
                                            'susceptibility. Y143H is a less-common mutation '
                                            'at this position that is likely a transitional '
                                            'mutation between the wildtype Y amino acid and '
                                            'the mutant R (which differs from the wildtype Y '
                                            'by two nucleotides).'},
                            'IN143KGSA': {'1': 'Y143K/G/S/A are rare mutations that cause '
                                            'intermediate reductions in RAL '
                                            'susceptibility. When they occur in '
                                            'combination with other accessory mutations, '
                                            'they may also cause low-level reductions in '
                                            'EVG susceptibility.'},
                            'IN145ACDEFGHIKLMNQRTVWY_-': {'1': 'P145S is a rare '
                                                            'non-polymorphic mutation '
                                                            'selected in vitro by EVG and '
                                                            'rarely in patients receiving '
                                                            'EVG. It causes high-level '
                                                            'resistance to EVG but not to '
                                                            'RAL or DTG. '
                                                            '$listMutsIn{145(NOT S)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN145S': {'1': 'P145S is a rare non-polymorphic mutation '
                                            'selected in vitro by EVG and rarely in patients '
                                            'receiving EVG. It causes high-level resistance '
                                            'to EVG but not to RAL or DTG.'},
                            'IN146ACDEFGHIKLMNRSTVWY_-': {'1': 'Q146P is a rare '
                                                            'non-polymorphic mutation '
                                                            'selected in vitro by EVG '
                                                            'which causes '
                                                            'low-to-intermediate '
                                                            'reductions in EVG '
                                                            'susceptibility. '
                                                            '$listMutsIn{146(NOT P)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN146P': {'1': 'Q146P is a rare non-polymorphic mutation '
                                            'selected in vitro by EVG which causes '
                                            'low-to-intermediate reductions in EVG '
                                            'susceptibility.'},
                            'IN147ACDEFHIKLMNPQRTVWY_-': {'1': 'S147G is a non-polymorphic '
                                                            'mutation selected primarily '
                                                            'in patients receiving EVG. It '
                                                            'moderately reduces EVG '
                                                            'susceptibility. It does not '
                                                            'reduce RAL, DTG, or BIC '
                                                            'susceptibility. '
                                                            '$listMutsIn{147(NOT G)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN147G': {'1': 'S147G is a non-polymorphic mutation selected '
                                            'primarily in patients receiving EVG. It '
                                            'moderately reduces EVG susceptibility. It does '
                                            'not reduce RAL, DTG, or BIC susceptibility.'},
                            'IN148ACDEFGILMPSTVWY_-': {'1': 'Q148H/K/R are non-polymorphic '
                                                            'mutations selected by RAL, EVG, '
                                                            'and rarely DTG. Q148H/R/K are '
                                                            'associated with high-level '
                                                            'reductions in RAL and EVG '
                                                            'susceptibility particularly when '
                                                            'they occur In combination with '
                                                            'E138 or G140 mutations. Alone, '
                                                            'Q148H/K/R have minimal effects '
                                                            'on DTG and BIC susceptibility. '
                                                            'But in combination with E138 and '
                                                            'G140 mutations they cause '
                                                            'moderate and occasionally '
                                                            'high-level reductions in DTG and '
                                                            'BIC susceptibility. '
                                                            '$listMutsIn{148(NOT HKRN)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN148HKR': {'1': 'Q148H/K/R are non-polymorphic mutations '
                                            'selected by RAL, EVG, and rarely DTG. '
                                            'Q148H/R/K are associated with high-level '
                                            'reductions in RAL and EVG susceptibility '
                                            'particularly when they occur In combination '
                                            'with E138 or G140 mutations. Alone, Q148H/K/R '
                                            'have minimal effects on DTG and BIC '
                                            'susceptibility. But in combination with E138 '
                                            'and G140 mutations they cause moderate and '
                                            'occasionally high-level reductions in DTG and '
                                            'BIC susceptibility.'},
                            'IN148N': {'1': 'Q148H/K/R are non-polymorphic mutations selected '
                                            'by RAL, EVG, and rarely DTG. Q148H/R/K are '
                                            'associated with high-level reductions in RAL and '
                                            'EVG susceptibility particularly when they occur '
                                            'In combination with E138 or G140 mutations. '
                                            'Alone, Q148H/K/R have minimal effects on DTG and '
                                            'BIC susceptibility. But in combination with E138 '
                                            'and G140 mutations they cause moderate and '
                                            'occasionally high-level reductions in DTG and '
                                            'BIC susceptibility. Q148N is a rare '
                                            'INSTI-selected mutation that causes ~3-fold '
                                            'reduced EVG susceptibility and may represent a '
                                            'reversion from Q148H or Q148K.'},
                            'IN149A': {'1': 'G149A selected in vivo by DTG in RAL experienced '
                                            'patients. It appears to have no effect by itself '
                                            'but in combination with mutations at positions '
                                            '140 and 148, it reduces DTG susceptibility. '},
                            'IN151A': {'1': 'V151A is an extremely rare non-polymorphic '
                                            'mutation associated with minimally reduced '
                                            'susceptibility to RAL and EVG.'},
                            'IN151CDEFGHKMNPQRSTWY_-': {'1': 'V151I is an accessory INSTI '
                                                            'selected mutation that occurs '
                                                            'in 1% to 5% of viruses from '
                                                            'ARV-naive patients depending on '
                                                            'subtype. It appears to have '
                                                            'little or no effect on INSTI '
                                                            'susceptibility. V151L is an '
                                                            'extremely rare non-polymorphic '
                                                            'mutation that confers '
                                                            'intermediate / high-level '
                                                            'reduced susceptibility to RAL '
                                                            'and EVG and low-level reduced '
                                                            'susceptibility to DTG. V151A is '
                                                            'an extremely rare '
                                                            'non-polymorphic mutation '
                                                            'associated with minimally '
                                                            'reduced susceptibility to RAL '
                                                            'and EVG. $listMutsIn{151(NOT '
                                                            'IAL)} is an unusual mutation at '
                                                            'this position.'},
                            'IN151I': {'1': 'V151I is an accessory INSTI selected mutation '
                                            'that occurs in 1% to 5% of viruses from '
                                            'ARV-naive patients depending on subtype. Alone, '
                                            'it appears to have little or no effect on INSTI '
                                            'susceptibility.'},
                            'IN151L': {'1': 'V151L is an extremely rare non-polymorphic '
                                            'mutation that confers intermediate / high-level '
                                            'reduced susceptibility to RAL and EVG and '
                                            'low-level reduced susceptibility to DTG.'},
                            'IN153ACDEGHIKLMNPQRTVW_-': {'1': 'S153Y/F are rare '
                                                            'non-polymorphic mutations '
                                                            'selected in vitro by EVG, DTG, '
                                                            'and BIC. S153Y/F reduce RAL, '
                                                            'DTG, and possibly BIC '
                                                            'susceptibility about 2-fold '
                                                            'and EVG susceptibility about '
                                                            '4-fold. $listMutsIn{153(NOT '
                                                            'YF)} is an unusual mutation at '
                                                            'this position.'},
                            'IN153YF': {'1': 'S153Y/F are rare non-polymorphic mutations '
                                            'selected in vitro by EVG, DTG, and BIC. S153Y/F '
                                            'reduce RAL, DTG, and possibly BIC '
                                            'susceptibility about 2-fold and EVG '
                                            'susceptibility about 4-fold.'},
                            'IN155ACDEFGIKLMPQRVWY_-': {'1': 'N155H is a non-polymorphic '
                                                            'mutation selected in patients '
                                                            'receiving RAL, EVG, and rarely '
                                                            'DTG. It is associated with '
                                                            'high-level reductions in RAL '
                                                            'and EVG susceptibility. It '
                                                            'causes low-level reductions in '
                                                            'DTG susceptibility. N155S/T are '
                                                            'rare non-polymorphic mutations '
                                                            'selected in vitro by '
                                                            'investigational INSTIs. They '
                                                            'reduce RAL and EVG '
                                                            'susceptibility somewhat less '
                                                            'than does N155H. '
                                                            '$listMutsIn{155(NOT HST)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN155H': {'1': 'N155H is a non-polymorphic mutation selected in '
                                            'patients receiving RAL, EVG, and rarely DTG. It '
                                            'is associated with high-level reductions in RAL '
                                            'and EVG susceptibility. It causes low-level '
                                            'reductions in DTG susceptibility.'},
                            'IN155ST': {'1': 'N155H is a non-polymorphic mutation selected in '
                                            'patients receiving RAL, EVG, and rarely DTG. It '
                                            'is associated with high-level reductions in RAL '
                                            'and EVG susceptibility. It causes low-level '
                                            'reductions in DTG susceptibility. N155S/T are '
                                            'extremely rare non-polymorphic mutations '
                                            'selected in vitro by investigational INSTIs. '
                                            'They reduce RAL and EVG susceptibility somewhat '
                                            'less than does N155H.'},
                            'IN157Q': {'1': 'E157Q is a polymorphic mutation selected in '
                                            'patients receiving RAL and EVG. It appears to '
                                            'have little effect on INSTI susceptibility.'},
                            'IN163RK': {'1': 'G163R/K are polymorphic in subtype F viruses '
                                            'from ARV-naive patients but are otherwise '
                                            'non-polymorphic. They are common INSTI-selected '
                                            'mutations. Alone, they have little, if any, '
                                            'effect on INSTI susceptibility.'},
                            'IN230N': {'1': 'S230N is a polymorphism that is not associated '
                                            'with reduced INSTI susceptibility.'},
                            'IN230R': {'1': 'S230R is a non-polymorphic mutation selected by '
                                            'RAL, EVG, and DTG. It causes low-level '
                                            'reductions in DTG susceptibility.'},
                            'IN232N': {'1': 'D232N is a common nonpolymorphic accessory '
                                            'mutation selected in patients receiving RAL and '
                                            'EVG.'},
                            'IN263ACDEFGHILMNPQSTVWY_-': {'1': 'R263K is selected in vitro by '
                                                            'EVG, DTG, and BIC, and in '
                                                            'patients receiving DTG. It '
                                                            'reduces DTG and BIC '
                                                            'susceptibility about 2-fold '
                                                            'and EVG susceptibility '
                                                            'somewhat more. '
                                                            '$listMutsIn{263(NOT K)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN263K': {'1': 'R263K is selected in vitro by EVG, DTG, and BIC, '
                                            'and in patients receiving DTG. It reduces DTG '
                                            'and BIC susceptibility about 2-fold and EVG '
                                            'susceptibility somewhat more. '},
                            'IN50I': {'1': 'M50I is a polymorphic mutation selected in vitro '
                                        'by DTG and BIC in combination with R263K. It '
                                        'appears to contribute to reduced DTG '
                                        'susceptibility in combination with R263K. '},
                            'IN51ACDEFGIKLMNPQRSTVW_-': {'1': 'H51Y is a rare non-polymorphic '
                                                            'accessory mutation selected in '
                                                            'patients receiving RAL and EVG '
                                                            'and in vitro by DTG. H51Y '
                                                            'reduces EVG susceptibility 2 '
                                                            'to 3-fold. It does not reduce '
                                                            'RAL or DTG susceptibility. '
                                                            '$listMutsIn{51(NOT Y)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN51Y': {'1': 'H51Y is a rare non-polymorphic accessory mutation '
                                        'selected in patients receiving RAL and EVG and in '
                                        'vitro by DTG. H51Y minimally reduces EVG '
                                        'susceptibility (~2 to 3-fold). It does not reduce '
                                        'RAL or DTG susceptibility.'},
                            'IN66A': {'1': 'T66A is a non-polymorphic mutation selected in '
                                        'patients receiving EVG and RAL, usually in '
                                        'combination with other INSTI-resistance '
                                        'mutations. It causes a moderate reduction in EVG '
                                        'susceptibility but does not appear to reduce RAL, '
                                        'DTG, or BIC susceptibility.'},
                            'IN66CDEFGHLMNPQRSVWY_-': {'1': 'T66A/I/K are non-polymorphic '
                                                            'mutations selected primarily by '
                                                            'EVG. T66K is associated with '
                                                            'high-level EVG resistance, '
                                                            'intermediate/high-level RAL '
                                                            'resistance, and low-level DTG '
                                                            'resistance. Its effect on BIC is '
                                                            'not known. $listMutsIn{66(NOT '
                                                            'IAK)} is an unusual mutation at '
                                                            'this position.'},
                            'IN66I': {'1': 'T66I is a non-polymorphic mutation selected in '
                                        'patients receiving EVG, RAL, and DTG. It reduces '
                                        'EVG susceptibility about 10-fold but does not '
                                        'reduce RAL, DTG, or BIC susceptibility.'},
                            'IN66K': {'1': 'T66K is a non-polymorphic mutation selected in '
                                        'patients receiving EVG. It is associated with '
                                        'high-level EVG resistance, '
                                        'intermediate/high-level RAL resistance, and '
                                        'low-level DTG resistance. Its effect on BIC is '
                                        'not known.'},
                            'IN74MIF': {'1': 'L74M/I are polymorphic accessory mutations '
                                            'commonly selected by each of the INSTIs. In '
                                            'ARV-naive patients, L74M occurs in 0.5% to 10% '
                                            'of patients and L74I occurs in 3% to 20% of '
                                            'patients depending on subtype. Alone, L74M/I '
                                            'have minimal, if any, effect on INSTI '
                                            'susceptibility. However, they contribute '
                                            'reduced susceptibility to each of the INSTIs '
                                            'when they occur with major INSTI-resistance '
                                            'mutations. L74F is a rare nonpolymorphic '
                                            'mutation which also contributes reduced '
                                            'susceptibility when it occurs with other '
                                            'INSTI-resistance mutations.'},
                            'IN92ACDFHIKLMNPRSTWY_-': {'1': 'E92Q is a common non-polymorphic '
                                                            'mutation selected in patients '
                                                            'receiving RAL and EVG. It '
                                                            'reduces RAL susceptibility 5 to '
                                                            '10-fold and EVG susceptibility '
                                                            '~30-fold. It is selected in '
                                                            'vitro by DTG and reduces DTG '
                                                            'susceptibility ~1.5-fold. '
                                                            '$listMutsIn{92(NOT GQV)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'IN92G': {'1': 'E92G is a rare non-polymorphic mutation that has '
                                        'been selected in patients receiving EVG. It '
                                        'moderately reduces EVG susceptibility but does '
                                        'not reduce susceptibility to RAL, DTG, or BIC.'},
                            'IN92Q': {'1': 'E92Q is a common non-polymorphic mutation '
                                        'selected in patients receiving RAL and EVG. It '
                                        'reduces RAL susceptibility 5 to 10-fold and EVG '
                                        'susceptibility ~30-fold. It is selected in vitro '
                                        'by DTG and reduces DTG susceptibility ~1.5-fold. '
                                        'It does not appear to reduce BIC susceptibility.'},
                            'IN92V': {'1': 'E92V is a rare non-polymorphic mutation selected '
                                        'in vitro by an investigational INSTI. It causes '
                                        'intermediate and high-level reductions in EVG and '
                                        'RAL susceptibility, respectively.'},
                            'IN95K': {'1': 'Q95K is a non-polymorphic INSTI-selected '
                                        'mutation. Alone, it has little if any effect on '
                                        'INSTI susceptibility.'},
                            'IN97A': {'1': 'T97A is a polymorphic INSTI-selected mutation '
                                        'that, depending on subtype, occurs in 1% to 5% of '
                                        'viruses from untreated persons. Alone, it has '
                                        'minimal effects on INSTI susceptibility but in '
                                        'combination with other major resistance '
                                        'mutations, it synergistically reduces '
                                        'susceptibility to EVG, RAL, DTG, and possibly '
                                        'BIC.'},
                            'PR10ACDEGHKMNPQSTW_-': {'1': 'L10F/I/V/R/Y are PI-selected '
                                                        'accessory mutations. '
                                                        '$listMutsIn{10(NOT FIVRY)} is a '
                                                        'highly unusual mutation at this '
                                                        'position.'},
                            'PR10F': {'1': 'L10F is a common non-polymorphic, PI-selected '
                                        'accessory mutation associated with reduced '
                                        'susceptibility to DRV, FPV, IDV, LPV, and NFV.'},
                            'PR10IV': {'1': 'L10I/V are polymorphic, PI-selected accessory '
                                            'mutations that increase the replication of '
                                            'viruses with other PI-resistance mutations.'},
                            'PR10RY': {'1': 'L10R/Y are rare, non-polymorphic PI-selected '
                                            'mutations. Their effects on PI susceptibility '
                                            'have not been well studied.'},
                            'PR11IL': {'1': 'V11I is a relatively non-polymorphic accessory '
                                            'mutation selected in patients receiving DRV. It '
                                            'is included in the Tibotec DRV genotypic '
                                            'susceptibility score. V11L is a nonpolymorphic '
                                            'PI-selected mutation associated with reduced DRV '
                                            'and FPV susceptibility when it occurs in '
                                            'combination with other PI-resistance '
                                            'mutations.\n'},
                            'PR20I': {'1': 'K20I is the consensus amino acid in subtype G and '
                                        'CRF02_AG. In subtypes B and C, K20I is a '
                                        'PI-selected accessory mutation that reduces NFV '
                                        'susceptibility.'},
                            'PR20MV': {'1': 'K20M/V are rare, relatively non-polymorphic '
                                            'PI-selected mutations that have not been well '
                                            'studied.'},
                            'PR20R': {'1': 'K20R is a highly polymorphic PI-selected '
                                        'accessory mutation.'},
                            'PR20T': {'1': 'K20T is a non-polymorphic accessory PI-selected '
                                        'mutation associated with reduced susceptibility '
                                        'to each of the PIs except DRV and TPV.'},
                            'PR23I': {'1': 'L23I is an uncommon non-polymorphic mutation '
                                        'selected primarily by NFV. It causes low-level '
                                        'NFV resistance.'},
                            'PR24ACDEGHKNPQRSTVWY_-': {'1': 'L24I is a non-polymorphic '
                                                            'mutation selected by IDV and '
                                                            'LPV. It contributes reduced '
                                                            'susceptibility to each PI except '
                                                            'DRV and TPV. $listMutsIn{24(NOT '
                                                            'FIM)} is a highly unusual '
                                                            'mutation at this position.'},
                            'PR24FM': {'1': 'L24I is a non-polymorphic mutation selected by '
                                            'IDV and LPV. It contributes reduced '
                                            'susceptibility to each PI except DRV and TPV. '
                                            'L24F/M are uncommon non-polymorphic PI-selected '
                                            'mutations. L24F has a susceptibility profile '
                                            'similar to L24I.'},
                            'PR24I': {'1': 'L24I is a non-polymorphic mutation selected by '
                                        'IDV and LPV. It contributes reduced '
                                        'susceptibility to each PI except DRV and TPV.'},
                            'PR30ACEFGHIKLMPQRSTVWY_-': {'1': 'D30N is a non-polymorphic '
                                                            'mutation that causes '
                                                            'high-level resistance to NFV. '
                                                            '$listMutsIn{30(NOT N)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'PR30N': {'1': 'D30N is a non-polymorphic mutation that causes '
                                        'high-level resistance to NFV.'},
                            'PR32ACDEFGHKLMNPQRSTWY_-': {'1': 'V32I is a non-polymorphic '
                                                            'PI-selected mutation '
                                                            'associated with reduced '
                                                            'susceptibility to each of the '
                                                            'PIs except SQV. It is included '
                                                            'in the Tibotec DRV genotypic '
                                                            'susceptibility score. '
                                                            '$listMutsIn{32(NOT I)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'PR32I': {'1': 'V32I is a non-polymorphic PI-selected mutation '
                                        'associated with reduced susceptibility to each of '
                                        'the PIs except SQV. It is included in the Tibotec '
                                        'DRV genotypic susceptibility score.'},
                            'PR33F': {'1': 'L33F is a relatively non-polymorphic accessory '
                                        'mutation selected by each of the PIs. In '
                                        'combination with other PI-resistance mutations, '
                                        'it is associated with reduced susceptibility to '
                                        'each of the PIs. It is included in the Tibotec '
                                        'DRV genotypic susceptibility score.'},
                            'PR33_': {'1': 'Insertions between positions 33 and 41 do not '
                                        'appear to be selected by PIs or to reduce PI '
                                        'susceptibility.'},
                            'PR34_': {'1': 'Insertions between positions 33 and 41 do not '
                                        'appear to be selected by PIs or to reduce PI '
                                        'susceptibility.'},
                            'PR35_': {'1': 'Insertions between positions 33 and 41 do not '
                                        'appear to be selected by PIs or to reduce PI '
                                        'susceptibility.'},
                            'PR36_': {'1': 'Insertions between positions 33 and 41 do not '
                                        'appear to be selected by PIs or to reduce PI '
                                        'susceptibility.'},
                            'PR37_': {'1': 'Insertions between positions 33 and 41 do not '
                                        'appear to be selected by PIs or to reduce PI '
                                        'susceptibility.'},
                            'PR38_': {'1': 'Insertions between positions 33 and 41 do not '
                                        'appear to be selected by PIs or to reduce PI '
                                        'susceptibility.'},
                            'PR39_': {'1': 'Insertions between positions 33 and 41 do not '
                                        'appear to be selected by PIs or to reduce PI '
                                        'susceptibility.'},
                            'PR40_': {'1': 'Insertions between positions 33 and 41 do not '
                                        'appear to be selected by PIs or to reduce PI '
                                        'susceptibility.'},
                            'PR41_': {'1': 'Insertions between positions 33 and 41 do not '
                                        'appear to be selected by PIs or to reduce PI '
                                        'susceptibility.'},
                            'PR43T': {'1': 'K43T is a non-polymorphic PI-selected accessory '
                                        'mutation. K43T is included in the '
                                        'Boehringer-Ingelheim TPV genotypic susceptibility '
                                        'score.'},
                            'PR46ACDEFGHKNPQRSTWY_-': {'1': 'M46I/L are relatively '
                                                            'non-polymorphic PI-selected '
                                                            'mutations. In combination with '
                                                            'other PI-resistance mutations, '
                                                            'they are associated with reduced '
                                                            'susceptibility to each of the '
                                                            'PIs except DRV. '
                                                            '$listMutsIn{46(NOT ILV)} is a '
                                                            'highly unusual mutation at this '
                                                            'position.'},
                            'PR46IL': {'1': 'M46I/L are relatively non-polymorphic '
                                            'PI-selected mutations. In combination with other '
                                            'PI-resistance mutations, they are associated '
                                            'with reduced susceptibility to each of the PIs '
                                            'except DRV.'},
                            'PR46V': {'1': 'M46I/L are relatively non-polymorphic PI-selected '
                                        'mutations. In combination with other '
                                        'PI-resistance mutations, they are associated with '
                                        'reduced susceptibility to each of the PIs except '
                                        'DRV. M46V is a rare non-polymorphic PI-selected '
                                        'mutation that has not been well studied.\n'},
                            'PR47A': {'1': 'I47A is a non-polymorphic mutation selected by '
                                        'LPV. It usually occurs in combination with V32I '
                                        'and in this context it confers high-level '
                                        'resistance to LPV and FPV and '
                                        'low/intermediate-resistance to the remaining PIs '
                                        'except ATV and SQV.'},
                            'PR47CDEFGHKLMNPQRSTWY_-': {'1': 'I47V is a non-polymorphic '
                                                            'PI-selected mutation associated '
                                                            'with reduced susceptibility to '
                                                            'each of the PIs except SQV and '
                                                            'ATV. I47A is a non-polymorphic '
                                                            'mutation selected by LPV. It '
                                                            'usually occurs in combination '
                                                            'with V32I and in this context '
                                                            'it confers high-level '
                                                            'resistance to LPV and FPV and '
                                                            'low/intermediate-resistance to '
                                                            'the remaining PIs except ATV '
                                                            'and SQV. $listMutsIn{47(NOT '
                                                            'AV)} is a highly unusual '
                                                            'mutation at this position.'},
                            'PR47V': {'1': 'I47V is a non-polymorphic PI-selected mutation '
                                        'associated with reduced susceptibility to each of '
                                        'the PIs except SQV and ATV. I47V is included in '
                                        'the Tibotec DRV genotypic susceptibility score.'},
                            'PR48ASTQL': {'1': 'G48V is a non-polymorphic mutation selected '
                                            'by SQV and, less often, by IDV and LPV. It '
                                            'confers high-level resistance to SQV, '
                                            'intermediate resistance to ATV, and low-level '
                                            'resistance to NFV, IDV and LPV. G48M is a '
                                            'non-polymorphic mutation selected in viruses '
                                            'with multiple PI-resistance mutations. Its '
                                            'effects appear to be similar to those of '
                                            'G48V. G48A/S/T/Q are rare non-polymorphic '
                                            'mutations selected in viruses with multiple '
                                            'PI-resistance mutations.'},
                            'PR48CDEFHIKNPRWY_-': {'1': 'G48V is a non-polymorphic mutation '
                                                        'selected by SQV and, less often, by '
                                                        'IDV and LPV. It confers high-level '
                                                        'resistance to SQV, '
                                                        'intermediate-level resistance to '
                                                        'ATV, and low-level resistance to '
                                                        'NFV, IDV and LPV. G48M is a '
                                                        'non-polymorphic mutation selected in '
                                                        'viruses with multiple PI-resistance '
                                                        'mutations. Its effects appear to be '
                                                        'similar to those of G48V. G48A/S/T/Q '
                                                        'are rare non-polymorphic mutations '
                                                        'selected in viruses with multiple '
                                                        'PI-resistance mutations. '
                                                        '$listMutsIn{48(NOT VMALSTQL)} is a '
                                                        'highly unusual mutation at this '
                                                        'position.'},
                            'PR48M': {'1': 'G48V is a non-polymorphic mutation selected by '
                                        'SQV and, less often, by IDV and LPV. It confers '
                                        'high-level resistance to SQV, intermediate '
                                        'resistance to ATV, and low-level resistance to '
                                        'NFV, IDV and LPV. G48M is a non-polymorphic '
                                        'mutation selected in viruses with multiple '
                                        'PI-resistance mutations. Its effects appear to be '
                                        'similar to those of G48V.'},
                            'PR48V': {'1': 'G48V is a non-polymorphic mutation selected by '
                                        'SQV and, less often, by IDV and LPV. It confers '
                                        'high-level resistance to SQV, intermediate '
                                        'resistance to ATV, and low-level resistance to '
                                        'NFV, IDV and LPV.'},
                            'PR50ACDEFGHKMNPQRSTWY_-': {'1': 'I50V is a non-polymorphic '
                                                            'mutation selected by DRV, LPV '
                                                            'and FPV. It causes low-level, '
                                                            'intermediate, and high-level '
                                                            'resistance to DRV, LPV, and '
                                                            'FPV, respectively. I50L is a '
                                                            'non-polymorphic mutation '
                                                            'selected by ATV. It causes '
                                                            'high-level resistance to ATV '
                                                            'and increases susceptibility to '
                                                            'the remaining PIs. '
                                                            '$listMutsIn{50(NOT LV)} is a '
                                                            'highly unusual mutation at this '
                                                            'position.'},
                            'PR50L': {'1': 'I50L is a non-polymorphic mutation selected by '
                                        'ATV. It causes high-level resistance to ATV and '
                                        'increases susceptibility to most of the remaining '
                                        'PIs.'},
                            'PR50V': {'1': 'I50V is a non-polymorphic mutation selected by '
                                        'DRV, LPV and FPV. It causes low-level, '
                                        'intermediate, and high-level resistance to DRV, '
                                        'LPV, and FPV, respectively. It is included in the '
                                        'Tibotec DRV genotypic susceptibility score.'},
                            'PR53L': {'1': 'F53L is a non-polymorphic accessory PI-selected '
                                        'mutation that reduces susceptibility primarily to '
                                        'ATV, SQV, and NFV.'},
                            'PR53Y': {'1': 'F53L is a non-polymorphic accessory PI-selected '
                                        'mutation which reduces susceptibility primarily '
                                        'to ATV, SQV, and NFV. F53Y is a rare '
                                        'non-polymorphic PI-selected mutation that has not '
                                        'been well studied. $listMutsIn{53(NOT LY)} is a '
                                        'highly unusual mutation at this position.'},
                            'PR54ATS': {'1': 'I54A/T/S are non-polymorphic PI-selected '
                                            'mutations that occur almost exclusively in '
                                            'viruses with multiple PI-resistance mutations. '
                                            'I54A/T/S are associated with reduced '
                                            'susceptibility to each of the PIs except DRV.'},
                            'PR54CDEFGHKNPQRWY_-': {'1': 'I54V/A/T/S reduce susceptibility to '
                                                        'each of the PIs except DRV. I54M '
                                                        'reduces susceptibility to each of '
                                                        'the PIs. I54L reduces '
                                                        'susceptibility to each of the PIs '
                                                        'except TPV. $listMutsIn{54(NOT '
                                                        'ALMSTV)} is a highly unusual '
                                                        'mutation at this position.'},
                            'PR54LM': {'1': 'I54M/L are non-polymorphic mutations selected '
                                            'primarily by DRV and FPV. I54M reduces '
                                            'susceptibility to each of the PIs. I54L reduces '
                                            'susceptibility to each of the PIs except TPV. '
                                            'I54M/L are each included in the Tibotec DRV '
                                            'genotypic susceptibility score.'},
                            'PR54V': {'1': 'I54V is a non-polymorphic PI-selected mutation '
                                        'that contributes reduced susceptibility to each '
                                        'of the PIs except DRV.'},
                            'PR58E': {'1': 'Q58E is a non-polymorphic accessory PI-selected '
                                        'mutation associated with reduced susceptibility '
                                        'to TPV and possibly other PIs.'},
                            'PR71IL': {'1': 'A71I/L are non-polymorphic, PI-selected '
                                            'accessory mutations that appear to increase the '
                                            'replication of viruses with other PI-resistance '
                                            'mutations.'},
                            'PR71TV': {'1': 'A71V/T are polymorphic, PI-selected accessory '
                                            'mutations that increase the replication of '
                                            'viruses with other PI-resistance mutations.'},
                            'PR73EFHIKLMNPQRWY_-': {'1': 'G73S/T/C/A/D/V are non-polymorphic '
                                                        'accessory PI-selected mutations. '
                                                        'They are associated primarily with '
                                                        'reduced susceptibility to ATV, SQV, '
                                                        'IDV, and NFV. $listMutsIn{73(NOT '
                                                        'ACDSTV)} is an unusual mutation at '
                                                        'this position that has not been '
                                                        'well studied.'},
                            'PR73STCADV': {'1': 'G73S/T/C/A are non-polymorphic accessory '
                                                'PI-selected mutations. They are associated '
                                                'primarily with reduced susceptibility to '
                                                'ATV, SQV, FPV, IDV, and NFV. G73V/D are rare '
                                                'non-polymorphic PI-selected mutations.'},
                            'PR74P': {'1': 'T74P is a non-polymorphic PI-selected accessory '
                                        'mutation that occurs primarily in viruses from '
                                        'patients who have received multiple PIs. It is '
                                        'associated with reduced susceptibility to each of '
                                        'the PIs. It is included in the '
                                        'Boehringer-Ingelheim TPV and Tibotec DRV '
                                        'genotypic susceptibility scores.'},
                            'PR74S': {'1': 'T74S is a PI-selected accessory mutation that is '
                                        'polymorphic in most non-B subtypes.'},
                            'PR76ACDEFGHIKMNPQRSTWY_-': {'1': 'L76V is a non-polymorphic '
                                                            'mutation selected by IDV, LPV '
                                                            'and DRV. It reduces '
                                                            'susceptibility to these PIs '
                                                            'and to FPV and NFV. It '
                                                            'increases susceptibility to '
                                                            'ATV, SQV and TPV. L76V is '
                                                            'included in the Tibotec DRV '
                                                            'genotypic susceptibility '
                                                            'score. $listMutsIn{76(NOT V)} '
                                                            'is a highly unusual mutation '
                                                            'at this position.'},
                            'PR76V': {'1': 'L76V is a non-polymorphic mutation selected by '
                                        'IDV, LPV and DRV. It reduces susceptibility to '
                                        'these PIs and to FPV and NFV. It increases '
                                        'susceptibility to ATV, SQV and TPV. L76V is '
                                        'included in the Tibotec DRV genotypic '
                                        'susceptibility score.'},
                            'PR82A': {'1': 'V82A is a non-polymorphic mutation selected '
                                        'primarily by IDV and LPV. It reduces '
                                        'susceptibility to these PIs and contributes '
                                        'cross-resistance to each of the remaining PIs '
                                        'except DRV and TPV.'},
                            'PR82C': {'1': 'V82C is an uncommon non-polymorphic mutation that '
                                        'occurs primarily in viruses with multiple other '
                                        'PI-resistance mutations. Its effects on PI '
                                        'susceptibility have not been well studied.'},
                            'PR82DEGHKNPQRWY_-': {'1': 'V82A/T/S/F/C/M are non-polymorphic '
                                                    'mutations selected primarily in '
                                                    'patients who have received IDV, LPV, '
                                                    'or multiple PIs. V82A/S/T are '
                                                    'associated with reduced '
                                                    'susceptibility to each of the PIs '
                                                    'except DRV. V82F is associated with '
                                                    'reduced susceptibility to each of the '
                                                    'PIs. V82C/M usually occur in viruses '
                                                    'with multiple additional '
                                                    'PI-resistance mutations and have not '
                                                    'been well studied. $listMutsIn{82(NOT '
                                                    'ATIFLMSC)} is a highly unusual '
                                                    'mutation at this position.'},
                            'PR82F': {'1': 'V82F is a non-polymorphic mutation selected '
                                        'primarily by IDV and LPV. It reduces '
                                        'susceptibility to these PIs and contributes '
                                        'cross-resistance to each of the remaining PIs.'},
                            'PR82I': {'1': 'V82I is a highly polymorphic mutation that is not '
                                        'selected by PIs. It is the consensus amino acid '
                                        'in subtype G viruses.'},
                            'PR82L': {'1': 'V82L is an uncommon non-polymorphic mutation '
                                        'selected by TPV. It reduces TPV susceptibility '
                                        'but its effects on other PIs have not been well '
                                        'studied.'},
                            'PR82M': {'1': 'In most subtypes, V82M is a 2-base-pair mutation '
                                        'that develops in viruses with multiple other '
                                        'PI-resistance mutations. In subtype G, V82M is a '
                                        '1-base-pair mutation. V82M reduces susceptibility '
                                        'to IDV, LPV and possibly other PIs.'},
                            'PR82TS': {'1': 'V82T/S are non-polymorphic PI-selected '
                                            'mutations. They are associated with reduced '
                                            'susceptibility to each of the PIs except DRV. '
                                            'V82T is included in the Boehringer-Ingelheim TPV '
                                            'genotypic susceptibility score.'},
                            'PR83D': {'1': 'N83D is a non-polymorphic mutation selected '
                                        'primarily in patients who have received multiple '
                                        'PIs. It is included in the Boehringer-Ingelheim '
                                        'genotypic susceptibility score for TPV.'},
                            'PR84AC': {'1': 'I84A is an extremely rare non-polymorphic '
                                            'PI-selected substrate-cleft mutation associated '
                                            'with high-level resistance to each of the PIs '
                                            'except possibly DRV. I84C is an extremely rare '
                                            'non-polymorphic PI-selected mutation associated '
                                            'with varying degrees of reduced susceptibility '
                                            'to each of the PIs.'},
                            'PR84DEFGHKLMNPQRSTWY_-': {'1': 'I84V/A/C are non-polymorphic '
                                                            'PI-selected substrate-cleft '
                                                            'mutations associated with '
                                                            'varying degrees of reduced '
                                                            'susceptibility to each of the '
                                                            'PIs. $listMutsIn{84(NOT ACV)} is '
                                                            'a highly unusual mutation at '
                                                            'this position.'},
                            'PR84V': {'1': 'I84V is a non-polymorphic mutation selected by '
                                        'each of the PIs. It causes high-level resistance '
                                        'to ATV, FPV, IDV, NFV and SQV, intermediate '
                                        'resistance to LPV and TPV, and low-level '
                                        'resistance to DRV.'},
                            'PR85V': {'1': 'I85V is a non-polymorphic PI-selected mutation. '
                                        'It has minimal, if any, effects on PI '
                                        'susceptibility.'},
                            'PR88ACEFHIKLMPQRVWY_-': {'1': 'N88S is a non-polymorphic '
                                                        'mutation usually selected by NFV, '
                                                        'ATV, or IDV which causes '
                                                        'high-level resistance to NFV and '
                                                        'ATV and low-level resistance to '
                                                        'IDV and SQV. It increases '
                                                        'susceptibility to FPV. N88D is a '
                                                        'non-polymorphic mutation selected '
                                                        'by NFV usually in combination '
                                                        'with D30N. $listMutsIn{88(NOT '
                                                        'DGST)} is a highly unusual '
                                                        'mutation at this position.'},
                            'PR88D': {'1': 'N88D is a non-polymorphic mutation selected by '
                                        'NFV usually in combination with D30N. It reduces '
                                        'NFV susceptibility and may cause low-level '
                                        'cross-resistance to ATV and SQV.'},
                            'PR88S': {'1': 'N88S is a non-polymorphic mutation usually '
                                        'selected by NFV, ATV, or IDV. It causes '
                                        'high-level resistance to NFV and ATV and '
                                        'low-level resistance to IDV and SQV. It increases '
                                        'susceptibility to DRV and FPV.'},
                            'PR88TG': {'1': 'N88G/T are extremely rare non-polymorphic '
                                            'PI-selected mutations that reduce susceptibility '
                                            'to NFV and ATV.'},
                            'PR89VT': {'1': 'L89V is a non-polymorphic PI-selected accessory '
                                            'mutation that contributes reduced susceptibility '
                                            'to FPV, DRV, NFV, and IDV. L89V is included in '
                                            'the Tibotec DRV genotypic susceptibility score. '
                                            'L89T is a rare non-polymorphic PI-selected '
                                            'mutation that has not been well studied.'},
                            'PR90ACDEFGHIKNPQRSTVWY_-': {'1': 'L90M is a non-polymorphic '
                                                            'PI-selected mutation which '
                                                            'reduces susceptibility to each '
                                                            'of the PIs except TPV and DRV. '
                                                            '$listMutsIn{90(NOT M)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'PR90M': {'1': 'L90M is a non-polymorphic PI-selected mutation '
                                        'that reduces susceptibility to each of the PIs '
                                        'except TPV and DRV.'},
                            'RT100ACDEFGHKMNPQRSTWY_-': {'1': 'L100I is a non-polymorphic '
                                                            'mutation that usually occurs '
                                                            'in combination with K103N. In '
                                                            'this setting it causes '
                                                            'high-level resistance to NVP, '
                                                            'EFV and RPV and intermediate '
                                                            'resistance to ETR. It has a '
                                                            'weight of 2.5 in the Tibotec '
                                                            'ETR genotypic susceptibility '
                                                            'score. $listMutsIn{100(NOT '
                                                            'IV)} is a highly unusual '
                                                            'mutation at this position.'},
                            'RT100I': {'1': 'L100I is a non-polymorphic mutation that usually '
                                            'occurs in combination with K103N. In this '
                                            'setting it causes high-level resistance to NVP, '
                                            'EFV, and RPV and intermediate resistance to ETR. '
                                            'It has a weight of 2.5 in the Tibotec ETR '
                                            'genotypic susceptibility score. Preliminary data '
                                            'suggests it is contributes low to intermediate '
                                            'reductions in DOR susceptibility when it occurs '
                                            'in combination with K103N.'},
                            'RT100V': {'1': 'L100I is a non-polymorphic mutation that usually '
                                            'occurs with K103N. In this setting it causes '
                                            'high-level resistance to NVP, EFV, and RPV and '
                                            'intermediate resistance to ETR. L100V is a rare '
                                            'non-polymorphic mutation with low / intermediate '
                                            'resistance to EFV, NVP, and RPV.'},
                            'RT101CDFGILMSVWY_-': {'1': 'K101E/H/P are NNRTI-resistance '
                                                        'mutations. $listMutsIn{101(NOT '
                                                        'AEHNPQRT)} is a highly unusual '
                                                        'mutation at this position.'},
                            'RT101E': {'1': 'K101E is a non-polymorphic primarily accessory '
                                            'mutation that causes intermediate resistance to '
                                            'NVP and RPV, low-level resistance to EFV, and '
                                            'potentially low-level resistance to ETR. It has '
                                            'a weight of 1.0 in the Tibotec ETR genotypic '
                                            'susceptibility score. It is associated with '
                                            'low-level reductions in DOR susceptibility.'},
                            'RT101H': {'1': 'K101H is a non-polymorphic accessory mutation '
                                            'selected by NVP, EFV and ETR. When present with '
                                            'other NNRTI-resistance mutations, K101H further '
                                            'reduces susceptibility to these NNRTIs. It has a '
                                            'weight of 1.0 in the Tibotec ETR genotypic '
                                            'susceptibility score. Its effect on DOR '
                                            'susceptibility is not known.'},
                            'RT101NAT': {'1': 'K101N/A/T are uncommon non-polymorphic '
                                            'NNRTI-selected mutation of uncertain '
                                            'phenotypic and clinical significance.'},
                            'RT101P': {'1': 'K101P is a non-polymorphic mutation that causes '
                                            'high-level resistance to each of the NNRTIs. It '
                                            'has a weight of 2.5 in the Tibotec ETR genotypic '
                                            'susceptibility score. Its does not appear to '
                                            'reduce DOR susceptibility.'},
                            'RT101Q': {'1': 'K101Q is a relatively non-polymorphic mutation '
                                            'that is weakly selected in patients receiving '
                                            'NVP and EFV. It is of uncertain phenotypic and '
                                            'clinical significance.'},
                            'RT103ACDFGILMPVWY_-': {'1': 'K103N/S/T/H are NNRTI-resistance '
                                                        'mutations. K103R/E/Q are variants '
                                                        'that do not reduce NNRTI '
                                                        'susceptibility. $listMutsIn{103(NOT '
                                                        'EHNQRST)} is a highly unusual '
                                                        'mutation at this position.'},
                            'RT103EQ': {'1': 'K103E/Q are rare mutations that have not been '
                                            'associated with reduced NNRTI susceptibility.'},
                            'RT103H': {'1': 'K103H is a rare non-polymorphic mutation that '
                                            'causes high-level resistance to NVP and EFV.'},
                            'RT103N': {'1': 'K103N is a non-polymorphic mutation that causes '
                                            'high-level reductions in NVP and EFV '
                                            'susceptibility.'},
                            'RT103R': {'1': 'K103R is a polymorphic mutation that alone has '
                                            'no effect on NNRTI susceptibility. However, in '
                                            'combination with V179D, it reduces NVP and EFV '
                                            'susceptibility about 15-fold.'},
                            'RT103S': {'1': 'K103S is a non-polymorphic mutation that causes '
                                            'high-level reductions in NVP susceptibility but '
                                            'intermediate reductions in EFV susceptibility. '
                                            'Because K103S is a 2-bp change from the wildtype '
                                            'K and a 1-bp change from K103N, persons with '
                                            'K103S may be likely to have once had K103N.'},
                            'RT103T': {'1': 'K103T is an extremely rare non-polymorphic '
                                            'mutation that appears to cause '
                                            'intermediate/high-level resistance to NVP. It '
                                            'has little if any effect on EFV susceptibility.'},
                            'RT106A': {'1': 'V106A is a non-polymorphic mutation that causes '
                                            'high-level resistance to NVP and intermediate '
                                            'resistance to EFV. It is selected in vitro and '
                                            'in vivo by DOR and alone it causes intermediate '
                                            'reductions in DOR susceptibility. In combination '
                                            'with other DOR-associated DRMs it is associated '
                                            'with high-level DOR resistance.'},
                            'RT106CDEFGHKLNPQRSTWY_-': {'1': 'V106A/M are NNRTI-resistance '
                                                            'mutations. $listMutsIn{106(NOT '
                                                            'AIM)} is a highly unusual '
                                                            'mutation at this position.'},
                            'RT106I': {'1': 'V106I is occurs in 1% to 2% of viruses from '
                                            'untreated persons. It contributes to reduced '
                                            'NNRTI susceptibility in combination with other '
                                            'mutations. It has a weight of 1.5 in the Tibotec '
                                            'ETR genotypic susceptibility score despite not '
                                            'contributing much to reduced ETR susceptibility. '
                                            'It likely plays a greater role in reducing DOR '
                                            'susceptibility, particularly in combination with '
                                            'other NNRTI-associated DRMs. '},
                            'RT106M': {'1': 'V106M is a non-polymorphic mutation that causes '
                                            'high-level resistance to NVP and EFV. It is '
                                            'selected in vitro and in vivo by DOR and '
                                            'preliminary data suggests it is associated with '
                                            'low/intermediate reductions in DOR '
                                            'susceptibility.'},
                            'RT108ACDEFGHKLMNPQRSTWY_-': {'1': 'V108I is a relatively '
                                                            'non-polymorphic accessory '
                                                            'mutation selected in vitro '
                                                            'and/or in vivo with each of '
                                                            'the NNRTIs. '
                                                            '$listMutsIn{108(NOT I)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'RT108I': {'1': 'V108I is a relatively non-polymorphic accessory '
                                            'mutation selected in vitro and/or in vivo with '
                                            'each of the NNRTIs. It causes low-level '
                                            'reductions in susceptibility to NVP and DOR. '
                                            'Alone, it does not appear to reduce '
                                            'susceptibility to EFV, ETR, or RPV.'},
                            'RT115ACDEGHIKLMNPQRSTVW_-': {'1': 'Y115F causes intermediate '
                                                            'resistance to ABC and '
                                                            'low-level resistance to TDF. '
                                                            '$listMutsIn{115(NOT F)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'RT115F': {'1': 'Y115F causes intermediate resistance to ABC and '
                                            'low-level resistance to TDF.'},
                            'RT116ACDEGHIKLMNPQRSTVW_-': {'1': 'F116Y usually occurs in '
                                                            'combination with the '
                                                            'multi-NRTI resistance '
                                                            'mutation Q151M. '
                                                            '$listMutsIn{116(NOT Y)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'RT116Y': {'1': 'F116Y usually occurs in combination with the '
                                            'multi-NRTI resistance mutation Q151M.'},
                            'RT118I': {'1': 'V118I is a polymorphic accessory NRTI-resistance '
                                            'mutation that often occurs in combination with '
                                            'multiple TAMs.'},
                            'RT132ML': {'1': 'I132M is an extremely rare non-polymorphic '
                                            'mutation associated with uncertain amount of '
                                            'reduced NVP and EFV susceptibility. I132L is a '
                                            'more common, non-polymorphic NNRTI-selected '
                                            'mutation that has not been well studied.'},
                            'RT138A': {'1': 'E138A is a common polymorphic accessory mutation '
                                            'weakly selected in patients receiving ETR and '
                                            'RPV. It reduces ETR and RPV susceptibility '
                                            '~2-fold. It has a weight of 1.5 in the Tibotec '
                                            'ETR genotypic susceptibility score.'},
                            'RT138CDFHILMNPSTVWY_-': {'1': 'E138K is a non-polymorphic '
                                                        'RPV-resistance mutation. '
                                                        'E138A/Q/G/R are less well studied '
                                                        'mutations associated with reduced '
                                                        'ETR and RPV susceptibility. '
                                                        '$listMutsIn{138(NOT AGKQR)} is an '
                                                        'unusual mutation at this '
                                                        'position.'},
                            'RT138K': {'1': 'E138K is a non-polymorphic mutation selected in '
                                            'a high proportion of patients receiving RPV. It '
                                            'reduces RPV susceptibility by 2 to 3-fold and in '
                                            'combination with K101E or the NRTI-resistance '
                                            'mutation M184I, it is sufficient to cause '
                                            'virological failure on a first-line '
                                            'RPV-containing regimen. E138K causes low-level '
                                            'cross-resistance to ETR but no cross-resistance '
                                            'to NVP, EFV, or DOR.'},
                            'RT138QG': {'1': 'E138Q/G are non-polymorphic accessory mutations '
                                            'frequently selected in patients receiving ETR '
                                            'and RPV and occasionally in patients receiving '
                                            'NVP and EFV. In most studies, they cause '
                                            'low-level reductions in susceptibility to NVP, '
                                            'RPV, and ETR.'},
                            'RT138R': {'1': 'E138R is a rare non-polymorphic accessory '
                                            'mutation selected in vitro by RPV. It is '
                                            'associated with 2 to 3-fold reduced '
                                            'susceptibility to ETR and RPV.'},
                            'RT151ACDEFGHIKNPRSTVWY_-': {'1': 'Q151M causes '
                                                            'intermediate/high-level '
                                                            'resistance to AZT, ddI, d4T '
                                                            'and ABC and low-level '
                                                            'resistance to TDF, 3TC and '
                                                            'FTC. In combination with '
                                                            'accessory mutations at '
                                                            'positions 62, 75, 77, and 116, '
                                                            'Q151M confers high-level '
                                                            'resistance to AZT, ddI, d4T '
                                                            'and ABC and intermediate '
                                                            'resistance to TDF, 3TC and '
                                                            'FTC. Q151L is an extremely '
                                                            'rare transitional mutation '
                                                            'that may precede the emergence '
                                                            'of the Q151M. '
                                                            '$listMutsIn{151(NOT ML)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'RT151L': {'1': 'Q151M causes intermediate/high-level resistance '
                                            'to AZT, ddI, d4T and ABC and low-level '
                                            'resistance to TDF, 3TC and FTC. In combination '
                                            'with accessory mutations at positions 62, 75, '
                                            '77, and 116, Q151M confers high-level resistance '
                                            'to AZT, ddI, d4T and ABC and intermediate '
                                            'resistance to TDF, 3TC and FTC. Q151L is an '
                                            'extremely rare transitional mutation that may '
                                            'precede the emergence of the Q151M.'},
                            'RT151M': {'1': 'Q151M causes intermediate/high-level resistance '
                                            'to AZT, ddI, d4T and ABC and low-level '
                                            'resistance to TDF, 3TC and FTC. In combination '
                                            'with accessory mutations at positions 62, 75, '
                                            '77, and 116, Q151M confers high-level resistance '
                                            'to AZT, ddI, d4T and ABC and intermediate '
                                            'resistance to TDF, 3TC and FTC.'},
                            'RT179ACGHKMNPQRSWY_-': {'1': 'V179D/E/F/T/L are accessory '
                                                        'NNRTI-resistance mutations. '
                                                        '$listMutsIn{179(NOT DEFILT)} is an '
                                                        'unusual mutation at this '
                                                        'position.'},
                            'RT179DE': {'1': 'V179D is a polymorphic accessory NNRTI-selected '
                                            'mutation. It contributes low-levels reductions '
                                            'in susceptibility to each of the NNRTIs. The '
                                            'combination of V179D and K103R act '
                                            'synergistically to reduce NVP and EFV '
                                            'susceptibility. V179D has a weight of 1.0 in '
                                            'the Tibotec ETR genotypic susceptibility score. '
                                            'V179E is a non-polymorphic mutation '
                                            'occasionally selected by NVP and EFV. V179E '
                                            'appears similar to V179D in its effects on '
                                            'NNRTIs. V179D/E do not appear to reduce the '
                                            'virological response to a first-line '
                                            'EFV-containing regimen.'},
                            'RT179F': {'1': 'V179F is a non-polymorphic mutation frequently '
                                            'selected in patients receiving ETR. It nearly '
                                            'always occurs in combination with Y181C. Alone '
                                            'V179F has little effect on NNRTI susceptibility. '
                                            'In combination with Y181C, however, it is '
                                            'associated with high-level ETR and RPV '
                                            'resistance. It has a weight of 1.5 in the '
                                            'Tibotec ETR genotypic susceptibility score.'},
                            'RT179I': {'1': 'V179I is a polymorphic mutation that is '
                                            'frequently selected in patients receiving ETR '
                                            'and RPV. But It has little, if any, direct '
                                            'effect on NNRTI susceptibility.'},
                            'RT179L': {'1': 'V179L is a rare non-polymorphic mutation '
                                            'occasionally selected in patients receiving '
                                            'NNRTIs. Its effects on NNRTI susceptibility have '
                                            'not been well studied. It is listed as an '
                                            'RPV-associated resistance mutation in the RPV '
                                            'package-insert.'},
                            'RT179T': {'1': 'V179T is a relatively rare non-polymorphic '
                                            'mutation occasionally selected in patients '
                                            'receiving NNRTIs. It is associated with minimal, '
                                            'if any, reduction in ETR and RPV susceptibility. '
                                            'It has a weight of 1.0 in the Tibotec ETR '
                                            'genotypic susceptibility score.'},
                            'RT181ADEHKLMNPQRTW_-': {'1': 'Y181C/I/V are associated with '
                                                        'intermediate or high-level '
                                                        'resistance to NVP, ETR, and RPV. '
                                                        '$listMutsIn{181(NOT CIVFSG)} is a '
                                                        'highly unusual mutation at this '
                                                        'position.'},
                            'RT181C': {'1': 'Y181C is a non-polymorphic mutation selected in '
                                            'patients receiving NVP, ETR and RPV. It reduces '
                                            'susceptibility to NVP, ETR, RPV, and EFV by '
                                            '>50-fold, 5-fold, 3-fold, and 2-fold, '
                                            'respectively. Although Y181C itself reduces EFV '
                                            'susceptibility by only 2-fold, it has been '
                                            'associated with a reduced response to an '
                                            'EFV-containing regimen in NNRTI-experienced '
                                            'patients. Y181C has a weight of 2.5 in the '
                                            'Tibotec ETR genotypic susceptibility score. '
                                            'Alone, it does not appear to reduce DOR '
                                            'susceptibility.'},
                            'RT181FSG': {'1': 'Y181F/S/G are rare non-polymorphic '
                                            'NNRTI-associated mutations that are usually '
                                            'present as part of an electrophoretic mixture. '
                                            'They are likely to represent transitional '
                                            'mutations between Y and I or V.'},
                            'RT181IV': {'1': 'Y181I/V are 2-base pair non-polymorphic '
                                            'mutations selected by NVP and ETR. They cause '
                                            'high-level resistance to NVP, ETR, and RPV. '
                                            'They each have a weight of 3.0 in the Tibotec '
                                            'ETR genotypic susceptibility score. Their '
                                            'effects on DOR have not been '
                                            'well-characterized.'},
                            'RT184ACDEFGHKLNPQRSTWY_-': {'1': 'M184V/I cause high-level in '
                                                            'vitro resistance to 3TC and '
                                                            'FTC and low-level resistance '
                                                            'to ddI and ABC. However, '
                                                            'M184V/I are not '
                                                            'contraindications to continued '
                                                            'treatment with 3TC or FTC '
                                                            'because they increase '
                                                            'susceptibility to AZT, TDF, '
                                                            'and d4T and are associated '
                                                            'with clinically significant '
                                                            'reductions in HIV-1 '
                                                            'replication. '
                                                            '$listMutsIn{184(NOT VI)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'RT184VI': {'1': 'M184V/I cause high-level in vitro resistance to '
                                            '3TC and FTC and low-level resistance to ddI and '
                                            'ABC. However, M184V/I are not contraindications '
                                            'to continued treatment with 3TC or FTC because '
                                            'they increase susceptibility to AZT, TDF and '
                                            'd4T and are associated with clinically '
                                            'significant reductions in HIV-1 replication.'},
                            'RT188ADEGIKMNPQRSTVW_-': {'1': 'Y188L/H/C are NNRTI-resistance '
                                                            'mutations. $listMutsIn{188(NOT '
                                                            'CFHL)} is a highly unusual '
                                                            'mutation at this position.'},
                            'RT188C': {'1': 'Y188C is a non-polymorphic mutation selected in '
                                            'patients receiving NVP and EFV. It confers '
                                            'high-level resistance to NVP and EFV.'},
                            'RT188F': {'1': 'Y188F is a rare non-polymorphic NNRTI-associated '
                                            'mutation that is usually present as part of an '
                                            'electrophoretic mixture. It appears to represent '
                                            'a transitional mutation between Y and L.'},
                            'RT188H': {'1': 'Y188H is a non-polymorphic mutation selected in '
                                            'patients receiving NVP and EFV. It causes about '
                                            '5 to 10-fold reduced susceptibility to NVP and '
                                            'EFV.'},
                            'RT188L': {'1': 'Y188L is a non-polymorphic mutation that causes '
                                            'high-level resistance to NVP, EFV, RPV, and DOR, '
                                            'and potentially low-level resistance to ETR.'},
                            'RT190A': {'1': 'G190A is a non-polymorphic mutation that causes '
                                            'high-level resistance to NVP and intermediate '
                                            'resistance to EFV. It has a weight of 1.0 in the '
                                            'Tibotec ETR genotypic susceptibility score but '
                                            'does not appear to be selected by ETR or RPV or '
                                            'to reduce their in vitro susceptibility in the '
                                            'absence of other NNRTI-resistance mutations. It '
                                            'also does not appear to reduce DOR '
                                            'susceptibility.'},
                            'RT190CTV': {'1': 'G190C/T/V are rare non-polymorphic mutations '
                                            'that cause high-level resistance to NVP and '
                                            'EFV. Their effects on ETR, RPV, and DOR '
                                            'susceptibility are not known.'},
                            'RT190DFHIKLMNPWY_-': {'1': 'G190A/S/E/Q/C/T/V are each '
                                                        'associated with reduced '
                                                        'susceptibility to one or more '
                                                        'NNRTIs. $listMutsIn{190(NOT '
                                                        'ACEQRSTV)} is a highly unusual '
                                                        'mutation at this position.'},
                            'RT190EQ': {'1': 'G190E is a non-polymorphic mutation that causes '
                                            'high-level resistance to each of the NNRTIs '
                                            'including DOR. G190Q is a less common '
                                            'non-polymorphic NNRTI-selected mutation that is '
                                            'associated with high-level NVP and EFV '
                                            'resistance. Its effects on RPV, ETR, and DOR '
                                            'susceptibility is not known.'},
                            'RT190R': {'1': 'G190R is strongly suspicious for being an '
                                            'artifact of APOBEC-mediated G-to-A '
                                            'hypermutation. It has not been associated with '
                                            'reduced NNRTI susceptibility.'},
                            'RT190S': {'1': 'G190S is a non-polymorphic mutation that causes '
                                            'high-level resistance to NVP and EFV. It has a '
                                            'weight of 1.5 in the Tibotec ETR genotypic '
                                            'susceptibility score but does not appear to be '
                                            'selected by ETR or RPV or to reduce their in '
                                            'vitro susceptibility in the absence of other '
                                            'NNRTI-resistance mutations. Preliminary data '
                                            'suggests it is associated with low/intermediate '
                                            'reductions in DOR susceptibility.'},
                            'RT210W': {'1': 'L210W is a TAM that usually occurs in '
                                            'combination with M41L and T215Y. The combination '
                                            'of M41, L210W and T215Y causes high-level '
                                            'resistance to AZT and d4T and intermediate to '
                                            'high-level resistance to ddI, ABC and TDF.'},
                            'RT215F': {'1': 'T215F is a TAM that causes '
                                            'intermediate/high-level resistance to AZT and '
                                            'd4T, low-level resistance to ddI, and '
                                            'potentially low-level resistance to ABC and '
                                            'TDF.'},
                            'RT215GHKMPQRW_-': {'1': 'T215Y/F cause intermediate/high-level '
                                                    'resistance to AZT and d4T and low-level '
                                                    'resistance to ABC, ddI and TDF. '
                                                    'T215S/C/D/E/I/V/N/A/L do not reduce '
                                                    'NRTI susceptibility but arise from '
                                                    'viruses that once contained T215Y/F. '
                                                    'The presence of one of these revertant '
                                                    'mutations may suggest that the patient '
                                                    'may have once had a majority virus '
                                                    'population with T215Y/F. '
                                                    '$listMutsIn{215(NOT YFSCDEIVALN)} is a '
                                                    'highly unusual mutation at this '
                                                    'position.'},
                            'RT215SCDEIVALN': {'1': 'T215Y/F cause intermediate/high-level '
                                                    'resistance to AZT and d4T, low-level '
                                                    'resistance to ddI, and potentially '
                                                    'low-level resistance to ABC and TDF. '
                                                    'T215S/C/D/E/I/V/N/A/L do not reduce NRTI '
                                                    'susceptibility but arise from viruses '
                                                    'that once contained T215Y/F. The '
                                                    'presence of one of these revertant '
                                                    'mutations suggests that the patient may '
                                                    'have once had a majority virus '
                                                    'population with T215Y/F.'},
                            'RT215Y': {'1': 'T215Y is a TAM that causes '
                                            'intermediate/high-level resistance to AZT and '
                                            'd4T, low-level resistance to ddI, and '
                                            'potentially low-level resistance to ABC and '
                                            'TDF.'},
                            'RT219ACDFGHILMPSTVY_-': {'1': 'K219Q/E/N/R/W are accessory TAMS '
                                                        'associated with reduced '
                                                        'susceptibility to AZT and '
                                                        'possibly d4T. $listMutsIn{219(NOT '
                                                        'QENRW)} is an unusual mutation at '
                                                        'this position.'},
                            'RT219NR': {'1': 'K219N/R are accessory TAMS that usually occur '
                                            'in combination with multiple other TAMs.'},
                            'RT219QE': {'1': 'K219Q/E are accessory TAMS associated with '
                                            'reduced susceptibility to AZT and possibly '
                                            'd4T.'},
                            'RT219W': {'1': 'K219W is an uncommon NRTI-selected mutation.'},
                            'RT221Y': {'1': 'H221Y is a non-polymorphic accessory mutation '
                                            'selected primarily by NVP and RPV. It frequently '
                                            'occurs in combination with Y181C.'},
                            'RT225ACDEFGIKLMNQRSTVWY_-': {'1': 'P225H is a non-polymorphic '
                                                            'EFV-selected mutation that '
                                                            'usually occurs in combination '
                                                            'with K103N. The combination '
                                                            'of P225H and K103N '
                                                            'synergistically reduces NVP, '
                                                            'EFV and DOR susceptibility. '
                                                            '$listMutsIn{225(NOT H)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'RT225H': {'1': 'P225H is a non-polymorphic EFV-selected mutation '
                                            'that usually occurs in combination with K103N. '
                                            'The combination of P225H and K103N '
                                            'synergistically reduces NVP, EFV and DOR '
                                            'susceptibility.'},
                            'RT227ADEGHKMNPQRSTWY_-': {'1': 'F227L usually occurs together '
                                                            'with V106A. This combination is '
                                                            'associated with high-level '
                                                            'resistance to NVP and EFV. F227C '
                                                            'is a rare non-polymorphic '
                                                            'mutation that usually occurs in '
                                                            'combination with other '
                                                            'NNRTI-resistance mutations. In '
                                                            'this context it causes '
                                                            'high-level resistance to each of '
                                                            'the NNRTIs. $listMutsIn{227(NOT '
                                                            'CILV)} is a highly unusual '
                                                            'mutation at this position.'},
                            'RT227C': {'1': 'F227C is a nonpolymorphic mutation selected in '
                                            'persons receiving DOR and rarely in persons '
                                            'receiving ETR and RPV. It usually occurs in '
                                            'combination with other DRMs and in this setting '
                                            'has consistently been associated with the '
                                            'highest possible levels of DOR resistance. It is '
                                            'also usually associated with intermediate or '
                                            'higher reductions in susceptibility to NVP, EFV, '
                                            'ETR, and RPV. '},
                            'RT227ILV': {'1': 'F227L is a non-polymorphic mutation that '
                                            'usually occurs in combination with V106A. It '
                                            'is selected in vivo and in vitro with both NVP '
                                            'and DOR. In this context it is associated with '
                                            'high-level reductions in NVP and DOR '
                                            'susceptibility and intermediate reductions in '
                                            'EFV susceptibility. F227I/V are extremely rare '
                                            'mutations that have been selected in vitro by '
                                            'DOR.'},
                            'RT230ACDEFGHKNPQRSTVWY_-': {'1': 'M230L causes intermediate to '
                                                            'high-level resistance to each '
                                                            'of the NNRTIs. M230I is an '
                                                            'extremely rare mutation '
                                                            'selected in vitro by RPV. '
                                                            '$listMutsIn{230(NOT IL)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'RT230I': {'1': 'M230I is an extremely rare mutation selected by '
                                            'RPV. Its effects on NNRTI susceptibility have '
                                            'not been well studied. It also often occurs as a '
                                            'result of APOBEC-mediated G-to-A hypermutation '
                                            'resulting in viruses that are likely to be '
                                            'noninfectious.'},
                            'RT230L': {'1': 'M230L is an uncommon non-polymorphic mutation '
                                            'selected in persons receiving EFV, NVP, and RPV. '
                                            'It causes intermediate to high-level resistance '
                                            'to each of the NNRTIs.'},
                            'RT234ACDEFGHKLMNPQRSTWY_-': {'1': 'L234I is a nonpolymorphic '
                                                            'mutation selected in persons '
                                                            'receiving NVP and EFV. It is '
                                                            'also selected in vitro by ETR '
                                                            'and DOR. In combination with '
                                                            'V106A, it is associated with '
                                                            'high-level DOR resistance. '
                                                            'Its effect on susceptibility '
                                                            'when it occurs alone has not '
                                                            'been studied. '
                                                            '$listMutsIn{108(NOT I)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'RT234I': {'1': 'L234I is a nonpolymorphic mutation selected in '
                                            'persons receiving NVP and EFV. It is also '
                                            'selected in vitro by ETR and DOR. In combination '
                                            'with V106A, it is associated with high-level DOR '
                                            'resistance. Its effect on susceptibility when it '
                                            'occurs alone has not been studied.'},
                            'RT236ACDEFGHIKMNQRSTVWY_-': {'1': 'P236L is a non-polymorphic '
                                                            'mutation that causes '
                                                            'high-level DLV resistance but '
                                                            'does not reduce '
                                                            'susceptibility to any other '
                                                            'NNRTIs. $listMutsIn{236(NOT '
                                                            'L)} is a highly unusual '
                                                            'mutation at this position.'},
                            'RT236L': {'1': 'P236L is a non-polymorphic mutation that causes '
                                            'high-level DLV resistance but does not reduce '
                                            'susceptibility to any other NNRTIs.'},
                            'RT238ACDEFGHILMPQSVWY_-': {'1': 'K238T reduces NVP and EFV '
                                                            'susceptibility by about 5-fold. '
                                                            'It may also reduce '
                                                            'susceptibility to ETR and RPV. '
                                                            '$listMutsIn{238(NOT NRT)} is a '
                                                            'highly unusual mutation at this '
                                                            'position.'},
                            'RT238R': {'1': 'K238R is a common polymorphism that does not '
                                            'reduce NNRTI susceptibility.'},
                            'RT238TN': {'1': 'K238T is a non-polymorphic mutation selected in '
                                            'patients receiving NVP and EFV. It usually '
                                            'occurs in combination with K103N. It reduces '
                                            'susceptibility to NVP and EFV by about 5-fold. '
                                            'It may also reduce susceptibility to ETR and '
                                            'RPV. K238N is a non-polymorphic accessory '
                                            'mutation that is also selected by NVP and EFV. '
                                            'It appears to have minimal, if any, effects on '
                                            'NNRTI susceptibility.'},
                            'RT318ACDEGHIKLMNPQRSTVW_-': {'1': 'Y318F is an uncommon mutation '
                                                            'that causes intermediate NVP '
                                                            'resistance and potentially '
                                                            'low-level EFV resistance. '
                                                            '$listMutsIn{318(NOT F)} is a '
                                                            'highly unusual mutation at '
                                                            'this position.'},
                            'RT318F': {'1': 'Y318F is an uncommon mutation that causes '
                                            'intermediate NVP resistance and potentially '
                                            'low-level EFV resistance. It has been selected '
                                            'in vivo by DOR but its effect on DOR '
                                            'susceptibility has not been studied.'},
                            'RT348ACDEFGHKLMPQRSTVWY_-': {'1': 'N348I is a non-polymorphic '
                                                            'accessory mutation selected '
                                                            'by the NRTIs AZT and d4T and '
                                                            'by NVP and EFV. Alone it '
                                                            'reduces AZT and NVP '
                                                            'susceptibility by about '
                                                            '3-fold and EFV susceptibility '
                                                            'by 2-fold. '
                                                            '$listMutsIn{348(NOT I)} is an '
                                                            'unusual mutation at this '
                                                            'position.'},
                            'RT348I': {'1': 'N348I is a non-polymorphic accessory mutation '
                                            'selected by NVP and EFV and the NRTIs AZT and '
                                            'D4T. Alone it reduces AZT and NVP susceptibility '
                                            'by about 3-fold and EFV susceptibility by '
                                            '2-fold.'},
                            'RT40F': {'1': 'E40F is a non-polymorphic accessory mutation '
                                        'selected by AZT and d4T. It usually occurs in '
                                        'combination with M41L, L210W and T215Y. In this '
                                        'context it is associated with reduced '
                                        'susceptibility to each of the NRTIs.'},
                            'RT41I': {'1': 'M41I is usually an artifact resulting from '
                                        'APOBEC3G-mediated hypermutation.'},
                            'RT41L': {'1': 'M41L is a TAM that usually occurs with T215Y. In '
                                        'combination, M41L plus T215Y confer intermediate '
                                        '/ high-level resistance to AZT and d4T and '
                                        'contribute to reduced ddI, ABC and TDF '
                                        'susceptibility.'},
                            'RT44AD': {'1': 'E44D is a relatively non-polymorphic accessory '
                                            'mutation and E44A is a nonpolymorphic accessory '
                                            'mutation. Each usually occurs with multiple '
                                            'TAMs.'},
                            'RT62V': {'1': 'A62V is an accessory mutation that often occurs '
                                        'in combination with the multi-NRTI resistance '
                                        'mutations K65R or Q151M. A62V is widespread in '
                                        'subtype A viruses in former Soviet Union '
                                        'countries but A62 is otherwise non-polymorphic.'},
                            'RT65ACDFGHILMPQSTVWY_-': {'1': 'K65R causes '
                                                            'intermediate/high-level '
                                                            'resistance to TDF, ddI, ABC and '
                                                            'd4T and low/intermediate '
                                                            'resistance to 3TC and FTC. K65R '
                                                            'increases susceptibility to AZT. '
                                                            '$listMutsIn{65(NOT NRE)} is a '
                                                            'highly unusual mutation at this '
                                                            'position.'},
                            'RT65E': {'1': 'K65R causes intermediate/high-level resistance to '
                                        'TDF, ddI, ABC and d4T and low/intermediate '
                                        'resistance to 3TC and FTC. K65R increases '
                                        'susceptibility to AZT. K65E is an extremely rare '
                                        'NRTI-selected mutation with markedly reduced '
                                        'replication fitness.'},
                            'RT65N': {'1': 'K65R causes intermediate/high-level resistance to '
                                        'TDF, ddI, ABC and d4T and low/intermediate '
                                        'resistance to 3TC and FTC. K65R increases '
                                        'susceptibility to AZT. K65N is a rare mutation '
                                        'with similar but less pronounced effects on NRTI '
                                        'susceptibility than K65R.'},
                            'RT65R': {'1': 'K65R causes intermediate/high-level resistance to '
                                        'TDF, ddI, ABC and d4T and low/intermediate '
                                        'resistance to 3TC and FTC. K65R increases '
                                        'susceptibility to AZT.'},
                            'RT67-': {'1': 'Amino acid deletions between codons 67 and 70 are '
                                        'rare and usually occur in combination with '
                                        'multiple TAMs, K65R, or the Q151M mutation '
                                        'complex. Deletions at position 67 are more often '
                                        'associated with multiple TAMs. Deletions at '
                                        'positions 69 and 70 are more often associated '
                                        'with K65R or the Q151M mutation complex.'},
                            'RT67ACFIKLMPQRVWY_': {'1': 'D67N is a non-polymorphic TAM '
                                                        'associated with low-level resistance '
                                                        'to AZT and d4T. When present with '
                                                        'other TAMs, it contributes reduced '
                                                        'susceptibility to ABC, ddI, and TDF. '
                                                        'D67G/E/S/T/H are non-polymorphic '
                                                        'NRTI-selected mutations that also '
                                                        'generally occur in viruses with '
                                                        'multiple TAMs. $listMutsIn{67(NOT '
                                                        'NGESTHd)} is a highly unusual '
                                                        'mutation at this position.'},
                            'RT67GESTH': {'1': 'D67N is a non-polymorphic TAM associated with '
                                            'low-level resistance to AZT and d4T. When '
                                            'present with other TAMs, it contributes '
                                            'reduced susceptibility to ABC, ddI, and TDF. '
                                            'D67G/E/S/T/H are non-polymorphic '
                                            'NRTI-selected mutations that also generally '
                                            'occur in viruses with multiple TAMs.'},
                            'RT67N': {'1': 'D67N is a non-polymorphic TAM associated with '
                                        'low-level resistance to AZT and d4T. When present '
                                        'with other TAMs, it contributes reduced '
                                        'susceptibility to ABC, ddI, and TDF.'},
                            'RT68-': {'1': 'Amino acid deletions between codons 67 and 70 are '
                                        'rare and usually occur in combination with '
                                        'multiple TAMs, K65R, or the Q151M mutation '
                                        'complex. Deletions at position 67 are more often '
                                        'associated with multiple TAMs. Deletions at '
                                        'positions 69 and 70 are more often associated '
                                        'with K65R or the Q151M mutation complex. '
                                        'Deletions at codon 68 are extremely rare and less '
                                        'well characterized.'},
                            'RT68G': {'1': 'S68G is a polymorphic mutation that is often '
                                        'selected in combination with K65R. It partially '
                                        'restores the replication defect associated with '
                                        'K65R.'},
                            'RT69-': {'1': 'Amino acid deletions between codons 67 and 70 are '
                                        'rare and usually occur in combination with '
                                        'multiple TAMs, K65R, or the Q151M mutation '
                                        'complex. Deletions at position 67 are more often '
                                        'associated with multiple TAMs. Deletions at '
                                        'positions 69 and 70 are more often associated '
                                        'with K65R or the Q151M mutation complex.'},
                            'RT69CFHKLMPQRVWY': {'1': 'T69D/N/G/S/A/I/E are NRTI-selected '
                                                    'mutations. $listMutsIn{69(NOT '
                                                    'DNGSAIEid)} is a highly unusual '
                                                    'mutation at this position.'},
                            'RT69D': {'1': 'T69D is a non-polymorphic mutation that reduces '
                                        'susceptibility to ddI and possibly d4T.'},
                            'RT69G': {'1': 'T69G is a rare non-polymorphic mutation that '
                                        'usually occurs in viruses with a deletion at '
                                        'codon 67 and multiple other NRTI-resistance '
                                        'mutations.'},
                            'RT69N': {'1': 'T69N is a relatively non-polymorphic mutation '
                                        'weakly selected in patients receiving NRTIs. In '
                                        'combination with TAMs, it may contribute '
                                        'minimally reduced susceptibility to ddI, d4T, and '
                                        'AZT.'},
                            'RT69SAIE': {'1': 'T69S/A/I/E are relatively polymorphic '
                                            'mutations weakly selected in patients '
                                            'receiving NRTIs. Their effects on NRTI '
                                            'susceptibility have not been well studied.'},
                            'RT69_': {'1': 'Amino acid insertions between codons 67 and 70 '
                                        'are by convention assigned to codon 69. Together '
                                        'with TAMs, they are associated with high-level '
                                        'resistance to AZT, d4T, ddI, ABC and TDF and '
                                        'intermediate to 3TC and FTC.'},
                            'RT70-': {'1': 'Amino acid deletions between codons 67 and 70 are '
                                        'rare and usually occur in combination with '
                                        'multiple TAMs, K65R, or the Q151M mutation '
                                        'complex. Deletions at position 67 are more often '
                                        'associated with multiple TAMs. Deletions at '
                                        'positions 69 and 70 are more often associated '
                                        'with K65R or the Q151M mutation complex.'},
                            'RT70ACDFHILMPVWY_': {'1': 'K70R causes intermediate resistance '
                                                    'to AZT and possibly low-level '
                                                    'resistance to d4T and TDF. K70E/G '
                                                    'cause low-level resistance to TDF, '
                                                    'ABC, DDI and possibly 3TC and FTC. '
                                                    'K70E increases susceptibility to AZT. '
                                                    '$listMutsIn{70(NOT REGQNSTd)} is an '
                                                    'unusual mutation at this position.'},
                            'RT70EG': {'1': 'K70E/G cause low-level resistance to TDF, ABC, '
                                            'DDI and possibly 3TC and FTC. K70E/G increase '
                                            'susceptibility to AZT.'},
                            'RT70QNST': {'1': 'K70E/G cause low-level resistance to TDF, ABC, '
                                            'DDI and possibly 3TC and FTC. K70E/G increase '
                                            'susceptibility to AZT. K70Q/N/S/T are rare '
                                            'non-polymorphic NRTI-selected mutations that '
                                            'appear to have resistance profiles similar to '
                                            'K70E/G.'},
                            'RT70R': {'1': 'K70R causes intermediate resistance to AZT and '
                                        'possibly low-level resistance to D4T, DDI, ABC '
                                        'and TDF.'},
                            'RT74ACDEFGHKMNPQRSTWY_-': {'1': 'L74V/I cause high-level '
                                                            'resistance to ddI and '
                                                            'intermediate resistance to ABC. '
                                                            '$listMutsIn{74(NOT VI)} is a '
                                                            'highly unusual mutation at this '
                                                            'position.'},
                            'RT74VI': {'1': 'L74V/I cause high-level resistance to ddI and '
                                            'intermediate resistance to ABC.'},
                            'RT75CDEFGHKNPQRWY_-': {'1': 'V75M/T/I/S/A are NRTI-selected '
                                                        'mutations. $listMutsIn{75(NOT '
                                                        'MTISAL)} is a highly unusual '
                                                        'mutation at this position.'},
                            'RT75I': {'1': 'V75I is a relatively non-polymorphic accessory '
                                        'mutation that often occurs in combination with '
                                        'the multi-NRTI resistance mutation Q151M. When '
                                        'V75I occurs alone its clinical significance is '
                                        'uncertain.'},
                            'RT75M': {'1': 'V75M causes intermediate d4T resistance, '
                                        'low-level ddI resistance, and potentially '
                                        'low-level AZT resistance.'},
                            'RT75SA': {'1': 'V75S/A are non-polymorphic mutations that appear '
                                            'to reduce susceptibility to d4T and ddI.'},
                            'RT75T': {'1': 'V75T causes high-level d4T resistance and '
                                        'intermediate ddI resistance.'},
                            'RT77ACDEGHIKMNPQRSTVWY_-': {'1': 'F77L usually occurs in '
                                                            'combination with the '
                                                            'multi-NRTI resistance mutation '
                                                            'Q151M. When it occurs alone, '
                                                            'its clinical significance is '
                                                            'uncertain. $listMutsIn{77(NOT '
                                                            'L)} is a highly unusual '
                                                            'mutation at this position.'},
                            'RT77L': {'1': 'F77L usually occurs in combination with the '
                                        'multi-NRTI resistance mutation Q151M. When it '
                                        'occurs alone, its clinical significance is '
                                        'uncertain.'},
                            'RT90I': {'1': 'V90I is a polymorphic accessory mutation weakly '
                                        'selected by each of the NNRTIs. It has a weight '
                                        'of 1.0 in the Tibotec ETR genotypic '
                                        'susceptibility score but is associated with '
                                        'minimal, if any, detectable reduction in NNRTI '
                                        'susceptibility.'},
                            'RT98G': {'1': 'A98G is a non-polymorphic accessory mutation '
                                        'associated with low-level reduced susceptibility '
                                        'to each of the NNRTIs.'}},
                'drugclass': {'INSTI': ['BIC', 'DTG', 'EVG', 'RAL'],
                            'NNRTI': ['DOR', 'EFV', 'ETR', 'NVP', 'RPV'],
                            'NRTI': ['ABC', 'AZT', 'D4T', 'DDI', 'FTC', '3TC', 'TDF'],
                            'PI': ['ATV/r',
                                    'DRV/r',
                                    'FPV/r',
                                    'IDV/r',
                                    'LPV/r',
                                    'NFV',
                                    'SQV/r',
                                    'TPV/r']},
                'gene': {'IN': ['INSTI'], 'PR': ['PI'], 'RT': ['NRTI', 'NNRTI']},
                'globalrange': {1: {'max': 9, 'min': -inf},
                                2: {'max': 14, 'min': 10},
                                3: {'max': 29, 'min': 15},
                                4: {'max': 59, 'min': 30},
                                5: {'max': inf, 'min': 60}},
                'level': {'1': {'Susceptible': 'S'},
                        '2': {'Potential Low-Level Resistance': 'S'},
                        '3': {'Low-Level Resistance': 'I'},
                        '4': {'Intermediate Resistance': 'I'},
                        '5': {'High-Level Resistance': 'R'}}}
        res_def_dict = self.algorithm.parse_definitions(self.root)

        self.assertEqual(exp_def_dict, res_def_dict)

    def testParseGlobalrange(self):
        exp_grange_dict = {1: {'min': -inf, 'max': 9},
                    2: {'min': 10, 'max': 14},
                    3: {'min': 15, 'max': 29},
                    4: {'min': 30, 'max': 59},
                    5: {'min': 60, 'max': inf}}

        # Setting params
        globalrange_element = self.root.find('.//GLOBALRANGE')
        scorerange = globalrange_element.text.split(',')
        res_grange_dict = self.algorithm.parse_globalrange({}, scorerange)

        self.assertEqual(res_grange_dict, exp_grange_dict)

    def testParseDrugs(self):
        self.definitions = self.algorithm.parse_definitions(self.root)

        exp_drugs = \
            {'3TC': ([{'group': [(62, 'V')], 'value': 5},
                      [{'group': [(65, 'N')], 'value': 15},
                       {'group': [(65, 'R')], 'value': 30}],
                      {'group': [(67, 'd')], 'value': 15},
                      {'group': [(68, 'd')], 'value': 15},
                      [{'group': [(69, 'i')], 'value': 30},
                       {'group': [(69, 'd')], 'value': 15}],
                      [{'group': [(70, 'E')], 'value': 10},
                       {'group': [(70, 'G')], 'value': 10},
                       {'group': [(70, 'N')], 'value': 10},
                       {'group': [(70, 'Q')], 'value': 10},
                       {'group': [(70, 'S')], 'value': 10},
                       {'group': [(70, 'T')], 'value': 10},
                       {'group': [(70, 'd')], 'value': 15}],
                      {'group': [(75, 'I')], 'value': 5},
                      {'group': [(77, 'L')], 'value': 5},
                      {'group': [(116, 'Y')], 'value': 5},
                      [{'group': [(151, 'L')], 'value': 10},
                       {'group': [(151, 'M')], 'value': 15}],
                      [{'group': [(184, 'I')], 'value': 60},
                       {'group': [(184, 'V')], 'value': 60}],
                      {'group': [(210, 'W'), (215, 'FY')], 'value': 5},
                      {'group': [(41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
                      {'group': [(41, 'L'), (215, 'FY')], 'value': 5},
                      {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                      {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                      {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 15}],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'ABC': ([{'group': [(41, 'L')], 'value': 5},
                      {'group': [(62, 'V')], 'value': 5},
                      [{'group': [(65, 'E')], 'value': 10},
                       {'group': [(65, 'N')], 'value': 30},
                       {'group': [(65, 'R')], 'value': 45}],
                      [{'group': [(67, 'E')], 'value': 5},
                       {'group': [(67, 'G')], 'value': 5},
                       {'group': [(67, 'H')], 'value': 5},
                       {'group': [(67, 'N')], 'value': 5},
                       {'group': [(67, 'S')], 'value': 5},
                       {'group': [(67, 'T')], 'value': 5},
                       {'group': [(67, 'd')], 'value': 30}],
                      {'group': [(68, 'd')], 'value': 15},
                      [{'group': [(69, 'G')], 'value': 10},
                       {'group': [(69, 'i')], 'value': 60},
                       {'group': [(69, 'd')], 'value': 15}],
                      [{'group': [(70, 'E')], 'value': 15},
                       {'group': [(70, 'G')], 'value': 15},
                       {'group': [(70, 'N')], 'value': 15},
                       {'group': [(70, 'Q')], 'value': 15},
                       {'group': [(70, 'R')], 'value': 5},
                       {'group': [(70, 'S')], 'value': 15},
                       {'group': [(70, 'T')], 'value': 15},
                       {'group': [(70, 'd')], 'value': 15}],
                      [{'group': [(74, 'I')], 'value': 30},
                       {'group': [(74, 'V')], 'value': 30}],
                      {'group': [(75, 'I')], 'value': 5},
                      {'group': [(77, 'L')], 'value': 5},
                      {'group': [(115, 'F')], 'value': 60},
                      {'group': [(116, 'Y')], 'value': 10},
                      [{'group': [(151, 'L')], 'value': 30},
                       {'group': [(151, 'M')], 'value': 60}],
                      [{'group': [(184, 'I')], 'value': 15},
                       {'group': [(184, 'V')], 'value': 15}],
                      {'group': [(210, 'W')], 'value': 5},
                      [{'group': [(215, 'A')], 'value': 5},
                       {'group': [(215, 'C')], 'value': 5},
                       {'group': [(215, 'D')], 'value': 5},
                       {'group': [(215, 'E')], 'value': 5},
                       {'group': [(215, 'F')], 'value': 10},
                       {'group': [(215, 'I')], 'value': 5},
                       {'group': [(215, 'L')], 'value': 5},
                       {'group': [(215, 'N')], 'value': 5},
                       {'group': [(215, 'S')], 'value': 5},
                       {'group': [(215, 'V')], 'value': 5},
                       {'group': [(215, 'Y')], 'value': 10}],
                      [{'group': [(219, 'E')], 'value': 5},
                       {'group': [(219, 'N')], 'value': 5},
                       {'group': [(219, 'Q')], 'value': 5},
                       {'group': [(219, 'R')], 'value': 5}],
                      {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                      {'group': [(41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
                      {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                      {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')], 'value': 5},
                      {'group': [(67, 'EGN'), (70, 'R'), (184, 'IV'), (219, 'ENQR')],
                       'value': 20},
                      {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                      {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                      {'group': [(74, 'IV'), (184, 'IV')], 'value': 15},
                      {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 10},
                      [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                      [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(41, 'L'), (215, 'FY')], 'value': 15}]],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'ATV/r': ([{'group': [(20, 'T')], 'value': 5},
                        [{'group': [(24, 'F')], 'value': 5},
                         {'group': [(24, 'I')], 'value': 10},
                         {'group': [(24, 'M')], 'value': 5}],
                        {'group': [(32, 'I')], 'value': 15},
                        {'group': [(33, 'F')], 'value': 5},
                        [{'group': [(46, 'I')], 'value': 10},
                         {'group': [(46, 'L')], 'value': 10},
                         {'group': [(46, 'V')], 'value': 10}],
                        {'group': [(47, 'V')], 'value': 10},
                        [{'group': [(48, 'A')], 'value': 10},
                         {'group': [(48, 'L')], 'value': 10},
                         {'group': [(48, 'M')], 'value': 30},
                         {'group': [(48, 'Q')], 'value': 10},
                         {'group': [(48, 'S')], 'value': 10},
                         {'group': [(48, 'T')], 'value': 10},
                         {'group': [(48, 'V')], 'value': 30}],
                        {'group': [(50, 'L')], 'value': 60},
                        {'group': [(53, 'L')], 'value': 10},
                        [{'group': [(54, 'A')], 'value': 15},
                         {'group': [(54, 'L')], 'value': 15},
                         {'group': [(54, 'M')], 'value': 15},
                         {'group': [(54, 'S')], 'value': 15},
                         {'group': [(54, 'T')], 'value': 15},
                         {'group': [(54, 'V')], 'value': 15}],
                        [{'group': [(73, 'A')], 'value': 10},
                         {'group': [(73, 'C')], 'value': 10},
                         {'group': [(73, 'D')], 'value': 5},
                         {'group': [(73, 'S')], 'value': 10},
                         {'group': [(73, 'T')], 'value': 10},
                         {'group': [(73, 'V')], 'value': 5}],
                        {'group': [(74, 'P')], 'value': 10},
                        [{'group': [(82, 'A')], 'value': 15},
                         {'group': [(82, 'C')], 'value': 15},
                         {'group': [(82, 'F')], 'value': 15},
                         {'group': [(82, 'L')], 'value': 10},
                         {'group': [(82, 'M')], 'value': 10},
                         {'group': [(82, 'S')], 'value': 30},
                         {'group': [(82, 'T')], 'value': 30}],
                        {'group': [(83, 'D')], 'value': 10},
                        [{'group': [(84, 'A')], 'value': 60},
                         {'group': [(84, 'C')], 'value': 60},
                         {'group': [(84, 'V')], 'value': 60}],
                        [{'group': [(88, 'D')], 'value': 10},
                         {'group': [(88, 'G')], 'value': 15},
                         {'group': [(88, 'S')], 'value': 60},
                         {'group': [(88, 'T')], 'value': 15}],
                        {'group': [(90, 'M')], 'value': 25},
                        {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                        {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                        {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                        {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(46, 'ILV'), (90, 'M')], 'value': 10},
                        {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                        {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                        {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                        {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                        {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                       {1: {'max': 9, 'min': -inf},
                        2: {'max': 14, 'min': 10},
                        3: {'max': 29, 'min': 15},
                        4: {'max': 59, 'min': 30},
                        5: {'max': inf, 'min': 60}}),
             'AZT': ([{'group': [(41, 'L')], 'value': 15},
                      {'group': [(62, 'V')], 'value': 5},
                      [{'group': [(65, 'N')], 'value': -10},
                       {'group': [(65, 'R')], 'value': -15}],
                      [{'group': [(67, 'E')], 'value': 10},
                       {'group': [(67, 'G')], 'value': 10},
                       {'group': [(67, 'H')], 'value': 10},
                       {'group': [(67, 'N')], 'value': 15},
                       {'group': [(67, 'S')], 'value': 10},
                       {'group': [(67, 'T')], 'value': 10},
                       {'group': [(67, 'd')], 'value': 30}],
                      [{'group': [(69, 'G')], 'value': 5},
                       {'group': [(69, 'i')], 'value': 60}],
                      [{'group': [(70, 'E')], 'value': -10},
                       {'group': [(70, 'G')], 'value': -10},
                       {'group': [(70, 'R')], 'value': 30}],
                      [{'group': [(75, 'I')], 'value': 5},
                       {'group': [(75, 'M')], 'value': 10}],
                      {'group': [(77, 'L')], 'value': 10},
                      {'group': [(116, 'Y')], 'value': 10},
                      [{'group': [(151, 'L')], 'value': 30},
                       {'group': [(151, 'M')], 'value': 60}],
                      [{'group': [(184, 'I')], 'value': -10},
                       {'group': [(184, 'V')], 'value': -10}],
                      {'group': [(210, 'W')], 'value': 15},
                      [{'group': [(215, 'A')], 'value': 20},
                       {'group': [(215, 'C')], 'value': 20},
                       {'group': [(215, 'D')], 'value': 20},
                       {'group': [(215, 'E')], 'value': 20},
                       {'group': [(215, 'F')], 'value': 40},
                       {'group': [(215, 'I')], 'value': 20},
                       {'group': [(215, 'L')], 'value': 20},
                       {'group': [(215, 'N')], 'value': 20},
                       {'group': [(215, 'S')], 'value': 20},
                       {'group': [(215, 'V')], 'value': 20},
                       {'group': [(215, 'Y')], 'value': 40}],
                      [{'group': [(219, 'E')], 'value': 10},
                       {'group': [(219, 'N')], 'value': 10},
                       {'group': [(219, 'Q')], 'value': 10},
                       {'group': [(219, 'R')], 'value': 10},
                       {'group': [(219, 'W')], 'value': 10}],
                      {'group': [(151, 'M'), (184, 'IV')], 'value': 10},
                      {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                      {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                      {'group': [(65, 'R'), (151, 'M')], 'value': 10},
                      {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')], 'value': 5},
                      {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                      {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                      {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 10},
                      [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                      [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(41, 'L'), (215, 'FY')], 'value': 10}]],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'BIC': ([{'group': [(51, 'Y')], 'value': 10},
                      [{'group': [(66, 'I')], 'value': 5},
                       {'group': [(66, 'K')], 'value': 15}],
                      {'group': [(92, 'Q')], 'value': 10},
                      {'group': [(118, 'R')], 'value': 30},
                      {'group': [(121, 'Y')], 'value': 10},
                      [{'group': [(138, 'A')], 'value': 10},
                       {'group': [(138, 'K')], 'value': 10},
                       {'group': [(138, 'T')], 'value': 10}],
                      [{'group': [(140, 'A')], 'value': 10},
                       {'group': [(140, 'C')], 'value': 10},
                       {'group': [(140, 'S')], 'value': 10}],
                      [{'group': [(143, 'A')], 'value': 5},
                       {'group': [(143, 'C')], 'value': 5},
                       {'group': [(143, 'G')], 'value': 5},
                       {'group': [(143, 'H')], 'value': 5},
                       {'group': [(143, 'K')], 'value': 5},
                       {'group': [(143, 'R')], 'value': 5},
                       {'group': [(143, 'S')], 'value': 5}],
                      [{'group': [(148, 'H')], 'value': 25},
                       {'group': [(148, 'K')], 'value': 30},
                       {'group': [(148, 'R')], 'value': 25}],
                      {'group': [(151, 'L')], 'value': 15},
                      [{'group': [(153, 'F')], 'value': 15},
                       {'group': [(153, 'Y')], 'value': 15}],
                      {'group': [(155, 'H')], 'value': 10},
                      {'group': [(230, 'R')], 'value': 10},
                      {'group': [(263, 'K')], 'value': 25},
                      {'group': [(138, 'AKT'), (118, 'R')], 'value': 10},
                      {'group': [(138, 'AKT'), (140, 'ACS')], 'value': 10},
                      {'group': [(138, 'AKT'), (148, 'HKR')], 'value': 10},
                      {'group': [(140, 'ACS'), (148, 'HKR')], 'value': 10},
                      {'group': [(140, 'ACS'), (148, 'HKR'), (149, 'A')], 'value': 10},
                      {'group': [(143, 'ACGHRS'), (163, 'R')], 'value': 5},
                      {'group': [(143, 'ACGHRS'), (230, 'R')], 'value': 5},
                      {'group': [(148, 'HKR'), (155, 'H')], 'value': 10},
                      {'group': [(148, 'HKR'), (163, 'KR')], 'value': 5},
                      {'group': [(157, 'Q'), (263, 'K')], 'value': 10},
                      {'group': [(51, 'Y'), (263, 'K')], 'value': 10},
                      {'group': [(74, 'FIM'), (118, 'R')], 'value': 10},
                      {'group': [(74, 'FIM'), (143, 'ACGHRS')], 'value': 5},
                      {'group': [(74, 'FIM'), (148, 'HKR')], 'value': 10},
                      {'group': [(92, 'Q'), (155, 'H')], 'value': 5},
                      {'group': [(97, 'A'), (118, 'R')], 'value': 10},
                      {'group': [(97, 'A'), (148, 'HKR')], 'value': 15}],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'D4T': ([{'group': [(41, 'L')], 'value': 15},
                      {'group': [(62, 'V')], 'value': 5},
                      [{'group': [(65, 'E')], 'value': 10},
                       {'group': [(65, 'N')], 'value': 30},
                       {'group': [(65, 'R')], 'value': 60}],
                      [{'group': [(67, 'E')], 'value': 10},
                       {'group': [(67, 'G')], 'value': 10},
                       {'group': [(67, 'H')], 'value': 10},
                       {'group': [(67, 'N')], 'value': 15},
                       {'group': [(67, 'S')], 'value': 10},
                       {'group': [(67, 'T')], 'value': 10},
                       {'group': [(67, 'd')], 'value': 30}],
                      {'group': [(68, 'd')], 'value': 30},
                      [{'group': [(69, 'D')], 'value': 10},
                       {'group': [(69, 'G')], 'value': 10},
                       {'group': [(69, 'i')], 'value': 60},
                       {'group': [(69, 'd')], 'value': 30}],
                      [{'group': [(70, 'E')], 'value': 15},
                       {'group': [(70, 'G')], 'value': 15},
                       {'group': [(70, 'N')], 'value': 15},
                       {'group': [(70, 'Q')], 'value': 15},
                       {'group': [(70, 'R')], 'value': 15},
                       {'group': [(70, 'S')], 'value': 15},
                       {'group': [(70, 'T')], 'value': 15},
                       {'group': [(70, 'd')], 'value': 30}],
                      [{'group': [(75, 'A')], 'value': 30},
                       {'group': [(75, 'I')], 'value': 5},
                       {'group': [(75, 'M')], 'value': 30},
                       {'group': [(75, 'S')], 'value': 30},
                       {'group': [(75, 'T')], 'value': 60}],
                      {'group': [(77, 'L')], 'value': 10},
                      {'group': [(116, 'Y')], 'value': 10},
                      [{'group': [(151, 'L')], 'value': 30},
                       {'group': [(151, 'M')], 'value': 60}],
                      [{'group': [(184, 'I')], 'value': -10},
                       {'group': [(184, 'V')], 'value': -10}],
                      {'group': [(210, 'W')], 'value': 15},
                      [{'group': [(215, 'A')], 'value': 20},
                       {'group': [(215, 'C')], 'value': 20},
                       {'group': [(215, 'D')], 'value': 20},
                       {'group': [(215, 'E')], 'value': 20},
                       {'group': [(215, 'F')], 'value': 40},
                       {'group': [(215, 'I')], 'value': 20},
                       {'group': [(215, 'L')], 'value': 20},
                       {'group': [(215, 'N')], 'value': 20},
                       {'group': [(215, 'S')], 'value': 20},
                       {'group': [(215, 'V')], 'value': 20},
                       {'group': [(215, 'Y')], 'value': 40}],
                      [{'group': [(219, 'E')], 'value': 10},
                       {'group': [(219, 'N')], 'value': 10},
                       {'group': [(219, 'Q')], 'value': 10},
                       {'group': [(219, 'R')], 'value': 10},
                       {'group': [(219, 'W')], 'value': 10}],
                      {'group': [(151, 'M'), (184, 'IV')], 'value': 10},
                      {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                      {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                      {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')], 'value': 5},
                      {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                      {'group': [(70, 'EGNQST'), (184, 'IV')], 'value': 10},
                      {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                      {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 10},
                      [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                      [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(41, 'L'), (215, 'FY')], 'value': 10}]],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'DDI': ([{'group': [(41, 'L')], 'value': 10},
                      {'group': [(62, 'V')], 'value': 5},
                      [{'group': [(65, 'E')], 'value': 10},
                       {'group': [(65, 'N')], 'value': 30},
                       {'group': [(65, 'R')], 'value': 60}],
                      [{'group': [(67, 'E')], 'value': 5},
                       {'group': [(67, 'G')], 'value': 5},
                       {'group': [(67, 'H')], 'value': 5},
                       {'group': [(67, 'N')], 'value': 5},
                       {'group': [(67, 'S')], 'value': 5},
                       {'group': [(67, 'T')], 'value': 5},
                       {'group': [(67, 'd')], 'value': 30}],
                      {'group': [(68, 'd')], 'value': 30},
                      [{'group': [(69, 'D')], 'value': 30},
                       {'group': [(69, 'G')], 'value': 10},
                       {'group': [(69, 'i')], 'value': 60},
                       {'group': [(69, 'd')], 'value': 30}],
                      [{'group': [(70, 'E')], 'value': 15},
                       {'group': [(70, 'G')], 'value': 15},
                       {'group': [(70, 'N')], 'value': 15},
                       {'group': [(70, 'Q')], 'value': 15},
                       {'group': [(70, 'R')], 'value': 10},
                       {'group': [(70, 'S')], 'value': 15},
                       {'group': [(70, 'T')], 'value': 15},
                       {'group': [(70, 'd')], 'value': 30}],
                      [{'group': [(74, 'I')], 'value': 60},
                       {'group': [(74, 'V')], 'value': 60}],
                      [{'group': [(75, 'A')], 'value': 15},
                       {'group': [(75, 'I')], 'value': 5},
                       {'group': [(75, 'M')], 'value': 15},
                       {'group': [(75, 'S')], 'value': 15},
                       {'group': [(75, 'T')], 'value': 30}],
                      {'group': [(77, 'L')], 'value': 10},
                      {'group': [(116, 'Y')], 'value': 10},
                      [{'group': [(151, 'L')], 'value': 30},
                       {'group': [(151, 'M')], 'value': 60}],
                      [{'group': [(184, 'I')], 'value': 10},
                       {'group': [(184, 'V')], 'value': 10}],
                      {'group': [(210, 'W')], 'value': 10},
                      [{'group': [(215, 'A')], 'value': 10},
                       {'group': [(215, 'C')], 'value': 10},
                       {'group': [(215, 'D')], 'value': 10},
                       {'group': [(215, 'E')], 'value': 10},
                       {'group': [(215, 'F')], 'value': 15},
                       {'group': [(215, 'I')], 'value': 10},
                       {'group': [(215, 'L')], 'value': 10},
                       {'group': [(215, 'N')], 'value': 10},
                       {'group': [(215, 'S')], 'value': 10},
                       {'group': [(215, 'V')], 'value': 10},
                       {'group': [(215, 'Y')], 'value': 15}],
                      [{'group': [(219, 'E')], 'value': 5},
                       {'group': [(219, 'N')], 'value': 5},
                       {'group': [(219, 'Q')], 'value': 5},
                       {'group': [(219, 'R')], 'value': 5},
                       {'group': [(219, 'W')], 'value': 5}],
                      {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                      {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                      {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')], 'value': 5},
                      {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                      {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                      {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 10},
                      [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                      [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(41, 'L'), (215, 'FY')], 'value': 10}]],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'DOR': ([{'group': [(98, 'G')], 'value': 15},
                      [{'group': [(100, 'I')], 'value': 15},
                       {'group': [(100, 'V')], 'value': 10}],
                      [{'group': [(101, 'E')], 'value': 15},
                       {'group': [(101, 'P')], 'value': 10}],
                      [{'group': [(106, 'A')], 'value': 60},
                       {'group': [(106, 'I')], 'value': 15},
                       {'group': [(106, 'M')], 'value': 30}],
                      {'group': [(108, 'I')], 'value': 10},
                      [{'group': [(181, 'C')], 'value': 10},
                       {'group': [(181, 'I')], 'value': 20},
                       {'group': [(181, 'V')], 'value': 20}],
                      [{'group': [(188, 'C')], 'value': 10},
                       {'group': [(188, 'F')], 'value': 30},
                       {'group': [(188, 'H')], 'value': 10},
                       {'group': [(188, 'L')], 'value': 60}],
                      [{'group': [(190, 'C')], 'value': 10},
                       {'group': [(190, 'E')], 'value': 60},
                       {'group': [(190, 'Q')], 'value': 60},
                       {'group': [(190, 'S')], 'value': 30},
                       {'group': [(190, 'T')], 'value': 10},
                       {'group': [(190, 'V')], 'value': 10}],
                      {'group': [(221, 'Y')], 'value': 10},
                      {'group': [(225, 'H')], 'value': 30},
                      [{'group': [(227, 'C')], 'value': 60},
                       {'group': [(227, 'I')], 'value': 30},
                       {'group': [(227, 'L')], 'value': 30},
                       {'group': [(227, 'V')], 'value': 30}],
                      [{'group': [(230, 'I')], 'value': 15},
                       {'group': [(230, 'L')], 'value': 60}],
                      {'group': [(234, 'I')], 'value': 30},
                      {'group': [(318, 'F')], 'value': 30},
                      {'group': [(100, 'I'), (103, 'N')], 'value': 15},
                      {'group': [(103, 'N'), (181, 'C')], 'value': 10},
                      {'group': [(106, 'A'), (227, 'L')], 'value': 15},
                      {'group': [(108, 'I'), (181, 'C')], 'value': 10},
                      {'group': [(108, 'I'), (234, 'I')], 'value': 15},
                      {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                      {'group': [(98, 'G'), (227, 'C')], 'value': 15},
                      [{'group': [(101, 'E'), (190, 'A')], 'value': 5},
                       {'group': [(101, 'E'), (190, 'S')], 'value': 5}],
                      [{'group': [(181, 'C'), (190, 'A')], 'value': 20},
                       {'group': [(181, 'C'), (190, 'CSTV')], 'value': 10}]],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'DRV/r': ([{'group': [(10, 'F')], 'value': 5},
                        {'group': [(32, 'I')], 'value': 15},
                        {'group': [(33, 'F')], 'value': 5},
                        [{'group': [(47, 'A')], 'value': 10},
                         {'group': [(47, 'V')], 'value': 10}],
                        [{'group': [(50, 'L')], 'value': -10},
                         {'group': [(50, 'V')], 'value': 20}],
                        [{'group': [(54, 'L')], 'value': 20},
                         {'group': [(54, 'M')], 'value': 20}],
                        {'group': [(74, 'P')], 'value': 5},
                        {'group': [(76, 'V')], 'value': 20},
                        {'group': [(82, 'F')], 'value': 15},
                        [{'group': [(84, 'A')], 'value': 30},
                         {'group': [(84, 'C')], 'value': 15},
                         {'group': [(84, 'V')], 'value': 15}],
                        {'group': [(88, 'S')], 'value': -5},
                        {'group': [(89, 'V')], 'value': 5},
                        {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                        {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                        {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                        {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                        {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                        {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                        {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                        {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                        {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                        {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                        {'group': [(54, 'LM'), (89, 'V')], 'value': 5}],
                       {1: {'max': 9, 'min': -inf},
                        2: {'max': 14, 'min': 10},
                        3: {'max': 29, 'min': 15},
                        4: {'max': 59, 'min': 30},
                        5: {'max': inf, 'min': 60}}),
             'DTG': ([{'group': [(51, 'Y')], 'value': 10},
                      [{'group': [(66, 'I')], 'value': 5},
                       {'group': [(66, 'K')], 'value': 15}],
                      {'group': [(92, 'Q')], 'value': 10},
                      {'group': [(118, 'R')], 'value': 30},
                      {'group': [(121, 'Y')], 'value': 10},
                      [{'group': [(138, 'A')], 'value': 10},
                       {'group': [(138, 'K')], 'value': 10},
                       {'group': [(138, 'T')], 'value': 10}],
                      [{'group': [(140, 'A')], 'value': 10},
                       {'group': [(140, 'C')], 'value': 10},
                       {'group': [(140, 'S')], 'value': 10}],
                      [{'group': [(143, 'A')], 'value': 5},
                       {'group': [(143, 'C')], 'value': 5},
                       {'group': [(143, 'G')], 'value': 5},
                       {'group': [(143, 'H')], 'value': 5},
                       {'group': [(143, 'K')], 'value': 5},
                       {'group': [(143, 'R')], 'value': 5},
                       {'group': [(143, 'S')], 'value': 5}],
                      [{'group': [(148, 'H')], 'value': 25},
                       {'group': [(148, 'K')], 'value': 30},
                       {'group': [(148, 'R')], 'value': 25}],
                      {'group': [(151, 'L')], 'value': 15},
                      [{'group': [(153, 'F')], 'value': 15},
                       {'group': [(153, 'Y')], 'value': 15}],
                      {'group': [(155, 'H')], 'value': 10},
                      {'group': [(230, 'R')], 'value': 20},
                      {'group': [(263, 'K')], 'value': 30},
                      {'group': [(138, 'AKT'), (118, 'R')], 'value': 10},
                      {'group': [(138, 'AKT'), (140, 'ACS')], 'value': 10},
                      {'group': [(138, 'AKT'), (148, 'HKR')], 'value': 10},
                      {'group': [(140, 'ACS'), (148, 'HKR')], 'value': 10},
                      {'group': [(140, 'ACS'), (148, 'HKR'), (149, 'A')], 'value': 10},
                      {'group': [(143, 'ACGHRS'), (163, 'R')], 'value': 5},
                      {'group': [(143, 'ACGHRS'), (230, 'R')], 'value': 5},
                      {'group': [(148, 'HKR'), (155, 'H')], 'value': 10},
                      {'group': [(148, 'HKR'), (163, 'KR')], 'value': 5},
                      {'group': [(157, 'Q'), (263, 'K')], 'value': 10},
                      {'group': [(51, 'Y'), (263, 'K')], 'value': 10},
                      {'group': [(74, 'FIM'), (118, 'R')], 'value': 10},
                      {'group': [(74, 'FIM'), (143, 'ACGHRS')], 'value': 5},
                      {'group': [(74, 'FIM'), (148, 'HKR')], 'value': 15},
                      {'group': [(92, 'Q'), (155, 'H')], 'value': 5},
                      {'group': [(97, 'A'), (118, 'R')], 'value': 10},
                      {'group': [(97, 'A'), (148, 'HKR')], 'value': 15}],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'EFV': ([{'group': [(98, 'G')], 'value': 15},
                      [{'group': [(100, 'I')], 'value': 60},
                       {'group': [(100, 'V')], 'value': 30}],
                      [{'group': [(101, 'E')], 'value': 15},
                       {'group': [(101, 'H')], 'value': 10},
                       {'group': [(101, 'P')], 'value': 60}],
                      [{'group': [(103, 'H')], 'value': 60},
                       {'group': [(103, 'N')], 'value': 60},
                       {'group': [(103, 'S')], 'value': 45},
                       {'group': [(103, 'T')], 'value': 15}],
                      [{'group': [(106, 'A')], 'value': 45},
                       {'group': [(106, 'M')], 'value': 60}],
                      {'group': [(108, 'I')], 'value': 10},
                      [{'group': [(138, 'G')], 'value': 10},
                       {'group': [(138, 'K')], 'value': 10},
                       {'group': [(138, 'Q')], 'value': 10},
                       {'group': [(138, 'R')], 'value': 10}],
                      [{'group': [(179, 'D')], 'value': 10},
                       {'group': [(179, 'E')], 'value': 10},
                       {'group': [(179, 'F')], 'value': 10},
                       {'group': [(179, 'L')], 'value': 10}],
                      [{'group': [(181, 'C')], 'value': 30},
                       {'group': [(181, 'F')], 'value': 15},
                       {'group': [(181, 'G')], 'value': 15},
                       {'group': [(181, 'I')], 'value': 30},
                       {'group': [(181, 'S')], 'value': 15},
                       {'group': [(181, 'V')], 'value': 30}],
                      [{'group': [(188, 'C')], 'value': 60},
                       {'group': [(188, 'F')], 'value': 60},
                       {'group': [(188, 'H')], 'value': 30},
                       {'group': [(188, 'L')], 'value': 60}],
                      [{'group': [(190, 'A')], 'value': 45},
                       {'group': [(190, 'C')], 'value': 60},
                       {'group': [(190, 'E')], 'value': 60},
                       {'group': [(190, 'Q')], 'value': 60},
                       {'group': [(190, 'S')], 'value': 60},
                       {'group': [(190, 'T')], 'value': 60},
                       {'group': [(190, 'V')], 'value': 60}],
                      {'group': [(221, 'Y')], 'value': 10},
                      {'group': [(225, 'H')], 'value': 45},
                      [{'group': [(227, 'C')], 'value': 45},
                       {'group': [(227, 'I')], 'value': 10},
                       {'group': [(227, 'L')], 'value': 15},
                       {'group': [(227, 'V')], 'value': 10}],
                      [{'group': [(230, 'I')], 'value': 15},
                       {'group': [(230, 'L')], 'value': 45}],
                      [{'group': [(238, 'N')], 'value': 10},
                       {'group': [(238, 'T')], 'value': 30}],
                      {'group': [(318, 'F')], 'value': 10},
                      {'group': [(101, 'E'), (181, 'C')], 'value': 5},
                      {'group': [(103, 'R'), (179, 'D')], 'value': 20},
                      {'group': [(106, 'A'), (227, 'L')], 'value': 15},
                      {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                      {'group': [(98, 'G'), (227, 'C')], 'value': 15}],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'ETR': ([{'group': [(98, 'G')], 'value': 10},
                      [{'group': [(100, 'I')], 'value': 30},
                       {'group': [(100, 'V')], 'value': 10}],
                      [{'group': [(101, 'E')], 'value': 15},
                       {'group': [(101, 'H')], 'value': 10},
                       {'group': [(101, 'P')], 'value': 60}],
                      {'group': [(106, 'I')], 'value': 10},
                      [{'group': [(138, 'A')], 'value': 10},
                       {'group': [(138, 'G')], 'value': 10},
                       {'group': [(138, 'K')], 'value': 10},
                       {'group': [(138, 'Q')], 'value': 10},
                       {'group': [(138, 'R')], 'value': 10}],
                      [{'group': [(179, 'D')], 'value': 10},
                       {'group': [(179, 'E')], 'value': 10},
                       {'group': [(179, 'F')], 'value': 15},
                       {'group': [(179, 'L')], 'value': 10}],
                      [{'group': [(181, 'C')], 'value': 30},
                       {'group': [(181, 'F')], 'value': 15},
                       {'group': [(181, 'G')], 'value': 15},
                       {'group': [(181, 'I')], 'value': 60},
                       {'group': [(181, 'S')], 'value': 15},
                       {'group': [(181, 'V')], 'value': 60}],
                      {'group': [(188, 'L')], 'value': 10},
                      [{'group': [(190, 'A')], 'value': 10},
                       {'group': [(190, 'C')], 'value': 10},
                       {'group': [(190, 'E')], 'value': 45},
                       {'group': [(190, 'Q')], 'value': 45},
                       {'group': [(190, 'S')], 'value': 10},
                       {'group': [(190, 'T')], 'value': 10},
                       {'group': [(190, 'V')], 'value': 10}],
                      {'group': [(221, 'Y')], 'value': 10},
                      {'group': [(227, 'C')], 'value': 30},
                      [{'group': [(230, 'I')], 'value': 15},
                       {'group': [(230, 'L')], 'value': 30}],
                      {'group': [(101, 'E'), (181, 'C')], 'value': 5},
                      {'group': [(101, 'E'), (188, 'L')], 'value': 5},
                      {'group': [(181, 'C'), (190, 'ACSTV')], 'value': 10},
                      {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                      {'group': [(98, 'G'), (227, 'C')], 'value': 15},
                      [{'group': [(101, 'E'), (190, 'A')], 'value': 5},
                       {'group': [(101, 'E'), (190, 'S')], 'value': 5}],
                      [{'group': [(179, 'F'), (181, 'C')], 'value': 15},
                       {'group': [(179, 'T'), (181, 'C')], 'value': 10}]],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'EVG': ([{'group': [(51, 'Y')], 'value': 15},
                      [{'group': [(66, 'A')], 'value': 60},
                       {'group': [(66, 'I')], 'value': 60},
                       {'group': [(66, 'K')], 'value': 60}],
                      [{'group': [(92, 'G')], 'value': 30},
                       {'group': [(92, 'Q')], 'value': 60},
                       {'group': [(92, 'V')], 'value': 60}],
                      {'group': [(95, 'K')], 'value': 10},
                      {'group': [(97, 'A')], 'value': 10},
                      {'group': [(118, 'R')], 'value': 30},
                      {'group': [(121, 'Y')], 'value': 60},
                      [{'group': [(138, 'A')], 'value': 15},
                       {'group': [(138, 'K')], 'value': 15},
                       {'group': [(138, 'T')], 'value': 15}],
                      [{'group': [(140, 'A')], 'value': 30},
                       {'group': [(140, 'C')], 'value': 30},
                       {'group': [(140, 'S')], 'value': 30}],
                      [{'group': [(143, 'A')], 'value': 10},
                       {'group': [(143, 'C')], 'value': 10},
                       {'group': [(143, 'G')], 'value': 10},
                       {'group': [(143, 'H')], 'value': 10},
                       {'group': [(143, 'K')], 'value': 10},
                       {'group': [(143, 'R')], 'value': 10},
                       {'group': [(143, 'S')], 'value': 10}],
                      {'group': [(145, 'S')], 'value': 60},
                      {'group': [(146, 'P')], 'value': 60},
                      {'group': [(147, 'G')], 'value': 60},
                      [{'group': [(148, 'H')], 'value': 60},
                       {'group': [(148, 'K')], 'value': 60},
                       {'group': [(148, 'N')], 'value': 15},
                       {'group': [(148, 'R')], 'value': 60}],
                      [{'group': [(151, 'A')], 'value': 30},
                       {'group': [(151, 'L')], 'value': 60}],
                      [{'group': [(153, 'F')], 'value': 15},
                       {'group': [(153, 'Y')], 'value': 15}],
                      [{'group': [(155, 'H')], 'value': 60},
                       {'group': [(155, 'S')], 'value': 30},
                       {'group': [(155, 'T')], 'value': 30}],
                      {'group': [(157, 'Q')], 'value': 10},
                      [{'group': [(163, 'K')], 'value': 15},
                       {'group': [(163, 'R')], 'value': 15}],
                      {'group': [(230, 'R')], 'value': 20},
                      {'group': [(232, 'N')], 'value': 10},
                      {'group': [(263, 'K')], 'value': 30},
                      {'group': [(138, 'AKT'), (118, 'R')], 'value': 10},
                      {'group': [(138, 'AKT'), (140, 'ACS')], 'value': 15},
                      {'group': [(140, 'ACS'), (148, 'HKR'), (149, 'A')], 'value': 10},
                      {'group': [(143, 'ACGHRS'), (163, 'R')], 'value': 5},
                      {'group': [(143, 'ACGHRS'), (230, 'R')], 'value': 5},
                      {'group': [(51, 'Y'), (263, 'K')], 'value': 15},
                      {'group': [(74, 'FIM'), (118, 'R')], 'value': 10},
                      {'group': [(74, 'FIM'), (143, 'ACGHRS')], 'value': 5},
                      {'group': [(74, 'FIM'), (148, 'HKR')], 'value': 15},
                      {'group': [(97, 'A'), (118, 'R')], 'value': 10},
                      {'group': [(97, 'A'), (143, 'ACGHRS')], 'value': 5}],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'FPV/r': ([{'group': [(10, 'F')], 'value': 15},
                        {'group': [(20, 'T')], 'value': 5},
                        [{'group': [(24, 'F')], 'value': 5},
                         {'group': [(24, 'I')], 'value': 10},
                         {'group': [(24, 'M')], 'value': 5}],
                        {'group': [(32, 'I')], 'value': 30},
                        {'group': [(33, 'F')], 'value': 10},
                        [{'group': [(46, 'I')], 'value': 10},
                         {'group': [(46, 'L')], 'value': 10},
                         {'group': [(46, 'V')], 'value': 10}],
                        [{'group': [(47, 'A')], 'value': 60},
                         {'group': [(47, 'V')], 'value': 35}],
                        [{'group': [(50, 'L')], 'value': -5},
                         {'group': [(50, 'V')], 'value': 60}],
                        [{'group': [(54, 'A')], 'value': 10},
                         {'group': [(54, 'L')], 'value': 60},
                         {'group': [(54, 'M')], 'value': 60},
                         {'group': [(54, 'S')], 'value': 10},
                         {'group': [(54, 'T')], 'value': 10},
                         {'group': [(54, 'V')], 'value': 10}],
                        [{'group': [(73, 'A')], 'value': 10},
                         {'group': [(73, 'C')], 'value': 10},
                         {'group': [(73, 'D')], 'value': 5},
                         {'group': [(73, 'S')], 'value': 10},
                         {'group': [(73, 'T')], 'value': 10},
                         {'group': [(73, 'V')], 'value': 5}],
                        {'group': [(74, 'P')], 'value': 10},
                        {'group': [(76, 'V')], 'value': 60},
                        [{'group': [(82, 'A')], 'value': 15},
                         {'group': [(82, 'C')], 'value': 15},
                         {'group': [(82, 'F')], 'value': 30},
                         {'group': [(82, 'L')], 'value': 15},
                         {'group': [(82, 'M')], 'value': 15},
                         {'group': [(82, 'S')], 'value': 15},
                         {'group': [(82, 'T')], 'value': 15}],
                        [{'group': [(84, 'A')], 'value': 60},
                         {'group': [(84, 'C')], 'value': 60},
                         {'group': [(84, 'V')], 'value': 60}],
                        {'group': [(88, 'S')], 'value': -10},
                        {'group': [(89, 'V')], 'value': 10},
                        {'group': [(90, 'M')], 'value': 20},
                        {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                        {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                        {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                        {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                        {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                        {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                        {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                        {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                        {'group': [(46, 'ILV'), (76, 'V')], 'value': 10},
                        {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(46, 'ILV'), (90, 'M')], 'value': 10},
                        {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                        {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                        {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                        {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                        {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                        {'group': [(54, 'LM'), (89, 'V')], 'value': 5},
                        {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                        {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                       {1: {'max': 9, 'min': -inf},
                        2: {'max': 14, 'min': 10},
                        3: {'max': 29, 'min': 15},
                        4: {'max': 59, 'min': 30},
                        5: {'max': inf, 'min': 60}}),
             'FTC': ([{'group': [(62, 'V')], 'value': 5},
                      [{'group': [(65, 'N')], 'value': 15},
                       {'group': [(65, 'R')], 'value': 30}],
                      {'group': [(67, 'd')], 'value': 15},
                      {'group': [(68, 'd')], 'value': 15},
                      [{'group': [(69, 'i')], 'value': 30},
                       {'group': [(69, 'd')], 'value': 15}],
                      [{'group': [(70, 'E')], 'value': 10},
                       {'group': [(70, 'G')], 'value': 10},
                       {'group': [(70, 'N')], 'value': 10},
                       {'group': [(70, 'Q')], 'value': 10},
                       {'group': [(70, 'S')], 'value': 10},
                       {'group': [(70, 'T')], 'value': 10},
                       {'group': [(70, 'd')], 'value': 15}],
                      {'group': [(75, 'I')], 'value': 5},
                      {'group': [(77, 'L')], 'value': 5},
                      {'group': [(116, 'Y')], 'value': 5},
                      [{'group': [(151, 'L')], 'value': 10},
                       {'group': [(151, 'M')], 'value': 15}],
                      [{'group': [(184, 'I')], 'value': 60},
                       {'group': [(184, 'V')], 'value': 60}],
                      {'group': [(210, 'W'), (215, 'FY')], 'value': 5},
                      {'group': [(41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
                      {'group': [(41, 'L'), (215, 'FY')], 'value': 5},
                      {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                      {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                      {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 15}],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'IDV/r': ([{'group': [(10, 'F')], 'value': 10},
                        {'group': [(20, 'T')], 'value': 5},
                        [{'group': [(24, 'F')], 'value': 5},
                         {'group': [(24, 'I')], 'value': 15},
                         {'group': [(24, 'M')], 'value': 5}],
                        {'group': [(32, 'I')], 'value': 15},
                        {'group': [(33, 'F')], 'value': 5},
                        [{'group': [(46, 'I')], 'value': 10},
                         {'group': [(46, 'L')], 'value': 10},
                         {'group': [(46, 'V')], 'value': 10}],
                        [{'group': [(47, 'A')], 'value': 15},
                         {'group': [(47, 'V')], 'value': 15}],
                        [{'group': [(48, 'A')], 'value': 10},
                         {'group': [(48, 'L')], 'value': 10},
                         {'group': [(48, 'M')], 'value': 10},
                         {'group': [(48, 'Q')], 'value': 10},
                         {'group': [(48, 'S')], 'value': 10},
                         {'group': [(48, 'T')], 'value': 10},
                         {'group': [(48, 'V')], 'value': 10}],
                        {'group': [(50, 'L')], 'value': -5},
                        [{'group': [(54, 'A')], 'value': 15},
                         {'group': [(54, 'L')], 'value': 10},
                         {'group': [(54, 'M')], 'value': 15},
                         {'group': [(54, 'S')], 'value': 15},
                         {'group': [(54, 'T')], 'value': 15},
                         {'group': [(54, 'V')], 'value': 15}],
                        [{'group': [(73, 'A')], 'value': 15},
                         {'group': [(73, 'C')], 'value': 15},
                         {'group': [(73, 'D')], 'value': 5},
                         {'group': [(73, 'S')], 'value': 15},
                         {'group': [(73, 'T')], 'value': 15},
                         {'group': [(73, 'V')], 'value': 5}],
                        {'group': [(74, 'P')], 'value': 10},
                        {'group': [(76, 'V')], 'value': 30},
                        [{'group': [(82, 'A')], 'value': 30},
                         {'group': [(82, 'C')], 'value': 15},
                         {'group': [(82, 'F')], 'value': 30},
                         {'group': [(82, 'L')], 'value': 10},
                         {'group': [(82, 'M')], 'value': 30},
                         {'group': [(82, 'S')], 'value': 30},
                         {'group': [(82, 'T')], 'value': 30}],
                        {'group': [(83, 'D')], 'value': 10},
                        [{'group': [(84, 'A')], 'value': 60},
                         {'group': [(84, 'C')], 'value': 60},
                         {'group': [(84, 'V')], 'value': 60}],
                        {'group': [(88, 'S')], 'value': 15},
                        {'group': [(89, 'V')], 'value': 5},
                        {'group': [(90, 'M')], 'value': 30},
                        {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                        {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                        {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                        {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                        {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                        {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                        {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                        {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                        {'group': [(46, 'ILV'), (76, 'V')], 'value': 10},
                        {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(46, 'ILV'), (90, 'M')], 'value': 10},
                        {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                        {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                        {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                        {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                        {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                        {'group': [(54, 'LM'), (89, 'V')], 'value': 5},
                        {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                        {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                       {1: {'max': 9, 'min': -inf},
                        2: {'max': 14, 'min': 10},
                        3: {'max': 29, 'min': 15},
                        4: {'max': 59, 'min': 30},
                        5: {'max': inf, 'min': 60}}),
             'LPV/r': ([{'group': [(10, 'F')], 'value': 5},
                        [{'group': [(24, 'F')], 'value': 5},
                         {'group': [(24, 'I')], 'value': 10},
                         {'group': [(24, 'M')], 'value': 5}],
                        {'group': [(32, 'I')], 'value': 15},
                        {'group': [(33, 'F')], 'value': 5},
                        [{'group': [(46, 'I')], 'value': 10},
                         {'group': [(46, 'L')], 'value': 10},
                         {'group': [(46, 'V')], 'value': 5}],
                        [{'group': [(47, 'A')], 'value': 60},
                         {'group': [(47, 'V')], 'value': 15}],
                        [{'group': [(48, 'A')], 'value': 10},
                         {'group': [(48, 'L')], 'value': 10},
                         {'group': [(48, 'M')], 'value': 10},
                         {'group': [(48, 'Q')], 'value': 10},
                         {'group': [(48, 'S')], 'value': 10},
                         {'group': [(48, 'T')], 'value': 10},
                         {'group': [(48, 'V')], 'value': 10}],
                        [{'group': [(50, 'L')], 'value': -10},
                         {'group': [(50, 'V')], 'value': 30}],
                        [{'group': [(54, 'A')], 'value': 15},
                         {'group': [(54, 'L')], 'value': 20},
                         {'group': [(54, 'M')], 'value': 20},
                         {'group': [(54, 'S')], 'value': 15},
                         {'group': [(54, 'T')], 'value': 15},
                         {'group': [(54, 'V')], 'value': 15}],
                        [{'group': [(73, 'A')], 'value': 5},
                         {'group': [(73, 'C')], 'value': 5},
                         {'group': [(73, 'D')], 'value': 5},
                         {'group': [(73, 'S')], 'value': 5},
                         {'group': [(73, 'T')], 'value': 5},
                         {'group': [(73, 'V')], 'value': 5}],
                        {'group': [(74, 'P')], 'value': 5},
                        {'group': [(76, 'V')], 'value': 30},
                        [{'group': [(82, 'A')], 'value': 30},
                         {'group': [(82, 'C')], 'value': 15},
                         {'group': [(82, 'F')], 'value': 30},
                         {'group': [(82, 'L')], 'value': 10},
                         {'group': [(82, 'M')], 'value': 25},
                         {'group': [(82, 'S')], 'value': 30},
                         {'group': [(82, 'T')], 'value': 30}],
                        [{'group': [(84, 'A')], 'value': 60},
                         {'group': [(84, 'C')], 'value': 30},
                         {'group': [(84, 'V')], 'value': 30}],
                        {'group': [(90, 'M')], 'value': 15},
                        {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                        {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                        {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                        {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                        {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                        {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                        {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                        {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                        {'group': [(46, 'ILV'), (76, 'V')], 'value': 10},
                        {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                        {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                        {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 5},
                        {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                        {'group': [(54, 'LM'), (89, 'V')], 'value': 5},
                        {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 5}],
                       {1: {'max': 9, 'min': -inf},
                        2: {'max': 14, 'min': 10},
                        3: {'max': 29, 'min': 15},
                        4: {'max': 59, 'min': 30},
                        5: {'max': inf, 'min': 60}}),
             'NFV': ([{'group': [(10, 'F')], 'value': 15},
                      {'group': [(20, 'T')], 'value': 15},
                      {'group': [(23, 'I')], 'value': 15},
                      [{'group': [(24, 'F')], 'value': 10},
                       {'group': [(24, 'I')], 'value': 10},
                       {'group': [(24, 'M')], 'value': 10}],
                      {'group': [(30, 'N')], 'value': 60},
                      {'group': [(32, 'I')], 'value': 15},
                      {'group': [(33, 'F')], 'value': 10},
                      {'group': [(43, 'T')], 'value': 10},
                      [{'group': [(46, 'I')], 'value': 30},
                       {'group': [(46, 'L')], 'value': 20},
                       {'group': [(46, 'V')], 'value': 20}],
                      [{'group': [(47, 'A')], 'value': 30},
                       {'group': [(47, 'V')], 'value': 20}],
                      [{'group': [(48, 'A')], 'value': 30},
                       {'group': [(48, 'L')], 'value': 30},
                       {'group': [(48, 'M')], 'value': 30},
                       {'group': [(48, 'Q')], 'value': 30},
                       {'group': [(48, 'S')], 'value': 30},
                       {'group': [(48, 'T')], 'value': 30},
                       {'group': [(48, 'V')], 'value': 30}],
                      {'group': [(50, 'V')], 'value': 15},
                      {'group': [(53, 'L')], 'value': 10},
                      [{'group': [(54, 'A')], 'value': 20},
                       {'group': [(54, 'L')], 'value': 20},
                       {'group': [(54, 'M')], 'value': 20},
                       {'group': [(54, 'S')], 'value': 20},
                       {'group': [(54, 'T')], 'value': 20},
                       {'group': [(54, 'V')], 'value': 20}],
                      {'group': [(58, 'E')], 'value': 10},
                      [{'group': [(73, 'A')], 'value': 15},
                       {'group': [(73, 'C')], 'value': 15},
                       {'group': [(73, 'D')], 'value': 10},
                       {'group': [(73, 'S')], 'value': 15},
                       {'group': [(73, 'T')], 'value': 15},
                       {'group': [(73, 'V')], 'value': 10}],
                      {'group': [(74, 'P')], 'value': 20},
                      {'group': [(76, 'V')], 'value': 10},
                      [{'group': [(82, 'A')], 'value': 30},
                       {'group': [(82, 'C')], 'value': 30},
                       {'group': [(82, 'F')], 'value': 30},
                       {'group': [(82, 'L')], 'value': 10},
                       {'group': [(82, 'M')], 'value': 30},
                       {'group': [(82, 'S')], 'value': 30},
                       {'group': [(82, 'T')], 'value': 30}],
                      {'group': [(83, 'D')], 'value': 15},
                      [{'group': [(84, 'A')], 'value': 60},
                       {'group': [(84, 'C')], 'value': 60},
                       {'group': [(84, 'V')], 'value': 60}],
                      [{'group': [(88, 'D')], 'value': 60},
                       {'group': [(88, 'G')], 'value': 30},
                       {'group': [(88, 'S')], 'value': 60},
                       {'group': [(88, 'T')], 'value': 30}],
                      {'group': [(89, 'V')], 'value': 10},
                      {'group': [(90, 'M')], 'value': 60},
                      {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                      {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                      {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                      {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                      {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                      {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                      {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                      {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                      {'group': [(46, 'ILV'), (76, 'V')], 'value': 10},
                      {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                      {'group': [(46, 'ILV'), (90, 'M')], 'value': 10},
                      {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                      {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                      {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                      {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                      {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                      {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                      {'group': [(54, 'LM'), (89, 'V')], 'value': 5},
                      {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                      {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'NVP': ([{'group': [(98, 'G')], 'value': 30},
                      [{'group': [(100, 'I')], 'value': 60},
                       {'group': [(100, 'V')], 'value': 30}],
                      [{'group': [(101, 'E')], 'value': 30},
                       {'group': [(101, 'H')], 'value': 15},
                       {'group': [(101, 'P')], 'value': 60}],
                      [{'group': [(103, 'H')], 'value': 60},
                       {'group': [(103, 'N')], 'value': 60},
                       {'group': [(103, 'S')], 'value': 60},
                       {'group': [(103, 'T')], 'value': 60}],
                      [{'group': [(106, 'A')], 'value': 60},
                       {'group': [(106, 'I')], 'value': 10},
                       {'group': [(106, 'M')], 'value': 60}],
                      {'group': [(108, 'I')], 'value': 15},
                      [{'group': [(138, 'G')], 'value': 10},
                       {'group': [(138, 'K')], 'value': 10},
                       {'group': [(138, 'Q')], 'value': 10},
                       {'group': [(138, 'R')], 'value': 10}],
                      [{'group': [(179, 'D')], 'value': 10},
                       {'group': [(179, 'E')], 'value': 10},
                       {'group': [(179, 'F')], 'value': 15},
                       {'group': [(179, 'L')], 'value': 10}],
                      [{'group': [(181, 'C')], 'value': 60},
                       {'group': [(181, 'F')], 'value': 60},
                       {'group': [(181, 'G')], 'value': 60},
                       {'group': [(181, 'I')], 'value': 60},
                       {'group': [(181, 'S')], 'value': 60},
                       {'group': [(181, 'V')], 'value': 60}],
                      [{'group': [(188, 'C')], 'value': 60},
                       {'group': [(188, 'F')], 'value': 60},
                       {'group': [(188, 'H')], 'value': 60},
                       {'group': [(188, 'L')], 'value': 60}],
                      [{'group': [(190, 'A')], 'value': 60},
                       {'group': [(190, 'C')], 'value': 60},
                       {'group': [(190, 'E')], 'value': 60},
                       {'group': [(190, 'Q')], 'value': 60},
                       {'group': [(190, 'S')], 'value': 60},
                       {'group': [(190, 'T')], 'value': 60},
                       {'group': [(190, 'V')], 'value': 60}],
                      {'group': [(221, 'Y')], 'value': 15},
                      {'group': [(225, 'H')], 'value': 45},
                      [{'group': [(227, 'C')], 'value': 45},
                       {'group': [(227, 'I')], 'value': 30},
                       {'group': [(227, 'L')], 'value': 30},
                       {'group': [(227, 'V')], 'value': 30}],
                      [{'group': [(230, 'I')], 'value': 30},
                       {'group': [(230, 'L')], 'value': 60}],
                      [{'group': [(238, 'N')], 'value': 10},
                       {'group': [(238, 'T')], 'value': 30}],
                      {'group': [(318, 'F')], 'value': 30},
                      {'group': [(348, 'I')], 'value': 15},
                      {'group': [(101, 'E'), (181, 'C')], 'value': 5},
                      {'group': [(103, 'R'), (179, 'D')], 'value': 20},
                      {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                      {'group': [(98, 'G'), (227, 'C')], 'value': 15}],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'RAL': ([{'group': [(51, 'Y')], 'value': 15},
                      [{'group': [(66, 'A')], 'value': 15},
                       {'group': [(66, 'I')], 'value': 15},
                       {'group': [(66, 'K')], 'value': 60}],
                      [{'group': [(92, 'G')], 'value': 15},
                       {'group': [(92, 'Q')], 'value': 30},
                       {'group': [(92, 'V')], 'value': 30}],
                      {'group': [(95, 'K')], 'value': 10},
                      {'group': [(97, 'A')], 'value': 10},
                      {'group': [(118, 'R')], 'value': 30},
                      {'group': [(121, 'Y')], 'value': 60},
                      [{'group': [(138, 'A')], 'value': 15},
                       {'group': [(138, 'K')], 'value': 15},
                       {'group': [(138, 'T')], 'value': 15}],
                      [{'group': [(140, 'A')], 'value': 30},
                       {'group': [(140, 'C')], 'value': 30},
                       {'group': [(140, 'S')], 'value': 30}],
                      [{'group': [(143, 'A')], 'value': 60},
                       {'group': [(143, 'C')], 'value': 60},
                       {'group': [(143, 'G')], 'value': 60},
                       {'group': [(143, 'H')], 'value': 60},
                       {'group': [(143, 'K')], 'value': 60},
                       {'group': [(143, 'R')], 'value': 60},
                       {'group': [(143, 'S')], 'value': 60}],
                      [{'group': [(148, 'H')], 'value': 60},
                       {'group': [(148, 'K')], 'value': 60},
                       {'group': [(148, 'N')], 'value': 10},
                       {'group': [(148, 'R')], 'value': 60}],
                      [{'group': [(151, 'A')], 'value': 15},
                       {'group': [(151, 'L')], 'value': 30}],
                      [{'group': [(155, 'H')], 'value': 60},
                       {'group': [(155, 'S')], 'value': 30},
                       {'group': [(155, 'T')], 'value': 30}],
                      {'group': [(157, 'Q')], 'value': 10},
                      [{'group': [(163, 'K')], 'value': 15},
                       {'group': [(163, 'R')], 'value': 15}],
                      {'group': [(230, 'R')], 'value': 20},
                      {'group': [(232, 'N')], 'value': 10},
                      {'group': [(263, 'K')], 'value': 25},
                      {'group': [(138, 'AKT'), (118, 'R')], 'value': 10},
                      {'group': [(138, 'AKT'), (140, 'ACS')], 'value': 15},
                      {'group': [(140, 'ACS'), (148, 'HKR'), (149, 'A')], 'value': 10},
                      {'group': [(74, 'FIM'), (118, 'R')], 'value': 10},
                      {'group': [(74, 'FIM'), (148, 'HKR')], 'value': 15},
                      {'group': [(97, 'A'), (118, 'R')], 'value': 10}],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'RPV': ([{'group': [(98, 'G')], 'value': 15},
                      [{'group': [(100, 'I')], 'value': 60},
                       {'group': [(100, 'V')], 'value': 15}],
                      [{'group': [(101, 'E')], 'value': 45},
                       {'group': [(101, 'H')], 'value': 10},
                       {'group': [(101, 'P')], 'value': 60}],
                      {'group': [(106, 'I')], 'value': 10},
                      [{'group': [(138, 'A')], 'value': 15},
                       {'group': [(138, 'G')], 'value': 15},
                       {'group': [(138, 'K')], 'value': 45},
                       {'group': [(138, 'Q')], 'value': 15},
                       {'group': [(138, 'R')], 'value': 15}],
                      [{'group': [(179, 'D')], 'value': 10},
                       {'group': [(179, 'E')], 'value': 10},
                       {'group': [(179, 'F')], 'value': 15},
                       {'group': [(179, 'L')], 'value': 15}],
                      [{'group': [(181, 'C')], 'value': 45},
                       {'group': [(181, 'F')], 'value': 30},
                       {'group': [(181, 'G')], 'value': 30},
                       {'group': [(181, 'I')], 'value': 60},
                       {'group': [(181, 'S')], 'value': 30},
                       {'group': [(181, 'V')], 'value': 60}],
                      [{'group': [(188, 'F')], 'value': 30},
                       {'group': [(188, 'L')], 'value': 60}],
                      [{'group': [(190, 'A')], 'value': 15},
                       {'group': [(190, 'C')], 'value': 10},
                       {'group': [(190, 'E')], 'value': 60},
                       {'group': [(190, 'Q')], 'value': 45},
                       {'group': [(190, 'S')], 'value': 15},
                       {'group': [(190, 'T')], 'value': 10},
                       {'group': [(190, 'V')], 'value': 10}],
                      {'group': [(221, 'Y')], 'value': 15},
                      {'group': [(227, 'C')], 'value': 45},
                      [{'group': [(230, 'I')], 'value': 30},
                       {'group': [(230, 'L')], 'value': 60}],
                      {'group': [(101, 'E'), (184, 'I')], 'value': 15},
                      {'group': [(103, 'R'), (179, 'D')], 'value': 15},
                      {'group': [(138, 'K'), (184, 'I')], 'value': 15},
                      {'group': [(181, 'C'), (190, 'ACSTV')], 'value': 10},
                      {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                      {'group': [(98, 'G'), (227, 'C')], 'value': 15},
                      [{'group': [(179, 'F'), (181, 'C')], 'value': 15},
                       {'group': [(179, 'T'), (181, 'C')], 'value': 10}]],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'SQV/r': ([{'group': [(20, 'T')], 'value': 5},
                        [{'group': [(24, 'F')], 'value': 5},
                         {'group': [(24, 'I')], 'value': 10},
                         {'group': [(24, 'M')], 'value': 5}],
                        {'group': [(33, 'F')], 'value': 5},
                        [{'group': [(46, 'I')], 'value': 10},
                         {'group': [(46, 'L')], 'value': 10},
                         {'group': [(46, 'V')], 'value': 5}],
                        [{'group': [(48, 'A')], 'value': 60},
                         {'group': [(48, 'L')], 'value': 60},
                         {'group': [(48, 'M')], 'value': 60},
                         {'group': [(48, 'Q')], 'value': 60},
                         {'group': [(48, 'S')], 'value': 60},
                         {'group': [(48, 'T')], 'value': 60},
                         {'group': [(48, 'V')], 'value': 60}],
                        [{'group': [(50, 'L')], 'value': -5},
                         {'group': [(50, 'V')], 'value': 15}],
                        {'group': [(53, 'L')], 'value': 15},
                        [{'group': [(54, 'A')], 'value': 15},
                         {'group': [(54, 'L')], 'value': 15},
                         {'group': [(54, 'M')], 'value': 15},
                         {'group': [(54, 'S')], 'value': 15},
                         {'group': [(54, 'T')], 'value': 15},
                         {'group': [(54, 'V')], 'value': 15}],
                        [{'group': [(73, 'A')], 'value': 15},
                         {'group': [(73, 'C')], 'value': 15},
                         {'group': [(73, 'D')], 'value': 10},
                         {'group': [(73, 'S')], 'value': 15},
                         {'group': [(73, 'T')], 'value': 15},
                         {'group': [(73, 'V')], 'value': 10}],
                        {'group': [(74, 'P')], 'value': 10},
                        [{'group': [(82, 'A')], 'value': 15},
                         {'group': [(82, 'C')], 'value': 15},
                         {'group': [(82, 'F')], 'value': 10},
                         {'group': [(82, 'L')], 'value': 10},
                         {'group': [(82, 'M')], 'value': 15},
                         {'group': [(82, 'S')], 'value': 15},
                         {'group': [(82, 'T')], 'value': 15}],
                        {'group': [(83, 'D')], 'value': 10},
                        [{'group': [(84, 'A')], 'value': 60},
                         {'group': [(84, 'C')], 'value': 60},
                         {'group': [(84, 'V')], 'value': 60}],
                        [{'group': [(88, 'D')], 'value': 10},
                         {'group': [(88, 'S')], 'value': 15}],
                        {'group': [(90, 'M')], 'value': 45},
                        {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                        {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(46, 'ILV'), (90, 'M')], 'value': 5},
                        {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                        {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                        {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                        {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                        {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                       {1: {'max': 9, 'min': -inf},
                        2: {'max': 14, 'min': 10},
                        3: {'max': 29, 'min': 15},
                        4: {'max': 59, 'min': 30},
                        5: {'max': inf, 'min': 60}}),
             'TDF': ([{'group': [(41, 'L')], 'value': 5},
                      {'group': [(62, 'V')], 'value': 5},
                      [{'group': [(65, 'E')], 'value': 10},
                       {'group': [(65, 'N')], 'value': 30},
                       {'group': [(65, 'R')], 'value': 60}],
                      [{'group': [(67, 'E')], 'value': 5},
                       {'group': [(67, 'G')], 'value': 5},
                       {'group': [(67, 'H')], 'value': 5},
                       {'group': [(67, 'N')], 'value': 5},
                       {'group': [(67, 'S')], 'value': 5},
                       {'group': [(67, 'T')], 'value': 5},
                       {'group': [(67, 'd')], 'value': 30}],
                      {'group': [(68, 'd')], 'value': 15},
                      [{'group': [(69, 'G')], 'value': 5},
                       {'group': [(69, 'i')], 'value': 60},
                       {'group': [(69, 'd')], 'value': 15}],
                      [{'group': [(70, 'E')], 'value': 15},
                       {'group': [(70, 'G')], 'value': 15},
                       {'group': [(70, 'N')], 'value': 15},
                       {'group': [(70, 'Q')], 'value': 15},
                       {'group': [(70, 'R')], 'value': 5},
                       {'group': [(70, 'S')], 'value': 15},
                       {'group': [(70, 'T')], 'value': 15},
                       {'group': [(70, 'd')], 'value': 15}],
                      {'group': [(74, 'I')], 'value': 5},
                      {'group': [(75, 'I')], 'value': 5},
                      {'group': [(77, 'L')], 'value': 5},
                      {'group': [(115, 'F')], 'value': 15},
                      {'group': [(116, 'Y')], 'value': 5},
                      [{'group': [(151, 'L')], 'value': 10},
                       {'group': [(151, 'M')], 'value': 15}],
                      [{'group': [(184, 'I')], 'value': -10},
                       {'group': [(184, 'V')], 'value': -10}],
                      {'group': [(210, 'W')], 'value': 5},
                      [{'group': [(215, 'A')], 'value': 5},
                       {'group': [(215, 'C')], 'value': 5},
                       {'group': [(215, 'D')], 'value': 5},
                       {'group': [(215, 'E')], 'value': 5},
                       {'group': [(215, 'F')], 'value': 10},
                       {'group': [(215, 'I')], 'value': 5},
                       {'group': [(215, 'L')], 'value': 5},
                       {'group': [(215, 'N')], 'value': 5},
                       {'group': [(215, 'S')], 'value': 5},
                       {'group': [(215, 'V')], 'value': 5},
                       {'group': [(215, 'Y')], 'value': 10}],
                      [{'group': [(219, 'E')], 'value': 5},
                       {'group': [(219, 'N')], 'value': 5},
                       {'group': [(219, 'Q')], 'value': 5},
                       {'group': [(219, 'R')], 'value': 5}],
                      {'group': [(151, 'M'), (184, 'IV')], 'value': 10},
                      {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                      {'group': [(41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
                      {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                       'value': 5},
                      {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                      {'group': [(65, 'R'), (151, 'M')], 'value': 10},
                      {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')], 'value': 5},
                      {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                      {'group': [(70, 'EGNQST'), (184, 'IV')], 'value': 10},
                      {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                      {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 15},
                      [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                      [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                       {'group': [(41, 'L'), (215, 'FY')], 'value': 10}]],
                     {1: {'max': 9, 'min': -inf},
                      2: {'max': 14, 'min': 10},
                      3: {'max': 29, 'min': 15},
                      4: {'max': 59, 'min': 30},
                      5: {'max': inf, 'min': 60}}),
             'TPV/r': ([{'group': [(24, 'I')], 'value': -5},
                        {'group': [(32, 'I')], 'value': 5},
                        {'group': [(33, 'F')], 'value': 10},
                        {'group': [(43, 'T')], 'value': 10},
                        [{'group': [(46, 'I')], 'value': 5},
                         {'group': [(46, 'L')], 'value': 10},
                         {'group': [(46, 'V')], 'value': 5}],
                        [{'group': [(47, 'A')], 'value': 30},
                         {'group': [(47, 'V')], 'value': 30}],
                        [{'group': [(50, 'L')], 'value': -5},
                         {'group': [(50, 'V')], 'value': -5}],
                        [{'group': [(54, 'A')], 'value': 20},
                         {'group': [(54, 'L')], 'value': -10},
                         {'group': [(54, 'M')], 'value': 20},
                         {'group': [(54, 'S')], 'value': 20},
                         {'group': [(54, 'T')], 'value': 20},
                         {'group': [(54, 'V')], 'value': 20}],
                        {'group': [(58, 'E')], 'value': 15},
                        {'group': [(74, 'P')], 'value': 25},
                        {'group': [(76, 'V')], 'value': -5},
                        [{'group': [(82, 'C')], 'value': 10},
                         {'group': [(82, 'L')], 'value': 45},
                         {'group': [(82, 'M')], 'value': 10},
                         {'group': [(82, 'S')], 'value': 30},
                         {'group': [(82, 'T')], 'value': 45}],
                        {'group': [(83, 'D')], 'value': 25},
                        [{'group': [(84, 'A')], 'value': 60},
                         {'group': [(84, 'C')], 'value': 30},
                         {'group': [(84, 'V')], 'value': 30}],
                        {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5}],
                       {1: {'max': 9, 'min': -inf},
                        2: {'max': 14, 'min': 10},
                        3: {'max': 29, 'min': 15},
                        4: {'max': 59, 'min': 30},
                        5: {'max': inf, 'min': 60}}),
             'abacavir': ([{'group': [(41, 'L')], 'value': 5},
                           {'group': [(62, 'V')], 'value': 5},
                           [{'group': [(65, 'E')], 'value': 10},
                            {'group': [(65, 'N')], 'value': 30},
                            {'group': [(65, 'R')], 'value': 45}],
                           [{'group': [(67, 'E')], 'value': 5},
                            {'group': [(67, 'G')], 'value': 5},
                            {'group': [(67, 'H')], 'value': 5},
                            {'group': [(67, 'N')], 'value': 5},
                            {'group': [(67, 'S')], 'value': 5},
                            {'group': [(67, 'T')], 'value': 5},
                            {'group': [(67, 'd')], 'value': 30}],
                           {'group': [(68, 'd')], 'value': 15},
                           [{'group': [(69, 'G')], 'value': 10},
                            {'group': [(69, 'i')], 'value': 60},
                            {'group': [(69, 'd')], 'value': 15}],
                           [{'group': [(70, 'E')], 'value': 15},
                            {'group': [(70, 'G')], 'value': 15},
                            {'group': [(70, 'N')], 'value': 15},
                            {'group': [(70, 'Q')], 'value': 15},
                            {'group': [(70, 'R')], 'value': 5},
                            {'group': [(70, 'S')], 'value': 15},
                            {'group': [(70, 'T')], 'value': 15},
                            {'group': [(70, 'd')], 'value': 15}],
                           [{'group': [(74, 'I')], 'value': 30},
                            {'group': [(74, 'V')], 'value': 30}],
                           {'group': [(75, 'I')], 'value': 5},
                           {'group': [(77, 'L')], 'value': 5},
                           {'group': [(115, 'F')], 'value': 60},
                           {'group': [(116, 'Y')], 'value': 10},
                           [{'group': [(151, 'L')], 'value': 30},
                            {'group': [(151, 'M')], 'value': 60}],
                           [{'group': [(184, 'I')], 'value': 15},
                            {'group': [(184, 'V')], 'value': 15}],
                           {'group': [(210, 'W')], 'value': 5},
                           [{'group': [(215, 'A')], 'value': 5},
                            {'group': [(215, 'C')], 'value': 5},
                            {'group': [(215, 'D')], 'value': 5},
                            {'group': [(215, 'E')], 'value': 5},
                            {'group': [(215, 'F')], 'value': 10},
                            {'group': [(215, 'I')], 'value': 5},
                            {'group': [(215, 'L')], 'value': 5},
                            {'group': [(215, 'N')], 'value': 5},
                            {'group': [(215, 'S')], 'value': 5},
                            {'group': [(215, 'V')], 'value': 5},
                            {'group': [(215, 'Y')], 'value': 10}],
                           [{'group': [(219, 'E')], 'value': 5},
                            {'group': [(219, 'N')], 'value': 5},
                            {'group': [(219, 'Q')], 'value': 5},
                            {'group': [(219, 'R')], 'value': 5}],
                           {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                            'value': 5},
                           {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                           {'group': [(41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
                           {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                            'value': 5},
                           {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                           {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')], 'value': 5},
                           {'group': [(67, 'EGN'), (70, 'R'), (184, 'IV'), (219, 'ENQR')],
                            'value': 20},
                           {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                           {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                           {'group': [(74, 'IV'), (184, 'IV')], 'value': 15},
                           {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 10},
                           [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                            {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                           [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                            {'group': [(41, 'L'), (215, 'FY')], 'value': 15}]],
                          {1: {'max': 9, 'min': -inf},
                           2: {'max': 14, 'min': 10},
                           3: {'max': 29, 'min': 15},
                           4: {'max': 59, 'min': 30},
                           5: {'max': inf, 'min': 60}}),
             'atazanavir/r': ([{'group': [(20, 'T')], 'value': 5},
                               [{'group': [(24, 'F')], 'value': 5},
                                {'group': [(24, 'I')], 'value': 10},
                                {'group': [(24, 'M')], 'value': 5}],
                               {'group': [(32, 'I')], 'value': 15},
                               {'group': [(33, 'F')], 'value': 5},
                               [{'group': [(46, 'I')], 'value': 10},
                                {'group': [(46, 'L')], 'value': 10},
                                {'group': [(46, 'V')], 'value': 10}],
                               {'group': [(47, 'V')], 'value': 10},
                               [{'group': [(48, 'A')], 'value': 10},
                                {'group': [(48, 'L')], 'value': 10},
                                {'group': [(48, 'M')], 'value': 30},
                                {'group': [(48, 'Q')], 'value': 10},
                                {'group': [(48, 'S')], 'value': 10},
                                {'group': [(48, 'T')], 'value': 10},
                                {'group': [(48, 'V')], 'value': 30}],
                               {'group': [(50, 'L')], 'value': 60},
                               {'group': [(53, 'L')], 'value': 10},
                               [{'group': [(54, 'A')], 'value': 15},
                                {'group': [(54, 'L')], 'value': 15},
                                {'group': [(54, 'M')], 'value': 15},
                                {'group': [(54, 'S')], 'value': 15},
                                {'group': [(54, 'T')], 'value': 15},
                                {'group': [(54, 'V')], 'value': 15}],
                               [{'group': [(73, 'A')], 'value': 10},
                                {'group': [(73, 'C')], 'value': 10},
                                {'group': [(73, 'D')], 'value': 5},
                                {'group': [(73, 'S')], 'value': 10},
                                {'group': [(73, 'T')], 'value': 10},
                                {'group': [(73, 'V')], 'value': 5}],
                               {'group': [(74, 'P')], 'value': 10},
                               [{'group': [(82, 'A')], 'value': 15},
                                {'group': [(82, 'C')], 'value': 15},
                                {'group': [(82, 'F')], 'value': 15},
                                {'group': [(82, 'L')], 'value': 10},
                                {'group': [(82, 'M')], 'value': 10},
                                {'group': [(82, 'S')], 'value': 30},
                                {'group': [(82, 'T')], 'value': 30}],
                               {'group': [(83, 'D')], 'value': 10},
                               [{'group': [(84, 'A')], 'value': 60},
                                {'group': [(84, 'C')], 'value': 60},
                                {'group': [(84, 'V')], 'value': 60}],
                               [{'group': [(88, 'D')], 'value': 10},
                                {'group': [(88, 'G')], 'value': 15},
                                {'group': [(88, 'S')], 'value': 60},
                                {'group': [(88, 'T')], 'value': 15}],
                               {'group': [(90, 'M')], 'value': 25},
                               {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                               {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                               {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                               {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                               {'group': [(46, 'ILV'), (90, 'M')], 'value': 10},
                               {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                               {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                               {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                               {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                               {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                               {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                              {1: {'max': 9, 'min': -inf},
                               2: {'max': 14, 'min': 10},
                               3: {'max': 29, 'min': 15},
                               4: {'max': 59, 'min': 30},
                               5: {'max': inf, 'min': 60}}),
             'bictegravir': ([{'group': [(51, 'Y')], 'value': 10},
                              [{'group': [(66, 'I')], 'value': 5},
                               {'group': [(66, 'K')], 'value': 15}],
                              {'group': [(92, 'Q')], 'value': 10},
                              {'group': [(118, 'R')], 'value': 30},
                              {'group': [(121, 'Y')], 'value': 10},
                              [{'group': [(138, 'A')], 'value': 10},
                               {'group': [(138, 'K')], 'value': 10},
                               {'group': [(138, 'T')], 'value': 10}],
                              [{'group': [(140, 'A')], 'value': 10},
                               {'group': [(140, 'C')], 'value': 10},
                               {'group': [(140, 'S')], 'value': 10}],
                              [{'group': [(143, 'A')], 'value': 5},
                               {'group': [(143, 'C')], 'value': 5},
                               {'group': [(143, 'G')], 'value': 5},
                               {'group': [(143, 'H')], 'value': 5},
                               {'group': [(143, 'K')], 'value': 5},
                               {'group': [(143, 'R')], 'value': 5},
                               {'group': [(143, 'S')], 'value': 5}],
                              [{'group': [(148, 'H')], 'value': 25},
                               {'group': [(148, 'K')], 'value': 30},
                               {'group': [(148, 'R')], 'value': 25}],
                              {'group': [(151, 'L')], 'value': 15},
                              [{'group': [(153, 'F')], 'value': 15},
                               {'group': [(153, 'Y')], 'value': 15}],
                              {'group': [(155, 'H')], 'value': 10},
                              {'group': [(230, 'R')], 'value': 10},
                              {'group': [(263, 'K')], 'value': 25},
                              {'group': [(138, 'AKT'), (118, 'R')], 'value': 10},
                              {'group': [(138, 'AKT'), (140, 'ACS')], 'value': 10},
                              {'group': [(138, 'AKT'), (148, 'HKR')], 'value': 10},
                              {'group': [(140, 'ACS'), (148, 'HKR')], 'value': 10},
                              {'group': [(140, 'ACS'), (148, 'HKR'), (149, 'A')],
                               'value': 10},
                              {'group': [(143, 'ACGHRS'), (163, 'R')], 'value': 5},
                              {'group': [(143, 'ACGHRS'), (230, 'R')], 'value': 5},
                              {'group': [(148, 'HKR'), (155, 'H')], 'value': 10},
                              {'group': [(148, 'HKR'), (163, 'KR')], 'value': 5},
                              {'group': [(157, 'Q'), (263, 'K')], 'value': 10},
                              {'group': [(51, 'Y'), (263, 'K')], 'value': 10},
                              {'group': [(74, 'FIM'), (118, 'R')], 'value': 10},
                              {'group': [(74, 'FIM'), (143, 'ACGHRS')], 'value': 5},
                              {'group': [(74, 'FIM'), (148, 'HKR')], 'value': 10},
                              {'group': [(92, 'Q'), (155, 'H')], 'value': 5},
                              {'group': [(97, 'A'), (118, 'R')], 'value': 10},
                              {'group': [(97, 'A'), (148, 'HKR')], 'value': 15}],
                             {1: {'max': 9, 'min': -inf},
                              2: {'max': 14, 'min': 10},
                              3: {'max': 29, 'min': 15},
                              4: {'max': 59, 'min': 30},
                              5: {'max': inf, 'min': 60}}),
             'darunavir/r': ([{'group': [(10, 'F')], 'value': 5},
                              {'group': [(32, 'I')], 'value': 15},
                              {'group': [(33, 'F')], 'value': 5},
                              [{'group': [(47, 'A')], 'value': 10},
                               {'group': [(47, 'V')], 'value': 10}],
                              [{'group': [(50, 'L')], 'value': -10},
                               {'group': [(50, 'V')], 'value': 20}],
                              [{'group': [(54, 'L')], 'value': 20},
                               {'group': [(54, 'M')], 'value': 20}],
                              {'group': [(74, 'P')], 'value': 5},
                              {'group': [(76, 'V')], 'value': 20},
                              {'group': [(82, 'F')], 'value': 15},
                              [{'group': [(84, 'A')], 'value': 30},
                               {'group': [(84, 'C')], 'value': 15},
                               {'group': [(84, 'V')], 'value': 15}],
                              {'group': [(88, 'S')], 'value': -5},
                              {'group': [(89, 'V')], 'value': 5},
                              {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                              {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                              {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                              {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                              {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                              {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                              {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                              {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                              {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                              {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                              {'group': [(54, 'LM'), (89, 'V')], 'value': 5}],
                             {1: {'max': 9, 'min': -inf},
                              2: {'max': 14, 'min': 10},
                              3: {'max': 29, 'min': 15},
                              4: {'max': 59, 'min': 30},
                              5: {'max': inf, 'min': 60}}),
             'didanosine': ([{'group': [(41, 'L')], 'value': 10},
                             {'group': [(62, 'V')], 'value': 5},
                             [{'group': [(65, 'E')], 'value': 10},
                              {'group': [(65, 'N')], 'value': 30},
                              {'group': [(65, 'R')], 'value': 60}],
                             [{'group': [(67, 'E')], 'value': 5},
                              {'group': [(67, 'G')], 'value': 5},
                              {'group': [(67, 'H')], 'value': 5},
                              {'group': [(67, 'N')], 'value': 5},
                              {'group': [(67, 'S')], 'value': 5},
                              {'group': [(67, 'T')], 'value': 5},
                              {'group': [(67, 'd')], 'value': 30}],
                             {'group': [(68, 'd')], 'value': 30},
                             [{'group': [(69, 'D')], 'value': 30},
                              {'group': [(69, 'G')], 'value': 10},
                              {'group': [(69, 'i')], 'value': 60},
                              {'group': [(69, 'd')], 'value': 30}],
                             [{'group': [(70, 'E')], 'value': 15},
                              {'group': [(70, 'G')], 'value': 15},
                              {'group': [(70, 'N')], 'value': 15},
                              {'group': [(70, 'Q')], 'value': 15},
                              {'group': [(70, 'R')], 'value': 10},
                              {'group': [(70, 'S')], 'value': 15},
                              {'group': [(70, 'T')], 'value': 15},
                              {'group': [(70, 'd')], 'value': 30}],
                             [{'group': [(74, 'I')], 'value': 60},
                              {'group': [(74, 'V')], 'value': 60}],
                             [{'group': [(75, 'A')], 'value': 15},
                              {'group': [(75, 'I')], 'value': 5},
                              {'group': [(75, 'M')], 'value': 15},
                              {'group': [(75, 'S')], 'value': 15},
                              {'group': [(75, 'T')], 'value': 30}],
                             {'group': [(77, 'L')], 'value': 10},
                             {'group': [(116, 'Y')], 'value': 10},
                             [{'group': [(151, 'L')], 'value': 30},
                              {'group': [(151, 'M')], 'value': 60}],
                             [{'group': [(184, 'I')], 'value': 10},
                              {'group': [(184, 'V')], 'value': 10}],
                             {'group': [(210, 'W')], 'value': 10},
                             [{'group': [(215, 'A')], 'value': 10},
                              {'group': [(215, 'C')], 'value': 10},
                              {'group': [(215, 'D')], 'value': 10},
                              {'group': [(215, 'E')], 'value': 10},
                              {'group': [(215, 'F')], 'value': 15},
                              {'group': [(215, 'I')], 'value': 10},
                              {'group': [(215, 'L')], 'value': 10},
                              {'group': [(215, 'N')], 'value': 10},
                              {'group': [(215, 'S')], 'value': 10},
                              {'group': [(215, 'V')], 'value': 10},
                              {'group': [(215, 'Y')], 'value': 15}],
                             [{'group': [(219, 'E')], 'value': 5},
                              {'group': [(219, 'N')], 'value': 5},
                              {'group': [(219, 'Q')], 'value': 5},
                              {'group': [(219, 'R')], 'value': 5},
                              {'group': [(219, 'W')], 'value': 5}],
                             {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                              'value': 5},
                             {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                             {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                              'value': 5},
                             {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                             {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')],
                              'value': 5},
                             {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')],
                              'value': 10},
                             {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                             {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 10},
                             [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                              {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                             [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                              {'group': [(41, 'L'), (215, 'FY')], 'value': 10}]],
                            {1: {'max': 9, 'min': -inf},
                             2: {'max': 14, 'min': 10},
                             3: {'max': 29, 'min': 15},
                             4: {'max': 59, 'min': 30},
                             5: {'max': inf, 'min': 60}}),
             'dolutegravir': ([{'group': [(51, 'Y')], 'value': 10},
                               [{'group': [(66, 'I')], 'value': 5},
                                {'group': [(66, 'K')], 'value': 15}],
                               {'group': [(92, 'Q')], 'value': 10},
                               {'group': [(118, 'R')], 'value': 30},
                               {'group': [(121, 'Y')], 'value': 10},
                               [{'group': [(138, 'A')], 'value': 10},
                                {'group': [(138, 'K')], 'value': 10},
                                {'group': [(138, 'T')], 'value': 10}],
                               [{'group': [(140, 'A')], 'value': 10},
                                {'group': [(140, 'C')], 'value': 10},
                                {'group': [(140, 'S')], 'value': 10}],
                               [{'group': [(143, 'A')], 'value': 5},
                                {'group': [(143, 'C')], 'value': 5},
                                {'group': [(143, 'G')], 'value': 5},
                                {'group': [(143, 'H')], 'value': 5},
                                {'group': [(143, 'K')], 'value': 5},
                                {'group': [(143, 'R')], 'value': 5},
                                {'group': [(143, 'S')], 'value': 5}],
                               [{'group': [(148, 'H')], 'value': 25},
                                {'group': [(148, 'K')], 'value': 30},
                                {'group': [(148, 'R')], 'value': 25}],
                               {'group': [(151, 'L')], 'value': 15},
                               [{'group': [(153, 'F')], 'value': 15},
                                {'group': [(153, 'Y')], 'value': 15}],
                               {'group': [(155, 'H')], 'value': 10},
                               {'group': [(230, 'R')], 'value': 20},
                               {'group': [(263, 'K')], 'value': 30},
                               {'group': [(138, 'AKT'), (118, 'R')], 'value': 10},
                               {'group': [(138, 'AKT'), (140, 'ACS')], 'value': 10},
                               {'group': [(138, 'AKT'), (148, 'HKR')], 'value': 10},
                               {'group': [(140, 'ACS'), (148, 'HKR')], 'value': 10},
                               {'group': [(140, 'ACS'), (148, 'HKR'), (149, 'A')],
                                'value': 10},
                               {'group': [(143, 'ACGHRS'), (163, 'R')], 'value': 5},
                               {'group': [(143, 'ACGHRS'), (230, 'R')], 'value': 5},
                               {'group': [(148, 'HKR'), (155, 'H')], 'value': 10},
                               {'group': [(148, 'HKR'), (163, 'KR')], 'value': 5},
                               {'group': [(157, 'Q'), (263, 'K')], 'value': 10},
                               {'group': [(51, 'Y'), (263, 'K')], 'value': 10},
                               {'group': [(74, 'FIM'), (118, 'R')], 'value': 10},
                               {'group': [(74, 'FIM'), (143, 'ACGHRS')], 'value': 5},
                               {'group': [(74, 'FIM'), (148, 'HKR')], 'value': 15},
                               {'group': [(92, 'Q'), (155, 'H')], 'value': 5},
                               {'group': [(97, 'A'), (118, 'R')], 'value': 10},
                               {'group': [(97, 'A'), (148, 'HKR')], 'value': 15}],
                              {1: {'max': 9, 'min': -inf},
                               2: {'max': 14, 'min': 10},
                               3: {'max': 29, 'min': 15},
                               4: {'max': 59, 'min': 30},
                               5: {'max': inf, 'min': 60}}),
             'doravirine': ([{'group': [(98, 'G')], 'value': 15},
                             [{'group': [(100, 'I')], 'value': 15},
                              {'group': [(100, 'V')], 'value': 10}],
                             [{'group': [(101, 'E')], 'value': 15},
                              {'group': [(101, 'P')], 'value': 10}],
                             [{'group': [(106, 'A')], 'value': 60},
                              {'group': [(106, 'I')], 'value': 15},
                              {'group': [(106, 'M')], 'value': 30}],
                             {'group': [(108, 'I')], 'value': 10},
                             [{'group': [(181, 'C')], 'value': 10},
                              {'group': [(181, 'I')], 'value': 20},
                              {'group': [(181, 'V')], 'value': 20}],
                             [{'group': [(188, 'C')], 'value': 10},
                              {'group': [(188, 'F')], 'value': 30},
                              {'group': [(188, 'H')], 'value': 10},
                              {'group': [(188, 'L')], 'value': 60}],
                             [{'group': [(190, 'C')], 'value': 10},
                              {'group': [(190, 'E')], 'value': 60},
                              {'group': [(190, 'Q')], 'value': 60},
                              {'group': [(190, 'S')], 'value': 30},
                              {'group': [(190, 'T')], 'value': 10},
                              {'group': [(190, 'V')], 'value': 10}],
                             {'group': [(221, 'Y')], 'value': 10},
                             {'group': [(225, 'H')], 'value': 30},
                             [{'group': [(227, 'C')], 'value': 60},
                              {'group': [(227, 'I')], 'value': 30},
                              {'group': [(227, 'L')], 'value': 30},
                              {'group': [(227, 'V')], 'value': 30}],
                             [{'group': [(230, 'I')], 'value': 15},
                              {'group': [(230, 'L')], 'value': 60}],
                             {'group': [(234, 'I')], 'value': 30},
                             {'group': [(318, 'F')], 'value': 30},
                             {'group': [(100, 'I'), (103, 'N')], 'value': 15},
                             {'group': [(103, 'N'), (181, 'C')], 'value': 10},
                             {'group': [(106, 'A'), (227, 'L')], 'value': 15},
                             {'group': [(108, 'I'), (181, 'C')], 'value': 10},
                             {'group': [(108, 'I'), (234, 'I')], 'value': 15},
                             {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                             {'group': [(98, 'G'), (227, 'C')], 'value': 15},
                             [{'group': [(101, 'E'), (190, 'A')], 'value': 5},
                              {'group': [(101, 'E'), (190, 'S')], 'value': 5}],
                             [{'group': [(181, 'C'), (190, 'A')], 'value': 20},
                              {'group': [(181, 'C'), (190, 'CSTV')], 'value': 10}]],
                            {1: {'max': 9, 'min': -inf},
                             2: {'max': 14, 'min': 10},
                             3: {'max': 29, 'min': 15},
                             4: {'max': 59, 'min': 30},
                             5: {'max': inf, 'min': 60}}),
             'efavirenz': ([{'group': [(98, 'G')], 'value': 15},
                            [{'group': [(100, 'I')], 'value': 60},
                             {'group': [(100, 'V')], 'value': 30}],
                            [{'group': [(101, 'E')], 'value': 15},
                             {'group': [(101, 'H')], 'value': 10},
                             {'group': [(101, 'P')], 'value': 60}],
                            [{'group': [(103, 'H')], 'value': 60},
                             {'group': [(103, 'N')], 'value': 60},
                             {'group': [(103, 'S')], 'value': 45},
                             {'group': [(103, 'T')], 'value': 15}],
                            [{'group': [(106, 'A')], 'value': 45},
                             {'group': [(106, 'M')], 'value': 60}],
                            {'group': [(108, 'I')], 'value': 10},
                            [{'group': [(138, 'G')], 'value': 10},
                             {'group': [(138, 'K')], 'value': 10},
                             {'group': [(138, 'Q')], 'value': 10},
                             {'group': [(138, 'R')], 'value': 10}],
                            [{'group': [(179, 'D')], 'value': 10},
                             {'group': [(179, 'E')], 'value': 10},
                             {'group': [(179, 'F')], 'value': 10},
                             {'group': [(179, 'L')], 'value': 10}],
                            [{'group': [(181, 'C')], 'value': 30},
                             {'group': [(181, 'F')], 'value': 15},
                             {'group': [(181, 'G')], 'value': 15},
                             {'group': [(181, 'I')], 'value': 30},
                             {'group': [(181, 'S')], 'value': 15},
                             {'group': [(181, 'V')], 'value': 30}],
                            [{'group': [(188, 'C')], 'value': 60},
                             {'group': [(188, 'F')], 'value': 60},
                             {'group': [(188, 'H')], 'value': 30},
                             {'group': [(188, 'L')], 'value': 60}],
                            [{'group': [(190, 'A')], 'value': 45},
                             {'group': [(190, 'C')], 'value': 60},
                             {'group': [(190, 'E')], 'value': 60},
                             {'group': [(190, 'Q')], 'value': 60},
                             {'group': [(190, 'S')], 'value': 60},
                             {'group': [(190, 'T')], 'value': 60},
                             {'group': [(190, 'V')], 'value': 60}],
                            {'group': [(221, 'Y')], 'value': 10},
                            {'group': [(225, 'H')], 'value': 45},
                            [{'group': [(227, 'C')], 'value': 45},
                             {'group': [(227, 'I')], 'value': 10},
                             {'group': [(227, 'L')], 'value': 15},
                             {'group': [(227, 'V')], 'value': 10}],
                            [{'group': [(230, 'I')], 'value': 15},
                             {'group': [(230, 'L')], 'value': 45}],
                            [{'group': [(238, 'N')], 'value': 10},
                             {'group': [(238, 'T')], 'value': 30}],
                            {'group': [(318, 'F')], 'value': 10},
                            {'group': [(101, 'E'), (181, 'C')], 'value': 5},
                            {'group': [(103, 'R'), (179, 'D')], 'value': 20},
                            {'group': [(106, 'A'), (227, 'L')], 'value': 15},
                            {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                            {'group': [(98, 'G'), (227, 'C')], 'value': 15}],
                           {1: {'max': 9, 'min': -inf},
                            2: {'max': 14, 'min': 10},
                            3: {'max': 29, 'min': 15},
                            4: {'max': 59, 'min': 30},
                            5: {'max': inf, 'min': 60}}),
             'elvitegravir': ([{'group': [(51, 'Y')], 'value': 15},
                               [{'group': [(66, 'A')], 'value': 60},
                                {'group': [(66, 'I')], 'value': 60},
                                {'group': [(66, 'K')], 'value': 60}],
                               [{'group': [(92, 'G')], 'value': 30},
                                {'group': [(92, 'Q')], 'value': 60},
                                {'group': [(92, 'V')], 'value': 60}],
                               {'group': [(95, 'K')], 'value': 10},
                               {'group': [(97, 'A')], 'value': 10},
                               {'group': [(118, 'R')], 'value': 30},
                               {'group': [(121, 'Y')], 'value': 60},
                               [{'group': [(138, 'A')], 'value': 15},
                                {'group': [(138, 'K')], 'value': 15},
                                {'group': [(138, 'T')], 'value': 15}],
                               [{'group': [(140, 'A')], 'value': 30},
                                {'group': [(140, 'C')], 'value': 30},
                                {'group': [(140, 'S')], 'value': 30}],
                               [{'group': [(143, 'A')], 'value': 10},
                                {'group': [(143, 'C')], 'value': 10},
                                {'group': [(143, 'G')], 'value': 10},
                                {'group': [(143, 'H')], 'value': 10},
                                {'group': [(143, 'K')], 'value': 10},
                                {'group': [(143, 'R')], 'value': 10},
                                {'group': [(143, 'S')], 'value': 10}],
                               {'group': [(145, 'S')], 'value': 60},
                               {'group': [(146, 'P')], 'value': 60},
                               {'group': [(147, 'G')], 'value': 60},
                               [{'group': [(148, 'H')], 'value': 60},
                                {'group': [(148, 'K')], 'value': 60},
                                {'group': [(148, 'N')], 'value': 15},
                                {'group': [(148, 'R')], 'value': 60}],
                               [{'group': [(151, 'A')], 'value': 30},
                                {'group': [(151, 'L')], 'value': 60}],
                               [{'group': [(153, 'F')], 'value': 15},
                                {'group': [(153, 'Y')], 'value': 15}],
                               [{'group': [(155, 'H')], 'value': 60},
                                {'group': [(155, 'S')], 'value': 30},
                                {'group': [(155, 'T')], 'value': 30}],
                               {'group': [(157, 'Q')], 'value': 10},
                               [{'group': [(163, 'K')], 'value': 15},
                                {'group': [(163, 'R')], 'value': 15}],
                               {'group': [(230, 'R')], 'value': 20},
                               {'group': [(232, 'N')], 'value': 10},
                               {'group': [(263, 'K')], 'value': 30},
                               {'group': [(138, 'AKT'), (118, 'R')], 'value': 10},
                               {'group': [(138, 'AKT'), (140, 'ACS')], 'value': 15},
                               {'group': [(140, 'ACS'), (148, 'HKR'), (149, 'A')],
                                'value': 10},
                               {'group': [(143, 'ACGHRS'), (163, 'R')], 'value': 5},
                               {'group': [(143, 'ACGHRS'), (230, 'R')], 'value': 5},
                               {'group': [(51, 'Y'), (263, 'K')], 'value': 15},
                               {'group': [(74, 'FIM'), (118, 'R')], 'value': 10},
                               {'group': [(74, 'FIM'), (143, 'ACGHRS')], 'value': 5},
                               {'group': [(74, 'FIM'), (148, 'HKR')], 'value': 15},
                               {'group': [(97, 'A'), (118, 'R')], 'value': 10},
                               {'group': [(97, 'A'), (143, 'ACGHRS')], 'value': 5}],
                              {1: {'max': 9, 'min': -inf},
                               2: {'max': 14, 'min': 10},
                               3: {'max': 29, 'min': 15},
                               4: {'max': 59, 'min': 30},
                               5: {'max': inf, 'min': 60}}),
             'emtricitabine': ([{'group': [(62, 'V')], 'value': 5},
                                [{'group': [(65, 'N')], 'value': 15},
                                 {'group': [(65, 'R')], 'value': 30}],
                                {'group': [(67, 'd')], 'value': 15},
                                {'group': [(68, 'd')], 'value': 15},
                                [{'group': [(69, 'i')], 'value': 30},
                                 {'group': [(69, 'd')], 'value': 15}],
                                [{'group': [(70, 'E')], 'value': 10},
                                 {'group': [(70, 'G')], 'value': 10},
                                 {'group': [(70, 'N')], 'value': 10},
                                 {'group': [(70, 'Q')], 'value': 10},
                                 {'group': [(70, 'S')], 'value': 10},
                                 {'group': [(70, 'T')], 'value': 10},
                                 {'group': [(70, 'd')], 'value': 15}],
                                {'group': [(75, 'I')], 'value': 5},
                                {'group': [(77, 'L')], 'value': 5},
                                {'group': [(116, 'Y')], 'value': 5},
                                [{'group': [(151, 'L')], 'value': 10},
                                 {'group': [(151, 'M')], 'value': 15}],
                                [{'group': [(184, 'I')], 'value': 60},
                                 {'group': [(184, 'V')], 'value': 60}],
                                {'group': [(210, 'W'), (215, 'FY')], 'value': 5},
                                {'group': [(41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
                                {'group': [(41, 'L'), (215, 'FY')], 'value': 5},
                                {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')],
                                 'value': 5},
                                {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')],
                                 'value': 10},
                                {'group': [(77, 'L'), (116, 'Y'), (151, 'M')],
                                 'value': 15}],
                               {1: {'max': 9, 'min': -inf},
                                2: {'max': 14, 'min': 10},
                                3: {'max': 29, 'min': 15},
                                4: {'max': 59, 'min': 30},
                                5: {'max': inf, 'min': 60}}),
             'etravirine': ([{'group': [(98, 'G')], 'value': 10},
                             [{'group': [(100, 'I')], 'value': 30},
                              {'group': [(100, 'V')], 'value': 10}],
                             [{'group': [(101, 'E')], 'value': 15},
                              {'group': [(101, 'H')], 'value': 10},
                              {'group': [(101, 'P')], 'value': 60}],
                             {'group': [(106, 'I')], 'value': 10},
                             [{'group': [(138, 'A')], 'value': 10},
                              {'group': [(138, 'G')], 'value': 10},
                              {'group': [(138, 'K')], 'value': 10},
                              {'group': [(138, 'Q')], 'value': 10},
                              {'group': [(138, 'R')], 'value': 10}],
                             [{'group': [(179, 'D')], 'value': 10},
                              {'group': [(179, 'E')], 'value': 10},
                              {'group': [(179, 'F')], 'value': 15},
                              {'group': [(179, 'L')], 'value': 10}],
                             [{'group': [(181, 'C')], 'value': 30},
                              {'group': [(181, 'F')], 'value': 15},
                              {'group': [(181, 'G')], 'value': 15},
                              {'group': [(181, 'I')], 'value': 60},
                              {'group': [(181, 'S')], 'value': 15},
                              {'group': [(181, 'V')], 'value': 60}],
                             {'group': [(188, 'L')], 'value': 10},
                             [{'group': [(190, 'A')], 'value': 10},
                              {'group': [(190, 'C')], 'value': 10},
                              {'group': [(190, 'E')], 'value': 45},
                              {'group': [(190, 'Q')], 'value': 45},
                              {'group': [(190, 'S')], 'value': 10},
                              {'group': [(190, 'T')], 'value': 10},
                              {'group': [(190, 'V')], 'value': 10}],
                             {'group': [(221, 'Y')], 'value': 10},
                             {'group': [(227, 'C')], 'value': 30},
                             [{'group': [(230, 'I')], 'value': 15},
                              {'group': [(230, 'L')], 'value': 30}],
                             {'group': [(101, 'E'), (181, 'C')], 'value': 5},
                             {'group': [(101, 'E'), (188, 'L')], 'value': 5},
                             {'group': [(181, 'C'), (190, 'ACSTV')], 'value': 10},
                             {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                             {'group': [(98, 'G'), (227, 'C')], 'value': 15},
                             [{'group': [(101, 'E'), (190, 'A')], 'value': 5},
                              {'group': [(101, 'E'), (190, 'S')], 'value': 5}],
                             [{'group': [(179, 'F'), (181, 'C')], 'value': 15},
                              {'group': [(179, 'T'), (181, 'C')], 'value': 10}]],
                            {1: {'max': 9, 'min': -inf},
                             2: {'max': 14, 'min': 10},
                             3: {'max': 29, 'min': 15},
                             4: {'max': 59, 'min': 30},
                             5: {'max': inf, 'min': 60}}),
             'fosamprenavir/r': ([{'group': [(10, 'F')], 'value': 15},
                                  {'group': [(20, 'T')], 'value': 5},
                                  [{'group': [(24, 'F')], 'value': 5},
                                   {'group': [(24, 'I')], 'value': 10},
                                   {'group': [(24, 'M')], 'value': 5}],
                                  {'group': [(32, 'I')], 'value': 30},
                                  {'group': [(33, 'F')], 'value': 10},
                                  [{'group': [(46, 'I')], 'value': 10},
                                   {'group': [(46, 'L')], 'value': 10},
                                   {'group': [(46, 'V')], 'value': 10}],
                                  [{'group': [(47, 'A')], 'value': 60},
                                   {'group': [(47, 'V')], 'value': 35}],
                                  [{'group': [(50, 'L')], 'value': -5},
                                   {'group': [(50, 'V')], 'value': 60}],
                                  [{'group': [(54, 'A')], 'value': 10},
                                   {'group': [(54, 'L')], 'value': 60},
                                   {'group': [(54, 'M')], 'value': 60},
                                   {'group': [(54, 'S')], 'value': 10},
                                   {'group': [(54, 'T')], 'value': 10},
                                   {'group': [(54, 'V')], 'value': 10}],
                                  [{'group': [(73, 'A')], 'value': 10},
                                   {'group': [(73, 'C')], 'value': 10},
                                   {'group': [(73, 'D')], 'value': 5},
                                   {'group': [(73, 'S')], 'value': 10},
                                   {'group': [(73, 'T')], 'value': 10},
                                   {'group': [(73, 'V')], 'value': 5}],
                                  {'group': [(74, 'P')], 'value': 10},
                                  {'group': [(76, 'V')], 'value': 60},
                                  [{'group': [(82, 'A')], 'value': 15},
                                   {'group': [(82, 'C')], 'value': 15},
                                   {'group': [(82, 'F')], 'value': 30},
                                   {'group': [(82, 'L')], 'value': 15},
                                   {'group': [(82, 'M')], 'value': 15},
                                   {'group': [(82, 'S')], 'value': 15},
                                   {'group': [(82, 'T')], 'value': 15}],
                                  [{'group': [(84, 'A')], 'value': 60},
                                   {'group': [(84, 'C')], 'value': 60},
                                   {'group': [(84, 'V')], 'value': 60}],
                                  {'group': [(88, 'S')], 'value': -10},
                                  {'group': [(89, 'V')], 'value': 10},
                                  {'group': [(90, 'M')], 'value': 20},
                                  {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                                  {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                                  {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                                  {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                                  {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                                  {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                                  {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                                  {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                                  {'group': [(46, 'ILV'), (76, 'V')], 'value': 10},
                                  {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                                  {'group': [(46, 'ILV'), (90, 'M')], 'value': 10},
                                  {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                                  {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                                  {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                                  {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                                  {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                                  {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                                  {'group': [(54, 'LM'), (89, 'V')], 'value': 5},
                                  {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                                  {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                                 {1: {'max': 9, 'min': -inf},
                                  2: {'max': 14, 'min': 10},
                                  3: {'max': 29, 'min': 15},
                                  4: {'max': 59, 'min': 30},
                                  5: {'max': inf, 'min': 60}}),
             'indinavir/r': ([{'group': [(10, 'F')], 'value': 10},
                              {'group': [(20, 'T')], 'value': 5},
                              [{'group': [(24, 'F')], 'value': 5},
                               {'group': [(24, 'I')], 'value': 15},
                               {'group': [(24, 'M')], 'value': 5}],
                              {'group': [(32, 'I')], 'value': 15},
                              {'group': [(33, 'F')], 'value': 5},
                              [{'group': [(46, 'I')], 'value': 10},
                               {'group': [(46, 'L')], 'value': 10},
                               {'group': [(46, 'V')], 'value': 10}],
                              [{'group': [(47, 'A')], 'value': 15},
                               {'group': [(47, 'V')], 'value': 15}],
                              [{'group': [(48, 'A')], 'value': 10},
                               {'group': [(48, 'L')], 'value': 10},
                               {'group': [(48, 'M')], 'value': 10},
                               {'group': [(48, 'Q')], 'value': 10},
                               {'group': [(48, 'S')], 'value': 10},
                               {'group': [(48, 'T')], 'value': 10},
                               {'group': [(48, 'V')], 'value': 10}],
                              {'group': [(50, 'L')], 'value': -5},
                              [{'group': [(54, 'A')], 'value': 15},
                               {'group': [(54, 'L')], 'value': 10},
                               {'group': [(54, 'M')], 'value': 15},
                               {'group': [(54, 'S')], 'value': 15},
                               {'group': [(54, 'T')], 'value': 15},
                               {'group': [(54, 'V')], 'value': 15}],
                              [{'group': [(73, 'A')], 'value': 15},
                               {'group': [(73, 'C')], 'value': 15},
                               {'group': [(73, 'D')], 'value': 5},
                               {'group': [(73, 'S')], 'value': 15},
                               {'group': [(73, 'T')], 'value': 15},
                               {'group': [(73, 'V')], 'value': 5}],
                              {'group': [(74, 'P')], 'value': 10},
                              {'group': [(76, 'V')], 'value': 30},
                              [{'group': [(82, 'A')], 'value': 30},
                               {'group': [(82, 'C')], 'value': 15},
                               {'group': [(82, 'F')], 'value': 30},
                               {'group': [(82, 'L')], 'value': 10},
                               {'group': [(82, 'M')], 'value': 30},
                               {'group': [(82, 'S')], 'value': 30},
                               {'group': [(82, 'T')], 'value': 30}],
                              {'group': [(83, 'D')], 'value': 10},
                              [{'group': [(84, 'A')], 'value': 60},
                               {'group': [(84, 'C')], 'value': 60},
                               {'group': [(84, 'V')], 'value': 60}],
                              {'group': [(88, 'S')], 'value': 15},
                              {'group': [(89, 'V')], 'value': 5},
                              {'group': [(90, 'M')], 'value': 30},
                              {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                              {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                              {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                              {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                              {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                              {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                              {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                              {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                              {'group': [(46, 'ILV'), (76, 'V')], 'value': 10},
                              {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                              {'group': [(46, 'ILV'), (90, 'M')], 'value': 10},
                              {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                              {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                              {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                              {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                              {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                              {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                              {'group': [(54, 'LM'), (89, 'V')], 'value': 5},
                              {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                              {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                             {1: {'max': 9, 'min': -inf},
                              2: {'max': 14, 'min': 10},
                              3: {'max': 29, 'min': 15},
                              4: {'max': 59, 'min': 30},
                              5: {'max': inf, 'min': 60}}),
             'lamivudine': ([{'group': [(62, 'V')], 'value': 5},
                             [{'group': [(65, 'N')], 'value': 15},
                              {'group': [(65, 'R')], 'value': 30}],
                             {'group': [(67, 'd')], 'value': 15},
                             {'group': [(68, 'd')], 'value': 15},
                             [{'group': [(69, 'i')], 'value': 30},
                              {'group': [(69, 'd')], 'value': 15}],
                             [{'group': [(70, 'E')], 'value': 10},
                              {'group': [(70, 'G')], 'value': 10},
                              {'group': [(70, 'N')], 'value': 10},
                              {'group': [(70, 'Q')], 'value': 10},
                              {'group': [(70, 'S')], 'value': 10},
                              {'group': [(70, 'T')], 'value': 10},
                              {'group': [(70, 'd')], 'value': 15}],
                             {'group': [(75, 'I')], 'value': 5},
                             {'group': [(77, 'L')], 'value': 5},
                             {'group': [(116, 'Y')], 'value': 5},
                             [{'group': [(151, 'L')], 'value': 10},
                              {'group': [(151, 'M')], 'value': 15}],
                             [{'group': [(184, 'I')], 'value': 60},
                              {'group': [(184, 'V')], 'value': 60}],
                             {'group': [(210, 'W'), (215, 'FY')], 'value': 5},
                             {'group': [(41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
                             {'group': [(41, 'L'), (215, 'FY')], 'value': 5},
                             {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                             {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')],
                              'value': 10},
                             {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 15}],
                            {1: {'max': 9, 'min': -inf},
                             2: {'max': 14, 'min': 10},
                             3: {'max': 29, 'min': 15},
                             4: {'max': 59, 'min': 30},
                             5: {'max': inf, 'min': 60}}),
             'lopinavir/r': ([{'group': [(10, 'F')], 'value': 5},
                              [{'group': [(24, 'F')], 'value': 5},
                               {'group': [(24, 'I')], 'value': 10},
                               {'group': [(24, 'M')], 'value': 5}],
                              {'group': [(32, 'I')], 'value': 15},
                              {'group': [(33, 'F')], 'value': 5},
                              [{'group': [(46, 'I')], 'value': 10},
                               {'group': [(46, 'L')], 'value': 10},
                               {'group': [(46, 'V')], 'value': 5}],
                              [{'group': [(47, 'A')], 'value': 60},
                               {'group': [(47, 'V')], 'value': 15}],
                              [{'group': [(48, 'A')], 'value': 10},
                               {'group': [(48, 'L')], 'value': 10},
                               {'group': [(48, 'M')], 'value': 10},
                               {'group': [(48, 'Q')], 'value': 10},
                               {'group': [(48, 'S')], 'value': 10},
                               {'group': [(48, 'T')], 'value': 10},
                               {'group': [(48, 'V')], 'value': 10}],
                              [{'group': [(50, 'L')], 'value': -10},
                               {'group': [(50, 'V')], 'value': 30}],
                              [{'group': [(54, 'A')], 'value': 15},
                               {'group': [(54, 'L')], 'value': 20},
                               {'group': [(54, 'M')], 'value': 20},
                               {'group': [(54, 'S')], 'value': 15},
                               {'group': [(54, 'T')], 'value': 15},
                               {'group': [(54, 'V')], 'value': 15}],
                              [{'group': [(73, 'A')], 'value': 5},
                               {'group': [(73, 'C')], 'value': 5},
                               {'group': [(73, 'D')], 'value': 5},
                               {'group': [(73, 'S')], 'value': 5},
                               {'group': [(73, 'T')], 'value': 5},
                               {'group': [(73, 'V')], 'value': 5}],
                              {'group': [(74, 'P')], 'value': 5},
                              {'group': [(76, 'V')], 'value': 30},
                              [{'group': [(82, 'A')], 'value': 30},
                               {'group': [(82, 'C')], 'value': 15},
                               {'group': [(82, 'F')], 'value': 30},
                               {'group': [(82, 'L')], 'value': 10},
                               {'group': [(82, 'M')], 'value': 25},
                               {'group': [(82, 'S')], 'value': 30},
                               {'group': [(82, 'T')], 'value': 30}],
                              [{'group': [(84, 'A')], 'value': 60},
                               {'group': [(84, 'C')], 'value': 30},
                               {'group': [(84, 'V')], 'value': 30}],
                              {'group': [(90, 'M')], 'value': 15},
                              {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                              {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                              {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                              {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                              {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                              {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                              {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                              {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                              {'group': [(46, 'ILV'), (76, 'V')], 'value': 10},
                              {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                              {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                              {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                              {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                              {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 5},
                              {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                              {'group': [(54, 'LM'), (89, 'V')], 'value': 5},
                              {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 5}],
                             {1: {'max': 9, 'min': -inf},
                              2: {'max': 14, 'min': 10},
                              3: {'max': 29, 'min': 15},
                              4: {'max': 59, 'min': 30},
                              5: {'max': inf, 'min': 60}}),
             'nelfinavir': ([{'group': [(10, 'F')], 'value': 15},
                             {'group': [(20, 'T')], 'value': 15},
                             {'group': [(23, 'I')], 'value': 15},
                             [{'group': [(24, 'F')], 'value': 10},
                              {'group': [(24, 'I')], 'value': 10},
                              {'group': [(24, 'M')], 'value': 10}],
                             {'group': [(30, 'N')], 'value': 60},
                             {'group': [(32, 'I')], 'value': 15},
                             {'group': [(33, 'F')], 'value': 10},
                             {'group': [(43, 'T')], 'value': 10},
                             [{'group': [(46, 'I')], 'value': 30},
                              {'group': [(46, 'L')], 'value': 20},
                              {'group': [(46, 'V')], 'value': 20}],
                             [{'group': [(47, 'A')], 'value': 30},
                              {'group': [(47, 'V')], 'value': 20}],
                             [{'group': [(48, 'A')], 'value': 30},
                              {'group': [(48, 'L')], 'value': 30},
                              {'group': [(48, 'M')], 'value': 30},
                              {'group': [(48, 'Q')], 'value': 30},
                              {'group': [(48, 'S')], 'value': 30},
                              {'group': [(48, 'T')], 'value': 30},
                              {'group': [(48, 'V')], 'value': 30}],
                             {'group': [(50, 'V')], 'value': 15},
                             {'group': [(53, 'L')], 'value': 10},
                             [{'group': [(54, 'A')], 'value': 20},
                              {'group': [(54, 'L')], 'value': 20},
                              {'group': [(54, 'M')], 'value': 20},
                              {'group': [(54, 'S')], 'value': 20},
                              {'group': [(54, 'T')], 'value': 20},
                              {'group': [(54, 'V')], 'value': 20}],
                             {'group': [(58, 'E')], 'value': 10},
                             [{'group': [(73, 'A')], 'value': 15},
                              {'group': [(73, 'C')], 'value': 15},
                              {'group': [(73, 'D')], 'value': 10},
                              {'group': [(73, 'S')], 'value': 15},
                              {'group': [(73, 'T')], 'value': 15},
                              {'group': [(73, 'V')], 'value': 10}],
                             {'group': [(74, 'P')], 'value': 20},
                             {'group': [(76, 'V')], 'value': 10},
                             [{'group': [(82, 'A')], 'value': 30},
                              {'group': [(82, 'C')], 'value': 30},
                              {'group': [(82, 'F')], 'value': 30},
                              {'group': [(82, 'L')], 'value': 10},
                              {'group': [(82, 'M')], 'value': 30},
                              {'group': [(82, 'S')], 'value': 30},
                              {'group': [(82, 'T')], 'value': 30}],
                             {'group': [(83, 'D')], 'value': 15},
                             [{'group': [(84, 'A')], 'value': 60},
                              {'group': [(84, 'C')], 'value': 60},
                              {'group': [(84, 'V')], 'value': 60}],
                             [{'group': [(88, 'D')], 'value': 60},
                              {'group': [(88, 'G')], 'value': 30},
                              {'group': [(88, 'S')], 'value': 60},
                              {'group': [(88, 'T')], 'value': 30}],
                             {'group': [(89, 'V')], 'value': 10},
                             {'group': [(90, 'M')], 'value': 60},
                             {'group': [(11, 'IL'), (32, 'I')], 'value': 5},
                             {'group': [(11, 'IL'), (54, 'LM')], 'value': 5},
                             {'group': [(32, 'I'), (47, 'AV')], 'value': 5},
                             {'group': [(32, 'I'), (54, 'LM')], 'value': 5},
                             {'group': [(32, 'I'), (76, 'V')], 'value': 5},
                             {'group': [(32, 'I'), (84, 'V')], 'value': 5},
                             {'group': [(32, 'I'), (89, 'V')], 'value': 5},
                             {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                             {'group': [(46, 'ILV'), (76, 'V')], 'value': 10},
                             {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                             {'group': [(46, 'ILV'), (90, 'M')], 'value': 10},
                             {'group': [(47, 'AV'), (54, 'LM')], 'value': 5},
                             {'group': [(47, 'AV'), (84, 'V')], 'value': 5},
                             {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                             {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                             {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                             {'group': [(54, 'LM'), (84, 'V')], 'value': 5},
                             {'group': [(54, 'LM'), (89, 'V')], 'value': 5},
                             {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                             {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                            {1: {'max': 9, 'min': -inf},
                             2: {'max': 14, 'min': 10},
                             3: {'max': 29, 'min': 15},
                             4: {'max': 59, 'min': 30},
                             5: {'max': inf, 'min': 60}}),
             'nevirapine': ([{'group': [(98, 'G')], 'value': 30},
                             [{'group': [(100, 'I')], 'value': 60},
                              {'group': [(100, 'V')], 'value': 30}],
                             [{'group': [(101, 'E')], 'value': 30},
                              {'group': [(101, 'H')], 'value': 15},
                              {'group': [(101, 'P')], 'value': 60}],
                             [{'group': [(103, 'H')], 'value': 60},
                              {'group': [(103, 'N')], 'value': 60},
                              {'group': [(103, 'S')], 'value': 60},
                              {'group': [(103, 'T')], 'value': 60}],
                             [{'group': [(106, 'A')], 'value': 60},
                              {'group': [(106, 'I')], 'value': 10},
                              {'group': [(106, 'M')], 'value': 60}],
                             {'group': [(108, 'I')], 'value': 15},
                             [{'group': [(138, 'G')], 'value': 10},
                              {'group': [(138, 'K')], 'value': 10},
                              {'group': [(138, 'Q')], 'value': 10},
                              {'group': [(138, 'R')], 'value': 10}],
                             [{'group': [(179, 'D')], 'value': 10},
                              {'group': [(179, 'E')], 'value': 10},
                              {'group': [(179, 'F')], 'value': 15},
                              {'group': [(179, 'L')], 'value': 10}],
                             [{'group': [(181, 'C')], 'value': 60},
                              {'group': [(181, 'F')], 'value': 60},
                              {'group': [(181, 'G')], 'value': 60},
                              {'group': [(181, 'I')], 'value': 60},
                              {'group': [(181, 'S')], 'value': 60},
                              {'group': [(181, 'V')], 'value': 60}],
                             [{'group': [(188, 'C')], 'value': 60},
                              {'group': [(188, 'F')], 'value': 60},
                              {'group': [(188, 'H')], 'value': 60},
                              {'group': [(188, 'L')], 'value': 60}],
                             [{'group': [(190, 'A')], 'value': 60},
                              {'group': [(190, 'C')], 'value': 60},
                              {'group': [(190, 'E')], 'value': 60},
                              {'group': [(190, 'Q')], 'value': 60},
                              {'group': [(190, 'S')], 'value': 60},
                              {'group': [(190, 'T')], 'value': 60},
                              {'group': [(190, 'V')], 'value': 60}],
                             {'group': [(221, 'Y')], 'value': 15},
                             {'group': [(225, 'H')], 'value': 45},
                             [{'group': [(227, 'C')], 'value': 45},
                              {'group': [(227, 'I')], 'value': 30},
                              {'group': [(227, 'L')], 'value': 30},
                              {'group': [(227, 'V')], 'value': 30}],
                             [{'group': [(230, 'I')], 'value': 30},
                              {'group': [(230, 'L')], 'value': 60}],
                             [{'group': [(238, 'N')], 'value': 10},
                              {'group': [(238, 'T')], 'value': 30}],
                             {'group': [(318, 'F')], 'value': 30},
                             {'group': [(348, 'I')], 'value': 15},
                             {'group': [(101, 'E'), (181, 'C')], 'value': 5},
                             {'group': [(103, 'R'), (179, 'D')], 'value': 20},
                             {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                             {'group': [(98, 'G'), (227, 'C')], 'value': 15}],
                            {1: {'max': 9, 'min': -inf},
                             2: {'max': 14, 'min': 10},
                             3: {'max': 29, 'min': 15},
                             4: {'max': 59, 'min': 30},
                             5: {'max': inf, 'min': 60}}),
             'raltegravir': ([{'group': [(51, 'Y')], 'value': 15},
                              [{'group': [(66, 'A')], 'value': 15},
                               {'group': [(66, 'I')], 'value': 15},
                               {'group': [(66, 'K')], 'value': 60}],
                              [{'group': [(92, 'G')], 'value': 15},
                               {'group': [(92, 'Q')], 'value': 30},
                               {'group': [(92, 'V')], 'value': 30}],
                              {'group': [(95, 'K')], 'value': 10},
                              {'group': [(97, 'A')], 'value': 10},
                              {'group': [(118, 'R')], 'value': 30},
                              {'group': [(121, 'Y')], 'value': 60},
                              [{'group': [(138, 'A')], 'value': 15},
                               {'group': [(138, 'K')], 'value': 15},
                               {'group': [(138, 'T')], 'value': 15}],
                              [{'group': [(140, 'A')], 'value': 30},
                               {'group': [(140, 'C')], 'value': 30},
                               {'group': [(140, 'S')], 'value': 30}],
                              [{'group': [(143, 'A')], 'value': 60},
                               {'group': [(143, 'C')], 'value': 60},
                               {'group': [(143, 'G')], 'value': 60},
                               {'group': [(143, 'H')], 'value': 60},
                               {'group': [(143, 'K')], 'value': 60},
                               {'group': [(143, 'R')], 'value': 60},
                               {'group': [(143, 'S')], 'value': 60}],
                              [{'group': [(148, 'H')], 'value': 60},
                               {'group': [(148, 'K')], 'value': 60},
                               {'group': [(148, 'N')], 'value': 10},
                               {'group': [(148, 'R')], 'value': 60}],
                              [{'group': [(151, 'A')], 'value': 15},
                               {'group': [(151, 'L')], 'value': 30}],
                              [{'group': [(155, 'H')], 'value': 60},
                               {'group': [(155, 'S')], 'value': 30},
                               {'group': [(155, 'T')], 'value': 30}],
                              {'group': [(157, 'Q')], 'value': 10},
                              [{'group': [(163, 'K')], 'value': 15},
                               {'group': [(163, 'R')], 'value': 15}],
                              {'group': [(230, 'R')], 'value': 20},
                              {'group': [(232, 'N')], 'value': 10},
                              {'group': [(263, 'K')], 'value': 25},
                              {'group': [(138, 'AKT'), (118, 'R')], 'value': 10},
                              {'group': [(138, 'AKT'), (140, 'ACS')], 'value': 15},
                              {'group': [(140, 'ACS'), (148, 'HKR'), (149, 'A')],
                               'value': 10},
                              {'group': [(74, 'FIM'), (118, 'R')], 'value': 10},
                              {'group': [(74, 'FIM'), (148, 'HKR')], 'value': 15},
                              {'group': [(97, 'A'), (118, 'R')], 'value': 10}],
                             {1: {'max': 9, 'min': -inf},
                              2: {'max': 14, 'min': 10},
                              3: {'max': 29, 'min': 15},
                              4: {'max': 59, 'min': 30},
                              5: {'max': inf, 'min': 60}}),
             'rilpivirine': ([{'group': [(98, 'G')], 'value': 15},
                              [{'group': [(100, 'I')], 'value': 60},
                               {'group': [(100, 'V')], 'value': 15}],
                              [{'group': [(101, 'E')], 'value': 45},
                               {'group': [(101, 'H')], 'value': 10},
                               {'group': [(101, 'P')], 'value': 60}],
                              {'group': [(106, 'I')], 'value': 10},
                              [{'group': [(138, 'A')], 'value': 15},
                               {'group': [(138, 'G')], 'value': 15},
                               {'group': [(138, 'K')], 'value': 45},
                               {'group': [(138, 'Q')], 'value': 15},
                               {'group': [(138, 'R')], 'value': 15}],
                              [{'group': [(179, 'D')], 'value': 10},
                               {'group': [(179, 'E')], 'value': 10},
                               {'group': [(179, 'F')], 'value': 15},
                               {'group': [(179, 'L')], 'value': 15}],
                              [{'group': [(181, 'C')], 'value': 45},
                               {'group': [(181, 'F')], 'value': 30},
                               {'group': [(181, 'G')], 'value': 30},
                               {'group': [(181, 'I')], 'value': 60},
                               {'group': [(181, 'S')], 'value': 30},
                               {'group': [(181, 'V')], 'value': 60}],
                              [{'group': [(188, 'F')], 'value': 30},
                               {'group': [(188, 'L')], 'value': 60}],
                              [{'group': [(190, 'A')], 'value': 15},
                               {'group': [(190, 'C')], 'value': 10},
                               {'group': [(190, 'E')], 'value': 60},
                               {'group': [(190, 'Q')], 'value': 45},
                               {'group': [(190, 'S')], 'value': 15},
                               {'group': [(190, 'T')], 'value': 10},
                               {'group': [(190, 'V')], 'value': 10}],
                              {'group': [(221, 'Y')], 'value': 15},
                              {'group': [(227, 'C')], 'value': 45},
                              [{'group': [(230, 'I')], 'value': 30},
                               {'group': [(230, 'L')], 'value': 60}],
                              {'group': [(101, 'E'), (184, 'I')], 'value': 15},
                              {'group': [(103, 'R'), (179, 'D')], 'value': 15},
                              {'group': [(138, 'K'), (184, 'I')], 'value': 15},
                              {'group': [(181, 'C'), (190, 'ACSTV')], 'value': 10},
                              {'group': [(98, 'G'), (181, 'C')], 'value': 5},
                              {'group': [(98, 'G'), (227, 'C')], 'value': 15},
                              [{'group': [(179, 'F'), (181, 'C')], 'value': 15},
                               {'group': [(179, 'T'), (181, 'C')], 'value': 10}]],
                             {1: {'max': 9, 'min': -inf},
                              2: {'max': 14, 'min': 10},
                              3: {'max': 29, 'min': 15},
                              4: {'max': 59, 'min': 30},
                              5: {'max': inf, 'min': 60}}),
             'saquinavir/r': ([{'group': [(20, 'T')], 'value': 5},
                               [{'group': [(24, 'F')], 'value': 5},
                                {'group': [(24, 'I')], 'value': 10},
                                {'group': [(24, 'M')], 'value': 5}],
                               {'group': [(33, 'F')], 'value': 5},
                               [{'group': [(46, 'I')], 'value': 10},
                                {'group': [(46, 'L')], 'value': 10},
                                {'group': [(46, 'V')], 'value': 5}],
                               [{'group': [(48, 'A')], 'value': 60},
                                {'group': [(48, 'L')], 'value': 60},
                                {'group': [(48, 'M')], 'value': 60},
                                {'group': [(48, 'Q')], 'value': 60},
                                {'group': [(48, 'S')], 'value': 60},
                                {'group': [(48, 'T')], 'value': 60},
                                {'group': [(48, 'V')], 'value': 60}],
                               [{'group': [(50, 'L')], 'value': -5},
                                {'group': [(50, 'V')], 'value': 15}],
                               {'group': [(53, 'L')], 'value': 15},
                               [{'group': [(54, 'A')], 'value': 15},
                                {'group': [(54, 'L')], 'value': 15},
                                {'group': [(54, 'M')], 'value': 15},
                                {'group': [(54, 'S')], 'value': 15},
                                {'group': [(54, 'T')], 'value': 15},
                                {'group': [(54, 'V')], 'value': 15}],
                               [{'group': [(73, 'A')], 'value': 15},
                                {'group': [(73, 'C')], 'value': 15},
                                {'group': [(73, 'D')], 'value': 10},
                                {'group': [(73, 'S')], 'value': 15},
                                {'group': [(73, 'T')], 'value': 15},
                                {'group': [(73, 'V')], 'value': 10}],
                               {'group': [(74, 'P')], 'value': 10},
                               [{'group': [(82, 'A')], 'value': 15},
                                {'group': [(82, 'C')], 'value': 15},
                                {'group': [(82, 'F')], 'value': 10},
                                {'group': [(82, 'L')], 'value': 10},
                                {'group': [(82, 'M')], 'value': 15},
                                {'group': [(82, 'S')], 'value': 15},
                                {'group': [(82, 'T')], 'value': 15}],
                               {'group': [(83, 'D')], 'value': 10},
                               [{'group': [(84, 'A')], 'value': 60},
                                {'group': [(84, 'C')], 'value': 60},
                                {'group': [(84, 'V')], 'value': 60}],
                               [{'group': [(88, 'D')], 'value': 10},
                                {'group': [(88, 'S')], 'value': 15}],
                               {'group': [(90, 'M')], 'value': 45},
                               {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5},
                               {'group': [(46, 'ILV'), (82, 'ACFLMST')], 'value': 10},
                               {'group': [(46, 'ILV'), (90, 'M')], 'value': 5},
                               {'group': [(53, 'L'), (90, 'M')], 'value': 10},
                               {'group': [(54, 'ALMSTV'), (82, 'ACFLMST')], 'value': 10},
                               {'group': [(54, 'ALMSTV'), (90, 'M')], 'value': 10},
                               {'group': [(73, 'ACSTV'), (90, 'M')], 'value': 10},
                               {'group': [(82, 'ACFLMST'), (90, 'M')], 'value': 10}],
                              {1: {'max': 9, 'min': -inf},
                               2: {'max': 14, 'min': 10},
                               3: {'max': 29, 'min': 15},
                               4: {'max': 59, 'min': 30},
                               5: {'max': inf, 'min': 60}}),
             'stavudine': ([{'group': [(41, 'L')], 'value': 15},
                            {'group': [(62, 'V')], 'value': 5},
                            [{'group': [(65, 'E')], 'value': 10},
                             {'group': [(65, 'N')], 'value': 30},
                             {'group': [(65, 'R')], 'value': 60}],
                            [{'group': [(67, 'E')], 'value': 10},
                             {'group': [(67, 'G')], 'value': 10},
                             {'group': [(67, 'H')], 'value': 10},
                             {'group': [(67, 'N')], 'value': 15},
                             {'group': [(67, 'S')], 'value': 10},
                             {'group': [(67, 'T')], 'value': 10},
                             {'group': [(67, 'd')], 'value': 30}],
                            {'group': [(68, 'd')], 'value': 30},
                            [{'group': [(69, 'D')], 'value': 10},
                             {'group': [(69, 'G')], 'value': 10},
                             {'group': [(69, 'i')], 'value': 60},
                             {'group': [(69, 'd')], 'value': 30}],
                            [{'group': [(70, 'E')], 'value': 15},
                             {'group': [(70, 'G')], 'value': 15},
                             {'group': [(70, 'N')], 'value': 15},
                             {'group': [(70, 'Q')], 'value': 15},
                             {'group': [(70, 'R')], 'value': 15},
                             {'group': [(70, 'S')], 'value': 15},
                             {'group': [(70, 'T')], 'value': 15},
                             {'group': [(70, 'd')], 'value': 30}],
                            [{'group': [(75, 'A')], 'value': 30},
                             {'group': [(75, 'I')], 'value': 5},
                             {'group': [(75, 'M')], 'value': 30},
                             {'group': [(75, 'S')], 'value': 30},
                             {'group': [(75, 'T')], 'value': 60}],
                            {'group': [(77, 'L')], 'value': 10},
                            {'group': [(116, 'Y')], 'value': 10},
                            [{'group': [(151, 'L')], 'value': 30},
                             {'group': [(151, 'M')], 'value': 60}],
                            [{'group': [(184, 'I')], 'value': -10},
                             {'group': [(184, 'V')], 'value': -10}],
                            {'group': [(210, 'W')], 'value': 15},
                            [{'group': [(215, 'A')], 'value': 20},
                             {'group': [(215, 'C')], 'value': 20},
                             {'group': [(215, 'D')], 'value': 20},
                             {'group': [(215, 'E')], 'value': 20},
                             {'group': [(215, 'F')], 'value': 40},
                             {'group': [(215, 'I')], 'value': 20},
                             {'group': [(215, 'L')], 'value': 20},
                             {'group': [(215, 'N')], 'value': 20},
                             {'group': [(215, 'S')], 'value': 20},
                             {'group': [(215, 'V')], 'value': 20},
                             {'group': [(215, 'Y')], 'value': 40}],
                            [{'group': [(219, 'E')], 'value': 10},
                             {'group': [(219, 'N')], 'value': 10},
                             {'group': [(219, 'Q')], 'value': 10},
                             {'group': [(219, 'R')], 'value': 10},
                             {'group': [(219, 'W')], 'value': 10}],
                            {'group': [(151, 'M'), (184, 'IV')], 'value': 10},
                            {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                             'value': 5},
                            {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                            {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                             'value': 5},
                            {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                            {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')],
                             'value': 5},
                            {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                            {'group': [(70, 'EGNQST'), (184, 'IV')], 'value': 10},
                            {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                            {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 10},
                            [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                             {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                            [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                             {'group': [(41, 'L'), (215, 'FY')], 'value': 10}]],
                           {1: {'max': 9, 'min': -inf},
                            2: {'max': 14, 'min': 10},
                            3: {'max': 29, 'min': 15},
                            4: {'max': 59, 'min': 30},
                            5: {'max': inf, 'min': 60}}),
             'tenofovir': ([{'group': [(41, 'L')], 'value': 5},
                            {'group': [(62, 'V')], 'value': 5},
                            [{'group': [(65, 'E')], 'value': 10},
                             {'group': [(65, 'N')], 'value': 30},
                             {'group': [(65, 'R')], 'value': 60}],
                            [{'group': [(67, 'E')], 'value': 5},
                             {'group': [(67, 'G')], 'value': 5},
                             {'group': [(67, 'H')], 'value': 5},
                             {'group': [(67, 'N')], 'value': 5},
                             {'group': [(67, 'S')], 'value': 5},
                             {'group': [(67, 'T')], 'value': 5},
                             {'group': [(67, 'd')], 'value': 30}],
                            {'group': [(68, 'd')], 'value': 15},
                            [{'group': [(69, 'G')], 'value': 5},
                             {'group': [(69, 'i')], 'value': 60},
                             {'group': [(69, 'd')], 'value': 15}],
                            [{'group': [(70, 'E')], 'value': 15},
                             {'group': [(70, 'G')], 'value': 15},
                             {'group': [(70, 'N')], 'value': 15},
                             {'group': [(70, 'Q')], 'value': 15},
                             {'group': [(70, 'R')], 'value': 5},
                             {'group': [(70, 'S')], 'value': 15},
                             {'group': [(70, 'T')], 'value': 15},
                             {'group': [(70, 'd')], 'value': 15}],
                            {'group': [(74, 'I')], 'value': 5},
                            {'group': [(75, 'I')], 'value': 5},
                            {'group': [(77, 'L')], 'value': 5},
                            {'group': [(115, 'F')], 'value': 15},
                            {'group': [(116, 'Y')], 'value': 5},
                            [{'group': [(151, 'L')], 'value': 10},
                             {'group': [(151, 'M')], 'value': 15}],
                            [{'group': [(184, 'I')], 'value': -10},
                             {'group': [(184, 'V')], 'value': -10}],
                            {'group': [(210, 'W')], 'value': 5},
                            [{'group': [(215, 'A')], 'value': 5},
                             {'group': [(215, 'C')], 'value': 5},
                             {'group': [(215, 'D')], 'value': 5},
                             {'group': [(215, 'E')], 'value': 5},
                             {'group': [(215, 'F')], 'value': 10},
                             {'group': [(215, 'I')], 'value': 5},
                             {'group': [(215, 'L')], 'value': 5},
                             {'group': [(215, 'N')], 'value': 5},
                             {'group': [(215, 'S')], 'value': 5},
                             {'group': [(215, 'V')], 'value': 5},
                             {'group': [(215, 'Y')], 'value': 10}],
                            [{'group': [(219, 'E')], 'value': 5},
                             {'group': [(219, 'N')], 'value': 5},
                             {'group': [(219, 'Q')], 'value': 5},
                             {'group': [(219, 'R')], 'value': 5}],
                            {'group': [(151, 'M'), (184, 'IV')], 'value': 10},
                            {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                             'value': 5},
                            {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                            {'group': [(41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
                            {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                             'value': 5},
                            {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                            {'group': [(65, 'R'), (151, 'M')], 'value': 10},
                            {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')],
                             'value': 5},
                            {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
                            {'group': [(70, 'EGNQST'), (184, 'IV')], 'value': 10},
                            {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                            {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 15},
                            [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                             {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                            [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                             {'group': [(41, 'L'), (215, 'FY')], 'value': 10}]],
                           {1: {'max': 9, 'min': -inf},
                            2: {'max': 14, 'min': 10},
                            3: {'max': 29, 'min': 15},
                            4: {'max': 59, 'min': 30},
                            5: {'max': inf, 'min': 60}}),
             'tipranavir/r': ([{'group': [(24, 'I')], 'value': -5},
                               {'group': [(32, 'I')], 'value': 5},
                               {'group': [(33, 'F')], 'value': 10},
                               {'group': [(43, 'T')], 'value': 10},
                               [{'group': [(46, 'I')], 'value': 5},
                                {'group': [(46, 'L')], 'value': 10},
                                {'group': [(46, 'V')], 'value': 5}],
                               [{'group': [(47, 'A')], 'value': 30},
                                {'group': [(47, 'V')], 'value': 30}],
                               [{'group': [(50, 'L')], 'value': -5},
                                {'group': [(50, 'V')], 'value': -5}],
                               [{'group': [(54, 'A')], 'value': 20},
                                {'group': [(54, 'L')], 'value': -10},
                                {'group': [(54, 'M')], 'value': 20},
                                {'group': [(54, 'S')], 'value': 20},
                                {'group': [(54, 'T')], 'value': 20},
                                {'group': [(54, 'V')], 'value': 20}],
                               {'group': [(58, 'E')], 'value': 15},
                               {'group': [(74, 'P')], 'value': 25},
                               {'group': [(76, 'V')], 'value': -5},
                               [{'group': [(82, 'C')], 'value': 10},
                                {'group': [(82, 'L')], 'value': 45},
                                {'group': [(82, 'M')], 'value': 10},
                                {'group': [(82, 'S')], 'value': 30},
                                {'group': [(82, 'T')], 'value': 45}],
                               {'group': [(83, 'D')], 'value': 25},
                               [{'group': [(84, 'A')], 'value': 60},
                                {'group': [(84, 'C')], 'value': 30},
                                {'group': [(84, 'V')], 'value': 30}],
                               {'group': [(46, 'IL'), (84, 'V'), (90, 'M')], 'value': 5}],
                              {1: {'max': 9, 'min': -inf},
                               2: {'max': 14, 'min': 10},
                               3: {'max': 29, 'min': 15},
                               4: {'max': 59, 'min': 30},
                               5: {'max': inf, 'min': 60}}),
             'zidovudine': ([{'group': [(41, 'L')], 'value': 15},
                             {'group': [(62, 'V')], 'value': 5},
                             [{'group': [(65, 'N')], 'value': -10},
                              {'group': [(65, 'R')], 'value': -15}],
                             [{'group': [(67, 'E')], 'value': 10},
                              {'group': [(67, 'G')], 'value': 10},
                              {'group': [(67, 'H')], 'value': 10},
                              {'group': [(67, 'N')], 'value': 15},
                              {'group': [(67, 'S')], 'value': 10},
                              {'group': [(67, 'T')], 'value': 10},
                              {'group': [(67, 'd')], 'value': 30}],
                             [{'group': [(69, 'G')], 'value': 5},
                              {'group': [(69, 'i')], 'value': 60}],
                             [{'group': [(70, 'E')], 'value': -10},
                              {'group': [(70, 'G')], 'value': -10},
                              {'group': [(70, 'R')], 'value': 30}],
                             [{'group': [(75, 'I')], 'value': 5},
                              {'group': [(75, 'M')], 'value': 10}],
                             {'group': [(77, 'L')], 'value': 10},
                             {'group': [(116, 'Y')], 'value': 10},
                             [{'group': [(151, 'L')], 'value': 30},
                              {'group': [(151, 'M')], 'value': 60}],
                             [{'group': [(184, 'I')], 'value': -10},
                              {'group': [(184, 'V')], 'value': -10}],
                             {'group': [(210, 'W')], 'value': 15},
                             [{'group': [(215, 'A')], 'value': 20},
                              {'group': [(215, 'C')], 'value': 20},
                              {'group': [(215, 'D')], 'value': 20},
                              {'group': [(215, 'E')], 'value': 20},
                              {'group': [(215, 'F')], 'value': 40},
                              {'group': [(215, 'I')], 'value': 20},
                              {'group': [(215, 'L')], 'value': 20},
                              {'group': [(215, 'N')], 'value': 20},
                              {'group': [(215, 'S')], 'value': 20},
                              {'group': [(215, 'V')], 'value': 20},
                              {'group': [(215, 'Y')], 'value': 40}],
                             [{'group': [(219, 'E')], 'value': 10},
                              {'group': [(219, 'N')], 'value': 10},
                              {'group': [(219, 'Q')], 'value': 10},
                              {'group': [(219, 'R')], 'value': 10},
                              {'group': [(219, 'W')], 'value': 10}],
                             {'group': [(151, 'M'), (184, 'IV')], 'value': 10},
                             {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')],
                              'value': 5},
                             {'group': [(41, 'L'), (210, 'W')], 'value': 10},
                             {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')],
                              'value': 5},
                             {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
                             {'group': [(65, 'R'), (151, 'M')], 'value': 10},
                             {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')],
                              'value': 5},
                             {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')],
                              'value': 10},
                             {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
                             {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 10},
                             [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
                              {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
                             [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
                              {'group': [(41, 'L'), (215, 'FY')], 'value': 10}]],
                            {1: {'max': 9, 'min': -inf},
                             2: {'max': 14, 'min': 10},
                             3: {'max': 29, 'min': 15},
                             4: {'max': 59, 'min': 30},
                             5: {'max': inf, 'min': 60}})}
        
        res_drugs = self.algorithm.parse_drugs(self.root)

        self.assertEqual(exp_drugs, res_drugs)

    def testParseCondition(self):
        # Setting params
        condition = self.root.find('.//CONDITION').text
        exp_drms = \
            [{'group': [(41, 'L')], 'value': 5},
             {'group': [(62, 'V')], 'value': 5},
             [{'group': [(65, 'E')], 'value': 10},
              {'group': [(65, 'N')], 'value': 30},
              {'group': [(65, 'R')], 'value': 45}],
             [{'group': [(67, 'E')], 'value': 5},
              {'group': [(67, 'G')], 'value': 5},
              {'group': [(67, 'H')], 'value': 5},
              {'group': [(67, 'N')], 'value': 5},
              {'group': [(67, 'S')], 'value': 5},
              {'group': [(67, 'T')], 'value': 5},
              {'group': [(67, 'd')], 'value': 30}],
             {'group': [(68, 'd')], 'value': 15},
             [{'group': [(69, 'G')], 'value': 10},
              {'group': [(69, 'i')], 'value': 60},
              {'group': [(69, 'd')], 'value': 15}],
             [{'group': [(70, 'E')], 'value': 15},
              {'group': [(70, 'G')], 'value': 15},
              {'group': [(70, 'N')], 'value': 15},
              {'group': [(70, 'Q')], 'value': 15},
              {'group': [(70, 'R')], 'value': 5},
              {'group': [(70, 'S')], 'value': 15},
              {'group': [(70, 'T')], 'value': 15},
              {'group': [(70, 'd')], 'value': 15}],
             [{'group': [(74, 'I')], 'value': 30}, {'group': [(74, 'V')], 'value': 30}],
             {'group': [(75, 'I')], 'value': 5},
             {'group': [(77, 'L')], 'value': 5},
             {'group': [(115, 'F')], 'value': 60},
             {'group': [(116, 'Y')], 'value': 10},
             [{'group': [(151, 'L')], 'value': 30}, {'group': [(151, 'M')], 'value': 60}],
             [{'group': [(184, 'I')], 'value': 15}, {'group': [(184, 'V')], 'value': 15}],
             {'group': [(210, 'W')], 'value': 5},
             [{'group': [(215, 'A')], 'value': 5},
              {'group': [(215, 'C')], 'value': 5},
              {'group': [(215, 'D')], 'value': 5},
              {'group': [(215, 'E')], 'value': 5},
              {'group': [(215, 'F')], 'value': 10},
              {'group': [(215, 'I')], 'value': 5},
              {'group': [(215, 'L')], 'value': 5},
              {'group': [(215, 'N')], 'value': 5},
              {'group': [(215, 'S')], 'value': 5},
              {'group': [(215, 'V')], 'value': 5},
              {'group': [(215, 'Y')], 'value': 10}],
             [{'group': [(219, 'E')], 'value': 5},
              {'group': [(219, 'N')], 'value': 5},
              {'group': [(219, 'Q')], 'value': 5},
              {'group': [(219, 'R')], 'value': 5}],
             {'group': [(40, 'F'), (41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
             {'group': [(41, 'L'), (210, 'W')], 'value': 10},
             {'group': [(41, 'L'), (210, 'W'), (215, 'FY')], 'value': 5},
             {'group': [(41, 'L'), (44, 'AD'), (210, 'W'), (215, 'FY')], 'value': 5},
             {'group': [(41, 'L'), (67, 'EGN'), (215, 'FY')], 'value': 5},
             {'group': [(67, 'EGN'), (215, 'FY'), (219, 'ENQR')], 'value': 5},
             {'group': [(67, 'EGN'), (70, 'R'), (184, 'IV'), (219, 'ENQR')], 'value': 20},
             {'group': [(67, 'EGN'), (70, 'R'), (219, 'ENQR')], 'value': 10},
             {'group': [(70, 'R'), (215, 'FY')], 'value': 5},
             {'group': [(74, 'IV'), (184, 'IV')], 'value': 15},
             {'group': [(77, 'L'), (116, 'Y'), (151, 'M')], 'value': 10},
             [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5},
              {'group': [(210, 'W'), (215, 'FY')], 'value': 10}],
             [{'group': [(41, 'L'), (215, 'ACDEILNSV')], 'value': 5},
              {'group': [(41, 'L'), (215, 'FY')], 'value': 15}]]
        
        res_drms = self.algorithm.parse_condition(condition)

        self.assertEqual(exp_drms, res_drms)

    def testParseScores(self):
        # Setting params
        drm_lib = []
        drm = 'MAX ((210W AND 215ACDEILNSV) => 5, (210W AND 215FY) => 10)'
        chunk = '(210W AND 215ACDEILNSV) => 5'
        iter = 0

        exp_drm_lib = [{'group': [(210, 'W'), (215, 'ACDEILNSV')], 'value': 5}]
        self.algorithm._parse_scores(drm_lib, drm, chunk, iter)
        self.assertEqual(exp_drm_lib, drm_lib)

        # if iter parameter is set to anything greater than 1
        with self.assertRaises(IndexError) as err: 
            self.algorithm._parse_scores(drm_lib, drm, chunk, 3)

        with self.assertRaises(IndexError) as err: 
            self.algorithm._parse_scores(drm_lib, drm, chunk, 2)

        # Setting params
        drm_lib = []
        drm = '(67EGN AND 70R AND 184IV AND 219ENQR) => 20'
        chunk = '(67EGN AND 70R AND 184IV AND 219ENQR) => 20'
        iter = 0

        exp_drm_lib = [{'group': [(67, 'EGN'), (70, 'R'), (184, 'IV'), (219, 'ENQR')], 'value': 20}]
        self.algorithm._parse_scores(drm_lib, drm, chunk, iter)
        self.assertEqual(exp_drm_lib, drm_lib)

        # if iter parameter is set to anything greater than 0
        with self.assertRaises(IndexError) as err: 
            self.algorithm._parse_scores(drm_lib, drm, chunk, 1)

        with self.assertRaises(IndexError) as err: 
            self.algorithm._parse_scores(drm_lib, drm, chunk, 2)

    def testParseComments(self):
        exp_comments = \
            {'IN': {'118ACDEFHIKLMNPQSTVWYid': 'IN118ACDEFHIKLMNPQSTVWY_-',
                    '118R': 'IN118R',
                    '119R': 'IN119R',
                    '121ACDEGHIKLMNPQRSTVWid': 'IN121ACDEGHIKLMNPQRSTVW_-',
                    '121Y': 'IN121Y',
                    '128T': 'IN128T',
                    '138CFGHILMNPQRSVWYid': 'IN138CFGHILMNPQRSVWY_-',
                    '138D': 'IN138D',
                    '138KAT': 'IN138KAT',
                    '140DEFHIKLMNPQRTVWYid': 'IN140DEFHIKLMNPQRTVWY_-',
                    '140SAC': 'IN140SAC',
                    '142T': 'IN142T',
                    '143CR': 'IN143CR',
                    '143DEFILMNPQTVWid': 'IN143DEFILMNPQTVW_-',
                    '143H': 'IN143H',
                    '143KGSA': 'IN143KGSA',
                    '145ACDEFGHIKLMNQRTVWYid': 'IN145ACDEFGHIKLMNQRTVWY_-',
                    '145S': 'IN145S',
                    '146ACDEFGHIKLMNRSTVWYid': 'IN146ACDEFGHIKLMNRSTVWY_-',
                    '146P': 'IN146P',
                    '147ACDEFHIKLMNPQRTVWYid': 'IN147ACDEFHIKLMNPQRTVWY_-',
                    '147G': 'IN147G',
                    '148ACDEFGILMPSTVWYid': 'IN148ACDEFGILMPSTVWY_-',
                    '148HKR': 'IN148HKR',
                    '148N': 'IN148N',
                    '149A': 'IN149A',
                    '151A': 'IN151A',
                    '151CDEFGHKMNPQRSTWYid': 'IN151CDEFGHKMNPQRSTWY_-',
                    '151I': 'IN151I',
                    '151L': 'IN151L',
                    '153ACDEGHIKLMNPQRTVWid': 'IN153ACDEGHIKLMNPQRTVW_-',
                    '153YF': 'IN153YF',
                    '155ACDEFGIKLMPQRVWYid': 'IN155ACDEFGIKLMPQRVWY_-',
                    '155H': 'IN155H',
                    '155ST': 'IN155ST',
                    '157Q': 'IN157Q',
                    '163RK': 'IN163RK',
                    '230N': 'IN230N',
                    '230R': 'IN230R',
                    '232N': 'IN232N',
                    '263ACDEFGHILMNPQSTVWYid': 'IN263ACDEFGHILMNPQSTVWY_-',
                    '263K': 'IN263K',
                    '50I': 'IN50I',
                    '51ACDEFGIKLMNPQRSTVWid': 'IN51ACDEFGIKLMNPQRSTVW_-',
                    '51Y': 'IN51Y',
                    '66A': 'IN66A',
                    '66CDEFGHLMNPQRSVWYid': 'IN66CDEFGHLMNPQRSVWY_-',
                    '66I': 'IN66I',
                    '66K': 'IN66K',
                    '74MIF': 'IN74MIF',
                    '92ACDFHIKLMNPRSTWYid': 'IN92ACDFHIKLMNPRSTWY_-',
                    '92G': 'IN92G',
                    '92Q': 'IN92Q',
                    '92V': 'IN92V',
                    '95K': 'IN95K',
                    '97A': 'IN97A'},
             'PR': {'10ACDEGHKMNPQSTWid': 'PR10ACDEGHKMNPQSTW_-',
                    '10F': 'PR10F',
                    '10IV': 'PR10IV',
                    '10RY': 'PR10RY',
                    '11IL': 'PR11IL',
                    '20I': 'PR20I',
                    '20MV': 'PR20MV',
                    '20R': 'PR20R',
                    '20T': 'PR20T',
                    '23I': 'PR23I',
                    '24ACDEGHKNPQRSTVWYid': 'PR24ACDEGHKNPQRSTVWY_-',
                    '24FM': 'PR24FM',
                    '24I': 'PR24I',
                    '30ACEFGHIKLMPQRSTVWYid': 'PR30ACEFGHIKLMPQRSTVWY_-',
                    '30N': 'PR30N',
                    '32ACDEFGHKLMNPQRSTWYid': 'PR32ACDEFGHKLMNPQRSTWY_-',
                    '32I': 'PR32I',
                    '33F': 'PR33F',
                    '33i': 'PR33_',
                    '34i': 'PR34_',
                    '35i': 'PR35_',
                    '36i': 'PR36_',
                    '37i': 'PR37_',
                    '38i': 'PR38_',
                    '39i': 'PR39_',
                    '40i': 'PR40_',
                    '41i': 'PR41_',
                    '43T': 'PR43T',
                    '46ACDEFGHKNPQRSTWYid': 'PR46ACDEFGHKNPQRSTWY_-',
                    '46IL': 'PR46IL',
                    '46V': 'PR46V',
                    '47A': 'PR47A',
                    '47CDEFGHKLMNPQRSTWYid': 'PR47CDEFGHKLMNPQRSTWY_-',
                    '47V': 'PR47V',
                    '48ASTQL': 'PR48ASTQL',
                    '48CDEFHIKNPRWYid': 'PR48CDEFHIKNPRWY_-',
                    '48M': 'PR48M',
                    '48V': 'PR48V',
                    '50ACDEFGHKMNPQRSTWYid': 'PR50ACDEFGHKMNPQRSTWY_-',
                    '50L': 'PR50L',
                    '50V': 'PR50V',
                    '53L': 'PR53L',
                    '53Y': 'PR53Y',
                    '54ATS': 'PR54ATS',
                    '54CDEFGHKNPQRWYid': 'PR54CDEFGHKNPQRWY_-',
                    '54LM': 'PR54LM',
                    '54V': 'PR54V',
                    '58E': 'PR58E',
                    '71IL': 'PR71IL',
                    '71TV': 'PR71TV',
                    '73EFHIKLMNPQRWYid': 'PR73EFHIKLMNPQRWY_-',
                    '73STCADV': 'PR73STCADV',
                    '74P': 'PR74P',
                    '74S': 'PR74S',
                    '76ACDEFGHIKMNPQRSTWYid': 'PR76ACDEFGHIKMNPQRSTWY_-',
                    '76V': 'PR76V',
                    '82A': 'PR82A',
                    '82C': 'PR82C',
                    '82DEGHKNPQRWYid': 'PR82DEGHKNPQRWY_-',
                    '82F': 'PR82F',
                    '82I': 'PR82I',
                    '82L': 'PR82L',
                    '82M': 'PR82M',
                    '82TS': 'PR82TS',
                    '83D': 'PR83D',
                    '84AC': 'PR84AC',
                    '84DEFGHKLMNPQRSTWYid': 'PR84DEFGHKLMNPQRSTWY_-',
                    '84V': 'PR84V',
                    '85V': 'PR85V',
                    '88ACEFHIKLMPQRVWYid': 'PR88ACEFHIKLMPQRVWY_-',
                    '88D': 'PR88D',
                    '88S': 'PR88S',
                    '88TG': 'PR88TG',
                    '89VT': 'PR89VT',
                    '90ACDEFGHIKNPQRSTVWYid': 'PR90ACDEFGHIKNPQRSTVWY_-',
                    '90M': 'PR90M'},
             'RT': {'100ACDEFGHKMNPQRSTWYid': 'RT100ACDEFGHKMNPQRSTWY_-',
                    '100I': 'RT100I',
                    '100V': 'RT100V',
                    '101CDFGILMSVWYid': 'RT101CDFGILMSVWY_-',
                    '101E': 'RT101E',
                    '101H': 'RT101H',
                    '101NAT': 'RT101NAT',
                    '101P': 'RT101P',
                    '101Q': 'RT101Q',
                    '103ACDFGILMPVWYid': 'RT103ACDFGILMPVWY_-',
                    '103EQ': 'RT103EQ',
                    '103H': 'RT103H',
                    '103N': 'RT103N',
                    '103R': 'RT103R',
                    '103S': 'RT103S',
                    '103T': 'RT103T',
                    '106A': 'RT106A',
                    '106CDEFGHKLNPQRSTWYid': 'RT106CDEFGHKLNPQRSTWY_-',
                    '106I': 'RT106I',
                    '106M': 'RT106M',
                    '108ACDEFGHKLMNPQRSTWYid': 'RT108ACDEFGHKLMNPQRSTWY_-',
                    '108I': 'RT108I',
                    '115ACDEGHIKLMNPQRSTVWid': 'RT115ACDEGHIKLMNPQRSTVW_-',
                    '115F': 'RT115F',
                    '116ACDEGHIKLMNPQRSTVWid': 'RT116ACDEGHIKLMNPQRSTVW_-',
                    '116Y': 'RT116Y',
                    '118I': 'RT118I',
                    '132ML': 'RT132ML',
                    '138A': 'RT138A',
                    '138CDFHILMNPSTVWYid': 'RT138CDFHILMNPSTVWY_-',
                    '138K': 'RT138K',
                    '138QG': 'RT138QG',
                    '138R': 'RT138R',
                    '151ACDEFGHIKNPRSTVWYid': 'RT151ACDEFGHIKNPRSTVWY_-',
                    '151L': 'RT151L',
                    '151M': 'RT151M',
                    '179ACGHKMNPQRSWYid': 'RT179ACGHKMNPQRSWY_-',
                    '179DE': 'RT179DE',
                    '179F': 'RT179F',
                    '179I': 'RT179I',
                    '179L': 'RT179L',
                    '179T': 'RT179T',
                    '181ADEHKLMNPQRTWid': 'RT181ADEHKLMNPQRTW_-',
                    '181C': 'RT181C',
                    '181FSG': 'RT181FSG',
                    '181IV': 'RT181IV',
                    '184ACDEFGHKLNPQRSTWYid': 'RT184ACDEFGHKLNPQRSTWY_-',
                    '184VI': 'RT184VI',
                    '188ADEGIKMNPQRSTVWid': 'RT188ADEGIKMNPQRSTVW_-',
                    '188C': 'RT188C',
                    '188F': 'RT188F',
                    '188H': 'RT188H',
                    '188L': 'RT188L',
                    '190A': 'RT190A',
                    '190CTV': 'RT190CTV',
                    '190DFHIKLMNPWYid': 'RT190DFHIKLMNPWY_-',
                    '190EQ': 'RT190EQ',
                    '190R': 'RT190R',
                    '190S': 'RT190S',
                    '210W': 'RT210W',
                    '215F': 'RT215F',
                    '215GHKMPQRWid': 'RT215GHKMPQRW_-',
                    '215SCDEIVALN': 'RT215SCDEIVALN',
                    '215Y': 'RT215Y',
                    '219ACDFGHILMPSTVYid': 'RT219ACDFGHILMPSTVY_-',
                    '219NR': 'RT219NR',
                    '219QE': 'RT219QE',
                    '219W': 'RT219W',
                    '221Y': 'RT221Y',
                    '225ACDEFGIKLMNQRSTVWYid': 'RT225ACDEFGIKLMNQRSTVWY_-',
                    '225H': 'RT225H',
                    '227ADEGHKMNPQRSTWYid': 'RT227ADEGHKMNPQRSTWY_-',
                    '227C': 'RT227C',
                    '227ILV': 'RT227ILV',
                    '230ACDEFGHKNPQRSTVWYid': 'RT230ACDEFGHKNPQRSTVWY_-',
                    '230I': 'RT230I',
                    '230L': 'RT230L',
                    '234ACDEFGHKLMNPQRSTWYid': 'RT234ACDEFGHKLMNPQRSTWY_-',
                    '234I': 'RT234I',
                    '236ACDEFGHIKMNQRSTVWYid': 'RT236ACDEFGHIKMNQRSTVWY_-',
                    '236L': 'RT236L',
                    '238ACDEFGHILMPQSVWYid': 'RT238ACDEFGHILMPQSVWY_-',
                    '238R': 'RT238R',
                    '238TN': 'RT238TN',
                    '318ACDEGHIKLMNPQRSTVWid': 'RT318ACDEGHIKLMNPQRSTVW_-',
                    '318F': 'RT318F',
                    '348ACDEFGHKLMPQRSTVWYid': 'RT348ACDEFGHKLMPQRSTVWY_-',
                    '348I': 'RT348I',
                    '40F': 'RT40F',
                    '41I': 'RT41I',
                    '41L': 'RT41L',
                    '44AD': 'RT44AD',
                    '62V': 'RT62V',
                    '65ACDFGHILMPQSTVWYid': 'RT65ACDFGHILMPQSTVWY_-',
                    '65E': 'RT65E',
                    '65N': 'RT65N',
                    '65R': 'RT65R',
                    '67ACFIKLMPQRVWYi': 'RT67ACFIKLMPQRVWY_',
                    '67GESTH': 'RT67GESTH',
                    '67N': 'RT67N',
                    '67d': 'RT67-',
                    '68G': 'RT68G',
                    '68d': 'RT68-',
                    '69CFHKLMPQRVWY': 'RT69CFHKLMPQRVWY',
                    '69D': 'RT69D',
                    '69G': 'RT69G',
                    '69N': 'RT69N',
                    '69SAIE': 'RT69SAIE',
                    '69d': 'RT69-',
                    '69i': 'RT69_',
                    '70ACDFHILMPVWYi': 'RT70ACDFHILMPVWY_',
                    '70EG': 'RT70EG',
                    '70QNST': 'RT70QNST',
                    '70R': 'RT70R',
                    '70d': 'RT70-',
                    '74ACDEFGHKMNPQRSTWYid': 'RT74ACDEFGHKMNPQRSTWY_-',
                    '74VI': 'RT74VI',
                    '75CDEFGHKNPQRWYid': 'RT75CDEFGHKNPQRWY_-',
                    '75I': 'RT75I',
                    '75M': 'RT75M',
                    '75SA': 'RT75SA',
                    '75T': 'RT75T',
                    '77ACDEGHIKMNPQRSTVWYid': 'RT77ACDEGHIKMNPQRSTVWY_-',
                    '77L': 'RT77L',
                    '90I': 'RT90I',
                    '98G': 'RT98G'}}

        res_comments = self.algorithm.parse_comments(self.root)

        self.assertEqual(exp_comments, res_comments)


if __name__ == '__main__':
    unittest.main()