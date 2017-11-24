import unittest
import sys
sys.path.append("/home/jasper/git/sierra-local")
import sierralocal
import reference
import hxb2

class DummyTest(unittest.TestCase):
    def setUp(self):
        # add unit test fixtures
        datum = 123.4
    def testSimple(self):
        expected = 'foobar'
        # do something that should yield expected output
        result = 'foobar'
        self.assertEqual(expected, result)

    def test(self):
        self.assertEqual(sierralocal.align('TAGACTTACCAC', 'ACTTAGCAT'), '-A--CTTAGCAT')
        self.assertEqual(sierralocal.align('ACTTAGCAT', 'TAGACTTACCAC'), 'ACTTACCAC')

    # Check out this error later
    # def test_empty(self):
    #    self.assertEqual(align('-----------','ACGTACGTACGT'),'')

    def test_deletions(self):
        self.assertEqual(
            sierralocal.align('AGTACGCTCGTAGCAT', 'AGTACTAGCAT'), 'AGTAC-----TAGCAT')

    def test_insertions(self):
        self.assertEqual(
            sierralocal.align('AGTACTAGCAT', 'ACTACGCTCGTAGCAT'), 'ACTACTAGCAT')


    def testTranslate(self):
        pro = 'PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF'
        self.assertEqual(sierralocal.translate(hxb2.pro),pro)

if __name__ == '__main__':
    unittest.main()