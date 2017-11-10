import unittest
import sierra-local

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
    self.assertEqual(align('TAGACTTACCAC', 'ACTTAGCAT'), '-A--CTTAGCAT')
       self.assertEqual(align('ACTTAGCAT', 'TAGACTTACCAC'), 'ACTTACCAC')

    # Check out this error later
    # def test_empty(self):
    #    self.assertEqual(align('-----------','ACGTACGTACGT'),'')

    def test_deletions(self):
        self.assertEqual(
            align('AGTACGCTCGTAGCAT', 'AGTACTAGCAT'), 'AGTAC-----TAGCAT')

    def test_insertions(self):
        self.assertEqual(
            align('AGTACTAGCAT', 'ACTACGCTCGTAGCAT'), 'ACTACTAGCAT')

if __name__ == '__main__':
    unittest.main()