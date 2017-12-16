import unittest
import sys
sys.path.append("/home/jasper/git/sierra-local")
import sierralocal
from nucaminohook import NucAminoAligner
import os

class DummyTest(unittest.TestCase):

    def setUp(self):
        self.cwd = os.getcwd()

    def testSequence(self):
        self.aligner = NucAminoAligner('/home/jasper/git/sierra-local/AY030621.fasta')
        self.aligner.align_file()
        mutations = self.aligner.get_mutations()
        HIVdb_results = {'ATV/r':45,'DRV/r':0,'LPV/r':20}
        for res in HIVdb_results:
            self.assertEqual(mutations[res], HIVdb_results[res])


if __name__ == '__main__':
    unittest.main()