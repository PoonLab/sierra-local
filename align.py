from gotoh2.aligner import Aligner
import unittest

def align(reference, query):
    """
    'Procrustean' alignment of query nucleotide sequence to reference sequence

    :param query: Nucleotide sequence as a string
    :param reference: Reference sequence as a string
    :return: Aligned sequence fitted wrt relative insertions to reference
    """

    nucleotides = 'ACGT'

    #Use gotoh2 to align the query to reference sequence
    assert type(query) is str, 'Query must be a string.'
    g2 = Aligner()
    g2_results = g2.align(query, reference)

    aligned_reference = g2_results[1]
    aligned_query = g2_results[0]
    fitted_query = ''

    #Remove insertions in the sequence relative to reference sequence
    for idx,char in enumerate(aligned_reference):
        #Iff the query does not have insertion relative to reference, add
        #the nucleotide/relative deletion to the output string
        if char in nucleotides:
            fitted_query += aligned_query[idx]

    return fitted_query

class TestAlignments(unittest.TestCase):
    def test(self):
        self.assertEqual(align('TAGACTTACCAC', 'ACTTAGCAT'),'-A--CTTAGCAT')
        self.assertEqual(align('ACTTAGCAT', 'TAGACTTACCAC'),'ACTTACCAC')

    def test_empty(self):
        self.assertEqual(align('-----------','ACGTACGTACGT'),'')

    def test_deletions(self):
        self.assertEqual(align('AGTACGCTCGTAGCAT', 'AGTACTAGCAT'),'AGTAC-----TAGCAT')

    def test_insertions(self):
        self.assertEqual(align('AGTACTAGCAT', 'ACTACGCTCGTAGCAT'),'ACTACTAGCAT')


if __name__ == '__main__':
    unittest.main()