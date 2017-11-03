from gotoh2.aligner import Aligner
import unittest

g2 = Aligner()

def align(reference, query):
    """
    'Procrustean' alignment of query nucleotide sequence to reference sequence

    :param query: Nucleotide sequence as a string
    :param reference: Reference sequence as a string
    :return: Aligned sequence fitted wrt relative insertions to reference
    """

    nucleotides = 'ACGT'

    #Use gotoh2 to align the query to reference sequence
    query = query.strip()
    assert type(query) is str, 'Query must be a string.'
    aligned_query, aligned_reference, aligned_score = g2.align(query, reference)
    fitted_query = ''

    ''' TODO:   flag indels not in sets of three
                exceptions for indels in proximity summing to 3n
    '''

    #Remove insertions in the sequence relative to reference sequence
    for idx,char in enumerate(aligned_reference):
        #Iff the query does not have insertion relative to reference, add
        #the nucleotide/relative deletion to the output string
        if char != '-':
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
