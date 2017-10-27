def translate(sequence):
    codon_table = {
        'TTT' : 'F', 'CTT' : 'L', 'ATT' : 'I', 'GTT' : 'V',
        'TTC' : 'F', 'CTC' : 'L', 'ATC' : 'I', 'GTC' : 'V',
        'TTA' : 'L', 'CTA' : 'L', 'ATA' : 'I', 'GTA' : 'V',
        'TTG' : 'L', 'CTG' : 'L', 'ATG' : 'M', 'GTG' : 'V',
        'TCT' : 'S', 'CCT' : 'P', 'ACT' : 'T', 'GCT' : 'A',
        'TCC' : 'S', 'CCC' : 'P', 'ACC' : 'T', 'GCC' : 'A',
        'TCA' : 'S', 'CCA' : 'P', 'ACA' : 'T', 'GCA' : 'A',
        'TCG' : 'S', 'CCG' : 'P', 'ACG' : 'T', 'GCG' : 'A',
        'TAT' : 'Y', 'CAT' : 'H', 'AAT' : 'N', 'GAT' : 'D',
        'TAC' : 'Y', 'CAC' : 'H', 'AAC' : 'N', 'GAC' : 'D',
        'TAA' : '*', 'CAA' : 'Q', 'AAA' : 'K', 'GAA' : 'E',
        'TAG' : '*', 'CAG' : 'Q', 'AAG' : 'K', 'GAG' : 'E',
        'TGT' : 'C', 'CGT' : 'R', 'AGT' : 'S', 'GGT' : 'G',
        'TGC' : 'C', 'CGC' : 'R', 'AGC' : 'S', 'GGC' : 'G',
        'TGA' : '*', 'CGA' : 'R', 'AGA' : 'R', 'GGA' : 'G',
        'TGG' : 'W', 'CGG' : 'R', 'AGG' : 'R', 'GGG' : 'G'
    }   

    peptide_sequence = ''

    if len(sequence) >= 3:
        for codon in [sequence[i:i+3] for i in range(0,len(sequence),3)]:
            if not codon_table[codon] == '*':
                peptide_sequence += codon_table[codon]
            else:
                break

    return peptide_sequence

def main():
    #HXB2 Integrase
    peptide = translate("TTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAG")
    print peptide
    print len(peptide), "amino acids"

if __name__ == '__main__':
    main()