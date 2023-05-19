import tempfile
import os

def get_input_sequences(handle, return_dict=False):
    """
    Parse open file as FASTA, return a list of sequences.

    :param handle:  Open stream to FASTA file in read mode
    :param return_dict:  Option to return a dictionary of sequences keyed by header
    :return:  if return_dict is True, returns dictionary; else a list of sequences
    """
    headers = []
    sequences = []

    header = ''
    sequence = ''
    for line in handle:
        if line.startswith('$'):
            # skip comment
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                headers.append(header)
                sequences.append(sequence)
                sequence = ''  # reset
            header = line.strip('>#\n')
        else:
            sequence += line.strip('\n').upper()

    # handle last entry
    headers.append(header)
    sequences.append(sequence)

    return dict(zip(headers, sequences)) if return_dict else sequences
#
# def generateTable():
#     """
#     Generates a dictionary of codon to amino acid mappings, including ambiguous combinations.
#     @return tripletTable: codon to amino acid dictionary
#     """
#     codonToAminoAcidMap = {
#         "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
#         "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
#         "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
#         "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
#         "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
#         "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
#         "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
#         "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
#         "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
#         "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
#         "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
#         "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
#         "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
#         "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
#         "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
#         "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
#     }
#     nas = ["A", "C", "G", "T", "R", "Y", "M", "W", "S", "K", "B", "D", "H", "V", "N"]
#     tripletTable = dict()
#     for i in range(len(nas)):
#         for j in range(len(nas)):
#             for k in range(len(nas)):
#                 triplet = nas[i] + nas[j] + nas[k]
#                 codons = enumerateCodonPossibilities(triplet)
#                 uniqueAAs = []
#                 for codon in codons:
#                     c = codonToAminoAcidMap[codon]
#                     if c not in uniqueAAs:
#                         uniqueAAs.append(c)
#                 if len(uniqueAAs) > 4:
#                     aas = "X"
#                 else:
#                     aas = ''.join(uniqueAAs)
#                 tripletTable[triplet] = aas
#     return tripletTable
#
# def enumerateCodonPossibilities(triplet):
#     """
#     Converts a potentially ambiguous nucleotide triplet into standard ATCG codons.
#     @param triplet: nucleotide triplet as a string
#     @return codonPossibilities: list of possible ATCG codons encoded by the triplet
#     """
#     ambiguityMap = {
#         "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
#         "R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"],
#         "W": ["A", "T"], "S": ["C", "G"], "K": ["G", "T"],
#         "B": ["C", "G", "T"], "D": ["A", "G", "T"],
#         "H": ["A", "C", "T"], "V": ["A", "C", "G"],
#         "N": ["A", "C", "G", "T"]
#     }
#     codonPossibilities = []
#     pos1, pos2, pos3 = triplet
#     for p1 in ambiguityMap[pos1]:
#         for p2 in ambiguityMap[pos2]:
#             for p3 in ambiguityMap[pos3]:
#                 codonPossibilities.append(p1 + p2 + p3)
#     return codonPossibilities