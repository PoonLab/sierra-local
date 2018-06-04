import subprocess
import os
from pathlib import Path
import csv
import re

class NucAminoAligner():
    def __init__(self):
        '''
        Initialize NucAmino for a specific input fasta file
        '''
        self.cwd = os.path.curdir
        items = os.listdir(self.cwd)
        self.nucamino_binary = 'nucamino'
        for name in items:
            if 'nucamino' in name and not (name.endswith('py') or name.endswith('pyc')):
                self.nucamino_binary = name
                break
        print("Found NucAmino binary", self.nucamino_binary)

    def align_file(self, filename): 
        '''
        Using subprocess to call NucAmino, generates an output .tsv containing mutation data for each sequence in the FASTA file
        '''
        self.inputname = filename
        self.outputname = os.path.splitext(filename)[0] + '.tsv'

        args = [
            "./{}".format(self.nucamino_binary),
            "hiv1b",
            "-q",
            "-i", "{}".format(self.inputname),
            "-g=POL",
            "-o", "{}".format(self.outputname)
        ]
        popen = subprocess.Popen(args)
        popen.wait()

    def gene_map(self):
        pol_start = 2085
        pol_nuc_map = {
            'PR':(2253,2549),
            'RT':(2550,3869),
            'IN':(4230,5096)
        }
        convert = lambda x:(x-pol_start)/3
        pol_aa_map = {}
        for key, val in pol_nuc_map.items():
            pol_aa_map[key] = (convert(val[0]), convert(val[1]))
        return pol_aa_map

    def get_genes(self, firstAA):
        """
        Determines the first POL gene that is present in the query sequence, by virtue of gene breakpoints
        @param firstAA: the first amino acid position of the query sequence, as determined by NucAmino
        @return: list of length 1 with a string of the gene in the sequence
        """
        start = 2085
        gene_list = ['PR', 'RT', 'IN']
        firstNA = [int((x - start)/3) for x in [2253, 2550, 4230, 5096]]
        gene_present = [firstNA[i+1] > firstAA > x for i, x in enumerate(firstNA[:-1:])]
        if any(gene_present):
            return [gene_list[gene_present.index(True)]] #first TRUE only
        return []

    def get_mutations(self, genemap, fastaFileName):
        '''
        From the tsv output of NucAmino, parses and adjusts indices and returns as lists.
        :return: list of sequence names, list of sequence mutation dictionaries.
        '''
        file_mutations = []
        file_genes = []
        file_firstlastNA = []
        sequence_headers = []
        tripletTable = self.generateTable()
        
        #grab sequences from file
        with open(fastaFileName, 'r') as handle:
            # names = []
            sequence_list = []
            sequence = ''
            for i in handle:
                if i[0] == '$':  # skip h info
                    continue
                elif i[0] == '>' or i[0] == '#':
                    if len(sequence) > 0:
                        sequence_list.append(sequence)
                        sequence = ''
                else:
                    sequence += i.strip('\n').upper()
            sequence_list.append(sequence)
        
        with open(os.path.splitext(fastaFileName)[0] + '.tsv','r') as nucamino_alignment: # open the NucAmino output file
            tsvin = csv.reader(nucamino_alignment, delimiter='\t')
            next(tsvin) #bypass the header row
            for idx, row in enumerate(tsvin): #iterate over sequences (1 per row)
                sequence_headers.append(row[0])
                # Let's use the gene map to figure the "reference position" for each gene
                # Nucamino outputs wonky positions, so calculate the shift so we can get positions relative to the start of each individual gene
                shift = 0
                firstAA = int(row[1])
                firstNA = int(row[3])
                lastNA = int(row[4])
                pol_map = self.gene_map()
                for key in pol_map:
                    low, hi = pol_map[key]
                    if low <= firstAA <= hi:
                        shift = low
                        break

                # Edge case of NO mutations
                if len(row[5]) == 0:
                    position_list = []
                    consensus_list = []
                    codon_list = []
                    aminoacid_list = []
                    gene_muts = {}
                # most cases fall here... mutations are present
                else:
                    #generate data lists for each sequence
                    mutation_list = row[5].split(',') #split list into individual mutations
                    position_list = [int(int(re.findall(r'\d+',x.split(':')[0])[0])-shift) for x in mutation_list] #residue position list
                    consensus_list = [re.findall(r'\D+',x.split(':')[0][:2])[0] for x in mutation_list] #consensus list
                    codon_list = [sequence_list[idx][(((p+int(shift)-int(row[1])+1)-1)*3):((p+int(shift)-int(row[1])+1)*3)] for p in position_list] #obtain the mutation codon directly from the query sequence
                    aminoacid_list = ['-' if x == '~~~' or x == '...' else self.translateNATriplet(x, tripletTable) for x in codon_list]
                    gene_muts = dict(zip(position_list, zip(consensus_list,aminoacid_list)))
                # append everything to list of lists of the entire sequence set
                file_mutations.append(gene_muts)
                file_genes.append(self.get_genes(firstAA))
                file_firstlastNA.append((firstNA, lastNA))
        assert len(file_mutations) == len(sequence_headers), "length of mutations dicts is not the same as length of names"
        return sequence_headers, file_genes, file_mutations, file_firstlastNA

# BELOW is an implementation of sierra's Java algorithm for determining codon ambiguity

    def translateNATriplet(self, triplet, tripletTable):
        if len(triplet) != 3:
            return "X"
        if '~' in triplet:
            return "X"
        if triplet in tripletTable:
            return tripletTable[triplet]
        return "X"

    def generateTable(self):
        codonToAminoAcidMap = {
            "TTT" : "F", "TTC" : "F", "TTA" : "L", "TTG" : "L", "CTT" : "L", "CTC" : "L", "CTA" : "L", "CTG" : "L", "ATT" : "I", "ATC" : "I", "ATA" : "I", "ATG" : "M", "GTT" : "V", "GTC" : "V", "GTA" : "V", "GTG" : "V", "TCT" : "S", "TCC" : "S", "TCA" : "S", "TCG" : "S", "CCT" : "P", "CCC" : "P", "CCA" : "P", "CCG" : "P", "ACT" : "T", "ACC" : "T", "ACA" : "T", "ACG" : "T", "GCT" : "A", "GCC" : "A", "GCA" : "A", "GCG" : "A", "TAT" : "Y", "TAC" : "Y", "TAA" : "*", "TAG" : "*", "CAT" : "H", "CAC" : "H", "CAA" : "Q", "CAG" : "Q", "AAT" : "N", "AAC" : "N", "AAA" : "K", "AAG" : "K", "GAT" : "D", "GAC" : "D", "GAA" : "E", "GAG" : "E", "TGT" : "C", "TGC" : "C", "TGA" : "*", "TGG" : "W", "CGT" : "R", "CGC" : "R", "CGA" : "R", "CGG" : "R", "AGT" : "S", "AGC" : "S", "AGA" : "R", "AGG" : "R", "GGT" : "G", "GGC" : "G", "GGA" : "G", "GGG" : "G"
        }
        nas = ["A","C","G","T","R","Y","M","W","S","K","B","D","H","V","N"]
        tripletTable = dict()
        for i in range(len(nas)):
            for j in range(len(nas)):
                for k in range(len(nas)):
                    triplet = nas[i] + nas[j] + nas[k]
                    codons = self.enumerateCodonPossibilities(triplet)
                    uniqueAAs = []
                    for codon in codons:
                        c = codonToAminoAcidMap[codon]
                        if c not in uniqueAAs:
                            uniqueAAs.append(c)
                    if len(uniqueAAs) > 4:
                        aas = "X"
                    else:
                        aas = ''.join(uniqueAAs)
                    tripletTable[triplet] = aas
        return tripletTable

    def enumerateCodonPossibilities(self, triplet):
        ambiguityMap = {
            "A" : ["A"],
            "C" : ["C"],
            "G" : ["G"],
            "T" : ["T"],
            "R" : ["A","G"],
            "Y" : ["C","T"],
            "M" : ["A","C"],
            "W" : ["A","T"],
            "S" : ["C","G"],
            "K" : ["G","T"],
            "B" : ["C","G","T"],
            "D" : ["A","G","T"],
            "H" : ["A","C","T"],
            "V" : ["A","C","G"],
            "N" : ["A","C","G","T"]
        }
        codonPossibilities = []
        pos1 = triplet[0]
        pos2 = triplet[1]
        pos3 = triplet[2]
        for p1 in ambiguityMap[pos1]:
            for p2 in ambiguityMap[pos2]:
                for p3 in ambiguityMap[pos3]:
                    codonPossibilities.append(p1+p2+p3)
        return codonPossibilities

# trim low quality nucleotides

    def trimLowQualities(sequence, firstAA, lastAA, mutations, frameshifts):
        SEQUENCE_SHRINKAGE_CUTOFF_PCNT = 30
        SEQUENCE_SHRINKAGE_WINDOW = 15
        SEQUENCE_SHRINKAGE_BAD_QUALITY_MUT_PREVALENCE = 0.1

        badPcnt = 0
        trimLeft = 0
        trimRight = 0
        problemSites = 0
        sinceLastBadQuality = 0
        proteinSize = lastAA - firstAA + 1
        candidates = []
        invalidSites = [False for i in range(proteinSize)]

        for mut in mutations:
            idx = mut.position - firstAA
            if not mut.isSequenced() and (mut.getHighestMutPrevalence() < SEQUENCE_SHRINKAGE_BAD_QUALITY_MUT_PREVALANCE or mut.AAs == 'X' or mut.isApobecMutation() or mut.hasStopCodon()):
                invalidSites[idx] = True

        for fs in frameshifts:
            idx = fs.getPosition() - firstAA
            invalidSites[idx] = True

        # forward scan for trimming left
        for idx in range(0, proteinSize):
            if sinceLastBadQuality > SEQUENCE_SHRINKAGE_WINDOW:
                break
            elif invalidSites[idx]:
                problemSites += 1
                trimLeft = idx + 1
                badPcnt = problemSites * 100 / trimLeft if trimLeft > 0 else 0
                if badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_POINT:
                    candidates.append(trimLeft)
                sinceLastBadQuality = 0
            else:
                sinceLastBadQuality += 1
        trimLeft = candidates[-1] if len(candidates) > 0 else 0
        candidates = []

        #backward scan for trimming right
        problemSites = 0
        sinceLastBadQuality = 0
        for idx in range(proteinSize-1, -1, -1):
            if sinceLastBadQuality > SEQUENCE_SHRINKAGE_WINDOW:
                break
            elif invalidSites[idx]:
                problemSites += 1
                trimRight = proteinSize - idx
                badPcnt = problemSites * 100 / trimRight if trimRight > 0 else 0
                if badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_POINT:
                    candidates.append(trimRight)
                sinceLastBadQuality = 0
            else:
                sinceLastBadQuality += 1
        trimRight = candidates[-1] if len(candidates) > 0 else 0
        return (trimLeft, trimRight)



if __name__ == '__main__':
    test = NucAminoAligner()
    print(test.translateNATriplet("YTD", test.generateTable()))