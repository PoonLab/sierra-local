import subprocess
import os
from pathlib import Path
import csv
import re

class NucAminoAligner():
    '''
    Initialize NucAmino for a specific input fasta file
    '''
    def __init__(self):
        self.cwd = os.path.curdir
        items = os.listdir(Path(os.path.dirname(__file__)))
        self.nucamino_binary = 'nucamino'
        for name in items:
            if 'nucamino' in name and not (name.endswith('py') or name.endswith('pyc')):
                self.nucamino_binary = name
                break
        print("Found NucAmino binary", self.nucamino_binary)
        self.tripletTable = self.generateTable()
        with open(Path(os.path.dirname(__file__))/'data'/'apobec.tsv','r') as csvfile:
            self.ApobecDRMs = list(csv.reader(csvfile, delimiter='\t'))

        with open(Path(os.path.dirname(__file__))/'data'/'INIPrevalences.tsv','r') as INIfile:
            INI_prevalences = list(csv.reader(INIfile, delimiter='\t'))
        self.INI_dict = {}
        col_idx = ["Naive:%" in cell for cell in INI_prevalences[0]]
        for row in INI_prevalences[1:]:
            for i, col in enumerate(col_idx):
                if col:
                    self.INI_dict[str(row[0]+row[1]+row[2]+INI_prevalences[0][i].split(':')[0])] = float(row[i])

        with open(Path(os.path.dirname(__file__))/'data'/'PIPrevalences.tsv','r') as PIfile:
            PI_prevalences = list(csv.reader(PIfile, delimiter='\t'))
        self.PI_dict = {}
        col_idx = ["Naive:%" in cell for cell in PI_prevalences[0]]
        for row in PI_prevalences[1:]:
            for i, col in enumerate(col_idx):
                if col:
                    self.PI_dict[str(row[0]+row[1]+row[2]+PI_prevalences[0][i].split(':')[0])] = float(row[i])

        with open(Path(os.path.dirname(__file__))/'data'/'RTIPrevalences.tsv','r') as RTIfile:
            RTI_prevalences = list(csv.reader(RTIfile, delimiter='\t'))
        self.RTI_dict = {}
        col_idx = ["Naive:%" in cell for cell in RTI_prevalences[0]]
        for row in RTI_prevalences[1:]:
            for i, col in enumerate(col_idx):
                if col:
                    self.RTI_dict[str(row[0]+row[1]+row[2]+RTI_prevalences[0][i].split(':')[0])] = float(row[i])


    '''
    Using subprocess to call NucAmino, generates an output .tsv containing mutation data for each sequence in the FASTA file
    '''
    def align_file(self, filename): 
        self.inputname = filename
        self.outputname = os.path.splitext(filename)[0] + '.tsv'

        args = [
            Path(os.path.dirname(__file__))/self.nucamino_binary,
            "hiv1b",
            "-q",
            "-i", "{}".format(self.inputname),
            "-g=POL",
            "-o", "{}".format(self.outputname)
        ]
        popen = subprocess.Popen(args)
        popen.wait()

    """
    Returns a dictionary with the amino acid position bounds for each gene in Pol, based on the HXB2 reference annotations.
    """
    def gene_map(self):
        pol_start = 2085
        pol_nuc_map = {
            'PR':(2253,2549),
            'RT':(2550,3869),
            'IN':(4230,5096)
        }
        convert = lambda x:int((x-pol_start)/3)
        pol_aa_map = {}
        for key, val in pol_nuc_map.items():
            pol_aa_map[key] = (convert(val[0]), convert(val[1]))
        return pol_aa_map

    """
    Determines the first POL gene that is present in the query sequence, by virtue of gene breakpoints
    @param firstAA: the first amino acid position of the query sequence, as determined by NucAmino
    @return: list of length 1 with a string of the gene in the sequence
    """
    def get_genes(self, firstAA):
        start = 2085
        gene_list = ['PR', 'RT', 'IN']
        firstNA = [int((x - start)/3) for x in [2253, 2550, 4230, 5096]]
        gene_present = [firstNA[i+1] > firstAA > x for i, x in enumerate(firstNA[:-1:])]
        if any(gene_present):
            return [gene_list[gene_present.index(True)]] #first TRUE only
        return []

    '''
    From the tsv output of NucAmino, parses and adjusts indices and returns as lists.
    :return: list of sequence names, list of sequence mutation dictionaries.
    '''
    def get_mutations(self, genemap, fastaFileName):
        file_mutations = []
        file_genes = []
        file_firstlastNA = []
        file_trims = []
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
                lastAA = int(row[2])
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
                    aminoacid_list = ['-' if x == '~~~' or x == '...' else self.translateNATriplet(x) for x in codon_list]
                    gene_muts = dict(zip(position_list, zip(consensus_list,aminoacid_list)))

                # trim low quality leading and trailing nucleotides
                trimLeft, trimRight = self.trimLowQualities(codon_list, shift, firstAA, lastAA, gene_muts, [], self.get_genes(firstAA)[0], row[0].split('.')[-1]) #TODO: get subtype via alignment
                if trimLeft + trimRight > 0:
                    print(row[0], trimLeft, trimRight)
                trimmed_gene_muts = {k:v for (k,v) in gene_muts.items() if (k >= firstAA - shift + trimLeft) and (k <= lastAA - shift - trimRight)}

                # append everything to list of lists of the entire sequence set
                file_mutations.append(trimmed_gene_muts)
                file_genes.append(self.get_genes(firstAA))
                file_firstlastNA.append((firstNA, lastNA))
                file_trims.append((trimLeft, trimRight))
        assert len(file_mutations) == len(sequence_headers), "length of mutations dicts is not the same as length of names"
        return sequence_headers, file_genes, file_mutations, file_firstlastNA, file_trims

# BELOW is an implementation of sierra's Java algorithm for determining codon ambiguity
    
    """
    Translates a nucleotide triplet into its amino acid mixture or ambiguity.
    @param triplet: nucleotide sequence as a string
    @return: translation of the triplet as a string
    """
    def translateNATriplet(self, triplet):
        if len(triplet) != 3:
            return "X"
        if '~' in triplet:
            return "X"
        if triplet in self.tripletTable:
            return self.tripletTable[triplet]
        return "X"

    """
    Generates a dictionary of codon to amino acid mappings, including ambiguous combinations.
    @return tripletTable: codon to amino acid dictionary
    """
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

    """
    Converts a potentially ambiguous nucleotide triplet into standard ATCG codons.
    @param triplet: nucleotide triplet as a string
    @return codonPossibilities: list of possible ATCG codons encoded by the triplet
    """
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
    
    """
    Filters low-quality leading and trailing nucleotides from a query.
    Removes large (length > SEQUENCE_TRIM_SITES_CUTOFF) low quality pieces.
    Low quality is defined as:
    (1) unusual mutation; or
    (2) 'X' in amino acid list; or
    (3) has a stop codon

    @param sequence: sequence
    @param firstAA: aligned position of first amino acid in query
    @param lastAA: aligned position of last amino acid in query
    @param mutations: list of mutations
    @param frameshifts: list of frameshifts
    @return: tuple of how many leading and trailing nucleotides to trim
    """
    def trimLowQualities(self, codon_list, shift, firstAA, lastAA, mutations, frameshifts, gene, subtype):
        SEQUENCE_SHRINKAGE_CUTOFF_PCNT = 30
        SEQUENCE_SHRINKAGE_WINDOW = 15
        SEQUENCE_SHRINKAGE_BAD_QUALITY_MUT_PREVALENCE = 0.1

        # print(mutations)
        # print(codon_list)

        badPcnt = 0
        trimLeft = 0
        trimRight = 0
        problemSites = 0
        sinceLastBadQuality = 0
        proteinSize = lastAA - firstAA + 1
        candidates = []
        invalidSites = [False for i in range(proteinSize)]

        # account for invalid sites
        for j, position in enumerate(mutations):
            idx = position - firstAA + shift
            # print(idx)
            if not self.isUnsequenced(codon_list[j]) and (self.getHighestMutPrevalence((position, mutations[position]), gene, subtype) < SEQUENCE_SHRINKAGE_BAD_QUALITY_MUT_PREVALENCE or mutations[position][1] == 'X' or self.isApobecDRM(gene, mutations[position][0], position, mutations[position][1]) or self.isStopCodon(codon_list[j])):
                invalidSites[idx] = True

        # for fs in frameshifts:
        #     idx = fs.getPosition() - firstAA
        #     invalidSites[idx] = True

        # forward scan for trimming left
        for idx in range(0, proteinSize):
            if sinceLastBadQuality > SEQUENCE_SHRINKAGE_WINDOW:
                break
            elif invalidSites[idx]:
                problemSites += 1
                trimLeft = idx + 1
                badPcnt = problemSites * 100 / trimLeft if trimLeft > 0 else 0
                if badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT:
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
                if badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT:
                    candidates.append(trimRight)
                sinceLastBadQuality = 0
            else:
                sinceLastBadQuality += 1
        trimRight = candidates[-1] if len(candidates) > 0 else 0
        return (trimLeft, trimRight)

    """
    Determines whether a triplet is unsequenced (has more than one N or deletion).
    "NNN", "NN-", "NNG" should be considered as unsequenced region.
    """
    def isUnsequenced(self, triplet):
        return (triplet.replace("-", "N").count("N") > 1) #TODO: incorporate !isInsertion &&

    def isStopCodon(self, triplet):
        return ("*" in self.translateNATriplet(triplet))

    def isApobecDRM(self, gene, consensus, position, AA):
        ls = [row[0:3] for row in self.ApobecDRMs[1:]]
        if [gene, consensus, str(position)] in ls:
            i = ls.index([gene, consensus, str(position)])
            for aa in AA:
                if aa in self.ApobecDRMs[1:][i][3]:
                    return True
        return False

    def getHighestMutPrevalence(self, mutation, gene, subtype):
        """
        #TODO
        @param mutation: a tuple(?) representing a specific position in the amino acid sequence 
                         that may contain multiple amino acids (polymorphic)
        @param gene: PR, RT, or INT
        @param subtype: predicted from function()
        @return: prevalence of the most common amino acid encoded at this position within the 
                 subtype alignment
        """
        position, aaseq = mutation
        cons, aas = aaseq
        aas = aas.replace(cons, '')  # ignore consensus
        aas = aas.replace('*', '')  # remove stop codons
        
        prevalence = 0.
        for aa in aas:
            aaPrevalence = self.getMutPrevalence(position, cons, aa, gene, subtype)
            prevalence = max(prevalence, aaPrevalence)

        return prevalence

    def getMutPrevalence(self, position, cons, aa, gene, subtype):
        # key = str(position)+str(cons)+str(aa)+str(subtype)

        # if gene == 'IN' and key in self.INI_dict:
        #     return self.INI_dict[key]

        # if gene == 'PR' and key in self.PI_dict:
        #     return self.PI_dict[key]

        # if gene == 'RT' and key in self.RTI_dict:
        #     return self.RTI_dict[key]

        key2 = str(position)+str(cons)+str(aa)+"All"

        if gene == 'IN' and key2 in self.INI_dict:
            return self.INI_dict[key2]

        if gene == 'PR' and key2 in self.PI_dict:
            return self.PI_dict[key2]

        if gene == 'RT' and key2 in self.RTI_dict:
            return self.RTI_dict[key2]

        return 0.0


if __name__ == '__main__':
    test = NucAminoAligner()
    assert test.translateNATriplet("YTD") == "LF"
    assert test.isStopCodon("TAG") == True
    assert test.isStopCodon("TAA") == True
    assert test.isStopCodon("TGA") == True
    assert test.isStopCodon("NNN") == False

    assert test.isUnsequenced("NNN") == True
    assert test.isUnsequenced("NN-") == True
    assert test.isUnsequenced("NNG") == True
    assert test.isUnsequenced("NTG") == False

    print(test.getMutPrevalence(6, 'D', 'E', 'IN', "CRF01_AE"))
    print(test.getMutPrevalence(6, 'D', 'E', 'IN', "G"))
    print(test.getMutPrevalence(6, 'E', 'D', 'RT', "A"))
