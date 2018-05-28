import subprocess
import os
from pathlib import Path
import csv
import re

class NucAminoAligner():
    def __init__(self, filename):
        '''
        Initialize NucAmino for a specific input fasta file
        '''
        self.inputname = filename
        self.outputname = os.path.splitext(filename)[0] + '.tsv'
        self.cwd = os.path.curdir
        items = os.listdir(self.cwd)
        self.nucamino_binary = 'nucamino'
        for name in items:
            if 'nucamino' in name and not (name.endswith('py') or name.endswith('pyc')):
                self.nucamino_binary = name
                break
        print("Found NucAmino binary", self.nucamino_binary)

    def align_file(self): 
        '''
        Using subprocess to call NucAmino, generates an output .tsv containing mutation data for each sequence in the FASTA file
        '''
        args = ("./{} hiv1b -i {} -g=POL -o {}".format(
            self.nucamino_binary, self.inputname, self.outputname)).split()
        # print(args)
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
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
        muts = []
        genes = []
        names = []
        codon_types = []
        tripletTable = self.generateTable()
        
        #grab sequences from file
        with open(fastaFileName, 'r') as handle:
            # names = []
            sequences = []
            sequence = ''
            for i in handle:
                if i[0] == '$':  # skip h info
                    continue
                elif i[0] == '>' or i[0] == '#':
                    if len(sequence) > 0:
                        # names.append(h)
                        sequences.append(sequence)
                        sequence = ''
                        h = i.strip('\n')[1:]
                    else:
                        h = i.strip('\n')[1:]
                else:
                    sequence += i.strip('\n').upper()
            # names.append(h)
            sequences.append(sequence)

        codon_typer = lambda x: 'Deletion' if x == 'NNN' else 'Normal' #lambda function to classify codons based on missing values ('N')
        
        with open(self.outputname,'r') as tsvin: # open the NucAmino output file
            tsvin = csv.reader(tsvin, delimiter='\t')
            next(tsvin) #bypass the header row
            for idx, row in enumerate(tsvin): #iterate over sequences (1 per row)
                names.append(row[0])
                # Let's use the gene map to figure the "reference position" for each gene
                # Nucamino outputs wonky positions, so calculate the shift so we can get positions relative to the start of each individual gene
                shift = 0
                firstAA = int(row[1])
                pol_map = self.gene_map()
                for key in pol_map:
                    low, hi = pol_map[key]
                    if low <= firstAA <= hi:
                        shift = low
                        break

                # Edge case of NO mutations
                if len(row[5]) == 0:
                    pos = []
                    orig_res = []
                    sub_res = []
                    gene_muts = {}
                # most cases fall here... mutations are present
                else:
                    #generate data lists for each sequence
                    mutationpairs = row[5].split(',') #split list into individual mutations
                    pos = [int(int(re.findall(r'\d+',x.split(':')[0])[0])-shift) for x in mutationpairs] #residue position list
                    orig_res = [re.findall(r'\D+',x.split(':')[0][:2])[0] for x in mutationpairs] #consensus list
                    # sub_res = [''.join(sorted(re.findall(r'(?<=\d)\D+', x.split(':')[0])[0])) for x in mutationpairs] #mixture/mutation list
                    codon = [x.split(':')[1] for x in mutationpairs] #list of the nucleotide codons
                    query_codons = [sequences[idx][(((p+int(shift)-int(row[1])+1)-1)*3):((p+int(shift)-int(row[1])+1)*3)] for p in pos] #obtain the mutation codon directly from the query sequence
                    # if '~' in sequences[idx]:
                    #     print(row[0], shift)
                    #     print(codon)
                    #     print(query_codons)
                    #     print(pos)
                    
                    codon_type = [codon_typer(x) for x in codon] # figure out if the codons are normal, frameshifted, or deletions
                    sub_res = ['-' if x == '~~~' else self.translateNATriplet(x, tripletTable) for x in query_codons]
                    # for i, c in enumerate(codon_type): # if frameshifted or deletion, the mutation must be altered before further processing
                    #     if c == 'Deletion':
                    #         sub_res[i] = '-'
                    #         print('deletion in sequence {}'.format(row[0]))
                    gene_muts = dict(zip(pos, zip(orig_res,sub_res)))
                # append everything to list of lists of the entire sequence set
                muts.append(gene_muts)
                genelist = self.get_genes(firstAA)
                genes.append(genelist)
                codon_types.append(codon_type)
        assert len(muts) == len(names), "length of mutations dicts is not the same as length of names"
        return names, genes, muts, codon_types

    def translateNATriplet(self, triplet, tripletTable):
        if len(triplet) != 3:
            return "X"
        if '~' in triplet:
            aa = "X"
        elif triplet in tripletTable:
            aa = tripletTable[triplet]
        else:
            aa = "X"
        return aa

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
                        uniqueAAs.append(codonToAminoAcidMap[codon])
                    aas = ""
                    if len(uniqueAAs) > 4:
                        aas = "X"
                    else:
                        for uniqueAA in uniqueAAs:
                            aas += uniqueAA
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
