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
        start = 2085
        gene_list = ['PR', 'RT', 'IN']
        firstNA = [int((x - start)/3) for x in [2253, 2550, 4230, 5096]]
        gene_present = [firstNA[i+1] > firstAA > x for i, x in enumerate(firstNA[:-1:])]
        if any(gene_present):
            return [gene_list[gene_present.index(True)]]
        return []

    def get_mutations(self, genemap):
        '''
        From the tsv output of NucAmino, parses and adjusts indices and returns as lists.
        :return: list of sequence names, list of sequence mutation dictionaries.
        '''
        muts = []
        genes = []
        names = []
        
        with open(self.outputname,'r') as tsvin:
            tsvin = csv.reader(tsvin, delimiter='\t')
            next(tsvin)
            for row in tsvin:
                names.append(row[0])
                # Let's use the gene map to figure the "reference position" for each gene
                # Nucamino outputs wonky positions 
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
                    mutationpairs = row[5].split(',') #split list into individual mutations
                    pos = [int(re.findall(r'\d+',x.split(':')[0])[0])-shift for x in mutationpairs]
                    orig_res = [re.findall(r'\D+',x.split(':')[0][:2])[0] for x in mutationpairs]
                    sub_res = [''.join(sorted(re.findall(r'(?<=\d)\D+', x.split(':')[0])[0])) for x in mutationpairs]
                    codon = [x.split(':')[1] for x in mutationpairs]
                    #deletion mixture handling
                    # sub_res = [sub_res[i] if not x else 'X' for i, x in enumerate(codon)]
                    # print(sub_res)
                    with open('codons.txt', 'a') as codonfile:
                        codonfile.write(str(codon)+'\n')

                    gene_muts = dict(zip(pos, zip(orig_res,sub_res)))
                muts.append(gene_muts)
                genelist = self.get_genes(firstAA)
                genes.append(genelist)
        assert len(muts) == len(names), "length of mutations dicts is not the same as length of names"
        return names, genes, muts
