import subprocess
import os
import csv
import re

class NucAminoAligner():
    def __init__(self, filename):
        '''
        Initialize NucAmino for a specific input fasta file
        '''
        self.inputname = filename
        self.outputname = os.path.splitext(filename)[0] + '.tsv'
        self.cwd = os.getcwd()+'/'
        items = os.listdir(".")
        self.nucamino_binary = 'nucamino'
        for name in items:
            if 'nucamino' in name and not (name.endswith('py') or name.endswith('pyc')):
                self.nucamino_binary = name
                break
        print(self.nucamino_binary)

    def align_file(self): 
        '''
        Using subprocess to call NucAmino, generates an output .tsv containing mutation data for each sequence in the FASTA file
        '''
        args = ("./{} hiv1b -i {} -g=POL -o {}".format(
            self.nucamino_binary, self.inputname, self.outputname)).split()
        print(args)
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

    def get_mutations(self, genemap):
        '''
        From the tsv output of NucAmino, parses and adjusts indices and returns as two lists.
        :return: list of sequence names, list of sequence mutation dictionaries.
        '''
        muts = []
        genes = []
        names = []
        with open(self.cwd+self.outputname,'r') as tsvin:
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
                    sub_res = [re.findall(r'(?<=\d)\D+', x.split(':')[0])[0] for x in mutationpairs]
                    #Filter out polymorphic variants with the same residue as the reference
                    # for index, sub in enumerate(sub_res):
                    #     if len(sub) > 1:
                    #         #print(sub)
                    #         sub_res[index] = sub.replace(orig_res[index],'')
                    #         #print(sub_res[index])
                    gene_muts = dict(zip(pos, zip(orig_res,sub_res)))
                muts.append(gene_muts)
                genelist = []
                for p in pos:
                    position = p+shift
                    for key, r in genemap.items():
                        if r[0] <= position <= r[1]:
                            if key not in genelist:
                                genelist.append(key)
                if len(genelist) == 0:
                    genelist = ['RT', 'PR', 'IN']
                genes.append(genelist)
        return names, genes, muts
