import subprocess
import os
import csv
import re

class NucAminoAligner():
    def __init__(self, filename):
        self.inputname = filename
        self.outputname = filename.replace('.fasta','.tsv')
        self.cwd = os.getcwd()+'/'
        self.nucamino_dir = '/'.join(self.cwd.split('/')[:-2]+['NucAmino','build'])+'/'

    def align_file(self): 
        '''
        Using subprocess to call NucAmino, generates an output .tsv containing mutation data for each sequence in the FASTA file
        '''
        args = (self.nucamino_dir+"nucamino hiv1b -i {} -g=POL -o {}".format(self.inputname,self.outputname)).split()
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
        popen.wait()

    def get_mutations(self):
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
                shift = int(row[1]) - 1
                mutationpairs = row[5].split(',')
                pos = [int(re.findall(r'\d+',x.split(':')[0])[0])-shift for x in mutationpairs]
                orig_res = [re.findall(r'\D+',x.split(':')[0][:2])[0] for x in mutationpairs]
                sub_res = [re.findall(r'(?<=\d)\D', x)[0] for x in mutationpairs]
                muts.append(dict(zip(pos, zip(orig_res,sub_res))))
                outlist = []
                for index,p in enumerate(pos):
                    outlist.append(str(orig_res[index]) + str(p) + str(sub_res[index]))
        return names, genes, muts

if __name__ == '__main__':
    n = NucAminoAligner('testsequences.fasta')
    n.align_file()
    print n.get_mutations()