import subprocess
import os
import csv
import re
cwd = os.getcwd()

def nucamino_align(file): 
    filename = cwd+"/"+str(file.split('.fasta')[0])
    args = (cwd+"/nucamino-linux-amd64 hiv1b -i {}.fasta -g=POL -o {}.tsv".format(filename,filename)).split()
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()

def parse_nucamino_results(filename):
    outdict = {}
    with open(filename,'r') as tsvin:
        tsvin = csv.reader(tsvin,delimiter='\t')
        data = list(tsvin)
        shift = int(data[1][1]) - 1
        mutationpairs = data[1][5].split(',')
        #print mutationpairs
        pos = [int(re.findall(r'\d+',x.split(':')[0])[0])-shift for x in mutationpairs]
        res = [re.findall(r'\D+',x.split(':')[0][1:])[0] for x in mutationpairs]
        muts = dict(zip(pos, res))
    return muts

if __name__ == '__main__':
    nucamino_align('AY030621.fasta')
    print parse_nucamino_results(cwd+'/AY030621.tsv')