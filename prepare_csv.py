import csv

lower = 1000
upper = 2000

output = open('IN.fasta','w')
with open('IN-data.tsv','r') as file:
    csvreader = list(csv.reader(file, delimiter='\t'))
    for line in csvreader[lower:upper]:
        header = '>'+str(line[2])+'.'+str(line[3]) #accession + Country
        seq = str(line[296]) # nucleotide sequence data
        if not len(seq) == 0: # Exclude those with no sequence data
            output.write(header+'\n'+seq+'\n')
output.close()

# regenerate the IN-sequence .tsv from the text file oops