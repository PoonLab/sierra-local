import csv

lower = 1
upper = 2649
count = 0
gene = 'RT'
output = open('{}.fasta'.format(gene),'w')
with open('{}-data.tsv'.format(gene),'r') as file:
    csvreader = list(csv.reader(file, delimiter='\t'))
    width = len(csvreader[0])
    for line in csvreader[lower:upper]:
        if count == 1000:
            break
        header = '>'+str(line[2])+'.'+str(line[3]) #accession + Country
        seq = str(line[width-1]) # nucleotide sequence data
        if not len(seq) == 0: # Exclude those with no sequence data
            output.write(header+'\n'+seq+'\n')
            count += 1

print count
output.close()

# regenerate the IN-sequence .tsv from the text file oops