import csv

lower = 1
upper = 30

output = open('IN.fasta','w')
with open('IN-data.tsv','r') as file:
    csvreader = list(csv.reader(file, delimiter='\t'))
    for line in csvreader[lower:upper]:
        header = '>'+str(line[2])+'.'+str(line[3])
        seq = str(line[296])
        if not len(seq) == 0:
            output.write(header+'\n'+seq+'\n')
output.close()

# regenerate the IN-sequence .tsv from the text file oops