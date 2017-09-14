## need a data structure for recording the location and composition of sequence insertions.
## AGTCAT---TTTGGACTC  # reference
## AGTCATCCCTTTGGA---  # query
## format of {'rpos': 1, 'seq': 'CCC'}
## codon location 0 (AGT)

def store_indels(reference, query):
    insertions = []
    deletions = []

    ## NOTE: this makes the assumption that the aligned sequence will always start with 'AGT'
    reflist = [reference[i:(i+3)] for i in range(0, len(reference), 3)]
    querylist = [query[i:(i+3)] for i in range(0, len(query), 3)]


    for index in range(0, len(reflist)):
        refcodon = reflist[index]
        querycodon = querylist[index]
        # insertions in query means there are '---' present in reference
        if refcodon.find('-') != -1:
            insertions.append({'rpos': index, 'inserted_seq': querycodon})
        # deletions in query means there are '---' present in query
        elif querycodon.find('-') != -1:
            deletions.append({'rpos': index, 'deleted_seq': refcodon})


    # may have to change the data structure of the return value
    return [insertions, deletions]










def main():
    indels = store_indels('AGTCAT---TTTGGACTC', 'AGTCATCCCTTTGGA---')
    print(indels)

main()