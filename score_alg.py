from hivdb import HIVdb
import os


""" score_drugs function iterates through each drug in the HIV database,
    with a given sequence it will calculate a score and store it within a resulting dictionary
    
    @param HIVdb: the database
    @param sequence: the given sequence
    @return result_dict: a dictionary holding the score results of each drug for the given sequence
"""
def score_drugs(HIVdb, seq_mutations, sequence):
    result_dict = {}
    for index, drug in enumerate(HIVdb):
        score = score_single(HIVdb, drug, seq_mutations, sequence)
        result_dict[drug] = score
    return result_dict



""" score_single function first checks if the drug is in the HIVdb
    if found, calculates score with a given drug and sequence according to Stanford algorithm

    @param drugname: name of the drug you want the score for
    @param sequence: user provided sequence of type str (tolerates whitespace on either side, will strip it out later)
    @return score: calculated drm mutation score
"""
def score_single(HIVdb, drugname, seq_mutations, sequence):
    # print("seq_mutations")
    # print(seq_mutations)
    assert drugname in HIVdb.keys(), "Drugname: %s not found." % drugname
    totalScore = 0
    partialScores = []
    mutations = []
    for condition in HIVdb[drugname][0]: # condition = potential DRMs associated with each drug
        candidates = [0]        # list of potential scores
        values = []
        residueAAtuples = []

        # separating values from groups of tuples in each condition, both appended to lists
        for gv_pairs in condition:
            # 'AND' or 'single-drm' condition
            if isinstance(gv_pairs, str):
                if gv_pairs == 'value':
                    values.append(condition[gv_pairs])
                else:
                    residueAAtuples.append(condition[gv_pairs])
            # 'MAX' or 'MAXAND' condition
            else:
                for item in gv_pairs:
                    if item == 'value':
                        values.append(gv_pairs[item])
                    else:
                        residueAAtuples.append(gv_pairs[item])

        # print("residueAAtuples")
        # print(residueAAtuples)
        residuelist = []

        for index, residueAA in enumerate(residueAAtuples): #iterate thru conditions e.g. ([(41, 'L'), (215, 'FY')])
            conditionTrue = True
            sequence_residues = []
            partialresidues = []
            for mutationpair in residueAA: # iterate thru DRM tuples in each condition
                residue = mutationpair[0]  #e.g. 41
                aminoacidlist = mutationpair[1] #e.g. 'L'

                present = False #assume DRM fulfilled
                if residue in seq_mutations: #check if the residue in the DRM is present in the sequence mutation list
                    for possibility in seq_mutations[residue][1]:
                        if possibility in aminoacidlist:
                            if not residue in partialresidues and not [residue] in residuelist:
                                partialresidues.append(residue)
                                sequence_residues.append(str(seq_mutations[residue][0])+str(residue)+str(possibility))
                                present = True
                                
                    if residue == 143 and present:
                        print(mutationpair, seq_mutations[residue], partialresidues)

                if not present:
                    conditionTrue = False

            if conditionTrue and sequence_residues != []:
                residuelist.append(partialresidues)
                totalScore += values[index]
                partialScores.append(values[index])
                mutations.append(sequence_residues)
        # condition = True
        # sequence_residues = []
        # for index, residueAAtuple in enumerate(residueAAtuples): #iterate over subconditions within a condition
        #     # print(residueAAtuple)
        #     subcondition = True
        #     for tup in residueAAtuple: #all elements of a subcondition must be true
        #         pair = False
        #         residue = tup[0]
        #         aminoacidlist = tup[1]
        #         if residue in seq_mutations:
        #             for possibility in seq_mutations[residue][1]:
        #                 if possibility in aminoacidlist: #tuple is TRUE
        #                     sequence_residues.append(str(seq_mutations[residue][0])+str(residue)+str(seq_mutations[residue][1]))
        #                     pair = True
        #                     break
        #         if pair == False:
        #             subcondition == False
        #             break
        #     if subcondition == False:
        #         condition == False

        # if condition and sequence_residues != []:
        #     totalScore += values[index]
        #     partialScores.append(values[index])
        #     mutations.append(sequence_residues)





    #print mutations
    return totalScore, partialScores, mutations