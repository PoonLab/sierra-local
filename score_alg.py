from hivdb import HIVdb
import os
import nucaminohook

cwd = os.getcwd()


""" score_drugs function iterates through each drug in the HIV database,
    with a given sequence it will calculate a score and store it within a resulting dictionary
    
    @param HIVdb: the database
    @param sequence: the given sequence
    @return result_dict: a dictionary holding the score results of each drug for the given sequence
"""
def score_drugs(HIVdb, seq_mutations):
    result_dict = {}
    for drug in HIVdb.keys():
        score = score_single(HIVdb, drug, seq_mutations)
        result_dict[drug] = score
    return result_dict



""" score_single function first checks if the drug is in the HIVdb
    if found, calculates score with a given drug and sequence according to Stanford algorithm

    @param drugname: name of the drug you want the score for
    @param sequence: user provided sequence of type str (tolerates whitespace on either side, will strip it out later)
    @return score: calculated drm mutation score
"""
def score_single(HIVdb, drugname, seq_mutations):
    assert drugname in HIVdb.keys(), "Drugname: %s not found." % drugname
    totalScore = 0
    partialScores = []
    mutations = []
    for condition in HIVdb[drugname][0]:
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

        for index, residueAA in enumerate(residueAAtuples):
            conditionTrue = True
            sequence_residues = []
            for mutationpair in residueAA:
                residue = mutationpair[0]
                aminoacidlist = mutationpair[1]
                if residue in seq_mutations.keys():
                    for possibility in seq_mutations[residue][1]:
                        if not possibility in aminoacidlist:
                            conditionTrue = False
                        else:
                            sequence_residues.append(str(seq_mutations[residue][0])+str(residue)+str(aminoacidlist))
                else:
                    conditionTrue = False
            if conditionTrue:
                totalScore += values[index]
                partialScores.append(values[index])
                mutations.append(sequence_residues)
    return totalScore, partialScores, mutations