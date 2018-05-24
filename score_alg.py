from hivdb import HIVdb
import os


""" score_drugs function iterates through each drug in the HIV database,
    with a given sequence it will calculate a score and store it within a resulting dictionary
    
    @param HIVdb: the database
    @param sequence: the given sequence
    @param codon_type: list of strings describing each mutation: "Normal", "Frameshift", or "Deletion"
    @return result_dict: a dictionary holding the score results of each drug for the given sequence
"""
def score_drugs(HIVdb, seq_mutations, codon_type):
    result_dict = {}
    for index, drug in enumerate(HIVdb):
        score = score_single(HIVdb, drug, seq_mutations, codon_type)
        result_dict[drug] = score
    return result_dict



""" score_single function first checks if the drug is in the HIVdb
    if found, calculates score with a given drug and sequence according to Stanford algorithm

    @param drugname: name of the drug you want the score for
    @param sequence: user provided sequence of type str (tolerates whitespace on either side, will strip it out later)
    @return score: calculated drm mutation score
"""
def score_single(HIVdb, drugname, seq_mutations, codon_type):
    assert drugname in HIVdb.keys(), "Drugname: %s not found." % drugname
    rec = lambda x: sum(map(rec, x)) if isinstance(x, list) else x
    single_score = 0
    partial_scores = []
    sequence_mutations = []
    residuelist = []

    for condition in HIVdb[drugname][0]: # condition = potential DRMs associated with each drug
        # candidates = [0]        # list of potential scores
        penalties = []
        residueAAtuples = []

        # separating penalty values from groups of tuples in each condition, both appended to lists
        for gv_pairs in condition:
            # 'AND' or 'single-drm' condition
            if isinstance(gv_pairs, str):
                if gv_pairs == 'value':
                    penalties.append(condition[gv_pairs])
                else:
                    residueAAtuples.append(condition[gv_pairs])
            # 'MAX' or 'MAXAND' condition
            else:
                for item in gv_pairs:
                    if item == 'value':
                        penalties.append(gv_pairs[item])
                    else:
                        residueAAtuples.append(gv_pairs[item])


        for index, residueAAtuple in enumerate(residueAAtuples): #iterate thru conditions e.g. ([(41, 'L'), (215, 'FY')])
            conditionTrue = True
            combination_mutations = []
            combination_positions = []
            for j, mutationpair in enumerate(residueAAtuple): # iterate thru DRM tuples in each condition
                position = mutationpair[0]  #e.g. 41
                aalist = mutationpair[1] #e.g. 'L'

                present = False #assume DRM unfulfilled
                if position in seq_mutations: #check if the position in the DRM is present in the sequence mutation list
                    if position == 69 and ('X' in seq_mutations[position][1] or '-' in seq_mutations[position][1]):
                        print(aalist)
                    for seq_mutation in seq_mutations[position][1]:
                        if '-' in seq_mutations[position][1]:
                            seq_mutation = 'd'
                        if seq_mutation in aalist:
                            if not position in combination_positions: # and not combination_positions+[residue] in residuelist:
                                mut = str(seq_mutations[position][0])+str(position)+str(seq_mutation)
                                combination_positions.append(position)
                                combination_mutations.append(mut)
                                present = True

                                
                if not present:
                    conditionTrue = False

            if conditionTrue:
                residuelist.append(combination_positions)
                partial_scores.append(penalties[index])
                sequence_mutations.append(combination_mutations)

    if len(residuelist) > 1:
        max_mask = []
        for index, combination_positions in enumerate(residuelist):
            if residuelist.count(combination_positions) == 1: # the only occurrence
                max_mask.append(True)
            else: # multiple occurrences
                indices = [i for i, x in enumerate(residuelist) if x == combination_positions]
                mx = rec(partial_scores[indices[0]])
                max_index = indices[0]
                # print(str(drugname), str(residuelist), str(indices), mx)
                for j in indices:
                    if rec(partial_scores[j]) > mx:
                        mx = rec(partial_scores[j])
                        max_index = j
                if max_index == index:
                    max_mask.append(True)
                else:
                    max_mask.append(False)
        partial_scores = [x for i,x in enumerate(partial_scores) if max_mask[i]]
        sequence_mutations = [x for i,x in enumerate(sequence_mutations) if max_mask[i]]

    single_score = rec(partial_scores)

    return single_score, partial_scores, sequence_mutations