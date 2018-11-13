from sierralocal.hivdb import HIVdb
import os


def score_drugs(HIVdb, gene, seq_mutations):
    """ score_drugs function iterates through each drug in the HIV database,
        with a given sequence it will calculate a score and store it within a resulting dictionary

        @param HIVdb: the database
        @param sequence: the given sequence
        @return result_dict: a dictionary holding the score results of each drug for the given sequence
    """
    result_dict = {}
    for drug_class in HIVdb.definitions['gene'][gene]:
        for drug in HIVdb.definitions['drugclass'][drug_class]:
            result_dict.update({drug: score_single(HIVdb, drug, seq_mutations)})
    return result_dict


def score_single(HIVdb, drugname, seq_mutations):
    """ score_single function first checks if the drug is in the HIVdb
        if found, calculates score with a given drug and sequence according to Stanford algorithm

        @param drugname: name of the drug you want the score for
        @param sequence: user provided sequence of type str (tolerates whitespace on either side, will strip it out later)
        @return score: calculated drm mutation score
    """
    assert drugname in HIVdb.drugs.keys(), "Drugname: %s not found." % drugname
    rec = lambda x: sum(map(rec, x)) if isinstance(x, list) else x

    sequence_partial_scores = []
    sequence_DRMs = []
    sequence_DRM_positions = []

    for condition in HIVdb.drugs[drugname][0]: # condition = potential DRMs associated with each drug
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

        # iterate thru conditions e.g. ([(41, 'L'), (215, 'FY')])
        for i, residueAAtuple in enumerate(residueAAtuples):
            condition_present = True
            DRM_mutations = []
            DRM_positions = []

            # iterate thru DRM tuples in each condition
            for j, mutationpair in enumerate(residueAAtuple):
                position = mutationpair[0]  #e.g. 41
                aminoacidlist = mutationpair[1]  #e.g. 'L'

                DRM_present = False  # assume DRM unfulfilled
                if position in seq_mutations:
                    # check if the position in the DRM is DRM_present in the sequence mutation list
                    for seq_mutation in seq_mutations[position][1]:
                        if '-' in seq_mutations[position][1]:
                            seq_mutation = 'd'
                        if '_' in seq_mutations[position][1]:
                            seq_mutation = 'i'
                        if seq_mutation in aminoacidlist:
                            if not position in DRM_positions:
                                # and not DRM_positions+[residue] in sequence_DRM_positions:
                                mut = str(seq_mutations[position][0]) + str(position) + str(seq_mutation)
                                DRM_positions.append(position)
                                DRM_mutations.append(mut)
                                DRM_present = True

                if not DRM_present:
                    condition_present = False

            if condition_present:
                sequence_DRM_positions.append(DRM_positions)
                sequence_partial_scores.append(penalties[i])
                sequence_DRMs.append(DRM_mutations)

    # filter out overlapping DRM combinations that have sub-maximal scores
    if len(sequence_DRM_positions) > 1:
        max_mask = []
        for index, combination_positions in enumerate(sequence_DRM_positions):
            if sequence_DRM_positions.count(combination_positions) == 1: # the only occurrence
                max_mask.append(True)
            else: # multiple occurrences
                indices = [i for i, x in enumerate(sequence_DRM_positions) if x == combination_positions]
                mx = rec(sequence_partial_scores[indices[0]])
                max_index = indices[0]
                # print(str(drugname), str(sequence_DRM_positions), str(indices), mx)
                for j in indices:
                    if rec(sequence_partial_scores[j]) > mx:
                        mx = rec(sequence_partial_scores[j])
                        max_index = j
                if max_index == index:
                    max_mask.append(True)
                else:
                    max_mask.append(False)
        sequence_partial_scores = [x for i,x in enumerate(sequence_partial_scores) if max_mask[i]]
        sequence_DRMs = [x for i,x in enumerate(sequence_DRMs) if max_mask[i]]

    total_score = rec(sequence_partial_scores)

    return total_score, sequence_partial_scores, sequence_DRMs
