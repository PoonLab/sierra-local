import re

def score_drugs(algorithm, gene, seq_mutations):
    """
    score_drugs function iterates through each drug in the HIV database,
    with a given sequence it will calculate a score and store it within 
    a resulting dictionary
    @param algorithm: sierralocal.hivdb.HIVdb, the database
    @param gene: str, query gene
    @param seq_mutations: dict, residue and codon pairs, keyed by 
    amino acid position
    @return result_dict: dict a dictionary holding the score results of each drug
    for the given sequence
    """
    result_dict = {}
    for drug_class in algorithm.definitions['gene'][gene]:
        for drug in algorithm.definitions['drugclass'][drug_class]:
            result_dict.update({drug: score_single(algorithm, drug, seq_mutations)})
    return result_dict


def score_single(algorithm, drugname, seq_mutations):
    """
    score_single function first checks if the drug is in the HIVdb if found,
    calculates score with a given drug and sequence according to Stanford 
    algorithm
    @param algorithm: sierralocal.hivdb.HIVdb, algorithm object
    @param drugname: str, name of the drug you want the score for
    @param seq_mutations: dict, residue and codon pairs, keyed by 
    amino acid position
    @return total_score: int, calculated drm mutation score,
            sequence_partial_scores: list, list of sequence penalty values,
            sequence_drms: list, list of strings containing information about
            position in seq_mutations
    """
    assert drugname in algorithm.drugs.keys(), "Drugname: %s not found." % drugname
    rec = lambda x: sum(map(rec, x)) if isinstance(x, list) else x

    sequence_partial_scores = []
    sequence_drms = []
    sequence_drm_positions = []

    for condition in algorithm.drugs[drugname][0]: # condition = potential DRMs associated with each drug
        penalties = []
        residue_aa_tuples = []

        # separating penalty values from groups of tuples in each condition,
        # both appended to lists
        for gv_pairs in condition:
            # 'AND' or 'single-drm' condition
            if isinstance(gv_pairs, str):
                if gv_pairs == 'value':
                    penalties.append(condition[gv_pairs])
                else:
                    residue_aa_tuples.append(condition[gv_pairs])
            # 'MAX' or 'MAXAND' condition
            else:
                for item in gv_pairs:
                    if item == 'value':
                        penalties.append(gv_pairs[item])
                    else:
                        residue_aa_tuples.append(gv_pairs[item])

        # iterate thru conditions e.g. ([(41, 'L'), (215, 'FY')])
        for i, residue_aa_tuple in enumerate(residue_aa_tuples):
            condition_present = True
            drm_mutations = []
            drm_positions = []

            # iterate thru DRM tuples in each condition
            for j, mutationpair in enumerate(residue_aa_tuple):
                position = mutationpair[0]  #e.g. 41
                aminoacidlist = mutationpair[1]  #e.g. 'L'

                DRM_present = False  # assume DRM unfulfilled
                if position in seq_mutations:
                    # check if the position in the DRM is DRM_present in the
                    # sequence mutation list
                    for seq_mutation in seq_mutations[position][1]:
                        if '-' in seq_mutations[position][1]:
                            seq_mutation = 'd'
                        if '_' in seq_mutations[position][1]:
                            seq_mutation = 'i'
                        if seq_mutation in aminoacidlist:
                            if not position in drm_positions:
                                # and not drm_positions+[residue] in sequence_drm_positions:
                                mut = str(seq_mutations[position][0]) \
                                    + str(position) \
                                    + str(seq_mutation)
                                drm_positions.append(position)
                                drm_mutations.append(mut)
                                DRM_present = True

                if not DRM_present:
                    condition_present = False

            if condition_present:
                sequence_drm_positions.append(drm_positions)
                sequence_partial_scores.append(penalties[i])
                sequence_drms.append(drm_mutations)

    # filter out overlapping DRM combinations that have sub-maximal scores
    if len(sequence_drm_positions) > 1:
        max_mask = []
        for index, combination_positions in enumerate(sequence_drm_positions):
            if sequence_drm_positions.count(combination_positions) == 1: # the only occurrence
                max_mask.append(True)
            else: # multiple occurrences
                indices = [i for i, x in enumerate(sequence_drm_positions)
                           if x == combination_positions]
                mx = rec(sequence_partial_scores[indices[0]])
                max_index = indices[0]
                # print(str(drugname), str(sequence_drm_positions), str(indices), mx)
                for j in indices:
                    if rec(sequence_partial_scores[j]) > mx:
                        mx = rec(sequence_partial_scores[j])
                        max_index = j
                if max_index == index:
                    max_mask.append(True)
                else:
                    max_mask.append(False)
        
        sequence_drms = merge_drm_positions(sequence_drms, max_mask)
        sequence_partial_scores = [x for i,x in enumerate(sequence_partial_scores)
                                   if max_mask[i]]

    total_score = rec(sequence_partial_scores)
    
    return total_score, sequence_partial_scores, sequence_drms

def merge_drm_positions(sequence_drm_positions, max_mask):
    """
    Sierrapy reports all triggering AAs, this function merges the AAs of the ones we are dropping to the highest scoring one
    this function is not the most efficient, it triples over the input list

    """
    result = []
    merged = []

    # collect all True entries as is
    for idx, drm_list in enumerate(sequence_drm_positions):
        if max_mask[idx]:
            result.append(list(drm_list))  # copy to avoid mutation
            merged.append([d[:-1] for d in drm_list])  # track prefixes

    # for False entries, append suffixes to corresponding True entry
    for idx, drm_list in enumerate(sequence_drm_positions):
        if not max_mask[idx]:
            for drm in drm_list:
                prefix, suffix = drm[:-1], drm[-1]

                # Find a matching prefix in already added result
                for r_idx, prefix_list in enumerate(merged):
                    if prefix in prefix_list:
                        # Find the exact position in inner list to modify
                        for p_idx, pfx in enumerate(prefix_list):
                            if pfx == prefix:
                                existing = result[r_idx][p_idx]
                                if suffix not in existing:
                                    result[r_idx][p_idx] += suffix
                        break

    return result
