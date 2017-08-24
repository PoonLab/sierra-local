from hivdb import HIVdb
import os


""" score_drugs function iterates through each drug in the HIV database,
    with a given sequence it will calculate a score and store it within a resulting dictionary
    
    @param HIVdb: the database
    @param sequence: the given sequence
    @return result_dict: a dictionary holding the score results of each drug for the given sequence
"""
def score_drugs(HIVdb, sequence):
    result_dict = {}
    for drug in HIVdb.keys():
        score = score_single(HIVdb, drug, sequence)
        result_dict[drug] = score
    return result_dict



""" score_single function first checks if the drug is in the HIVdb
    if found, calculates score with a given drug and sequence according to Stanford algorithm

    @param drugname: name of the drug you want the score for
    @param sequence: user provided sequence of type str (tolerates whitespace on either side, will strip it out later)
    @return score: calculated drm mutation score
"""
def score_single(HIVdb, drugname, sequence):
    # assert drugname in HIVdb.keys(), "Drugname: %s not found." % drugname
    score = 0
    for condition in HIVdb[drugname]:

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

        iter = 0  # iter keeps track of the associated index in the values list
        for residueAA in residueAAtuples:
            count = 0  # count makes sure all the tuples conditions within a residueAAtuples group is satisfied
            for tuple in residueAA:
                # TODO: if IndexError, continue as well (don't throw error just because sequence isn't that long)
                # TODO: also might qualify for a FAIL status for this particular sequence (update data structure for a PASS/FAIL)
                if not sequence[tuple[0] - 1] in tuple[1]:
                    continue
                else:
                    count += 1
                if count == len(residueAA):
                    candidates.append(values[iter])
            iter += 1

        # take the max of what's in the list of potential scores (candidates) and update total score
        # doesn't matter for the single drm or combo condition because they only have one associated value anyways
        score += max(candidates)

    return score




""" testing """
def main():
    path = os.getcwd() + '/HIVDB.xml'
    algorithm = HIVdb(path)
    definitions = algorithm.parse_definitions(algorithm.root)
    #print(definitions['gene'])
    #print(definitions['level'])
    #print(definitions['drugclass'])
    #print(definitions['globalrange'])
    #print(definitions['comment'])
    database = algorithm.parse_drugs(algorithm.root)
    scores = score_drugs(database, 'DAAAAAGAAEFHKJDLSHJDFKSLDHJFKSLAFLAHSJDKFLAHKDFLAHSASDFKASDFLASDJFKALSDFRETORJTIETGOENRTIEROTNOERNTIENTIERTERTERJDKFLASKDJFHAKSIEKRJRNNFMFMMMFJHAKDHJFKHJFKHDKJSHLFKJHSKFHDSHFJDHSFFMMFMFMFMDFLAKJSDHFALKSDJHFASKDJHFALSDKJFHAJSDKFLAKSDJFHJSDKFLAKSJDFHALSKDJFHAKSDJFHALKSDJFHLKSDJHFLAKSDJHFLAKDJFHALKSJDHFALKSJDHFKAJSDHFJSKLAFKSJDHFJAKSLFJDHSFLKAJHFALSKDJFHSLADKJFAHSJDFKAJDHFKSLGAAAATCTQAAAAAAAAAA')
    #print(scores)
    comments = algorithm.parse_comments(algorithm.root)
    #print(comments)

main()

