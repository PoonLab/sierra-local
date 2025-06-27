from ast import Index
import os
import re
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as xml
from sierralocal.hivdb import HIVdb
import hashlib


class JSONWriter():
    def __init__(self, algorithm):
        # possible alternative drug abbrvs
        self.names = {'3TC': 'LMV'}

        # Set up algorithm data
        self.algorithm = algorithm
        self.algorithm.root = xml.parse(str(self.algorithm.xml_filename)).getroot()
        self.algorithm.version = self.algorithm.root.find('ALGVERSION').text
        self.algorithm.version_date = self.algorithm.root.find('ALGDATE').text
        self.definitions = self.algorithm.parse_definitions(self.algorithm.root)
        self.levels = self.definitions['level']
        self.globalrange = self.definitions['globalrange']
        self.database = self.algorithm.parse_drugs(self.algorithm.root)
        self.comments = self.algorithm.parse_comments(self.algorithm.root)

        # Load comments files stored locally. These are distributed in the repo for now.
        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'INSTI-comments.csv')
        with open(dest, 'r') as insti_file:
            self.insti_comments = dict(csv.reader(insti_file, delimiter='\t'))

        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'PI-comments.csv')
        with open(dest, 'r') as pi_file:
            self.pi_comments = dict(csv.reader(pi_file, delimiter='\t'))

        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'RT-comments.csv')
        with open(dest, 'r') as rt_file:
            self.rt_comments = dict(csv.reader(rt_file, delimiter='\t'))

        # make dictionary for isUnusual
        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'rx-all_subtype-all.csv')
        with open(dest, 'r', encoding='utf-8-sig') as is_unusual_file:
            is_unusual_file = csv.DictReader(is_unusual_file)
            self.is_unusual_dic = {}
            for row in is_unusual_file:
                gene = row['gene']
                pos = row['position']
                aa = row['aa']
                unusual = row['isUnusual']
                if gene not in self.is_unusual_dic:
                    self.is_unusual_dic.update({gene: {}})
                if pos not in self.is_unusual_dic[gene]:
                    self.is_unusual_dic[gene].update({pos: {}})
                self.is_unusual_dic[gene][pos].update({aa: unusual})

        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'sdrms_hiv1.csv')
        with open(dest, 'r', encoding='utf-8-sig') as sdrm_files:
            sdrm_files = csv.DictReader(sdrm_files)
            self.sdrm_dic = {}
            for row in sdrm_files:
                gene = row['gene']
                pos = row['position']
                aa = row['aa']
                if gene not in self.sdrm_dic:
                    self.sdrm_dic.update({gene: {}})
                if pos not in self.sdrm_dic[gene]:
                    self.sdrm_dic[gene].update({pos: aa})
                else:
                    self.sdrm_dic[gene][pos] += aa

        # make dictionary for APOBEC DRMS
        with open(algorithm.json_filename, 'r') as jsonfile:
            self.apobec_drm_dic = {}
            for entry in json.load(jsonfile):
                gene = entry['gene']
                position = entry['position']
                aa = entry['aa']

                if not gene in self.apobec_drm_dic:
                    self.apobec_drm_dic[gene] = {}

                if not position in self.apobec_drm_dic[gene]:
                    self.apobec_drm_dic[gene].update({position: aa})
                else:
                    self.apobec_drm_dic[gene][position] += aa

        # make dictionary for primary type
        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'mutation-type-pairs_hiv1.csv')
        with open(dest, 'r', encoding='utf-8-sig') as mut_type_pairs1_files:
            mut_type_pairs1_files = csv.DictReader(mut_type_pairs1_files)
            self.primary_type_dic = {}
            for row in mut_type_pairs1_files:
                gene = row['gene']
                pos = row['position']
                aa = row['aas']
                mut = row['mutationType']
                if gene not in self.primary_type_dic:
                    self.primary_type_dic.update({gene: {}})
                if pos not in self.primary_type_dic[gene]:
                    self.primary_type_dic[gene].update({pos: {}})
                self.primary_type_dic[gene][pos].update({aa: mut})

        # make dictionary for apobec mutations
        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'apobecs.csv')
        with open(dest, 'r', encoding='utf-8-sig') as apobec_mutations:
            apobec_mutations = csv.DictReader(apobec_mutations)
            self.apobec_mutations_dic = {}
            for row in apobec_mutations:
                gene = row['gene']
                pos = row['position']
                aa = row['aa']
                if gene not in self.apobec_mutations_dic:
                    self.apobec_mutations_dic.update({gene: {}})
                if pos not in self.apobec_mutations_dic[gene]:
                    self.apobec_mutations_dic[gene].update({pos: ""})
                if aa not in self.apobec_mutations_dic[gene][pos]:
                    self.apobec_mutations_dic[gene][pos] = self.apobec_mutations_dic[gene][pos] + aa

    def format_validation_results(self, validated):
        """
        Returns a list of dictionaries detailing validation results,
        meant for results output.
        @param validated: list, list of tuples from validate_sequence() function
        @return validation_results: list, list of dictionaries
        """
        validation_results = [{'level': v[0], 'message': v[1]} for v in validated]
        return validation_results

    def format_drug_resistance(self, scores, sequence_name, gene, ambiguous, names):
        """
        Returns formatted drug resistance and score breakdowns,
        meant for results output.
        @param scores: dict, results of one query from score_alg.score_drugs() {drug: (total score, [scores/SNP], [SNPS])}
        @param gene: str, gene found in the sequence
        @param seqeunce_name: str, sequence name of fasta input, use for ambiguous dict
        @param ambiguous: dict, {sequence name: position of NNN NAs}
        @param names: list, [sequence names]
        @return drug_resistance: dict, one dictionary encoding
        scores and descriptions
        """
        # remove all ambiguous NNN positions
        for drug, info in scores.items():
            muts_scores, drug_muts = info[1], info[2]
            new = [[] for i in drug_muts] # list of all the SNPS
            scores[drug] = list(scores[drug])
            inds = set() # holds index within score[1] list that are valid

            # not fastest way, but works around not mutating while iterating
            for index, position in enumerate(drug_muts):
                for ind2, mut in enumerate(position):
                    mut_pos = re.search(r'\d+', mut).group() # 184 in M184IV
                    if int(mut_pos) in ambiguous[sequence_name][gene]:
                        inds.add(index)
                        inds.add(index)
                    else:
                        new[index].append(mut)

            new1 = [i for i in new if i] # New list of mutations that doesn't result in X AA
            # if the drug score in info[1] index is in inds, it is amb position, so exclude
            new2 = [score for index, score in enumerate(muts_scores) if index not in inds]
                        
            scores[drug] = [sum(new2), new2, new1]
        drug_resistance = {}
        drug_resistance['version'] = {}
        drug_resistance['version']['text'] = self.algorithm.version
        drug_resistance['version']['publishDate'] = self.algorithm.version_date
        drug_resistance['gene'] = {'name': gene}
        drug_scores = []

        for drugclass in self.definitions['gene'][gene]:
            classlist = self.definitions['drugclass'][drugclass]
            for drug in classlist:
                drug_score = {}

                # Infer resistance level text from the score and globalrange
                resistance_level = -1
                for key in self.globalrange:
                    maximum = float(self.globalrange[key]['max'])
                    minimum = float(self.globalrange[key]['min'])
                    if minimum <= scores[drug][0] <= maximum:
                        resistance_level = str(key)
                        break
                resistance_text = self.levels[resistance_level]

                drug_score['drugClass'] = {'name': drugclass}  # e.g. NRTI
                drug_score['drug'] = {}
                if drug in self.names:
                    drug_score['drug']['name'] = self.names[drug]
                else:
                    drug_score['drug']['name'] = drug.replace('/r', '')
                drug_score['drug']['displayAbbr'] = drug
                drug_score['score'] = scores[drug][0]

                # infer level from score of this drug
                if scores[drug][0] < 10:
                    drug_score['level'] = 1
                elif scores[drug][0] < 15:
                    drug_score['level'] = 2
                elif scores[drug][0] < 30:
                    drug_score['level'] = 3
                elif scores[drug][0] < 60:
                    drug_score['level'] = 4
                else:
                    drug_score['level'] = 5

                # create partial score, for each mutation, datastructure
                drug_score['partialScores'] = []
                for index, pscore in enumerate(scores[drug][1]):

                    pscore = float(pscore)
                    if not pscore == 0.0:
                        pscoredict = {}
                        pscoredict['mutations'] = []

                        for combination in scores[drug][2][index]:
                            # find the mutation classification "type" based on the gene
                            type_ = drugclass
                            pos = re.findall(u'([0-9]+)', combination)[0]
                            muts = re.search(r'\d([A-Za-z]+)', combination).group(1)
                            if gene == 'IN':
                                for key in self.insti_comments:
                                    if pos in key and any(c in key for c in muts):
                                        type_ = self.insti_comments[key]
                                        break
                            elif gene == 'PR':
                                for key in self.pi_comments:
                                    if pos in key and any(c in key for c in muts):
                                        type_ = self.pi_comments[key]
                                        break
                            elif gene == 'RT':
                                for key in self.rt_comments:
                                    if pos in key and any(c in key for c in muts):
                                        type_ = self.rt_comments[key]
                                        break
                            mut = {}
                            mut['text'] = combination.replace('d', 'Deletion')
                            mut['primaryType'] = type_
                            mut['triggeredAAs'] = muts
                            mut['comments'] = [{
                                'type': type_,
                                'text': self.find_comment(gene,
                                                          combination,
                                                          self.comments,
                                                          self.definitions['comment'])
                            }]
                            pscoredict['mutations'].append(mut)
                        pscoredict['score'] = pscore
                        drug_score['partialScores'].append(pscoredict)
                drug_score['text'] = list(resistance_text)[0]  # resistance level
                drug_scores.append(drug_score)
        drug_resistance['drugScores'] = drug_scores
        return drug_resistance

    def format_aligned_gene_sequences(self, ordered_mutation_list,
                                      gene, first_last_aa):
        """
        Main function to format mutations into dictionary for
        results output
        @param ordered_mutation_list: list, ordered list of mutations
        in the query sequence relative to reference
        @param gene: str, genes found in the query sequence
        @param first_last_aa: tuple, positions tuple of the first and last
        nucleotide in the query sequence
        @return dic: dict, dictionary describing mutations in a single
        query sequence
        """
        dic = {}
        dic['firstAA'] = int(first_last_aa[0])
        dic['lastAA'] = int(first_last_aa[1])
        dic['gene'] = {'name': gene, 'length': None}  # TODO: output length
        dic['mutations'] = []
        dic['SDRMs'] = []
        mutation_line = []
        mutation_line.extend([" - " for i in range(int(first_last_aa[0]), int(first_last_aa[1]) + 1)])

        for idx, mutation in enumerate(ordered_mutation_list):
            mutdict = {}
            mutdict['consensus'] = mutation[2]
            mutdict['position'] = int(mutation[0])
            mutdict['AAs'] = "".join(mutation[1])
            mutdict['isInsertion'] =  mutation[1].startswith('_')
            mutdict['isDeletion'] = mutation[1] == '-'
            mutdict['isApobecMutation'] = self.is_apobec_mutation(gene,
                                                                  mutation[0],
                                                                  mutation[1])
            mutdict['isApobecDRM'] = self.is_apobec_drm(gene,
                                                        mutation[2],
                                                        mutation[0],
                                                        mutation[1])
            mutdict['isUnusual'] = self.is_unusual(gene,
                                                   mutation[0],
                                                   mutation[1],
                                                   mutation[3])
            mutdict['isSDRM'] = self.is_sdrm(gene,
                                             mutation[0],
                                             mutation[1])
            if self.is_sdrm(gene,
                            mutation[0],
                            mutation[1]):
                dic['SDRMs'].append({'text': mutation[2] + str(mutation[0]) + mutation[3]})
            mutdict['hasStop'] = self.has_stop(mutation, mutation[3])
            mutdict['primaryType'] = self.primary_type(gene,
                                                       mutation[0],
                                                       mutation[1])
            if mutdict['AAs'] == '-':
                mutdict['text'] = mutation[2] + str(mutation[0]) + 'del'
            else:
                mutdict['text'] = mutation[2] + str(mutation[0]) + mutation[3]
            if int(first_last_aa[0]) <= int(mutation[0]) <= int(first_last_aa[1]):
                mutation_line[int(mutation[0]) - int(first_last_aa[1]) - 1] = f"{''.join(sorted(mutation[1])):^3}"
            dic['mutations'].append(mutdict)
        return dic

    def format_input_sequence(self, header, sequence):
        out = {
            'header': header,
            'SHA512': hashlib.sha512(str.encode(sequence)).hexdigest()
        }
        return out

    def write_to_json(self, filename, file_headers, file_scores,
                      file_genes, file_mutation_lists, file_sequence_lengths,
                      file_trims, file_subtypes, na_sequence, ambiguous, names):
        """
        The main function to write passed result to a JSON file
        @param filename: str, the file path to write the JSON to
        @param file_headers: list, list of sequence header strings
        @param file_scores: list, list of single genes in queries
        list of sequence scores
        @param file_genes: list, list of single genes in queries
        @param file_mutation_lists: ordered list of mutations 
        in the query sequence relative to reference
        @param file_sequence_lengths: list, list of lists of ints denoting
        sequence lengths
        @param file_trims: list, list of lists of tuples of ints 
        @param file_subtypes: list, list of subtype strings
        @param na_sequence: dict, {sequence name: associated NA sequence}
        @param ambiguous: dict, {sequence name: positions of NNN triplet NA}
        @param names: list, [sequence names]
        """
        out = []
        for index, scores in enumerate(file_scores):
            genes = file_genes[index]

            data = {}
            data['inputSequence'] = self.format_input_sequence(file_headers[index], na_sequence[file_headers[index]])
            data['subtypeText'] = file_subtypes[index]

            validation = self.validate_sequence(genes,
                                                file_sequence_lengths[index],
                                                file_trims[index])
            data['validationResults'] = self.format_validation_results(validation)

            data['alignedGeneSequences'] = []
            data['drugResistance'] = []

            if not 'CRITICAL' in validation:
                for idx, gene_info in enumerate(genes):
                    gene, first_aa, last_aa, first_na, last_na = gene_info
                    omlist = file_mutation_lists[index][idx]
                    nalist = (first_aa, last_aa)
                    data['alignedGeneSequences'].append(
                        self.format_aligned_gene_sequences(omlist,
                                                           gene,
                                                           nalist)
                    )
                    data['drugResistance'].append(
                        self.format_drug_resistance(scores[idx], file_headers[index], gene, ambiguous, names)
                    )

            out.append(data)

        # dump data to JSON file
        with open(filename, 'w+') as outfile:
            json.dump(out, outfile, indent=2)
            print("Writing JSON to file {}".format(filename))

    def validate_sequence(self, genes, lengths, seq_trims):
        """
        Function to validate a sequence and return a
        list of validation results
        @param genes: list, list of single genes in queries
        @param lengths: list, list of lists of ints denoting
        sequence lengths
        @param seq_trims: list, list of lists of tuples of ints
        @return validation_results: list, lsit of either warning 
        or critical error messages
        """
        validation_results = []

        for index, gene in enumerate(genes):
            length = lengths[index]
            seq_trim = seq_trims[index]

            # Length validation
            if ('RT' in gene and length < 200) or ('PR' in gene and length < 80) or \
                    ('IN' in gene and length < 200):
                validation_results.append(
                    ('WARNING',
                     "The {} sequence contains just {} codons, which is not sufficient "
                     "for a comprehensive interpretation.".format(gene, int(length)))
                )
            elif ('RT' in gene and length < 150) or ('PR' in gene and length < 60) or \
                    ('IN' in gene and length < 100):
                validation_results.append(
                    ('SEVERE WARNING',
                     "The {} sequence contains just {} codons, which is not "
                     "sufficient for a comprehensive interpretation.".format(gene, int(length)))
                )

            # Gene validation
            if len(gene) == 0:
                validation_results.append(('CRITICAL', 'Unable to process sequence.'))

            if seq_trim[0] > 0:
                validation_results.append(
                    ('WARNING',
                     "The {} sequence had {} amino acid{} trimmed from its 5\u2032-end "
                     "due to poor quality.".format(gene, seq_trim[0], "s" if seq_trim[0] > 1 else ""))
                )

            if seq_trim[1] > 0:
                validation_results.append(
                    ('WARNING',
                     "The {} sequence had {} amino acid{} trimmed from its 3\u2032-end "
                     "due to poor quality.".format(gene, seq_trim[1], "s" if seq_trim[1] > 1 else ""))
                )

        return validation_results

    def find_comment(self, gene, mutation, comments, details):
        """
        @param gene: str, genes found in the query sequence
        @param mutation: str, TODO: incomplete
        @param comments: dict, value in comments attribute of HIVdb object
        @param details: dict, value of HIVdb object's definitions 
        attribute's "comment" key
        @return: str, TODO: incomplete
        """
        trunc_mut = re.findall(r'\d+\D', mutation)[0]  # 163K
        pos = re.findall(u'([0-9]+)', trunc_mut)[0]
        muts = re.findall(u'(?<=[0-9])([A-Za-z])+', trunc_mut)[0]
        for g, mutationdict in comments.items():
            for item in mutationdict.keys():
                if pos in item and muts in item:
                    full_mut = mutationdict[item]
                    if full_mut in details and g == gene:
                        return details[full_mut]['1']

    def is_apobec_drm(self, gene, consensus, position, AA):
        """
        see if specific amino acid mutation is an apobec drm through checking hivbd facts
        @param gene: str, RT, IN, PR
        @param consensus: str, consensus amino acid
        @param position: int, position of mutation relative to POL
        @param AA: new amino acid
        @return: bool
        """
        if gene in self.apobec_drm_dic:
            if position in self.apobec_drm_dic[gene]:
                for aa in AA:
                    if aa in self.apobec_drm_dic[gene][position]:
                        return True
        return False

    def is_sdrm(self, gene, position, AA):
        """
        see if specific amino acid mutation is a sdrm through checking hivbd facts
        @param gene: str, RT, IN, PR
        @param position: int, position of mutation relative to POL
        @param AA: new amino acid
        @return: bool
        """
        position = str(position)
        if gene in self.sdrm_dic:
            if position in self.sdrm_dic[gene]:
                for aa in AA:
                    if aa in self.sdrm_dic[gene][position]:
                        return True
        return False

    def has_stop(self, ordered_mut_list_index, text):
        """
        see if specific amino acid mutation results in nonsense mutation
        @param: ordered_mut_list_index, list, mutations of specific NA sequence
        @param: text, str, translated amino acid from codon of mutation
        @return: bool
        """
        if text == 'X':
            return False

        if "*" in ordered_mut_list_index[1]:
            return True
        return False

    def is_apobec_mutation(self, gene, position, AA):
        """
        see if specific amino acid mutation is an apobec mutation through checking hivbd facts
        @param gene: str, RT, IN, PR
        @param position: int, position of mutation relative to POL
        @param AA: new amino acid
        @return: bool
        """
        position = str(position)
        if gene in self.apobec_mutations_dic:
            if position in self.apobec_mutations_dic[gene]:
                for aa in AA:
                    if aa in self.apobec_mutations_dic[gene][position]:
                        return True
        return False

    def is_unusual(self, gene, position, AA, text):
        """
        see if specific amino acid mutation 'is unusual' through checking hivbd facts
        @param gene: str, RT, IN, PR
        @param position: int, position of mutation relative to POL
        @param AA: new amino acid
        @return: bool
        """
        position = str(position)
        # The AA == X is stated in hivdb sierra core java files
        if AA == 'X':
            return True
        if text == 'X':  # this just fixes most of the errors, can't find source
            return False
        if gene in self.is_unusual_dic:
            if position in self.is_unusual_dic[gene]:
                for aa in AA:
                    if aa in self.is_unusual_dic[gene][position]:
                        if self.is_unusual_dic[gene][position][aa].lower() == 'true':
                            return True
        return False

    def primary_type(self, gene, position, AA):
        """
        see if specific amino acid's primary type through checking hivbd facts
        @param gene: str, RT, IN, PR
        @param position: int, position of mutation relative to POL
        @param AA: new amino acid
        @return: bool
        """
        position = str(position)
        if gene in self.primary_type_dic:
            if position in self.primary_type_dic[gene]:
                keys = self.primary_type_dic[gene][position].keys()
                for aa in AA:
                    for key in keys:
                        if aa in key:
                            return self.primary_type_dic[gene][position][key]
        return "Other"


if __name__ == "__main__":
    writer = JSONWriter(HIVdb())
    assert (writer.is_apobec_drm("IN", "G", 163, "TRAG")) == True