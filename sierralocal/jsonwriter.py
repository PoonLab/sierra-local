import json
import xml.etree.ElementTree as xml
import re
from sierralocal.hivdb import HIVdb
import csv
from pathlib import Path
import sys, os
from sierralocal.utils import generateTable, enumerateCodonPossibilities
from sierralocal.nucaminohook import NucAminoAligner
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
        dest = algorithm.json_filename
        with open(dest,'r') as csvfile:
            self.ApobecDRMs = json.load(csvfile)

        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'INSTI-comments.csv')
        with open(dest, 'r') as INSTI_file:
            self.INSTI_comments = dict(csv.reader(INSTI_file, delimiter='\t'))

        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'PI-comments.csv')
        with open(dest, 'r') as PI_file:
            self.PI_comments = dict(csv.reader(PI_file, delimiter='\t'))

        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'RT-comments.csv')
        with open(dest, 'r') as RT_file:
            self.RT_comments = dict(csv.reader(RT_file, delimiter='\t'))

        # make dictionary for isUnusual
        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'hivfacts' / 'data' / 'aapcnt' / 'rx-all_subtype-all.csv')
        with open(dest, 'r', encoding='utf-8-sig') as isUnusual_file:
            isUnusual_file = csv.DictReader(isUnusual_file)
            self.rx_all_subtype_all = {}
            for row in isUnusual_file:
                gene = row['gene']
                pos = row['position']
                aa = row['aa']
                unusual = row['isUnusual']
                if gene not in self.rx_all_subtype_all:
                    self.rx_all_subtype_all.update({gene: {}})
                if pos not in self.rx_all_subtype_all[gene]:
                    self.rx_all_subtype_all[gene].update({pos: {}})
                self.rx_all_subtype_all[gene][pos].update({aa: unusual})

        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'hivfacts' / 'data' / 'sdrms_hiv1.csv')
        with open(dest, 'r', encoding='utf-8-sig') as SDRM1_files:
            SDRM1_files = csv.DictReader(SDRM1_files)
            self.SDRMs = {}
            for row in SDRM1_files:
                gene = row['gene']
                pos = row['position']
                aa = row['aa']
                if gene not in self.SDRMs:
                    self.SDRMs.update({gene: {}})
                if pos not in self.SDRMs[gene]:
                    self.SDRMs[gene].update({pos: aa})

        # make dictionary for APOBEC DRMS
        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'hivfacts' / 'data' / 'apobecs' / 'apobec_drms.csv')
        with open(dest, 'r', encoding='utf-8-sig') as APOBEC_file:
            APOBEC_file = csv.DictReader(APOBEC_file)
            self.isApobecDRMs = {}
            for row in APOBEC_file:
                gene = row['gene']
                pos = row['position']
                aa = row['aa']
                if gene not in self.isApobecDRMs:
                    self.isApobecDRMs.update({gene: {}})
                if pos not in self.isApobecDRMs[gene]:
                    self.isApobecDRMs[gene].update({pos: {aa}})

        # make dictionary for primary type
        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'hivfacts' / 'data' / 'mutation-type-pairs_hiv1.csv')
        with open(dest, 'r', encoding='utf-8-sig') as mutTypePairs1_files:
            mutTypePairs1_files = csv.DictReader(mutTypePairs1_files)
            self.mutTypePairs = {}
            for row in mutTypePairs1_files:
                gene = row['gene']
                pos = row['position']
                aa = row['aas']
                mut = row['mutationType']
                if gene not in self.mutTypePairs:
                    self.mutTypePairs.update({gene: {}})
                if pos not in self.mutTypePairs[gene]:
                    self.mutTypePairs[gene].update({pos: {}})
                self.mutTypePairs[gene][pos].update({aa: mut})

        # make dictionary for apobec mutations
        dest = str(Path(os.path.dirname(__file__)) / 'data' / 'hivfacts' / 'data' / 'apobecs' / 'apobecs.csv')
        with open(dest, 'r', encoding='utf-8-sig') as apobec_muts:
            apobec_muts = csv.DictReader(apobec_muts)
            self.Apobec_Mutations = {}
            for row in apobec_muts:
                gene = row['gene']
                pos = row['position']
                aa = row['aa']
                if gene not in self.Apobec_Mutations:
                    self.Apobec_Mutations.update({gene: {}})
                if pos not in self.Apobec_Mutations[gene]:
                    self.Apobec_Mutations[gene].update({pos: ""})
                if aa not in self.Apobec_Mutations[gene][pos]:
                    self.Apobec_Mutations[gene][pos] = self.Apobec_Mutations[gene][pos] + aa


    def formatValidationResults(self, validated):
        """
        Returns a list of dictionaries detailing validation results, meant for results output.
        @param validated: list of tuples from validateSequence() function
        @return validationResults: list of dictionaries
        """
        validationResults = [{'level': v[0], 'message': v[1]} for v in validated]
        return validationResults

    def formatDrugResistance(self, scores, gene):
        """
        Returns formatted drug resistance and score breakdowns, meant for results output.
        @param scores: results of one query from score_alg.score_drugs()
        @param genes: gene found in the sequence
        @return drugResistance: a list of one dictionary encoding scores and descriptions
        """
        drugResistance = {}
        drugResistance['version'] = {}
        drugResistance['version']['text'] = self.algorithm.version
        drugResistance['version']['publishDate'] = self.algorithm.version_date
        drugResistance['gene'] = {'name': gene}
        drugScores = []

        for drugclass in self.definitions['gene'][gene]:
            classlist = self.definitions['drugclass'][drugclass]
            for drug in classlist:
                drugScore = {}

                # Infer resistance level text from the score and globalrange
                resistancelevel = -1
                for key in self.globalrange:
                    maximum = float(self.globalrange[key]['max'])
                    minimum = float(self.globalrange[key]['min'])
                    if minimum <= scores[drug][0] <= maximum:
                        resistancelevel = str(key)
                        break
                resistance_text = self.levels[resistancelevel]

                drugScore['drugClass'] = {'name': drugclass}  # e.g. NRTI
                drugScore['drug'] = {}
                if drug in self.names:
                    drugScore['drug']['name'] = self.names[drug]
                else:
                    drugScore['drug']['name'] = drug.replace('/r', '')
                drugScore['drug']['displayAbbr'] = drug

                # infer level from score of this drug
                if scores[drug][0] < 10:
                    drugScore['level'] = 1
                elif scores[drug][0] < 15:
                    drugScore['level'] = 2
                elif scores[drug][0] < 30:
                    drugScore['level'] = 3
                elif scores[drug][0] < 60:
                    drugScore['level'] = 4
                else:
                    drugScore['level'] = 5

                # create partial score, for each mutation, datastructure
                drugScore['partialScores'] = []
                for index, pscore in enumerate(scores[drug][1]):
                    pscore = float(pscore)
                    if not pscore == 0.0:
                        pscoredict = {}
                        pscoredict['mutations'] = []
                        for combination in scores[drug][2][index]:
                            # find the mutation classification "type" based on the gene
                            type_ = drugclass
                            pos = re.findall(u'([0-9]+)', combination)[0]
                            muts = re.findall(u'(?<=[0-9])([A-Za-z])+', combination)[0]
                            # print(pos, muts)
                            if gene == 'IN':
                                for key in self.INSTI_comments:
                                    if pos in key and muts in key:
                                        type_ = self.INSTI_comments[key]
                                        break
                            elif gene == 'PR':
                                for key in self.PI_comments:
                                    if pos in key and muts in key:
                                        type_ = self.PI_comments[key]
                                        break
                            elif gene == 'RT':
                                for key in self.RT_comments:
                                    if pos in key and muts in key:
                                        type_ = self.RT_comments[key]
                                        break
                            mut = {}
                            mut['text'] = combination.replace('d', 'Deletion')
                            mut['primaryType'] = type_
                            mut['comments'] = [{
                                'type': type_,
                                'text': self.findComment(gene, combination, self.comments, self.definitions['comment'])
                            }]
                            pscoredict['mutations'].append(mut)
                        pscoredict['score'] = pscore
                        drugScore['partialScores'].append(pscoredict)
                drugScore['text'] = list(resistance_text)[0]  # resistance level
                drugScores.append(drugScore)
        drugResistance['drugScores'] = drugScores
        return drugResistance

    def formatAlignedGeneSequences(self, ordered_mutation_list, gene, firstlastAA, aligned_NA, sequence_pos):
        """
        Main function to format mutations into dictionary for results output
        @param ordered_mutation_list: ordered list of mutations in the query sequence relative to reference
        @param genes: genes found in the query sequence
        @param firstlastNA: positions tuple of the first and last nucleotide in the query sequence
        @return: list of one dictionary element describing mutations in a single query sequence
        """
        dic = {}
        dic['firstAA'] = int(firstlastAA[0])
        dic['lastAA'] = int(firstlastAA[1])
        dic['gene'] = {'name': gene, 'length': None}
        dic['mutations'] = []
        dic['SDRMs'] = []
        mutation_line = []
        mutation_line.extend([" - " for i in range(int(firstlastAA[0]), int(firstlastAA[1]) + 1)])

        for mutation in ordered_mutation_list:
            mutdict = {}
            mutdict['consensus'] = mutation[2]
            mutdict['position'] = int(mutation[0])
            mutdict['AAs'] = "".join(sorted(mutation[1]))
            mutdict['isInsertion'] = mutation[1] == '_'
            mutdict['isDeletion'] = mutation[1] == '-'
            mutdict['isApobecMutation'] = self.isApobecMutation(gene, mutation[0], mutation[1])
            mutdict['isApobecDRM'] = self.isApobecDRM(gene, mutation[2], mutation[0], mutation[1])
            mutdict['isUnusual'] = self.isUnusual(gene, mutation[0], mutation[1])
            mutdict['isSDRM'] = self.isSDRM(gene, mutation[0], mutation[1])
            if self.isSDRM(gene, mutation[0], mutation[1]):
                dic['SDRMs'].append({'text': mutation[2] + str(mutation[0]) + "".join(sorted(mutation[1]))})
            mutdict['hasStop'] = self.hasStop(mutation)
            mutdict['primaryType'] = self.primaryType(gene, mutation[0], mutation[1])
            mutdict['text'] = mutation[2] + str(mutation[0]) + "".join(sorted(mutation[1]))
            if int(firstlastAA[0]) <= int(mutation[0]) <= int(firstlastAA[1]):
                mutation_line[int(mutation[0]) - int(firstlastAA[1]) - 1] = f"{''.join(sorted(mutation[1])):^3}"
            dic['mutations'].append(mutdict)

        dic['alignedNAs'] = aligned_NA
        nuc = self.prettyNAsequence(aligned_NA)
        AAs, dic['alignedAAs'] = self.prettyAAsequence(aligned_NA)

        dic['prettyPairwise'] = {"positionLine": [f"{i:^3}" for i in range(int(firstlastAA[0]), int(firstlastAA[1]) + 1)],
                                 "refAALine": AAs,
                                 "alignedNAsLine": nuc,
                                 "mutationLine": mutation_line
                                }

        return dic

    def formatInputSequence(self, header, sequence):
        out = {
            'header': header,
            'SHA512': hashlib.sha512(str.encode(sequence)).hexdigest()
        }
        return out

    def write_to_json(self, filename, file_headers, file_scores, file_genes, file_mutation_lists,
                      file_sequence_lengths, file_trims, file_subtypes, alignedNA, sequence_pos):
        '''
        The main function to write passed result to a JSON file
        :param filename: the filename to write the JSON to
        :param names: list of sequence headers
        :param scores: list of sequence scores
        :param file_genes: list of single genes in queries
        :param ordered_mutation_list: ordered list of mutations in the query sequence relative to reference
        '''
        out = []
        for index, scores in enumerate(file_scores):
            genes = file_genes[index]

            data = {}
            data['inputSequence'] = self.formatInputSequence(file_headers[index], alignedNA[file_headers[index]])
            data['strain'] = {"name": None}  # TODO: add strain
            data['subtypeText'] = file_subtypes[index]

            validation = self.validateSequence(genes, file_sequence_lengths[index], file_trims[index])
            data['validationResults'] = self.formatValidationResults(validation)

            data['alignedGeneSequences'] = []
            data['drugResistance'] = []
            if not 'CRITICAL' in validation:
                for idx, gene_info in enumerate(genes):
                    gene, firstAA, lastAA, firstNA, lastNA = gene_info
                    omlist = file_mutation_lists[index][idx]
                    nalist = (firstAA, lastAA)
                    data['alignedGeneSequences'].append(
                        self.formatAlignedGeneSequences(omlist, gene, nalist, alignedNA[file_headers[index]],
                                                        sequence_pos)
                    )
                    data['drugResistance'].append(self.formatDrugResistance(scores[idx], gene))

            out.append(data)

        # dump data to JSON file
        with open(filename, 'w+') as outfile:
            json.dump(out, outfile, indent=2)
            print("Writing JSON to file {}".format(filename))

    def validateSequence(self, genes, lengths, seq_trims):
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

    def findComment(self, gene, mutation, comments, details):
        trunc_mut = re.findall(r'\d+\D', mutation)[0]  # 163K
        pos = re.findall(u'([0-9]+)', trunc_mut)[0]
        muts = re.findall(u'(?<=[0-9])([A-Za-z])+', trunc_mut)[0]
        for g, mutationdict in comments.items():
            for item in mutationdict.keys():
                if pos in item and muts in item:
                    full_mut = mutationdict[item]
                    if full_mut in details and g == gene:
                        return details[full_mut]['1']

    def isApobecDRM(self, gene, consensus, position, AA):
        ls = [[row['gene'], str(row['position'])] for row in self.ApobecDRMs]
        if [gene, str(position)] in ls:
            i = ls.index([gene, str(position)])
            for aa in AA:
                if aa in self.ApobecDRMs[i]['aa']:
                    return True
        return False

    def isSDRM(self, gene, position, AA):
        position = str(position)
        if gene in self.SDRMs:
            if position in self.SDRMs[gene]:
                for aa in AA:
                    if aa in self.SDRMs[gene][position]:
                        return True
        return False

    def hasStop(self, ordered_mut_list_index):
        if "*" in ordered_mut_list_index[1]:
            return True
        return False

    def isApobecMutation(self, gene, position, AA):
        position = str(position)
        if gene in self.Apobec_Mutations:
            if position in self.Apobec_Mutations[gene]:
                for aa in AA:
                    if aa in self.isApobecDRMs[gene][position]:
                        return True
        return False

    def isUnusual(self, gene, position, AA):
        position = str(position)
        if gene in self.rx_all_subtype_all:
            if position in self.rx_all_subtype_all[gene]:
                for aa in AA:
                    if aa in self.rx_all_subtype_all[gene][position]:
                        if self.rx_all_subtype_all[gene][position][aa].lower() in "true":
                            return True
        return False

    def primaryType(self, gene, position, AA):
        position = str(position)
        if gene in self.mutTypePairs:
            if position in self.mutTypePairs[gene]:
                keys = self.mutTypePairs[gene][position].keys()
                for aa in AA:
                    for key in keys:
                        if aa in key:
                            return self.mutTypePairs[gene][position][key]
        return "Other"

    def prettyNAsequence(self, sequence):
        """
        param: sequence, string of NA of specific sequence
        return: nuc, list of NA separated into codons
        """
        nuc = []
        for index in range(0, len(sequence), 3):
            nuc.append(f"{sequence[index:index + 3]:^3}")
        return nuc

    def prettyAAsequence(self, sequence):
        """
        takes NA, generates codon table, then converts NA to AA while accounting for deletions or insertions
        param: sequence, NA string of strain
        return List of AA and a string of the previous list as AAs joined
        """
        table = generateTable()
        aligned_NA = ''.join(e for e in sequence if e.isalnum())
        AAs = []
        for index in range(0, len(aligned_NA), 3):
            if len(aligned_NA[index:index + 3]) == 3:
                AAs.append(f"{table[aligned_NA[index:index + 3]]:^3}")
        return AAs, ''.join(AAs).replace(" ", "")


if __name__ == "__main__":
    writer = JSONWriter(HIVdb())
    assert (writer.isApobecDRM("IN", "G", 163, "TRAG")) == True
