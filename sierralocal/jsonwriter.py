import json
import xml.etree.ElementTree as xml
import re
from sierralocal.hivdb import HIVdb
import csv
from pathlib import Path
import sys, os

class JSONWriter():
    def __init__(self, algorithm):
        # possible alternative drug abbrvs
        self.names = {'3TC':'LMV'}

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

        dest = str(Path(os.path.dirname(__file__))/'data'/'INSTI-comments.csv')
        with open(dest, 'r') as INSTI_file:
            self.INSTI_comments = dict(csv.reader(INSTI_file, delimiter='\t'))

        dest = str(Path(os.path.dirname(__file__))/'data'/'PI-comments.csv')
        with open(dest, 'r') as PI_file:
            self.PI_comments = dict(csv.reader(PI_file, delimiter='\t'))

        dest = str(Path(os.path.dirname(__file__))/'data'/'RT-comments.csv')
        with open(dest, 'r') as RT_file:
            self.RT_comments = dict(csv.reader(RT_file, delimiter='\t'))


    def formatValidationResults(self, validated):
        """
        Returns a list of dictionaries detailing validation results, meant for results output.
        @param validated: list of tuples from validateSequence() function
        @return validationResults: list of dictionaries
        """
        validationResults = [{'level':v[0], 'message':v[1]} for v in validated]
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
        drugResistance['gene'] = {'name' : gene}
        drugScores = []

        for drugclass in self.definitions['gene'][gene]:
            classlist = self.definitions['drugclass'][drugclass]
            for drug in classlist:
                drugScore = {}

                #Infer resistance level text from the score and globalrange
                resistancelevel = -1
                for key in self.globalrange:
                    maximum = float(self.globalrange[key]['max'])
                    minimum = float(self.globalrange[key]['min'])
                    if minimum <= scores[drug][0] <= maximum:
                        resistancelevel = str(key)
                        break
                resistance_text = self.levels[resistancelevel]

                drugScore['drugClass'] = {'name': drugclass} #e.g. NRTI
                drugScore['drug'] = {}
                if drug in self.names:
                    drugScore['drug']['name'] = self.names[drug]
                else:
                    drugScore['drug']['name'] = drug.replace('/r', '')

                drugScore['drug']['displayAbbr'] = drug
                drugScore['score'] = float(scores[drug][0]) # score for this paritcular drug
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
                            pos = re.findall(u'([0-9]+)',combination)[0]
                            muts = re.findall(u'(?<=[0-9])([A-Za-z])+',combination)[0]
                            #print(pos, muts)
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
                                'type' : type_,
                                'text' : self.findComment(gene, combination, self.comments, self.definitions['comment'])
                            }]
                            pscoredict['mutations'].append(mut)                    
                        pscoredict['score'] = pscore
                        drugScore['partialScores'].append(pscoredict)
                drugScore['text'] = list(resistance_text)[0] #resistance level
                drugScores.append(drugScore)
        drugResistance['drugScores'] = drugScores
        return drugResistance

    def formatAlignedGeneSequences(self, ordered_mutation_list, gene, firstlastAA):
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
        dic['gene'] = {'name': gene, 'length': None}  # TODO: output length
        dic['mutations'] = []
        for mutation in ordered_mutation_list:
            mutdict = {}
            mutdict['consensus'] = mutation[2]
            mutdict['position'] = int(mutation[0])
            mutdict['AAs'] = mutation[1]
            mutdict['isInsertion'] = mutation[1] == '_'
            mutdict['isDeletion'] = mutation[1] == '-'
            mutdict['isApobecDRM'] = self.isApobecDRM(gene, mutation[2], mutation[0], mutation[1])
            dic['mutations'].append(mutdict)
        return dic

    def formatInputSequence(self, header):
        out = {
            'header' : header
            #TODO: SHA512 hash of the gene sequence
        }
        return out

    def write_to_json(self, filename, file_headers, file_scores, file_genes, file_mutation_lists,
                      file_sequence_lengths, file_trims, file_subtypes):
        '''
        The main function to write passed result to a JSON file
        :param filename: the filename to write the JSON to
        :param names: list of sequence headers
        :param scores: list of sequence scores
        :param genes: list of single genes in queries
        :param ordered_mutation_list: ordered list of mutations in the query sequence relative to reference
        '''
        out = []
        for index, scores in enumerate(file_scores):
            genes = file_genes[index]

            data = {}
            data['inputSequence'] = self.formatInputSequence(file_headers[index])
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
                        self.formatAlignedGeneSequences(omlist, gene, nalist)
                    )
                    data['drugResistance'].append(self.formatDrugResistance(scores[idx], gene))

            out.append(data)

        # dump data to JSON file
        with open(filename,'w+') as outfile:
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
        trunc_mut = re.findall(r'\d+\D',mutation)[0] #163K
        pos = re.findall(u'([0-9]+)',trunc_mut)[0]
        muts = re.findall(u'(?<=[0-9])([A-Za-z])+',trunc_mut)[0]
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

if __name__ == "__main__":
    writer = JSONWriter(HIVdb())
    assert (writer.isApobecDRM("IN", "G", 163, "TRAG")) == True