import subprocess
import os
from pathlib import Path
import csv
import re
from sierralocal.subtyper import Subtyper
import sys
import platform
from csv import DictReader
import json
import codecs
import tempfile


class NucAminoAligner():
    """
    Initialize NucAmino for a specific input fasta file
    """

    def __init__(self, algorithm, binary=None, program='post'):
        """

        :param binary:  Absolute path to nucamino binary
        """
        self.cwd = os.path.curdir
        self.reader = codecs.getreader('utf-8')  # for parsing byte strings from json in align_file()

        if program == 'post':
            print('Aligning using post-align')
            pass  # TODO include post align as a submodule
        else:  # get necessary binaries for nucAmino
            if binary is None:
                target = 'nucamino-{}-{}'.format(
                    platform.system().lower(),
                    'amd64' if platform.architecture()[0] == '64bit' else '386'
                )

                # autodetect nucamino binary
                bin_dir = Path(os.path.dirname(__file__)).joinpath('bin')
                items = os.listdir(str(bin_dir))
                self.nucamino_binary = None
                for name in items:
                    fn = os.path.splitext(name)[0]
                    if fn == target:
                        self.nucamino_binary = bin_dir.joinpath(name)
                        break

                if self.nucamino_binary is None:
                    sys.exit('Failed to locate expected NucAmino binary {}. '.format(target) +
                             'Please download binary from ' +
                             'http://github.com/hivdb/nucamino/releases')

                print("Found NucAmino binary", self.nucamino_binary)
            else:
                self.nucamino_binary = binary

        self.tripletTable = self.generateTable()

        # with open(str(Path(os.path.dirname(__file__))/'data'/'apobec.tsv'), 'r') as csvfile:
        with open(algorithm.tsv_filename) as csvfile:
            self.ApobecDRMs = list(csv.reader(csvfile, delimiter='\t'))

        self.PI_dict = self.prevalence_parser('PIPrevalences.tsv')
        self.RTI_dict = self.prevalence_parser('RTIPrevalences.tsv')
        self.INI_dict = self.prevalence_parser('INIPrevalences.tsv')

        # initialize gene map
        self.pol_start = 2085
        self.pol_nuc_map = {
            'PR': (2253, 2549),
            'RT': (2550, 4229),  # incorrectly includes RNAse, emulating sierrapy
            'IN': (4230, 5096)
        }
        self.gene_map = self.create_gene_map()

        # initialize Subtyper class
        self.typer = Subtyper()

    def prevalence_parser(self, filename):
        '''
        Abstracted method for reading ARV prevalence TSV and returning a dictionary of these data.
        There are two entries for each subtype for treatment-naive and -experienced populations,
        respectively, e.g., B:RTI_Naive:% and B:RTI:%
        We take the larger of the two.

        :param filename:  Name of TSV file to parse
        :return:  Dictionary of position-consensus-mutation-subtype keys to %prevalence in naive
                  populations
        '''
        handle = open(str(Path(os.path.dirname(__file__)) / 'data' / filename), 'r')
        table = csv.DictReader(handle, delimiter='\t')

        keys = [fn for fn in table.fieldnames if fn.endswith(':%')]
        result = {}
        for row in table:
            for key in keys:
                value = float(row[key])
                subtype = key.split(':')[0]
                label = row['Pos'] + row['Cons'] + row['Mut'] + subtype
                if label in result and result[label] < value:
                    # take the greater percentage of naive or experienced groups
                    result[label] = value
                else:
                    result.update({label: value})
        return result

    def get_aligned_seq(self, nuc, sites):
        """
        NucAmino does not return the aligned nucleotide sequence, but its JSON output provides
        sufficient information to reconstitute this sequence.
        :param nuc:  NucleicAcidsLine field from JSON record
        :param sites:  AlignedSites field from JSON record
        :return:  Aligned nt sequence
        """

        # FIXME: this isn't handling insertions and deletions properly
        # FIXME: the NucAmino coordinate system does not adapt to indels
        aligned = ''
        skip = 0
        pad = 0  # used to accommodate deletions upstream of PR start codon
        for si, site in enumerate(sites):
            if skip > 0:
                skip -= 1
                continue

            codon = nuc[3 * si:(3 * si + 3)]

            lengthNA = site['LengthNA']
            if lengthNA < 3:
                pad += 3 - lengthNA

            posAA = site['PosAA'] - pad // 3
            if posAA < 57:
                # codon is upstream of 5'-PR
                continue

            if lengthNA == 3:
                aligned += codon
            elif lengthNA < 3:
                # deletion
                aligned += codon.replace(' ', '-')
            else:
                # insertion
                skip = (lengthNA - 3) / 3
                aligned += codon

        return aligned

    def makeReferenceFASTA(self, fragmentName, refSeq):
        tempFasta = tempfile.NamedTemporaryFile('w', prefix='postalign-ref-', suffix='.fas', delete=False)
        tempFasta.write(">Ref_{}\n".format(fragmentName))
        tempFasta.write("{}\n".format(refSeq))
        tempFasta.close()
        return os.path.abspath(tempFasta.name)

    def getConfigField(self, config, field):
        resultmap = {}
        for entry in config['fragmentConfig']:
            if field == 'refSequence':
                if field in entry:
                    resultmap.update(
                        {entry['fragmentName']: self.makeReferenceFASTA(entry['fragmentName'], entry['refSequence'])})
            else:
                if field in entry:
                    if entry['fragmentName'] not in resultmap:
                        resultmap.update({entry['fragmentName']: [entry[field]]})
                    else:
                        resultmap[entry['fragmentName']].append(entry[field])

        return resultmap

    def align_file(self, filename, program='post'):
        '''
        Using subprocess to call NucAmino or postalign, generates an output .tsv containing mutation
        data for each sequence in the FASTA file.
        Reconstitute aligned codon sequence from NucAmino output.
        For each codon in NucleicAcidsLine:
        - if LengthNA < 3, the codon has a deletion
        - if LengthNA == 3+n where n>0, the following n bases are insertions to be removed

        @param filename:  Path to FASTA file to process
        program: 'post' indicates running postalign, anything else goes to nucamino
        '''
        # TODO: check that file is FASTA format

        # postalign doesn't have the aligned sequence in its output or the input sequence
        # making a dictionary to keep track of the needed sequence so we can align it like we did with nucAmino
        inputSequences = {}
        # remove illegal characters
        tf = tempfile.NamedTemporaryFile(mode='w', delete=False)
        with open(filename) as handle:
            for line in handle:
                if not line.startswith('>'):
                    line = line.replace('~', '').replace('-', '').replace('.', '')
                    inputSequences[name] = line
                else:
                    name = line[1:].strip()
                    inputSequences.update({name: ''})
                tf.write(line)
        tf.close()

        if program == 'post':
            # incorporate config file as apart of sierralocal
            with open(str(Path(os.path.dirname(__file__)) / 'data' / 'alignment-config_hiv1.json'), 'r') as f:
                config = json.load(f)

            MIN_MATCH_PCNT = self.getConfigField(config=config, field='minMatchPcnt')
            MIN_NUM_OF_AA = self.getConfigField(config=config, field='minNumOfAA')
            FROM_FRAGMENT = self.getConfigField(config=config, field='fromFragment')
            GENE = self.getConfigField(config=config, field='geneName')
            REF_RANGES = self.getConfigField(config=config, field='refRanges')
            REF_SEQUENCE = self.getConfigField(config=config, field='refSequence')
            POST_PROCESSORS = self.getConfigField(config=config, field='postProcessors')
            MINIMAP2_OPTS = self.getConfigField(config=config, field='minimap2Opts')

            # hold the output of postalign
            tfPostOut = tempfile.NamedTemporaryFile(mode='w', delete=False)
            for refFragmentName in REF_SEQUENCE:
                refSeqFile = REF_SEQUENCE[refFragmentName]
                cmd = [
                    'postalign',
                    '-i', tf.name,
                    '-o', tfPostOut.name,
                    '-f', 'MINIMAP2',
                    '-r', refSeqFile
                ]
                if refFragmentName in MINIMAP2_OPTS:
                    cmd.append('--minimap2-opts')
                    cmd.append(MINIMAP2_OPTS[refFragmentName][0])

                for op in POST_PROCESSORS[refFragmentName][0]:
                    cmd.append(op)

                cmd.append('save-json')
                for fragmentName in FROM_FRAGMENT.keys():
                    if FROM_FRAGMENT[fragmentName][0] == refFragmentName:
                        cmd.append(fragmentName)
                        for range in REF_RANGES[fragmentName][0]:
                            cmd.append(str(range[0]))
                            cmd.append(str(range[1]))
                _ = subprocess.check_call(cmd)
                for f in REF_SEQUENCE:
                    os.remove(REF_SEQUENCE[f])

            output = []
            with open(tfPostOut.name, 'r') as postOut:
                postOut = json.load(postOut)
                for sequence in postOut:
                    result = {'AlignedSites': [],
                              'FirstAA': None,
                              'FirstNA': None,
                              'FrameShifts': [],
                              'LastAA': None,
                              'LastNA': None,
                              'Mutations': [],
                              'Name': '',
                              'Sequence': ''
                              }

                    for field, data in sequence.items():
                        # postalign has different AlignedSites output
                        # missing fields lastNA and firstNA

                        if field == 'GeneReports':
                            # this field holds all the fields of that of nucAMINO
                            for protein in data:

                                # sequenced if the report field is not empty
                                # NucAmino is positioned on POL
                                if not protein['Error'] and protein['Report'] and 'pol' in protein['Gene']:
                                    for key, info in protein['Report'].items():
                                        if key == 'AlignedSites':
                                            for i in info:
                                                # PosNA shows all three positions, should only show first
                                                i['PosNA'] = i['PosNAs'][0]
                                                i.pop('PosNAs')
                                                i['PosAA'] = i['PosAA'] + 1

                                            result['AlignedSites'] += info

                                        elif key == 'Mutations':
                                            for mutation in info:
                                                mutation['ReferenceText'] = mutation['RefAminoAcidText']
                                                mutation.pop('RefAminoAcidText')
                                                mutation['Position'] += 1
                                            result['Mutations'] += info

                                        else:
                                            result.update({key: info})
                                        # print(result['AlignedSites'])

                        else:
                            if field in result:
                                result.update({field: data})
                    # manually add in the last and first NA based on aligned sequences
                    # in the case it finds no genes, it will report no AlignedSites
                    if result['AlignedSites']:
                        result['AlignedSites'] = sorted(result['AlignedSites'], key=lambda x: x['PosAA'])
                        result['FirstNA'] = result['AlignedSites'][0]['PosNA']
                        result['LastNA'] = result['AlignedSites'][-1]['PosNA']
                        result['FirstAA'] = result['AlignedSites'][0]['PosAA']
                        result['LastAA'] = result['AlignedSites'][-1]['PosAA']
                        result['Sequence'] = self.get_aligned_seq(inputSequences[result['Name']], result['AlignedSites'])
                    output.append(result)
            # os.remove(tfPostOut.name)
            os.remove(tf.name)

            with open('org2.json', 'w') as out:
                json.dump(result, out, indent=4, sort_keys=True)

            return output

        else:  # call nucAmino instead of post-align
            args = [
                '{}'.format(self.nucamino_binary),  # in case of byte-string
                "hiv1b",
                "-q",
                "-i", tf.name,
                "-g=POL",
                '--output-format', 'json',
            ]
            p = subprocess.Popen(args, stdout=subprocess.PIPE)  # , encoding='utf8')

            result = json.load(self.reader(p.stdout))
            records = []

            for record in result['POL']:
                try:
                    sites = record['Report']['AlignedSites']
                except:
                    print(record)
                    raise
                nuc = record['Report']['NucleicAcidsLine']

                records.append({
                    'Name': record['Name'],
                    'FirstAA': record['Report']['FirstAA'],  # relative to start of pol
                    'LastAA': record['Report']['LastAA'],
                    'FirstNA': record['Report']['FirstNA'],
                    'LastNA': record['Report']['LastNA'],
                    'Mutations': record['Report']['Mutations'],
                    'Frameshifts': record['Report']['FrameShifts'],
                    'AlignedSites': sites,
                    'Sequence': self.get_aligned_seq(nuc, sites)
                })
            os.remove(tf.name)
            return records

    def create_gene_map(self):
        """
        Returns a dictionary with the AMINO ACID position bounds for each gene in Pol,
        based on the HXB2 reference annotations.
        """
        # start and end nucleotide coordinates in HXB2 pol
        convert = lambda x: int((x - self.pol_start) / 3)
        pol_aa_map = {}
        for key, val in self.pol_nuc_map.items():
            pol_aa_map[key] = (convert(val[0]), convert(val[1]))
        return pol_aa_map

    def get_genes(self, polAlignedSites, polFirstAA, polLastAA):
        """
        Determines the first POL gene that is present in the query sequence, by virtue of gene breakpoints
        TODO: sierra uses different minimum numbers of sites per gene (40, 60 and 30 for PR, RT and IN)
        @param firstAA: the first amino acid position of the query sequence, as determined by NucAmino
        @return: list of length 1 with a string of the gene in the sequence
        """
        # good here
        min_overlap = {'PR': 40, 'RT': 60, 'IN': 30}
        genes = []
        for gene, bounds in self.gene_map.items():
            aaStart, aaEnd = bounds
            geneLength = aaEnd - aaStart + 1
            try:
                overlap = min(aaEnd, polLastAA) - max(aaStart, polFirstAA)
                if overlap < min_overlap[gene]:
                    # discard alignment of this gene, too short
                    continue
                alignedSites = filter(lambda x: x['PosAA'] >= aaStart and x['PosAA'] <= aaEnd, polAlignedSites)
                alignedSites = list(alignedSites)

                firstAA = max(polFirstAA - aaStart, 1)
                lastAA = min(polLastAA - aaStart + 1, geneLength)
                firstNA = alignedSites[0]['PosNA']
                lastNA = alignedSites[-1]['PosNA'] - 1 + alignedSites[-1]['LengthNA']

                genes.append((gene, firstAA, lastAA, firstNA, lastNA))
            except:
                # if post-align returns nothing, need to raise the message
                # "There were no PR, RT, and IN genes found, refuse to process."
                pass

        return genes

    def get_mutations(self, records, do_subtype=False):
        '''
        From the tsv output of NucAmino, parses and adjusts indices and returns as lists.

        TSV has mutations format I59V:GTC,N93S:AGT

        JSON has mutations format (NUCAMINO):
        [{'AminoAcidText': 'V', 'InsertedCodonsText': '', 'IsDeletion': False, 'Control': '...',
         'ReferenceText': 'I', 'InsertedAminoAcidsText': '', 'IsPartial': False, 'NAPosition': 7,
         'IsInsertion': False, 'Position': 59, 'CodonText': 'GTC'},
         {'AminoAcidText': 'S', 'InsertedCodonsText': '', 'IsDeletion': False, 'Control': '...',
         'ReferenceText': 'N', 'InsertedAminoAcidsText': '', 'IsPartial': False, 'NAPosition': 109,
         'IsInsertion': False, 'Position': 93, 'CodonText': 'AGT'}]

        :param fastaFileName:  FASTA input processed by NucAmino
        :return: list of sequence names, list of sequence mutation dictionaries.
        '''
        file_mutations = []
        file_genes = []
        file_trims = []
        sequence_headers = []
        subtypes = []
        for record in records:
            sequence_headers.append(record['Name'])

            polFirstAA = record['FirstAA']
            polLastAA = record['LastAA']

            # predict subtype
            subtype = ''
            if do_subtype:
                offset = (polFirstAA - 57) * 3
                if offset < 0:
                    offset = 0  # align_file() will have trimmed sequence preceding PR
                subtype = self.typer.getClosestSubtype(record['Sequence'], offset)

            ### good at genes
            genes = self.get_genes(record['AlignedSites'], polFirstAA, polLastAA)
            trimmed_gene_muts = []
            trims = []
            first_lastNAs = []
            just_genes = []

            for gene, firstAA, lastAA, firstNA, lastNA in genes:
                just_genes.append(gene)
                first_lastNAs.append((firstNA, lastNA))

                # {'PR': (56, 154), 'RT': (155, 714), 'IN': (715, 1003)}
                left, right = self.gene_map[gene]
                codon_list = []
                gene_muts = {}
                for mut in record['Mutations']:
                    position = mut['Position']
                    if position < left or right < position:
                        continue  # mutation in other gene
                    codon = mut['CodonText']
                    gene_muts.update(
                        {position - left: (mut['ReferenceText'], self.translateNATriplet(codon))}
                    )
                    codon_list.append(codon)

                # trim low quality leading and trailing nucleotides
                trimLeft, trimRight = self.trimLowQualities(
                    codon_list, left, firstAA, lastAA, mutations=gene_muts,
                    frameshifts=[], gene=gene, subtype=subtype
                )
                trims.append((trimLeft, trimRight))
                trimmed_gene_muts.append(
                    {k: v for k, v in gene_muts.items() if
                     (k >= firstAA + trimLeft) and
                     (k <= lastAA - trimRight)}
                )

            # update lists
            file_mutations.append(trimmed_gene_muts)
            file_genes.append(genes)
            file_trims.append(trims)
            subtypes.append(subtype)

        assert len(file_mutations) == len(sequence_headers), \
            "error: length of mutations dicts is not the same as length of names"
        return sequence_headers, file_genes, file_mutations, file_trims, subtypes

    # BELOW is an implementation of sierra's Java algorithm for determining codon ambiguity

    def translateNATriplet(self, triplet):
        """
        Translates a nucleotide triplet into its amino acid mixture or ambiguity.
        @param triplet: nucleotide sequence as a string
        @return: translation of the triplet as a string
        """
        if len(triplet) == 0:
            return '-'
        if len(triplet) != 3:
            return "X"
        if '~' in triplet:
            return "X"
        return self.tripletTable.get(triplet, 'X')

    def generateTable(self):
        """
        Generates a dictionary of codon to amino acid mappings, including ambiguous combinations.
        @return tripletTable: codon to amino acid dictionary
        """
        codonToAminoAcidMap = {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        }
        nas = ["A", "C", "G", "T", "R", "Y", "M", "W", "S", "K", "B", "D", "H", "V", "N"]
        tripletTable = dict()
        for i in range(len(nas)):
            for j in range(len(nas)):
                for k in range(len(nas)):
                    triplet = nas[i] + nas[j] + nas[k]
                    codons = self.enumerateCodonPossibilities(triplet)
                    uniqueAAs = []
                    for codon in codons:
                        c = codonToAminoAcidMap[codon]
                        if c not in uniqueAAs:
                            uniqueAAs.append(c)
                    if len(uniqueAAs) > 4:
                        aas = "X"
                    else:
                        aas = ''.join(uniqueAAs)
                    tripletTable[triplet] = aas
        return tripletTable

    def enumerateCodonPossibilities(self, triplet):
        """
        Converts a potentially ambiguous nucleotide triplet into standard ATCG codons.
        @param triplet: nucleotide triplet as a string
        @return codonPossibilities: list of possible ATCG codons encoded by the triplet
        """
        ambiguityMap = {
            "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
            "R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"],
            "W": ["A", "T"], "S": ["C", "G"], "K": ["G", "T"],
            "B": ["C", "G", "T"], "D": ["A", "G", "T"],
            "H": ["A", "C", "T"], "V": ["A", "C", "G"],
            "N": ["A", "C", "G", "T"]
        }
        codonPossibilities = []
        pos1, pos2, pos3 = triplet
        for p1 in ambiguityMap[pos1]:
            for p2 in ambiguityMap[pos2]:
                for p3 in ambiguityMap[pos3]:
                    codonPossibilities.append(p1 + p2 + p3)
        return codonPossibilities

    def trimLowQualities(self, codon_list, shift, firstAA, lastAA, mutations, frameshifts, gene, subtype):
        """
        Filters low-quality leading and trailing nucleotides from a query.
        Removes large (length > SEQUENCE_TRIM_SITES_CUTOFF) low quality pieces.
        Low quality is defined as:
        (1) unusual mutation; or
        (2) 'X' in amino acid list; or
        (3) has a stop codon

        @param sequence: sequence
        @param firstAA: aligned position of first amino acid in query
        @param lastAA: aligned position of last amino acid in query
        @param mutations: dictionary of <position>: (<wt>, <mutant>) pairs
        @param frameshifts: list of frameshifts
        @return: tuple of how many leading and trailing nucleotides to trim
        """
        SEQUENCE_SHRINKAGE_CUTOFF_PCNT = 30
        SEQUENCE_SHRINKAGE_WINDOW = 15
        SEQUENCE_SHRINKAGE_BAD_QUALITY_MUT_PREVALENCE = 0.1

        # print(mutations)
        # print(codon_list)

        badPcnt = 0
        problemSites = 0
        sinceLastBadQuality = 0
        proteinSize = lastAA - firstAA + 1

        candidates = []
        invalidSites = [False for i in range(proteinSize)]

        # account for invalid sites
        for j, position in enumerate(mutations):
            idx = position - firstAA  # + shift
            if not self.isUnsequenced(codon_list[j]):
                highest_prev = self.getHighestMutPrevalence((position, mutations[position]), gene, subtype)
                is_weird = (highest_prev < SEQUENCE_SHRINKAGE_BAD_QUALITY_MUT_PREVALENCE)
                is_ambig = (mutations[position][1] == 'X')
                is_apobec = self.isApobecDRM(gene, mutations[position][0], position,
                                             mutations[position][1])
                is_stop_codon = self.isStopCodon(codon_list[j])
                reasons = [is_weird, is_ambig, is_apobec, is_stop_codon]
                if any(reasons):
                    # print(idx, position, reasons, highest_prev, subtype, position, mutations[position], gene)
                    invalidSites[idx] = True

        # for fs in frameshifts:
        #     idx = fs.getPosition() - firstAA
        #     invalidSites[idx] = True

        # forward scan for trimming left
        for idx in range(0, proteinSize):
            if sinceLastBadQuality > SEQUENCE_SHRINKAGE_WINDOW:
                break

            if invalidSites[idx]:
                problemSites += 1
                trimLeft = idx + 1
                badPcnt = problemSites * 100 / trimLeft if trimLeft > 0 else 0
                if badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT:
                    candidates.append(trimLeft)
                sinceLastBadQuality = 0
            else:
                sinceLastBadQuality += 1

        trimLeft = candidates[-1] if len(candidates) > 0 else 0
        candidates = []

        # backward scan for trimming right
        problemSites = 0
        sinceLastBadQuality = 0
        for idx in range(proteinSize - 1, -1, -1):
            if sinceLastBadQuality > SEQUENCE_SHRINKAGE_WINDOW:
                break
            elif invalidSites[idx]:
                problemSites += 1
                trimRight = proteinSize - idx
                badPcnt = problemSites * 100 / trimRight if trimRight > 0 else 0
                if badPcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT:
                    candidates.append(trimRight)
                sinceLastBadQuality = 0
            else:
                sinceLastBadQuality += 1
        trimRight = candidates[-1] if len(candidates) > 0 else 0

        return (trimLeft, trimRight)

    """
    Determines whether a triplet is unsequenced (has more than one N or deletion).
    "NNN", "NN-", "NNG" should be considered as unsequenced region.
    """

    def isUnsequenced(self, triplet):
        return (triplet.replace("-", "N").count("N") > 1)  # TODO: incorporate !isInsertion &&

    def isStopCodon(self, triplet):
        return ("*" in self.translateNATriplet(triplet))

    def isApobecDRM(self, gene, consensus, position, AA):
        ls = [row[0:3] for row in self.ApobecDRMs[1:]]
        if [gene, consensus, str(position)] in ls:
            i = ls.index([gene, consensus, str(position)])
            for aa in AA:
                if aa in self.ApobecDRMs[1:][i][3]:
                    return True
        return False

    def getHighestMutPrevalence(self, mutation, gene, subtype):
        """
        #TODO
        @param mutation: a tuple(?) representing a specific position in the amino acid sequence 
                         that may contain multiple amino acids (polymorphic)
        @param gene: PR, RT, or INT
        @param subtype: predicted from function()
        @return: prevalence of the most common amino acid encoded at this position within the 
                 subtype alignment
        """
        position, aaseq = mutation
        cons, aas = aaseq
        aas = aas.replace(cons, '')  # ignore consensus
        aas = aas.replace('*', '')  # remove stop codons

        prevalence = 0.
        for aa in aas:
            aaPrevalence = self.getMutPrevalence(position, cons, aa, gene, subtype)
            prevalence = max(prevalence, aaPrevalence)

        return prevalence

    def getMutPrevalence(self, position, cons, aa, gene, subtype):
        key2 = str(position) + str(cons) + str(aa) + subtype
        # print(key2)

        if gene == 'IN' and key2 in self.INI_dict:
            return self.INI_dict[key2]

        if gene == 'PR' and key2 in self.PI_dict:
            return self.PI_dict[key2]

        if gene == 'RT' and key2 in self.RTI_dict:
            return self.RTI_dict[key2]

        return 100.0


if __name__ == '__main__':
    from sierralocal.hivdb import HIVdb

    algorithm = HIVdb()
    test = NucAminoAligner(algorithm)
    """
    assert test.translateNATriplet("YTD") == "LF"
    assert test.isStopCodon("TAG") == True
    assert test.isStopCodon("TAA") == True
    assert test.isStopCodon("TGA") == True
    assert test.isStopCodon("NNN") == False

    assert test.isUnsequenced("NNN") == True
    assert test.isUnsequenced("NN-") == True
    assert test.isUnsequenced("NNG") == True
    assert test.isUnsequenced("NTG") == False

    print(test.getMutPrevalence(6, 'D', 'E', 'IN', "CRF01_AE"))
    print(test.getMutPrevalence(6, 'D', 'E', 'IN', "G"))
    print(test.getMutPrevalence(6, 'E', 'D', 'RT', "A"))

    print(test.gene_map)
    print(test.get_gene(594))
    """
    fn = sys.argv[1]
    results = test.align_file(fn)
    res = test.get_mutations(results)
    print(res)
