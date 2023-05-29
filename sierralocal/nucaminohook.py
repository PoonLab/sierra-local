import os
import csv
import sys
import json
import codecs
import platform
import tempfile
import subprocess
from pathlib import Path
from sierralocal.subtyper import Subtyper


class NucAminoAligner():
    """
    Initialize NucAmino for a specific input fasta file
    """

    def __init__(self, algorithm, binary=None, program='post'):
        """
        @param binary: str, Absolute path to nucamino binary
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
        self.triplet_table = self.generate_table()

        # with open(str(Path(os.path.dirname(__file__))/'data'/'apobec.tsv'), 'r') as csvfile:
        with open(algorithm.json_filename) as jsonfile:
            # self.apobec_drms = list(csv.reader(csvfile, delimiter='\t'))
            self.apobec_drms = json.load(jsonfile)

        self.pi_dict = self.prevalence_parser('PIPrevalences.tsv')
        self.rti_dict = self.prevalence_parser('RTIPrevalences.tsv')
        self.ini_dict = self.prevalence_parser('INIPrevalences.tsv')

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
        """
        Abstracted method for reading ARV prevalence TSV and 
        returning a dictionary of these data. There are two 
        entries for each subtype for treatment-naive and 
        experienced populations,
        respectively, e.g., B:RTI_Naive:% and B:RTI:%
        We take the larger of the two.
        @param filename: str, Name of TSV file to parse
        @return: dict, Dictionary of position-consensus-mutation-subtype 
                keys to %prevalence in naive populations
        """
        with open(str(Path(os.path.dirname(__file__)) / 'data' / filename), 'r') as handle:
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
        NucAmino does not return the aligned nucleotide sequence,
        but its JSON output provides sufficient information to 
        reconstitute this sequence.
        @param nuc: str, NucleicAcidsLine field from JSON record
        @param sites: list, AlignedSites field from JSON record
        @return: aligned: str, Aligned nt sequence
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
        """
        Using subprocess to call NucAmino, generates an output .tsv
        containing mutation data for each sequence in the FASTA file.
        Reconstitute aligned codon sequence from NucAmino output.
        For each codon in NucleicAcidsLine:
        - if LengthNA < 3, the codon has a deletion
        - if LengthNA == 3+n where n>0, the following n bases are insertions 
        to be removed
        @param filename: str, Path to FASTA file to process
        @return: list, list of records
        """

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
                        result['Sequence'] = self.get_aligned_seq(inputSequences[result['Name']],
                                                                  result['AlignedSites'])
                    output.append(result)
            # os.remove(tfPostOut.name)
            os.remove(tf.name)

            return output

        else:  # call nucAmino instead of post-align
            args = [
                '{}'.format(self.nucamino_binary),  # in case of byte-string
                'align',
                "hiv1b",
                'pol',
                "-q",
                "-i", tf.name,
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
        Returns a dictionary with the AMINO ACID position bounds
        for each gene in Pol, based on the HXB2 reference annotations.
        @return pol_aa_map: dict
        """
        # start and end nucleotide coordinates in HXB2 pol
        convert = lambda x: int((x - self.pol_start) / 3)
        pol_aa_map = {}
        for key, val in self.pol_nuc_map.items():
            pol_aa_map[key] = (convert(val[0]), convert(val[1]))
        return pol_aa_map

    def get_genes(self, pol_aligned_sites, pol_first_aa, pol_last_aa):
        """
        Determines the first POL gene that is present in 
        the query sequence, by virtue of gene breakpoints
        TODO: sierra uses different minimum numbers of sites per gene
        (40, 60 and 30 for PR, RT and IN)
        @param pol_aligned_sites: list, sublist holds alignment program
        output aligned POL sites
        @param pol_first_aa: int, location of first amino acid in pol
        @param pol_last_aa: int, location of last amino acid in pol
        @return: list, list of length 1
        [gene, first aa position in pol, last aa position in pol
         first na position in pol, last na position in pol]
        """
        # good here
        min_overlap = {'PR': 40, 'RT': 60, 'IN': 30}
        genes = []
        for gene, bounds in self.gene_map.items():
            aa_start, aa_end = bounds
            geneLength = aa_end - aa_start + 1
            try:
                overlap = min(aa_end, pol_last_aa) - max(aa_start, pol_first_aa)
                if overlap < min_overlap[gene]:
                    # discard alignment of this gene, too short
                    continue

                aligned_sites = filter(lambda x: x['PosAA'] >= aa_start and x['PosAA'] <= aa_end,
                                       pol_aligned_sites)
                aligned_sites = list(aligned_sites)

                first_aa = max(pol_first_aa - aa_start, 1)
                last_aa = min(pol_last_aa - aa_start + 1, geneLength)
                first_na = aligned_sites[0]['PosNA']
                last_na = aligned_sites[-1]['PosNA'] - 1 + aligned_sites[-1]['LengthNA']

                genes.append((gene, first_aa, last_aa, first_na, last_na))
            except:
                # if post-align returns nothing, need to raise message
                # "There were no PR, RT, and IN genes found, refuse to process."
                pass
        return genes

    def get_mutations(self, records, do_subtype=False):
        """
        From the tsv output of NucAmino, parses and adjusts indices
        and returns as lists.
        TSV has mutations format I59V:GTC,N93S:AGT
        JSON has mutations format:
        [{'AminoAcidText': 'V', 'InsertedCodonsText': '', 'IsDeletion': False, 'Control': '...',
          'ReferenceText': 'I', 'InsertedAminoAcidsText': '', 'IsPartial': False, 'NAPosition': 7,
          'IsInsertion': False, 'Position': 59, 'CodonText': 'GTC'},
         {'AminoAcidText': 'S', 'InsertedCodonsText': '', 'IsDeletion': False, 'Control': '...',
          'ReferenceText': 'N', 'InsertedAminoAcidsText': '', 'IsPartial': False, 'NAPosition': 109,
          'IsInsertion': False, 'Position': 93, 'CodonText': 'AGT'}]
        @param records: list, output of align_file function
        @param do_subtype: bool, ???
        @return: sequence_header, list of gene names
        @return: file_genes, output of get_genes() per each input sequence in input file
        @return: file_mutations, list, holds list of amino acid mutations based on position of POL
        [[consensus, aa, translated aa (text field)]]
        @return: file_trims, ???
        @return: subtypes, ???
        """
        file_mutations = []
        file_genes = []
        file_trims = []
        sequence_headers = []
        subtypes = []
        for record in records:
            sequence_headers.append(record['Name'])

            pol_first_aa = record['FirstAA']
            pol_last_aa = record['LastAA']

            # predict subtype
            subtype = ''
            if do_subtype:
                offset = (pol_first_aa - 57) * 3
                if offset < 0:
                    offset = 0  # align_file() will have trimmed sequence preceding PR
                subtype = self.typer.get_closest_subtype(record['Sequence'],
                                                         offset)
            genes = self.get_genes(record['AlignedSites'],
                                   pol_first_aa,
                                   pol_last_aa)
            trimmed_gene_muts = []
            trims = []
            first_last_nas = []
            just_genes = []

            for gene, first_aa, last_aa, first_na, last_na in genes:
                just_genes.append(gene)
                first_last_nas.append((first_na, last_na))

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
                        {position - left: (mut['ReferenceText'],  # consensus
                                           mut['AminoAcidText'],  # AAs line in output
                                           self.translate_na_triplet(codon)  # Text line in the output
                                           )}
                    )
                    codon_list.append(codon)

                # trim low quality leading and trailing nucleotides
                trim_left, trim_right = self.trim_low_qualities(
                    codon_list, left, first_aa, last_aa, mutations=gene_muts,
                    frameshifts=[], gene=gene, subtype=subtype
                )
                trims.append((trim_left, trim_right))
                trimmed_gene_muts.append(
                    {k: v for k, v in gene_muts.items() if
                     (k >= first_aa + trim_left) and
                     (k <= last_aa - trim_right)}
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

    def translate_na_triplet(self, triplet):
        """
        Translates a nucleotide triplet into its amino 
        acid mixture or ambiguity.
        @param triplet: str, nucleotide sequence as a string
        @return: str, translation of the triplet as a string
        """
        if len(triplet) == 0:
            return '-'
        if len(triplet) != 3:
            return "X"
        if '~' in triplet:
            return "X"
        return self.triplet_table.get(triplet, 'X')

    def generate_table(self):
        """
        Generates a dictionary of codon to amino acid mappings,
        including ambiguous combinations.
        @return triplet_table: dict, codon to amino acid dictionary
        """
        codon_to_aminoacid_map = {
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
        nas = ["A", "C", "G", "T", "R", "Y", "M", "W",
               "S", "K", "B", "D", "H", "V", "N"]
        triplet_table = dict()
        for i in range(len(nas)):
            for j in range(len(nas)):
                for k in range(len(nas)):
                    triplet = nas[i] + nas[j] + nas[k]
                    codons = self.enumerate_codon_possibilities(triplet)
                    unique_aas = []
                    for codon in codons:
                        c = codon_to_aminoacid_map[codon]
                        if c not in unique_aas:
                            unique_aas.append(c)
                    if len(unique_aas) > 4:
                        aas = "X"
                    else:
                        aas = ''.join(unique_aas)
                    triplet_table[triplet] = aas
        return triplet_table

    def enumerate_codon_possibilities(self, triplet):
        """
        Converts a potentially ambiguous nucleotide triplet into
        standard ATCG codons.
        @param triplet: str, nucleotide triplet as a string
        @return codon_possibilities: list, list of possible ATCG codons
        encoded by the triplet
        """
        ambiguity_map = {
            "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
            "R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"],
            "W": ["A", "T"], "S": ["C", "G"], "K": ["G", "T"],
            "B": ["C", "G", "T"], "D": ["A", "G", "T"],
            "H": ["A", "C", "T"], "V": ["A", "C", "G"],
            "N": ["A", "C", "G", "T"]
        }
        codon_possibilities = []
        pos1, pos2, pos3 = triplet
        for p1 in ambiguity_map[pos1]:
            for p2 in ambiguity_map[pos2]:
                for p3 in ambiguity_map[pos3]:
                    codon_possibilities.append(p1 + p2 + p3)
        return codon_possibilities

    def trim_low_qualities(self, codon_list, shift, first_aa, last_aa, mutations, frameshifts, gene, subtype):
        """
        Filters low-quality leading and trailing nucleotides from a query.
        Removes large (length > SEQUENCE_TRIM_SITES_CUTOFF) low quality pieces.
        Low quality is defined as:
        (1) unusual mutation; or
        (2) 'X' in amino acid list; or
        (3) has a stop codon

        @param codon_list: list, ???
        @param shift: int, ???
        @param first_aa: int, aligned position of first amino acid in query
        @param last_aa: int, aligned position of last amino acid in query
        @param mutations: dict, dictionary of <position>: (<wt>, <mutant>) pairs
        @param frameshifts: list, list of frameshifts
        @param gene: str, list of frameshifts
        @param subtype: str, ???
        @return: tuple of how many leading and trailing nucleotides to trim
        """
        SEQUENCE_SHRINKAGE_CUTOFF_PCNT = 30
        SEQUENCE_SHRINKAGE_WINDOW = 15
        SEQUENCE_SHRINKAGE_BAD_QUALITY_MUT_PREVALENCE = 0.1

        bad_pcnt = 0
        problem_sites = 0
        since_last_bad_quality = 0
        protien_size = last_aa - first_aa + 1

        candidates = []
        invalid_sites = [False for i in range(protien_size)]

        # account for invalid sites
        for j, position in enumerate(mutations):
            idx = position - first_aa  # + shift
            if not self.is_unsequenced(codon_list[j]):
                highest_prev = self.get_highest_mut_prevalance(
                    (position, mutations[position]),
                    gene,
                    subtype
                )
                is_weird = (highest_prev < SEQUENCE_SHRINKAGE_BAD_QUALITY_MUT_PREVALENCE)
                is_ambig = (mutations[position][1] == 'X')
                is_apobec = self.is_apobec_drm(gene,
                                               mutations[position][0],
                                               position,
                                               mutations[position][1])
                is_stop_codon = self.is_stop_codon(codon_list[j])
                reasons = [is_weird, is_ambig, is_apobec, is_stop_codon]
                if any(reasons):
                    # print(idx, position, reasons, highest_prev, subtype, position, mutations[position], gene)
                    invalid_sites[idx] = True

        # for fs in frameshifts:
        #     idx = fs.getPosition() - firstAA
        #     invalidSites[idx] = True

        # forward scan for trimming left
        for idx in range(0, protien_size):
            if since_last_bad_quality > SEQUENCE_SHRINKAGE_WINDOW:
                break

            if invalid_sites[idx]:
                problem_sites += 1
                trim_left = idx + 1
                bad_pcnt = problem_sites * 100 / trim_left if trim_left > 0 else 0
                if bad_pcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT:
                    candidates.append(trim_left)
                since_last_bad_quality = 0
            else:
                since_last_bad_quality += 1

        trim_left = candidates[-1] if len(candidates) > 0 else 0
        candidates = []

        # backward scan for trimming right
        problem_sites = 0
        since_last_bad_quality = 0
        for idx in range(protien_size - 1, -1, -1):
            if since_last_bad_quality > SEQUENCE_SHRINKAGE_WINDOW:
                break
            elif invalid_sites[idx]:
                problem_sites += 1
                trim_right = protien_size - idx
                bad_pcnt = problem_sites * 100 / trim_right if trim_right > 0 else 0
                if bad_pcnt > SEQUENCE_SHRINKAGE_CUTOFF_PCNT:
                    candidates.append(trim_right)
                since_last_bad_quality = 0
            else:
                since_last_bad_quality += 1
        trim_right = candidates[-1] if len(candidates) > 0 else 0

        return (trim_left, trim_right)

    def is_unsequenced(self, triplet):
        """
        Determines whether a triplet is unsequenced (has more than one N or
        deletion). "NNN", "NN-", "NNG" should be considered as unsequenced region.
        """
        return (triplet.replace("-", "N").count("N") > 1)  # TODO: incorporate !isInsertion &&

    def is_stop_codon(self, triplet):
        return ("*" in self.translate_na_triplet(triplet))

    def is_apobec_drm(self, gene, consensus, position, AA):  # pragma: no cover
        ls = [[row['gene'], str(row['position'])] for row in self.apobec_drms]
        if [gene, str(position)] in ls:
            i = ls.index([gene, str(position)])
            for aa in AA:
                if aa in self.apobec_drms[i]['aa']:
                    return True
        return False

    def get_highest_mut_prevalance(self, mutation, gene, subtype):
        """
        #TODO
        @param mutation: tuple, a tuple representing a specific 
        position in the amino acid sequence that may contain 
        multiple amino acids (polymorphic)
        @param gene: str, PR, RT, or INT
        @param subtype: str, predicted from Subtyper.get_closest_subtype()
        @return: float, prevalence of the most common amino acid encoded 
        at this position within the subtype alignment
        """
        position, aa_seq = mutation
        cons, aas, text = aa_seq
        aas = aas.replace(cons, '')  # ignore consensus
        aas = aas.replace('*', '')  # remove stop codons

        prevalence = 0.
        for aa in aas:
            aa_prevalence = self.get_mut_prevalence(position, cons, aa, gene, subtype)
            prevalence = max(prevalence, aa_prevalence)

        return prevalence

    def get_mut_prevalence(self, position, cons, aa, gene, subtype):
        """
        ???
        """
        key2 = str(position) + str(cons) + str(aa) + subtype

        if gene == 'IN' and key2 in self.ini_dict:
            return self.ini_dict[key2]

        if gene == 'PR' and key2 in self.pi_dict:
            return self.pi_dict[key2]

        if gene == 'RT' and key2 in self.rti_dict:
            return self.rti_dict[key2]

        return 100.0


if __name__ == '__main__':
    from sierralocal.hivdb import HIVdb

    algorithm = HIVdb()
    test = NucAminoAligner(algorithm)
    """
    assert test.translate_na_triplet("YTD") == "LF"
    assert test.is_stop_codon("TAG") == True
    assert test.is_stop_codon("TAA") == True
    assert test.is_stop_codon("TGA") == True
    assert test.is_stop_codon("NNN") == False

    assert test.is_unsequenced("NNN") == True
    assert test.is_unsequenced("NN-") == True
    assert test.is_unsequenced("NNG") == True
    assert test.is_unsequenced("NTG") == False

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
