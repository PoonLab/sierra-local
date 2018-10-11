from pathlib import Path
import csv
import os
#import re
import sys
from sierralocal.utils import get_input_sequences
#from sierralocal.nucaminohook import NucAminoAligner
from time import time

class Subtyper():
    def __init__(self):
        self.ambiguities = {
            "A": "A", "C": "C", "G": "G", "T": "T",
            "R": "AG", "Y": "CT", "M": "AC", "W": "AT", "S": "CG", "K": "GT",
            "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
            "N": "ACGT", '~': '', '.': '', '-': 'ACGT'
        }
        self.offset = {'PR': 0, 'RT': 297, 'IN': 1977}
        self.simple_subtypes = {}
        self.subtype_references = self.getSubtypeReferences()
        self.properties = self.getGenotypeProperties()


    def getSubtypeReferences(self):
        """
        Parse FASTA file containing genotype reference sequences.
        Also extract simple subtypes from reference labels.
        :return: a dictionary of label:sequence pairs
        """
        # FIXME: this is a hard-coded data filename
        filepath = Path(os.path.dirname(__file__))/'data'/'genotype-references.9c610d61.fasta'
        handle = open(str(filepath))
        result = get_input_sequences(handle, return_dict=True)
        for label, seq in result.items():
            subtype = label.split('|')[1]
            if subtype not in self.simple_subtypes:
                self.simple_subtypes.update({subtype: []})
            self.simple_subtypes[subtype].append(label)
        return result

    def uncorrectedDistance(self, seq, ref):
        """
        p-distance for nucleotide sequences.
        Assuming that sequences are aligned, count the number of base
        differences while accommodating ambiguous base calls.

        :param seq:  first input sequence
        :param ref:  second input sequences
        :return:  proportion of base differences
        """

        # TODO: Step 1: mask SDRM in seq with ref nucleotide
        count = 0
        for idx, nuc in enumerate(seq):
            # print(idx, nuc)
            rnuc = ref[idx]
            if nuc == '-' or nuc == rnuc:
                continue
            if not ref[idx] in self.ambiguities[nuc]:
                count += 1
        return count / len(seq.strip('-'))

    def getDistances(self, sequence, offset=0):
        """
        Uncorrected pairwise distances: the number of identical bases divided by
        the sequence length, ignoring positions with gaps.

        :param sequence:  query nucleotide sequence
        :param offset:  in the case of RT or IN, query sequence is missing
                        upstream nucleotides so we shift the reference
        :return: a dictionary of (sequence label, distance) key-value pairs
        """
        dists = {}
        #TODO: generate concatenated sequences, with masked SDRMS
        for header, reference in self.subtype_references.items():
            dists[header] = self.uncorrectedDistance(sequence, reference[offset:])
        return dists


    def getGenotypeProperties(self):
        """
        Parse CSV file that summarizes the subtypes and CRFs in the HIVdb database
        :return:  A nested dictionary keyed by subtype/CRF name
        """
        # FIXME: this is a hard-coded path
        filepath = str(Path(os.path.dirname(__file__))/'data'/ 'genotype-properties.cc00f512.csv')
        props = {}
        with open(filepath) as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                rname = row['name']
                props.update({rname: {
                    'parent': row['parentSubtype'],
                    'dist_upper': float(row['distanceUpperLimit'].strip('%'))/100,
                    'is_simple_CRF': bool(row['isSimpleCRF']),
                    'class_level': row['classificationLevel'],
                    'breakpoints': row['breakPoints']
                }})
        return props


    def getClosestSubtype(self, gene, sequence, offset):
        dists = self.getDistances(sequence, offset)

        intermed = [(v, k) for k, v in dists.items()]
        intermed.sort(reverse=False)  # increasing order
        closestDist, closestMatch = intermed[0]

        # sort dist dict
        #sorted_dists = dict([(k, dists[k]) for k in sorted(dists, key=dists.get, reverse=False)])
        #closestMatch = next(iter(sorted_dists.keys()))
        #closestDist = sorted_dists[closestMatch]
        #closestSubtype = re.findall('\|(.*?)\||$', closestMatch)[0]

        # FIXME: where do we get this number from?
        if closestDist > 0.11:
            return "Unknown"

        closestSubtype = closestMatch.split('|')[1]  # subtype is second field

        # is closest match within tolerance for subtype?
        if closestDist < self.properties[closestSubtype]['dist_upper']:
            if self.properties[closestSubtype]['is_simple_CRF']:
                sufficientCoverage = True  # TODO: WIP, method with breakpoints
                if sufficientCoverage:
                    return closestSubtype
                return self.properties[closestSubtype]['parent']  # never reached
            else:
                return closestSubtype
        else:
            # distance exceeds upper limit
            if closestSubtype in ["CRF01_AE", "CRF02_AG"] or \
                    self.properties[closestSubtype]['is_simple_CRF']:

                if self.properties[closestSubtype]['parent'] == '-':
                    return "Unknown"

                # if simple recombinant, select most representative parent subtype
                parents = [x.strip() for x in self.properties[closestSubtype]['parent'].split("+")]

                # search through parent subtypes
                nextClosestParent = None
                min_dist = 1e6
                for parent in parents:
                    if parent == 'Unknown':
                        continue

                    for ref in self.simple_subtypes[parent]:
                        this_dist = dists[ref]
                        if this_dist < min_dist:
                            nextClosestParent = parent
                            min_dist = this_dist

                if abs(min_dist - closestDist) <= 0.01:
                    return nextClosestParent
                else:
                    return self.properties[closestSubtype]['parent']

            # not a recombinant
            if self.properties[closestSubtype]['class_level'] == "SUBTYPE":
                return closestSubtype

        return "Unknown"


def main():
    from sierralocal.nucaminohook import NucAminoAligner
    subtyper = Subtyper()
    # test out all reference sequences on this system
    #for h, ref in subtyper.subtype_references.items():
    #    subtype = re.findall('\|(.*?)\||$', h)[0]
    #    guess = subtyper.getClosestSubtype(ref)
    #    if not guess == subtype:
    #        print(h, guess, subtype)
    if len(sys.argv) < 2:
        print('Cmdline use: python3 subtyper.py <input FASTA>')
        sys.exit()

    input_file = sys.argv[1]

    # call nucamino
    aligner = NucAminoAligner()
    aligner.align_file(input_file)

    with open(input_file) as handle:
        sequence_list = get_input_sequences(handle, return_dict=False)

    # open the NucAmino output file
    with open(os.path.splitext(input_file)[0] + '.tsv', 'r') as nucamino_alignment:
        tsvin = csv.reader(nucamino_alignment, delimiter='\t')
        next(tsvin)  # bypass the header row

        for idx, row in enumerate(tsvin):  # iterate over sequences (1 per row)
            header, firstAA, _, firstNA, _, _ = row[:6]
            sequence = sequence_list[idx]  # NucAmino preserves input order

            firstAA = int(firstAA)  # relative to pol start
            firstNA = int(firstNA)  # this should be used to adjust offset

            gene = aligner.get_gene(firstAA)
            assert gene is not None, "Fatal error in get_mutations"

            closestSubtype = subtyper.getClosestSubtype(
                gene=gene, sequence=sequence, offset=((firstAA-1)-56)*3
            )
            print (header, closestSubtype)


if __name__ == '__main__':
    main()
