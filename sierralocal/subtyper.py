import os
import csv
import sys
from pathlib import Path
from time import time

from sierralocal.utils import get_input_sequences
#import re
#from sierralocal.nucaminohook import NucAminoAligner

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
        self.subtype_references = self.get_subtype_references()
        self.properties = self.get_genotype_properties()

    def get_subtype_references(self):
        """
        Parse FASTA file containing genotype reference sequences.
        Also extract simple subtypes from reference labels.
        @return result: dict, a dictionary of label:sequence pairs
        """
        # FIXME: this is a hard-coded data filename
        filepath = Path(os.path.dirname(__file__))/'data'/'genotype-references.9c610d61.fasta'
        with open(str(filepath)) as handle:
            result = get_input_sequences(handle, return_dict=True)
            for label, seq in result.items():
                subtype = label.split('|')[1]
                if subtype not in self.simple_subtypes:
                    self.simple_subtypes.update({subtype: []})
                self.simple_subtypes[subtype].append(label)
            return result

    def uncorrected_distance(self, seq, ref):
        """
        p-distance for nucleotide sequences.
        Assuming that sequences are aligned, count the number of base
        differences while accommodating ambiguous base calls.
        @param seq: str, first input sequence
        @param ref: str, second input sequences
        @return: float, proportion of base differences
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

    def get_distances(self, sequence, offset=0):
        """
        Uncorrected pairwise distances: the number of identical bases divided by
        the sequence length, ignoring positions with gaps.
        @param sequence: str, query nucleotide sequence
        @param offset: int, in the case of RT or IN, query sequence is missing
                        upstream nucleotides so we shift the reference
        @return: dict, a dictionary of (sequence label, distance) key-value pairs
        """
        dists = {}
        #TODO: generate concatenated sequences, with masked SDRMS
        for header, reference in self.subtype_references.items():
            dists[header] = self.uncorrected_distance(sequence, reference[offset:])
        return dists

    def get_genotype_properties(self):
        """
        Parse CSV file that summarizes the subtypes and CRFs in the HIVdb database
        @return props: dict, A nested dictionary keyed by subtype/CRF name
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

    def get_closest_subtype(self, sequence, offset):
        """
        Uses the output of get distances to get the closest subtype.
        @param sequence: str, query nucleotide sequence
        @param offset: int, in the case of RT or IN, query sequence is missing
                        upstream nucleotides so we shift the reference
        """
        dists = self.get_distances(sequence, offset)

        intermed = [(v, k) for k, v in dists.items()]
        intermed.sort(reverse=False)  # increasing order
        closest_dist, closest_match = intermed[0]

        # sort dist dict
        #sorted_dists = dict([(k, dists[k]) for k in sorted(dists, key=dists.get, reverse=False)])
        #closest_match = next(iter(sorted_dists.keys()))
        #closest_dist = sorted_dists[closest_match]
        #closest_subtype = re.findall('\|(.*?)\||$', closest_match)[0]

        # FIXME: where do we get this number from?
        if closest_dist > 0.11:
            return "Unknown"

        closest_subtype = closest_match.split('|')[1]  # subtype is second field

        # is closest match within tolerance for subtype?
        if closest_dist < self.properties[closest_subtype]['dist_upper']:
            if self.properties[closest_subtype]['is_simple_CRF']:
                sufficient_coverage = True  # TODO: WIP, method with breakpoints
                if sufficient_coverage:
                    return closest_subtype
                return self.properties[closest_subtype]['parent']  # never reached
            else:
                return closest_subtype
        else:
            # distance exceeds upper limit
            if closest_subtype in ["CRF01_AE", "CRF02_AG"] or \
                    self.properties[closest_subtype]['is_simple_CRF']:

                if self.properties[closest_subtype]['parent'] == '-':
                    return "Unknown"

                # if simple recombinant, select most representative parent subtype
                parents = [x.strip() for x in self.properties[closest_subtype]['parent'].split("+")]

                # search through parent subtypes
                next_closest_parent = None
                min_dist = 1e6
                for parent in parents:
                    if parent == 'Unknown':
                        continue

                    for ref in self.simple_subtypes[parent]:
                        this_dist = dists[ref]
                        if this_dist < min_dist:
                            next_closest_parent = parent
                            min_dist = this_dist

                if abs(min_dist - closest_dist) <= 0.01:
                    return next_closest_parent
                else:
                    return self.properties[closest_subtype]['parent']

            # not a recombinant
            if self.properties[closest_subtype]['class_level'] == "SUBTYPE":
                return closest_subtype

        return "Unknown"


def main():
    if len(sys.argv) < 2:
        print('Cmdline use: python3 subtyper.py <input FASTA>')
        sys.exit()

    # not running as package, so circular dependency is not an issue
    from sierralocal.hivdb import HIVdb
    from sierralocal.nucaminohook import NucAminoAligner
    algorithm = HIVdb()
    subtyper = Subtyper()

    # call nucamino
    input_file = sys.argv[1]
    aligner = NucAminoAligner(algorithm)
    records = aligner.align_file(input_file)

    for record in records:
        offset = max((record['FirstAA'] - 57) * 3, 0)
        subtype = subtyper.get_closest_subtype(record['Sequence'], offset)
        print('{}\t{}'.format(record['Name'], subtype))


if __name__ == '__main__':
    main()
