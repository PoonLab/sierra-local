from pathlib import Path
import csv
import os
import re
from sierralocal.utils import get_input_sequences

class Subtyper():
    def __init__(self):
        self.subtype_references = self.getSubtypeReferences()
        self.properties = self.getGenotypeProperties()
        self.simple_subtypes = {}

    def getSubtypeReferences(self):
        """
        Parse FASTA file containing genotype reference sequences.
        Also extract simple subtypes from reference labels.
        :return:
        """
        # FIXME: this is a hard-coded data filename
        filepath = Path(os.path.dirname(__file__))/'data'/'genotype-references.9c610d61.fasta'
        handle = open(filepath)
        result = get_input_sequences(handle, return_dict=True)
        for label, seq in result.items():
            subtype = label.split('|')[1]
            if subtype not in self.simple_subtypes:
                self.simple_subtypes.update({subtype: []})
            self.simple_subtypes[subtype].append(label)
        return result

    def uncorrectedDistance(self, seq, ref):
        possibilies = {
            "A" : ["A"],
            "C" : ["C"],
            "G" : ["G"],
            "T" : ["T"],
            "R" : ["A","G"],
            "Y" : ["C","T"],
            "M" : ["A","C"],
            "W" : ["A","T"],
            "S" : ["C","G"],
            "K" : ["G","T"],
            "B" : ["C","G","T"],
            "D" : ["A","G","T"],
            "H" : ["A","C","T"],
            "V" : ["A","C","G"],
            "N" : ["A","C","G","T"],
            '~' : [],
            '.' : [],
            '-' : []
        }
        # Assume sequence has been aligned already
        # TODO: Step 1: mask SDRM in seq with ref nucleotide
        count = 0
        for idx, nuc in enumerate(seq):
            # print(idx, nuc)
            if not ref[idx] in possibilies[nuc]:
                count += 1
        return count / len(seq)

    def getDistances(self, sequence):
        """
        Uncorrected pairwise distances: the number of identical bases divided by
        the sequence length, ignoring positions with gaps.
        :param sequence:
        :return: a dictionary of (sequence label, distance) key-value pairs
        """
        dists = {}
        #TODO: generate concatenated sequences, with masked SDRMS
        for header, reference in self.subtype_references.items():
            dists[header] = self.uncorrectedDistance(sequence, reference)
        return dists


    def getGenotypeProperties(self):
        """
        Parse CSV file that summarizes the subtypes and CRFs in the HIVdb database
        :return:  A nested dictionary keyed by subtype/CRF name
        """
        # FIXME: this is a hard-coded path
        filepath = Path(os.path.dirname(__file__))/'data'/ 'genotype-properties.cc00f512.csv'
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


    def getClosestSubtype(self, sequence):
        dists = self.getDistances(sequence)

        intermed = [(v, k) for k, v in dists.items()]
        intermed.sort(reverse=False)  # increasing order
        closestMatch, closestDist = intermed[0]

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
                # if simple recombinant, select most representative parent subtype
                parents = [x.strip() for x in self.properties[closestSubtype]['parent'].split("+")]

                # search through parent subtypes
                nextClosestParent = None
                min_dist = 1e6
                for parent in parents:
                    for ref in self.simple_subtypes[parent]:
                        this_dist = dists[ref]
                        if this_dist < min_dist:
                            nextClosestParent = parent
                            min_dist = this_dist

                if abs(min_dist - closestDist) <= 0.01:
                    return nextClosestParent
                else:
                    return self.properties[closestSubtype]['parent']

            if self.properties[closestSubtype]['class_level'] == "SUBTYPE":
                return closestSubtype

def main():
    subtyper = Subtyper()
    # test out all reference sequences on this system
    for h, ref in subtyper.subtype_references.items():
        subtype = re.findall('\|(.*?)\||$', h)[0]
        guess = subtyper.getClosestSubtype(ref)
        if not guess == subtype:
            print(h, guess, subtype)


if __name__ == '__main__':
    main()
