from pathlib import Path
import csv
import os
import re

class Subtyper():
	def __init__(self):
		self.subtype_references = self.getSubtypeReferences()
		self.name, self.parentSubtype, self.distanceUpperLimit, self.isSimpleCRF, self.classificationLevel, self.breakPoints = self.getGenotypeProperties()

	def getSubtypeReferences(self):
		filepath = Path(os.path.dirname(__file__))/'data'/'genotype-references.9c610d61.fasta'
		sequences = []
		names = []
		sequence = ''

		with open(filepath, 'r') as handle:
		    res = {}
		    sequence = ''
		    for i in handle:
		        if i[0] == '$': # skip h info
		            continue
		        elif i[0] == '>' or i[0] == '#':
		            if len(sequence) > 0:
		                res.update({h: sequence})
		                sequence = ''   # reset containers
		                h = i.strip('\n')[1:]
		            else:
		                h = i.strip('\n')[1:]
		        else:
		            sequence += i.strip('\n').upper()
		    res.update({h: sequence})
		return res

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
		# Step 1: mask SDRM in seq with ref nucleotide
		count = 0
		for idx, nuc in enumerate(seq):
			# print(idx, nuc)
			if not ref[idx] in possibilies[nuc]:
				count += 1
		return count / len(seq)

	def getDistances(self, sequence):
		dists = {}
		#uncorrected pairwise distances: the number of identical bases divided by the number of bases being compared, ignoring any positions with gaps
		#generate concatenated sequenced, with masked SDRMS
		# print(self.subtype_references)
		for header, reference in self.subtype_references.items():
			dists[header] = self.uncorrectedDistance(sequence, reference)
		return dists

	def getGenotypeProperties(self):
		filepath = Path(os.path.dirname(__file__))/'data'/ 'genotype-properties.cc00f512.csv'
		name = {}
		parentSubtype = {}
		distanceUpperLimit = {}
		isSimpleCRF = {}
		classificationLevel = {}
		breakpoints = {}

		with open(filepath, 'r') as file:
			reader = csv.reader(file, delimiter='\t')
			next(reader, None)
			for row in reader:
				name[row[0]] = row[0]
				parentSubtype[row[0]] = row[1]
				distanceUpperLimit[row[0]] = float(row[2].replace('%',''))/100
				isSimpleCRF[row[0]] = bool(row[3])
				classificationLevel[row[0]] = row[4]
				breakpoints[row[0]] = row[5]
		return name, parentSubtype, distanceUpperLimit, isSimpleCRF, classificationLevel, breakpoints

	def getClosestSubtype(self, sequence):
		dists = self.getDistances(sequence)
		# sort dist dict
		sorted_dists = dict([(k, dists[k]) for k in sorted(dists, key=dists.get, reverse=False)])


		closestMatch = next(iter(sorted_dists.keys()))
		closestSubtype = re.findall('\|(.*?)\||$', closestMatch)[0]

		if sorted_dists[closestMatch] > 0.11:
			return "Unknown"

		if sorted_dists[closestMatch] < self.distanceUpperLimit[closestSubtype]:
			if self.isSimpleCRF[closestSubtype] == False:
				return closestSubtype

			if self.isSimpleCRF[closestSubtype] == True:
				sufficientCoverage = True #TODO method with breakpoints
				if sufficientCoverage:
					return closestSubtype
				return self.parentSubtype[closestSubtype]

		if sorted_dists[closestMatch] > self.distanceUpperLimit[closestSubtype]:
			if closestSubtype in ["CRF01_AE", "CRF02_AG"] or self.isSimpleCRF[closestSubtype] == True:
				parents = [x.trim() for x in self.parentSubtype[closestSubtype].split("+")]
				nextClosestParent = parents[0]
				for parent in parents:
					if sorted_dists[parent] > sorted_dists[nextClosestParent]:
						nextClosestParent = parent
				if abs(sorted_dists[nextClosestParent] - sorted_dists[closestSubtype]) <= 0.01:
					return nextClosestParent
				else:
					return self.parentSubtype[closestSubtype]
			if self.classificationLevel[closestSubtype] == "SUBTYPE":
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