from pathlib import Path
import csv

class Subtyper():
	def __init__(self):
		self.subtype_references = self.getSubtypeReferences()
		self.name, self.parentSubtype, self.distanceUpperLimit, self.isSimpleCRF, self.classificationLevel, self.breakPoints = self.getGenotypeProperties()

	def getSubtypeReferences(self):
		filepath = Path('.') / 'data' / 'genotype-references.9c610d61.fasta'
		return {}

	def getDistances(self, sequence):
		dists = {}
		#uncorrected pairwise distances: the number of identical bases divided by the number of bases being compared, ignoring any positions with gaps
		#generate concatenated sequenced, with masked SDRMS
		for reference in self.subtype_references:
			dist = 0.00 #calculate distance
			dists[reference] = dist
		return dists

	def getGenotypeProperties(self):
		filepath = Path('.') / 'data' / 'genotype-properties.cc00f512.csv'
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
				distanceUpperLimit[row[0]] = row[2]
				isSimpleCRF[row[0]] = row[3]
				classificationLevel[row[0]] = row[4]
				breakpoints[row[0]] = row[5]
		print(classificationLevel)
		return name, parentSubtype, distanceUpperLimit, isSimpleCRF, classificationLevel, breakpoints

	def getClosestSubtype(self, sequence):
		dists = self.getDistances(sequence)

		# sort dist dict
		closestMatch = dists[0]

		if dists[closestMatch] > 0.11:
			return "Unknown"

		if dists[closestMatch] < self.distanceUpperLimit[closestMatch]:
			if self.isSimpleCRF[closestMatch] == "False":
				return closestMatch

			if self.isSimpleCRF[closestMatch] == "True":
				sufficientCoverage = True #TODO method with breakpoints
				if sufficientCoverage:
					return closestMatch
				return self.parentSubtype[closestMatch]

		if dists[closestMatch] > self.distanceUpperLimit[closestMatch]:
			if closestMatch in ["CRF01_AE", "CRF02_AG"] or self.isSimpleCRF[closestMatch] == "True":
				parents = [x.trim() for x in self.parentSubtype[closestMatch].split("+")]
				nextClosestParent = parents[0]
				for parent in parents:
					if dists[parent] > dists[nextClosestParent]:
						nextClosestParent = parent
				if abs(dists[nextClosestParent] - dists[closestMatch]) <= 0.01:
					return nextClosestParent
				else:
					return self.parentSubtype[closestMatch]
			if self.classificationLevel[closestMatch] == "SUBTYPE":
				return closestMatch

def main():
	subtyper = Subtyper()

if __name__ == '__main__':
	main()