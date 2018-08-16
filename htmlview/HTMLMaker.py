import json
from pathlib import Path
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help='Input file.')
    args = parser.parse_args()
    return args

def main():
	args = parse_args()

	html_template = [
		"<html>",
		"<head>",
		# '<link rel="stylesheet" href="css/styles.css">',
		'<link rel="stylesheet" href="css/bootstrap.css" crossorigin="anonymous">',
		"</head>",
		"<body>",
	]

	with open(Path(args.file), "r") as file:
		js = json.load(file)

	#create resistance info
	html_template += ["<div class='col-md-8'>"]
	for idx, entry in enumerate(js):
		inputSequence = entry['inputSequence']
		subtypeText = entry['subtypeText']
		validationResults = entry['validationResults']
		alignedGeneSequences = entry['alignedGeneSequences'][0]
		firstAA = alignedGeneSequences['firstAA']
		lastAA = alignedGeneSequences['lastAA']
		gene = alignedGeneSequences['gene']
		mutations = alignedGeneSequences['mutations'][0]
		drugResistance = entry['drugResistance'][0]
		drugScores = drugResistance['drugScores']

		#scoring information
		html_template += ["<div class='page-header' tabindex={}><h2>{}</h2>".format(idx+1, inputSequence['header'])]
		html_template += ["<table class='table table-striped'>"]
		html_template += ['<thead><tr><th>Drug</th><th>Score</th><th>Breakdown</th></tr></thead>']
		html_template += ['<tbody>']
		#populate table
		for drugScore in drugScores:
			html_template += ["<tr><td>{}</td>".format(drugScore['drug']['displayAbbr'])]
			html_template += ["<td>{} ({})</td>".format(str(drugScore['score']), drugScore['text'])]
			pscoretext = ''
			for partialScore in drugScore['partialScores']:
				pscoretext += str(partialScore['score']) + ' '
				pscoretext += '+'.join([str(m['text']) for m in partialScore['mutations']])
				pscoretext += ' '
			html_template += ['<td>{}</td>'.format(pscoretext)]
			html_template +=  ["</tr>"]
		html_template += ['</tbody></table>']
		html_template += ['</div>']

	html_template.append("</div>")

	html_template.append("<br>")

	html_template += ["</body>","</html>"]

	with open("results.html", "w") as htmlout:
		htmlout.write("\n".join(html_template))

if __name__ == "__main__":
	main()