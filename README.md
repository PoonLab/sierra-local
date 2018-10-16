
[![PyPI](https://img.shields.io/pypi/v/sierralocal.svg)](https://pypi.org/project/sierralocal/)

# sierra-local
sierra-local is a local implementation of Stanford University's HIVdb web interface Sierra. It generates drug resistance predictions given HIV-1 sequence data, and outputs this in a standard format.

## Rationale

The Stanford HIVdb algortihm is a widely used method for predicting the drug resistance phenotype of an HIV-1 infection based on its genetic sequence, specifically the complete or partial sequence of the genomic regions encoding the primary targets of modern antiretroviral therapy.  Prediction of HIV-1 drug resistance is an important component in the routine clinical management of HIV-1 infection, being faster and more cost-effective than the direct measurement of drug resistance by culturing virus isolates in the laboratory.  The HIVdb algorithm is essentially rules-based classifier that is actively maintained and released to the public domain in the ASI ([Algorithm Specification Interface](http://jcm.asm.org/content/41/6/2792.short)) exchange format, demonstrating a laudable commitment by the HIVdb developers to open-source research and clinical practice.

The HIVdb algorithm is usually accessed through a web service hosted at Stanford University ([Sierra](https://hivdb.stanford.edu/hivdb)).  While this is a convenient format for many clinical laboratories, it requires a network connection and the transmission of potentially sensitive patient-derived data to a remote server.  Transmitting sequence data over the web may present a bottleneck for laboratories located at sites that are geographically distant from the host server, or where network traffic is prone to service disruptions.  Furthermore, the use of HIV-1 sequence data in criminal cases raises significant issues around data privacy.

Our objective was to build a lightweight, open-source Python implementation of the HIVdb algorithm for processing data on a local computer without sending any data over the network.  During the development of sierra-local, the maintainers of [Sierra](https://github.com/hivdb/sierra) released the [source code](https://github.com/hivdb/sierra) for their web service under a permissive free software license (GPL v3.0).  We were thrilled that the [HIVdb developers](https://github.com/hivdb) elected to release their server code, but we remained committed to complete sierra-local so the HIV research and clinical communities can process their own data without needing to install and maintain an Apache server, build an SQL database, or to install a sizeable number of software [dependencies](https://github.com/hivdb/sierra#dependency-lists).


## Dependencies
We tried to minimize dependencies:
- Python 3 (tested on [Python 3.6.5](https://www.python.org/downloads/release/python-365/))
- [NucAmino](https://github.com/hivdb/nucamino) `v0.1.3` or later
- An internet connection to set up or update HIVdb files

## Installation

### For Use
Currently Linux only:
```bash
pip3 install sierralocal
```

### For Development
1. Clone this repository.
    ```
    git clone https://github.com/PoonLab/sierra-local.git
    ```
2. Download the correct [pre-compiled NucAmino](https://github.com/hivdb/nucamino) binary for your platform, *e.g.*:
   ```
   # for Linux
   wget https://github.com/hivdb/nucamino/releases/download/v0.1.3/nucamino-linux-amd64
   ```
    and place it in the `sierralocal` directory (not the root directory `sierra-local`).  
    You might also need to modify the user permissions for the binary file; for example: `chmod 755 nucamino`.


## Using sierra-local

### Command-line interface (CLI)

We provide shell access to the program:
```
sierralocal SEQUENCES.fasta -o OUTPUT.json
```

Optional flags:
```
-xml		Path to HIVDB ASI XML file
-skipalign	Skip NucAmino alignment if TSV file already present.
-cleanup	Deletes NucAmino alignment file after processing.
```
### Graphical user interface (GUI)
From the root directory of the project:
```
python3 gui.py
```

### As a Python module
If you did not install via `pip` or `pip3`, must install as a package locally. Navigate to the project root directory:
```bash
pip3 install -e .
```
```python
import sierralocal
sierralocal.score(filename)
#TODO: make sierralocal.score() return the JSON text. Currently it just runs the program and the JSON is stored as a file.
```

## About Us
This project was developed at the Poon Lab under the Department of Pathology and Laboratory Medicine, Schulich School of Medicine and Dentistry, Western University, London, Ontario.
