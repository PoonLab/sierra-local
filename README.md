
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
- [NucAmino](https://github.com/hivdb/nucamino) `v0.1.3` or later (included with package).

## Installation

On a Linux system, you can install *sierra-local* as follows:
```
git clone http://github.com/PoonLab/sierra-local
cd sierra-local
sudo python3 setup.py install
```
Note that you need super-user privileges to install the package by this method.  For more detailed instrucitons, please refer to the document [INSTALL.md](INSTALL.md) that should be located in the root directory of this Python package.

## Using sierra-local

### Command-line interface (CLI)

To run a quick example, use the following sequence of commands:
```console
art@Jesry:~/git/sierra-local$ python3 scripts/retrieve_hivdb_data.py RT RT.fa
art@Jesry:~/git/sierra-local$ sierralocal RT.fa
searching path /usr/local/lib/python3.6/dist-packages/sierralocal/data/HIVDB*.xml
Error: could not find local copy of HIVDB XML, attempting download...
Updated HIVDB XML from https://hivdb.stanford.edu/assets/media/HIVDB_8.7.9e470b87.xml into /usr/local/lib/python3.6/dist-packages/sierralocal/data/HIVDB_8.7.9e470b87.xml
/usr/local/lib/python3.6/dist-packages/sierralocal/data/apobec.tsv
Error: could not retrieve APOBEC DRM data
Updated APOBEC DRMs from https://hivdb.stanford.edu/assets/media/apobec-drms.5b7e1215.tsv into /usr/local/lib/python3.6/dist-packages/sierralocal/data/apobec.tsv
HIVdb version 8.7
Found NucAmino binary /usr/local/lib/python3.6/dist-packages/sierralocal/bin/nucamino-linux-amd64
Aligned RT.fa
100 sequences found in file RT.fa.
Writing JSON to file RT_results.json
Time elapsed: 7.476 seconds (17.803 it/s)
```
`retrieve_hivdb_data.py` is a Python script that we provided to download small samples of HIV-1 sequence data from the Stanford HIVdb database.  In this case, we have retrieved 100 reverse transcriptase (RT) sequences and processsed them with the *sierra-local* pipeline.  By default, the results are written to the file `[FASTA basename]_results.json`:
```console
art@Jesry:~/git/sierra-local$ head RT_results.json 
[
  {
    "inputSequence": {
      "header": "U54771.CM240.CRF01_AE.0"
    },
    "subtypeText": "CRF01_AE",
    "validationResults": [],
    "alignedGeneSequences": [
      {
        "firstAA": 1,
```


### As a Python module
If you have downloaded the package source to your computer, you can also run *sierra-local* as a Python module from the root directory of the package.  In the following example, we are calling the main function of *sierra-local* from an interactive Python session:
```console
art@Jesry:~/git/sierra-local$ git clone http://github.com/PoonLab/sierra-local
art@Jesry:~/git/sierra-local$ cd sierra-local
art@Jesry:~/git/sierra-local$ python3
Python 3.6.6 (default, Sep 12 2018, 18:26:19) 
[GCC 8.0.1 20180414 (experimental) [trunk revision 259383]] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from sierralocal.main import sierralocal
>>> sierralocal('RT.fa', 'RT.json')
searching path /home/art/git/sierra-local/sierralocal/data/HIVDB*.xml
Error: could not find local copy of HIVDB XML, attempting download...
Updated HIVDB XML from https://hivdb.stanford.edu/assets/media/HIVDB_8.7.9e470b87.xml into /home/art/git/sierra-local/sierralocal/data/HIVDB_8.7.9e470b87.xml
/home/art/git/sierra-local/sierralocal/data/apobec.tsv
Error: could not retrieve APOBEC DRM data
Updated APOBEC DRMs from https://hivdb.stanford.edu/assets/media/apobec-drms.5b7e1215.tsv into /home/art/git/sierra-local/sierralocal/data/apobec.tsv
HIVdb version 8.7
Found NucAmino binary /home/art/git/sierra-local/sierralocal/bin/nucamino-linux-amd64
Aligned RT.fa
100 sequences found in file RT.fa.
Writing JSON to file RT.json
(100, 1.9369409084320068)
```
Note that this doesn't require any `sudo` privileges.

### Graphical user interface (GUI)
We have a GUI in development that uses the Python Tkinter framework.  To use this script, you need to check out the *gui* branch of this repository and then call the script from the root directory of the project:
```
python3 gui.py
```


## About Us
This project was developed at the Poon Lab under the Department of Pathology and Laboratory Medicine, Schulich School of Medicine and Dentistry, Western University, London, Ontario.
