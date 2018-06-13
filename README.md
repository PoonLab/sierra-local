# sierra-local
sierra-local is a local implementation of Stanford University's HIVdb web interface Sierra. It generates drug resistance predictions given HIV-1 sequence data, and outputs this in a standard format.

## Rationale

The Stanford HIVdb algortihm is a widely used method for predicting the drug resistance phenotype of an HIV-1 infection based on its genetic sequence, specifically the complete or partial sequence of the genomic regions encoding the primary targets of modern antiretroviral therapy.  Prediction of HIV-1 drug resistance is an important component in the routine clinical management of HIV-1 infection, being faster and more cost-effective than the direct measurement of drug resistance by culturing virus isolates in the laboratory.  The HIVdb algorithm is essentially rules-based classifier that is actively maintained and released to the public domain in the ASI ([Algorithm Specification Interface](http://jcm.asm.org/content/41/6/2792.short)) exchange format, demonstrating a laudable commitment by the HIVdb developers to open-source research and clinical practice.

The HIVdb algorithm is usually accessed through a web service hosted at Stanford University ([Sierra](https://hivdb.stanford.edu/hivdb)).  While this is a convenient format for many clinical laboratories, it requires a network connection and the transmission of potentially sensitive patient-derived data to a remote server.  Transmitting sequence data over the web may present a bottleneck for laboratories located at sites that are geographically distant from the host server, or where network traffic is prone to service disruptions.  Furthermore, the use of HIV-1 sequence data in criminal cases raises significant issues around data privacy.

Our objective was to build a lightweight, open-source Python implementation of the HIVdb algorithm for processing data on a local computer without sending any data over the network.  During the development of sierra-local, the maintainers of [Sierra](https://github.com/hivdb/sierra) released the [source code](https://github.com/hivdb/sierra) for their web service under a permissive free software license (GPL v3.0).  We were thrilled that the [HIVdb developers](https://github.com/hivdb) elected to release their server code3, but we remained committed to complete sierra-local so the HIV research and clinical communities can process their own data without needing to install and maintain an Apache server, build an SQL database, or to install a sizeable number of software [dependencies](https://github.com/hivdb/sierra#dependency-lists).


## Dependencies
- Python 3
- [NucAmino](https://github.com/hivdb/nucamino)

## Installation
1. Clone this repository.
    ```bash
    git clone https://github.com/PoonLab/sierra-local.git
    ```
2. Download the correct [pre-compiled NucAmino binary](https://github.com/hivdb/nucamino#download-binaries) for your platform and place it in the `sierra-local` root directory.
3. Pull algorithm and essential data using `update_HIVDB.py`. This needs to be done before first use.
    ```bash
    python3 ./scripts/update_HIVDB.py
    ```

## Using sierra-local
### Example use case:

Navigate to the root directory of the project.
```bash
python3 sierralocal.py SEQUENCES.fasta -o OUTPUT.json
```
### Optional flags:
```bash
-xml		Path to HIVDB ASI XML file
-skipalign	Skip NucAmino alignment if TSV file already present.
-cleanup	Deletes NucAmino alignment file after processing.
```

## About Us
This project was developed at the Poon Lab under the Department of Pathology and Laboratory Medicine, Schulich School of Medicine and Dentistry, Western University, London, Ontario.
