
[![PyPI](https://img.shields.io/pypi/v/sierralocal.svg)](https://pypi.org/project/sierralocal/)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01186/status.svg)](https://doi.org/10.21105/joss.01186)

# sierra-local
sierra-local is a Python 3 implementation of the [Stanford University HIV Drug Resistance Database](https://hivdb.stanford.edu/) (HIVdb) [Sierra web service](https://hivdb.stanford.edu/page/webservice/) for generating drug resistance predictions from HIV-1 sequence data. This Python package enables laboratories to run this prediction algorithm without needing to transmit patient data over the network, and confers full control over [data provenance](https://en.wikipedia.org/wiki/Data_lineage#Data_Provenance) and security.

## Rationale

The Stanford HIVdb algortihm is a widely used method for predicting the drug resistance phenotype of an HIV-1 infection based on its genetic sequence, specifically the complete or partial sequence of the genomic regions encoding the primary targets of modern antiretroviral therapy.  Prediction of HIV-1 drug resistance is an important component in the routine clinical management of HIV-1 infection, being faster and more cost-effective than the direct measurement of drug resistance by culturing virus isolates in the laboratory.  The HIVdb algorithm is essentially rules-based classifier that is actively maintained and released to the public domain in the ASI ([Algorithm Specification Interface](http://jcm.asm.org/content/41/6/2792.short)) exchange format, demonstrating a laudable commitment by the HIVdb developers to open-source research and clinical practice.

The HIVdb algorithm is usually accessed through a web service hosted at Stanford University ([Sierra](https://hivdb.stanford.edu/hivdb)).  While this is a convenient format for many clinical laboratories, it requires a network connection and the transmission of potentially sensitive patient-derived data to a remote server.  Transmitting sequence data over the web may present a bottleneck for laboratories located at sites that are geographically distant from the host server, or where network traffic is prone to service disruptions.  Furthermore, the use of HIV-1 sequence data in criminal cases raises significant issues around data privacy.

Our objective was to build a lightweight, open-source Python implementation of the HIVdb algorithm for processing data on a local computer without sending any data over the network.  During the development of sierra-local, the maintainers of [Sierra](https://hivdb.stanford.edu/page/webservice/) released the [source code](https://github.com/hivdb/sierra) for their web service under a permissive free software license (GPL v3.0).  We were thrilled that the [HIVdb developers](https://github.com/hivdb) elected to release their server code, but we remained committed to complete sierra-local so the HIV research and clinical communities can process their own data without needing to install and maintain an Apache server, build an SQL database, or to install a sizeable number of software [dependencies](https://github.com/hivdb/sierra#dependency-lists).


## Dependencies
We tried to minimize dependencies:
- Python 3 (tested on [Python 3.6.5](https://www.python.org/downloads/release/python-365/))
- Python module [requests](https://pypi.org/project/requests/)
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
searching path /usr/local/lib/python3.6/dist-packages/sierralocal/data/apobec*.tsv
HIVdb version 8.8
Found NucAmino binary /usr/local/lib/python3.6/dist-packages/sierralocal/bin/nucamino-linux-amd64
Aligned RT.fa
100 sequences found in file RT.fa.
Writing JSON to file RT_results.json
Time elapsed: 5.2538 seconds (19.149 it/s)
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

We can also specify a different ASI (XML) file representing an earlier version of the HIVdb algorithm to reprocess the same data, and save the output to another file:
```console
art@Jesry:~/git/sierra-local$ sierralocal -xml sierralocal/data/HIVDB_8.5.d926dfff.xml RT.fa -o RT-v8.5.json
/usr/local/lib/python3.6/dist-packages/sierralocal/data/apobec.tsv
HIVdb version 8.5
Found NucAmino binary /usr/local/lib/python3.6/dist-packages/sierralocal/bin/nucamino-linux-amd64
Aligned RT.fa
100 sequences found in file RT.fa.
Writing JSON to file RT-v8.5.json
Time elapsed: 5.5014 seconds (18.337 it/s)
```

We find that switching versions of the algorithm from 8.5 to 8.7 results in substantial changes in resistance scores for these data with the introduction of a new drug [doravirine](https://aidsinfo.nih.gov/drugs/546/doravirine/0/patient) (DOR).  In addition, two of 100 cases were scored differently:
```console
art@Jesry:~/git/sierra-local$ python3 scripts/json2csv.py RT_results.json RT-v8.7.json.csv
art@Jesry:~/git/sierra-local$ python3 scripts/json2csv.py RT-v8.5.json RT-v8.5.json.csv
art@Jesry:~/git/sierra-local$ R
```
```R
> v5 <- read.csv('RT-v8.5.json.csv')
> v7 <- read.csv('RT-v8.7.json.csv')
> v7 <- v7[,-which(names(v7)=='DOR')]
> temp <- sapply(1:nrow(v5), function(i) any(v5[i,] != v7[i,]))
> which(temp)
[1] 23 63
> v5[23,]
                name subtype ABC AZT D4T DDI FTC LMV TDF EFV ETR NVP RPV
23 Y14503.BCF13.O.22 Group O  15 -10 -10  10  60  60 -10  50  45  95  65
> v7[23,]
                name subtype ABC AZT D4T DDI FTC LMV TDF EFV ETR NVP RPV
23 Y14503.BCF13.O.22 Group O  15 -10 -10  10  60  60 -10  50  55 105  75
> v5[63,]
                name subtype ABC AZT D4T DDI FTC LMV TDF EFV ETR NVP RPV
63 AF102332.A11.B.62       B  90 115 115  90  80  80  60   0   0   0   0
> v7[63,]
                name subtype ABC AZT D4T DDI FTC LMV TDF EFV ETR NVP RPV
63 AF102332.A11.B.62       B  90 115 115  90  80  80  60   0  10  10  10
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
searching path /home/art/git/sierra-local/sierralocal/data/apobec*.tsv
HIVdb version 8.8
Found NucAmino binary /home/art/git/sierra-local/sierralocal/bin/nucamino-linux-amd64
Aligned RT.fa
100 sequences found in file RT.fa.
Writing JSON to file RT.json
(100, 0.04676532745361328)
```
Note that this doesn't require any `sudo` privileges.

### Graphical user interface (GUI)
We have a GUI in development that uses the Python Tkinter framework.  To use this script, you need to check out the *gui* branch of this repository and then call the script from the root directory of the project:
```
python3 gui.py
```


## Updating the algorithm

The Stanford HIVdb database regularly updates its resistance genotyping algorithm and publishes the associated ASI2 XML file on their website.  In previous versions of *sierra-local*, we used Python to automatically query this website and download the newest version if it was not already present on the user's computer.  Subsequent changes to the Stanford HIVdb website, however, meant that users would have to install several additional dependencies in order for Python to locate the required files.  As a result, we decided to make the `updater.py` script an optional step of the pipeline.

To run `updater.py`, you need to install the following requirements:
* [Google Chrome](https://www.google.com/chrome/)
* The appropriate [chromedriver](https://chromedriver.chromium.org/) for your operating system *and* [version of Google Chrome](https://chromedriver.chromium.org/downloads/version-selection).
* The Python module [Selenium](https://pypi.org/project/selenium/)

You'll also need to be working with a local clone or download of this code repository, because you will probably need to modify a line of the `updater.py` script.

Next, you need to locate your `chromedriver` binary.  For the purpose of demonstration, I happened to leave this binary in my `~/Downloads` folder.  Use this information to modify the following lines in `updater.py`:
```python
# this needs to be modified to point Python to your local chromedriver binary
mod_path = Path(os.path.dirname(__file__))
driver_path = mod_path.parent / 'bin/chromedriver'
```
Otherwise, Python will be unable to locate the Chrome browser binary and throws a `WebDriverException`:
```console
art@orolo:~/git/sierra-local$ python3 sierralocal/updater.py
...
Traceback (most recent call last):
  File "sierralocal/updater.py", line 19, in <module>
    browser = webdriver.Chrome(executable_path=str(driver_path), chrome_options=options)
  File "/home/art/.local/lib/python3.6/site-packages/selenium/webdriver/chrome/webdriver.py", line 73, in __init__
    self.service.start()
  File "/home/art/.local/lib/python3.6/site-packages/selenium/webdriver/common/service.py", line 83, in start
    os.path.basename(self.path), self.start_error_message)
selenium.common.exceptions.WebDriverException: Message: 'chromedriver' executable needs to be in PATH. Please see https://sites.google.com/a/chromium.org/chromedriver/home
```

For my system running Ubuntu, I modified the `updater.py` script as follows:
```python
#driver_path = mod_path.parent / 'bin/chromedriver'
driver_path = Path('/home/art/Downloads/chromedriver')
```

Manually running the script enabled me to grab the most recent versions of the ASI2 and APOBEC files from the HIVdb webserver:
```console
art@orolo:~/git/sierra-local$ sudo python3 sierralocal/updater.py 
Updated HIVDB XML from https://hivdb.stanford.edu/assets/media/HIVDB_8.8.a126e04c.xml into sierralocal/data/HIVDB_8.8.a126e04c.xml
Updated APOBEC DRMs from https://hivdb.stanford.edu/assets/media/apobec-drms.221b0330.tsv into sierralocal/data/apobec.tsv
```

Now of course, it would be much simpler to manually point your browser to the Stanford HIVdb website and download these files yourself, but in some applications there may be a benefit to automating this step.


## About Us
This project was developed at the [Poon lab](http://github.com/PoonLab) within the [Department of Pathology and Laboratory Medicine](https://www.schulich.uwo.ca/pathol/), Schulich School of Medicine and Dentistry, [Western University](http://uwo.ca), London, Ontario.  Development of *sierra-local* was supported in part by a grant from the [Canadian Institutes of Health Research](http://www.cihr-irsc.gc.ca/e/193.html) (PJT-156178).

If you use *sierra-local* for your work, please cite the following paper:
* **sierra-local: A lightweight standalone application for drug resistance prediction.**
  Jasper C Ho, Garway T Ng, Mathias Renaud, Art FY Poon, (2019). 
  *Journal of Open Source Software*, 4(33), 1186, https://doi.org/10.21105/joss.01186

If you want to reference the validation of *sierra-local* on HIV-1 *pol* data, please cite the following preprint:
* **sierra-local: A lightweight standalone application for secure HIV-1 drug resistance prediction.**
  Jasper C Ho, Garway T Ng, Mathias Renaud, Art FY Poon.
  *bioRxiv* 393207; doi: https://doi.org/10.1101/393207
