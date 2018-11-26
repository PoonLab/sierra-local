---
title: 'sierra-local: A lightweight standalone application for drug resistance prediction'
tags:
  - Python
  - bioinformatics
  - HIV/AIDS
  - drug resistance
  - sequence analysis
  - clinical virology
authors:
  - name: Jasper C. Ho
    orcid: 0000-0003-0693-8940
    affiliation: 1
  - name: Garway T. Ng
    affiliation: 1
  - name: Mathias Renaud
    affiliation: 1
  - name: Art F. Y. Poon
    orcid: 0000-0003-3779-154X
    affiliation: "1, 2, 3"
affiliations:
  - name: Department of Pathology and Laboratory Medicine, Western University, London, ON, Canada
    index: 1
  - name: Department of Microbiology and Immunology, Western University, London, ON, Canada
    index: 2
  - name: Department of Applied Mathematics, Western University, London, ON, Canada
    index: 3
date: November 23, 2018
bibliography: main.bib
---

# Summary

Genotypic resistance interpretation systems for the prediction and
interpretation of HIV-1 antiretroviral resistance are an important part
of the clinical management of HIV-1 infection. Current interpretation
systems are generally hosted on remote webservers that enable clinical
laboratories to generate resistance predictions easily and quickly from
patient HIV-1 sequences encoding the primary targets of modern
antiretroviral therapy. However they also potentially compromise a
health provider’s ethical, professional, and legal obligations to data
security, patient information confidentiality, and data provenance.
Furthermore, reliance on web-based algorithms makes the clinical
management of HIV-1 dependent on a network connection. Here, we describe
the development and validation of *sierra-local*, an open-source
implementation of the Stanford HIVdb genotypic resistance interpretation
system for local execution, which aims to resolve the ethical, legal,
and infrastructure issues associated with remote computing. This package
reproduces the HIV-1 resistance scoring by the web-based Stanford HIVdb
algorithm with a high degree of concordance (99.997%) and a higher level
of performance than current methods of accessing HIVdb programmatically.


# Background and Rationale

Genotype-based prediction of human immunodeficiency virus type 1 (HIV-1)
drug resistance is an important component for the routine clinical
management of HIV-1 infection [@tural2002clinical; @Gunthard2014].
Detecting the presence of viruses carrying mutations that confer drug
resistance enables physicians to select an optimal drug combination for
that patient’s treatment regimen. Furthermore, genotyping by bulk
sequencing is a cost-effective alternative to the direct measurement of
drug resistance from culturing virus isolates in a laboratory
[@Mayer2001]. Provided access to affordable bulk sequencing at an
accredited laboratory for clinical microbiology, the interpretation of
HIV-1 sequence variation is the primary obstacle to utilizing resistance
genotyping for HIV-1 care. 

Fortunately, there are several HIV-1 drug
resistance interpretation algorithms that can be accessed at no cost
through web applications or services hosted by remote network servers,
such as the Standard University HIV Drug Resistance Database (HIVdb)
[@shafer2006], Agence Nationale de Recherche sur le SIDA (ANRS) AC11
[@Meynard2002], and Rega Institute [@VanLaethem2002] algorithms. The
Stanford HIVdb interpretation system can be accessed either through a
web browser at <http://hivdb.stanford.edu/hivdb> or programmatically
through its Sierra Web Service [@tang2012hivdb], which requires the
transmission of an HIV-1 sequence from a local computer over the network
to the remote server. This is a convenient arrangement for clinical
laboratories because there is no need to install any specialized
software, web browsers are ubiquitous and most users are familiar with
submitting web forms.  Alternatively, a laboratory may host an
instance of the HIVdb Sierra Web Service itself, which was recently made
possible with the release of the Sierra source code. This approach,
however, requires the configuration of a web server, the Apache Tomcat
web container, and a large number of Java libraries.

There are a number of disadvantages to accessing
interpretation systems over a network connection. First, HIV-1 sequences
are sensitive patient information, not only because infection with HIV-1
remains a highly stigmatized condition, but also because sequence data
have been used as evidence in the criminal prosecution of individuals
for engaging in sexual intercourse without disclosing their infection
status, leading to virus transmission [@bernard2007hiv]. Once sequence
data have been transmitted to a remote server, one cedes all control
over data security. Preventing the onward distribution of the data and
deleting the data once the analysis is complete, for instance, is
entirely the responsibility of the system administrators of the host
server. Furthermore, unless the host server employs a secure transfer
protocol, the unencrypted data are transmitted in the clear between a
number of intermediary web servers, exposing these data to a
‘man-in-the-middle’ attack [@patil2014big].

Second, the algorithm hosted on the server is effectively a black box —
one has no insight into how resistance predictions are generated. Even
if a version of the algorithm has been released into the public domain,
one cannot be certain that the exact same algorithm was applied to their
transmitted data. Importantly, different versions of a given algorithm
can output significantly different resistance predictions, with the
general trend being an increase in both resistance scores and predicted
resistance levels [@Hart2018]. In addition to contributing to
inconsistencies in algorithm outputs, this makes it difficult to track
data provenance, *i.e.*, the historical record of data processing, that
has become recognized as a critical gap in the workflows of clinical
laboratories. For instance, the College of American Pathologists
recently issued new accreditation requirements stipulating that clinical
laboratories must track the specific version of software programs used
to process patient data [@aziz2014college]. Thus, a reliance on
web-based systems creates significant issues for the reproducibility and
quality assurance of clinical workflows. The Stanford HIVdb web service
(Sierra [@tang2012hivdb]), for instance, automatically utilizes the most
recent version of the HIVdb algorithm. While this constraint ensures
that users employ the most up-to-date algorithm, it also introduces
hidden changes to clinical pipelines, which may have been locally
validated on older versions of the algorithm.

Third, dependence on a web resource may cause problems when the
laboratory cannot access the host server, either due to local or
regional network outages, or because the host server is malfunctioning
or offline. In our experience, the web servers hosting the more popular
HIV drug resistance interpretation algorithms such as the Stanford HIVdb
database are reliable and well-maintained. However, it is not unusual
for other web-based algorithms to be relocated or go offline when the
developers move to other institutions or lack the resources to maintain
the service.

One of the important features of the Stanford HIVdb algorithm is that it
is regularly updated and released into the public domain in a
standardized XML-based interchange format — the Algorithm Specification
Interface version 2 (ASI2) format [@betts2003asi] — that was formulated
and published by the same developers in conjunction with the Frontier
Science Foundation. Here, we describe the implementation and validation
of *sierra-local*, an open-source Python package for local execution of
the HIVdb algorithm in the ASI2 format. This package utilizes, but does
not require, a network connection to synchronize its local ASI2 file and
reference data with the latest releases on the Stanford HIVdb web
server. Our objective was to release a lightweight alternative to
transmitting HIV-1 sequences to the HIVdb web server that minimizes the
number of software dependencies, and that produces the exact same
interpretations as the Sierra web service for all available HIV-1
sequences in the Stanford database.


# Validation

We obtained the entirety of the genotype-treatment correlation datasets
available through the Stanford HIV Drug Resistance Database (HIVdb
[@shafer2006rationale]) on May 7 2018.  After screening for invalid data, 
the resulting dataset contained 103,711 HIV-1 protease, 110,222 reverse 
transcriptase and 11,769 integrase entries.  We scored these data with 
both *sierra-local* and SierraPy (version 0.2.1) using the HIVdb version 
8.5 algorithm on both platforms. Because the algorithm was updated to 
version 8.6.1 during the validation experiments, we used the newer 
version for the HIV-1 integrase data sets since the update mostly 
affected the interpretation of mutations within this region. Out of 
the total 226,702 sequence records, the predicted resistance 
scores were completely identical in 226,696 (99.997%). 

In addition, we retrieved 7 population-based HIV-1 *pol* datasets from 
Genbank using the NCBI PopSet interface (<http://www.ncbi.nlm.nih.gov/popset>). 
These datasets were selected from the most recent uploads of substantial 
numbers of HIV-1 sequences covering the regions encoding both PR and RT, 
and representing a diversity of HIV-1 subtypes and sampling locations 
around the world.  All resistance scores for all 1,006 sequences were 
completely concordant between the pipelines.

[lrL[2in]{}rL[1in]{}r]{} Country/region & Sample & Subtypes & Sequence &
Accession & Ref.\
& size & & length (nt) & numbers &\
Brazil & 103 & B (100%) & 1262.0 & MF545238 – MF545340 &\
Ethiopia & 67 & C (97.0%), B (1.5%), A (1.5%) & 1042.1 & MH324937 –
MH325003 & [@arimide2018hiv]\
Guinea-Bissau & 54 & CRF02\_AG (88.9%), A (5.6%), CRF06\_cpx (1.8%) &
1035.0 & MH605452 – MH605505 & [@wilhelmson2018prevalence]\
Hong Kong & 284 & C (36.0%), CRF07\_BC (36.0%), CRF02\_AG (8.8%) &
1157.8 & MH757122 – MH757405 &\
South Africa & 212 & C (100%) & 1195.0 & MH920641 – MH920852 &
[@rasmussen2017external]\
Tajikistan & 146 & A (97.3%), CRF02\_AG (2.0%), CRF63\_02A1 (0.7%) &
1351.1 & MH543115 – MH543260 &\
Uganda & 140 & D (99.3%), CRF10\_CD (0.7%) & 864.0 & MH925538 – MH925677
&\


# Performance

Performance was evaluated on a workstation running Ubuntu
18.04 LTS with an Intel Xeon E5-1650 v4 hexa-core CPU at 3.60 GHz and 16
GB of DDR4-2400 RAM with a gigabit network connection. *sierra-local*
achieved mean \[range\] processing speeds of 47.08 sequences/second
(seq/s) \[45.07, 48.49\] for PR, 16.20 seq/s \[14.01, 19.97\] for RT,
and 14.99 seq/s \[14.79, 15.56\] for IN. A substantial fraction of
processing time was consumed by subtyping. SierraPy, with the same
dataset as previously described, yielded mean processing speeds of 16.01
seq/s \[12.88, 17.60\] for PR, 6.12 seq/s \[4.83, 7.54\] for RT, and
5.19 seq/s \[5.05, 5.47\] for IN. Although the size of sequence batches
used in this performance comparison likely is a factor in the results by
virtue of file writing and reading being done once per batch, the large
batch size used minimizes the effect of these I/O processes on the
overall runtime. Overall, *sierra-local* is able to process and return
results for submitted query HIV-1 *pol* sequences roughly 3 times faster
than SierraPy, depending on the nature of the sequences and the type of
local computing resources available. 

# Concluding remarks

The distribution of the HIVdb resistance genotyping algorithm in a
standardized format (ASI [@betts2003asi]) is an important resource for
HIV-1 research and clinical management, and an exemplary case of open
science. *sierra-local* provides a convenient framework to generate HIV
drug resistance predictions from ASI releases in a secure environment
and confers full control over data provenance. The ability to apply
ASI-encoded algorithms locally (offline) also makes this part of the
laboratory workflow robust to network availabilty may be particularly
important for laboratories situated in resource-limited settings. We 
hope this lightweight, open-source implementation of the HIVdb ASI 
will further democratize HIV drug resistance genotyping across providers of 
HIV care.


# Acknowledgements

We thank Philip Tzou for bringing NucAmino to our attention, and for his
contributions to open science in the release of the Stanford HIVdb
resistance program source code. This work was supported in part by the
Government of Canada through Genome Canada and the Ontario Genomics
Institute (OGI-131) and by a grant from the Canadian Institutes of
Health Research (PJT-156178). The funders had no role in study design,
data collection and interpretation, or the decision to submit the work
for publication.


# References

