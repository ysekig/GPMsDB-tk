# GPMsDB-tk

GPMsDB-tk v1.0.1 was released on March 7, 2023. 

GPMsDB-tk/GPMsDB-dbtk are software toolkits for assigning taxonomic identification to user-provided MALDI-TOF mass spectrometry profiles obtained from bacterial and archaeal cultured isolates. They take advantages of a newly developed database of protein mass profiles predicted from ~200,000 bacterial and archaeal genome sequences. This toolkit is also designed to work with customized databases, allowing microbial identification based on user-provided genome/metagenome-assembled genome (MAG) sequences. The GPMsDB-tk is open source and released under the GNU General Public License (Version 3). 

Please post questions and issues related to GPMsDB-tk on the Issues section of the GitHub repository.

## Installing and using GPMsDB-tk

Prerequisites
* Python (version 3.7 or higher)
* Cython (version 0.29.1 or higher)
* matplotlib (version 3.5.0 or higher)

In the source directory, the following command will compile and install the software in your python environment.
```bash
git clone https://github.com/ysekig/GPMsDB-tk
cd GPMsDB-tk
python setup.py install
```
During the installation, you may see some deprecation warnings like “easy_install command is deprecated” but this will not cause any issues for GPMsDB-tk.

GPMsDB-tk requires ~7 GB of external data that needs to be downloaded from Open Science Framework (OSF) (DOI: [10.17605/OSF.IO/YA5HF](https://doi.org/10.17605/OSF.IO/YA5HF)) or Zenodo (DOI: [10.5281/zenodo.7703482](https://zenodo.org/record/7703483#.ZAbPNS_3Jf0)) and unarchived:

```bash
tar xvzf R01-RS95.tar.gz
```

GPMsDB-tk requires an environment variable named GPMsDB_PATH to be set to the directory containing the unarchived reference data.
```bash
export GPMsDB_PATH=/path/to/release/package/
```

If you are interested in customizing the database with user-provided genomes/metagenome-assembled genomes (MAGs), [GPMsDB-dbtk](https://github.com/ysekig/GPMsDB-dbtk) should also be installed.

## Features

* Peak-list characterization:
  * inspect       -> Inspection of a peak-list and generate peak plots
  * adjust        -> m/z adjustment for given peak-list
* Strain identification based on peak-list(s)
  * identify      -> Search for the best-matching genome(s) without m/z adjustment
  * identify_wf   -> Full identification workflow (option "-aa" should be set for m/z adjustment)
  * identify_bwf   -> Full identification workflow for a batch of files (option "-aa" should be set for m/z adjustment)
* Peak annotation
  * peak          -> Annotate protein names and Tigrfam/Pfam markers genes
  * peak_wf       -> Full peak-list characterization workflow (option "-aa" should be set for m/z adjustment)
  * peak_bwf      -> Full peak-list characterization workflow for a batch of files (option "-aa" should be set for m/z adjustment)

## Output values for the option "identify"

* protein_hit: number of hits with all proteins predicted for reference genome
* ribosomal_hit: number of hits with ribosomal proteins predicted for reference genome
* score: matching value calculated for reference genome (higher better matching)
* probability value: the frequency of appearance of a given score inferred based on scores from 100 randomely selected reference genomes
* likelihood(%): likelihood of correct identification (empirical, based on ribosomal_hit and probability.
* ncbi_name: NCBI oraganism name (genus/species) for reference genome
* ncbi_strain: NCBI strain name for reference genome
* taxonomy_gtdb: GTDB taxonomy string for reference genome

## Bug Reports

Please report bugs through the GitHub issues system, or contact Yuji Sekiguchi (y.sekiguchi@aist.go.jp)

## Copyright

Copyright (C) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)

This package is under the conditions of the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License. See LICENSE for further details.
