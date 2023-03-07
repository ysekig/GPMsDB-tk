# GPMsDB-tk

GPMsDB-tk v1.0.1 was released on March 7, 2023. 

GPMsDB-tk/GPMsDB-dbtk are software toolkits for assigning taxonomic identification to user-provided MALDI-TOF mass spectrometry profiles obtained from bacterial and archaeal cultured isolates. They take advantages of a newly developed database of protein mass profiles predicted from ~200,000 bacterial and archaeal genome sequences. This toolkit is also designed to work with customized databases, allowing microbial identification based on user-provided genome/metagenome-assembled genome (MAG) sequences. The ms-identification-tk is open source and released under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License. 

Please post questions and issues related to ms-identification-tk on the Issues section of the GitHub repository.

## Installing and using GPMsDB-tk

Prerequisites
* Python (version 3.7 or higher)
* Cython (version 0.29.1 or higher)
* matplotlib (version 3.5.0 or higher)

In the source directory, the following command will compile and install the software in your python environment.
```bash
python setup.py install
```

GPMsDB-tk requires ~60 GB of external data that needs to be downloaded and unarchived:
```bash
wget https://....data.tar.gz
tar xvzf ....data.tar.gz
```

GPMsDB-tk requires an environment variable named GPMsDB_PATH to be set to the directory containing the unarchived reference data.
```bash
export GPMsDB_PATH=/path/to/release/package/
```

### Features

* Peak-list characterization:
  * inspect       -> Inspection of a peak-list and generate peak plots
  * adjust        -> m/z adjustment for given peak-list
* Strain identification based on peak-list(s)
  * identify      -> Search for the best-matching genome(s) without m/z adjustment
  * identify_wf   -> Full identification workflow
	* identify_bwf   -> Full identification workflow for a batch of files
* Peak annotation
  * peak          -> Annotate protein names and Tigrfam/Pfam markers genes
  * peak_wf       -> Full peak-list characterization workflow
  * peak_bwf      -> Full peak-list characterization workflow for a batch of files
			
## Bug Reports

Please report bugs through the GitHub issues system, or contact Yuji Sekiguchi (y.sekiguchi@aist.go.jp)

## Copyright

Copyright (C) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)

This package is under the conditions of the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License. See LICENSE for further details.
