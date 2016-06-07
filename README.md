# cbio-proteomics
This is a repository for all of the data processing scripts for the transfer of [CPTAC](https://cptac-data-portal.georgetown.edu/cptacPublic/) data into [cBioPortal](https://github.com/cBioPortal/cbioportal) as part of 2016's [Google Summer of Code](https://developers.google.com/open-source/gsoc/). The purpose of this is not only to produce flat text files for import into the cBioPortal database, but it's also to do data exploration and cross-dataset normalization. Eventually this will all be incorporated into the [cBioPortal visualization interface](http://www.cbioportal.org/), and the necessary files will be merged into the main repository.

## Usage
Cloning the repo and running the commands in `pipeline.sh` into a Bash shell should basically reproduce everything in my current working directory, including all subfolders, static files, and plots that are not shown.

## Acknowledgements
Thanks to my PI David Fenyo and the GSoC mentors at MSKCC, JJ Gao and Zack Heins, for guidance.
