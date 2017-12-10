# CPTAC-proteomics-pipeline
This is a repository for all of the data processing scripts for the transfer of [CPTAC](https://cptac-data-portal.georgetown.edu/cptacPublic/) data into [cBioPortal](https://github.com/cBioPortal/cbioportal) as part of 2016's [Google Summer of Code](https://developers.google.com/open-source/gsoc/). The purpose of this is not only to produce flat text files for import into the cBioPortal database, but it's also to do data exploration and cross-dataset normalization. This has been incorporated into the [cBioPortal visualization interface](http://www.cbioportal.org/).

## Usage
This is a pretty specific package, so we designed it so that it was easy to use on-the-fly. First, clone the repo and `cd` in:

    git clone https://github.com/cBioPortal/CPTAC-proteomics-pipeline.git
    cd CPTAC-proteomics-pipeline

If you would like to have all the CPTAC files we used, please run the `wget` script:

    ./wget.sh

Please visit the [tutorial](https://github.com/cBioPortal/CPTAC-proteomics-pipeline/blob/master/tutorial/cbio_import_tutorial.ipynb), which goes through all the elements of the API.

NOTE: As shown in the tutorial, to import the classes, just add the relative location of the `ms2cbioportal.py` script to your current working directory. For example, since the tutorial is nested inside the repo:

    import sys
    sys.path.append('../')

## Acknowledgements
Thanks to my PI David Fenyo and the GSoC mentors at MSKCC, JJ Gao and Zack Heins, for guidance.
