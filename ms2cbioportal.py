import sys
import abc
import argparse
import re
import numpy as np
import pandas as pd


class MSTable(object):
    
    __metaclass__  = abc.ABCMeta
    
    def __init__(self, filename, sample_regex, ptm_type=None):
        self.df = None
        
        # read in data from csv file
        self.load_data(filename)
        
        # perform tests for data columns and data length
        self.check_data()
        
        # clean and format the data
        self.format_genes(ptm_type)
        self.clean_rows()
        self.clean_columns(sample_regex)
        self.average_duplicates()
        self.clean_values()
    
    def __str__(self):
        if self.df:
            return str(self.df)
        else:
            return str(pd.DataFrame())
    
    def load_data(self, filename):
        try:
            self.df = pd.read_csv(filename, sep='\t', header=0)
        except OSError:
            print('File {0} not found. Could not be converted into a Pandas dataframe.').format(filename)
            raise
        return
        
    def clean_columns(self, sample_regex):
        sample_columns = [col for col in self.df.columns if re.match(sample_regex, col)]
        if not sample_columns:
            raise ValueError('`sample_regex` was not able to find any columns.')
        self.df = self.df[sample_columns]
        return
    
    @abc.abstractmethod
    def check_data(self):
        return
    
    @abc.abstractmethod
    def clean_rows(self):
        return
    
    @abc.abstractmethod
    def format_genes(self):
        return
    
    def clean_values(self):
        self.df.replace([np.inf, -np.inf], np.nan, inplace=True)
        return
    
    def average_duplicates(self):
        nondup_df = self.df[self.df.columns[~self.df.columns.duplicated(keep=False)]]
        uniq_dup_cols = list(set(self.df.columns[self.df.columns.duplicated(keep=False)]))
        dedup_df = pd.DataFrame(index=nondup_df.index, columns=uniq_dup_cols)
        for col in uniq_dup_cols:
            dedup_df[col] = self.df[col].mean(axis=1)
        self.df = pd.concat((nondup_df, dedup_df), axis=1)
        return
    
    def rename_columns(self, renamer):
        self.df.columns = [renamer(col) for col in self.df.columns]
        return
    
    def vertical_concat(self, table):
        self.df = pd.concat((self.df, table.df), axis=0)
    
    def write_csv(self, filename):
        self.df.to_csv(filename, sep='\t', header=True, index=True)
        return


class CDAPTable(MSTable):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def check_data(self):
        assert len(self.df.index) > 0
        assert 'Gene' in self.df.columns
        return
    
    def clean_rows(self):
        self.df = self.df[~self.df['Gene'].str.contains('Mean|Median|StdDev', na=True)]
        return
    
    def format_genes(self, ptm_type):
        if 'Phosphosite' in self.df.columns:
            anno = []
            for g, ph in zip(self.df['Gene'], self.df['Phosphosite']):
                if not np.isnan(ph):
                    ptm_site = ph.split(':')[1]
                    ptm_build = re.sub('[A-Z]+', lambda m: '_'+m.group(0), ptm_site.upper())[1:]
                    anno.append('{0}|{0}_{1}{2}'.format(g, ptm_type, ptm_build))
                else:
                    anno.append('{0}|{0}'.format(g))
        else:
            anno = ['{0}|{0}'.format(g) for g in self.df['Gene']]
        self.df.index = anno
        self.df.index.names = ['Hugo_Symbol']
        return


class CDAPiTraqTable(CDAPTable):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class CDAPPrecursorAreaTable(CDAPTable):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def clean_values(self):
        self.df = np.log2(self.df)
        self.df.replace([np.inf, -np.inf], np.nan, inplace=True)
        return


class MaxQuantProteomeTable(MSTable):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def check_data(self):
        assert len(self.df.index) > 0
        assert all(col in self.df.columns for col in ('Protein IDs', 
                                                      'Gene names', 
                                                      'Q-value'))
        return
    
    def clean_rows(self):
        self.df = self.df[~self.df['Protein IDs'].str.contains('REV|CON', na=True)]
        self.df = self.df[pd.notnull(self.df['Gene names'])]
        self.df = self.df[self.df['Q-value'] < 0.05]
        return
    
    def format_genes(self, ptm_type):
        anno = [g.split(';')[0] for g in self.df['Gene names']]
        self.df.index = anno
        self.df.index.names = ['Hugo_Symbol']
        return


class MaxQuantPTMTable(MSTable):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def check_data(self):
        assert len(self.df.index) > 0
        assert all(col in self.df.columns for col in ('Protein', 
                                                      'Amino acid', 
                                                      'Positions within proteins', 
                                                      'Gene names', 
                                                      'Localization prob'))
        return
    
    def clean_rows(self):
        self.df = self.df[~self.df['Protein'].str.contains('REV|CON', na=True)]
        self.df = self.df[pd.notnull(self.df['Gene names'])]
        self.df = self.df[self.df['Localization prob'] > 0.75]
        return
    
    def format_genes(self, ptm_type):
        prot_pos = [p.split(';')[0] for p in self.df['Positions within proteins']]
        aa_col = self.df['Amino acid'].tolist()
        genes = [g.split(';')[0] for g in self.df['Gene names']]
        anno = []
        for g, aa, pos in zip(genes, aa_col, prot_pos):
            anno.append('{0}|{0}_{1}{2}{3}'.format(g, ptm_type, aa, pos))
        self.df.index = anno
        self.df.index.names = ['Hugo_Symbol']
        return


META_TEXT = """cancer_study_identifier: {0}
genetic_alteration_type: PROTEIN_LEVEL
datatype: Z-SCORE
stable_id: protein_quantification
profile_description: Protein Quantification (Mass Spec)
show_profile_in_analysis_tab: true
profile_name: Protein levels (mass spectrometry by CPTAC)
data_filename: {1}"""


class MSMeta(object):
    def __init__(self, cancer_id, prot_file):
        self.text = META_TEXT.format(cancer_id, prot_file)
    
    def __str__(self):
        return self.text
    
    def write(self, filename):
        with open(filename, 'w') as f:
            f.write(self.text)
        return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This script converts proteome abundance data and PTM data from CPTAC or MaxQuant into importable text files for cBioPortal. It will also output a meta file for use with cbioportalImporter.py. Use for one tissue type (i.e. breast cancer) at a time.')
    
    parser.add_argument('--proteome-files', 
        metavar='<filename_1,filename_2,filename_n>', 
        help='comma-separated list of protein abundance filenames', 
        required=True)
    parser.add_argument('--proteome-pipeline', 
        choices=('itraq', 'precursor_area'), 
        help='processing pipeline (i.e. "itraq" for Breast and Ovarian, "precursor_area for Colon")', 
        required=True)
    
    args = parser.parse_args()
