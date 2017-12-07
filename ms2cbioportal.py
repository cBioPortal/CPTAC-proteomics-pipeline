import sys
import argparse
import re
import numpy as np
import pandas as pd


class MSTable(object):
    def __init__(self, filenames, sample_regex, exp_type):
        self.df = None
        self.load_data(filenames)
        self.clean_rows(exp_type)
        self.clean_columns(exp_type, sample_regex)
        self.format_genes(exp_type)
        self.clean_values(exp_type)
    
    def __str__(self):
        if self.df:
            return str(self.df)
        else:
            return str(pd.DataFrame())
    
    def load_data(self, filenames):
        dfs = []
        for fname in filenames:
            try:
                df = pd.read_csv(fname, sep='\t', header=0)
                dfs.append(df)
            except OSError:
                print('File {0} not found. Could not be converted into a Pandas dataframe.').format(fname)
                raise
        self.df = pd.concat(dfs, axis=0)
        return
    
    def clean_columns(self, exp_type, sample_regex):
        sample_columns = [col for col in self.df.columns if re.match(sample_regex, col)]
        if not sample_columns:
            raise ValueError('`sample_regex` was not able to find any columns.')
        self.df = self.df[sample_columns]
        return
    
    def clean_rows(self, exp_type):
        if exp_type in ('cdap_itraq', 'cdap_pecursor_area'):
            self.df = self.df[~self.df['Gene'].str.contains('Mean|Median|StdDev', na=True)]
        elif exp_type == 'maxquant_proteome':
            self.df = self.df[~self.df['Protein IDs'].str.contains('REV|CON', na=True)]
            self.df = self.df[pd.notnull(self.df['Gene names'])]
            self.df = self.df[self.df['Q-value'] < 0.05]
        elif exp_type == 'maxquant_ptm':
            self.df = self.df[~self.df['Protein'].str.contains('REV|CON', na=True)]
            self.df = self.df[pd.notnull(self.df['Gene names'])]
            self.df = self.df[self.df['Localization prob'] > 0.75]
        return
    
    def format_genes(self, exp_type):
        if exp_type =='cdap_itraq':
            pass
        elif exp_type == 'cdap_pecursor_area':
            pass
        elif exp_type == 'maxquant_proteome':
            pass
        elif exp_type == 'maxquant_ptm':
            pass
    
    def clean_values(self, exp_type):
        self.df.replace([np.inf, -np.inf], np.nan, inplace=True)
        if exp_type == 'cdap_precursor_area':
            self.df = np.log2(self.df)
    
    def rename_columns(self):
        return







