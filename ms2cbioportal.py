import abc
import argparse
import re
import pickle
import numpy as np
import pandas as pd


BASE_ATTR = """
    
    Parameters
    ----------
    filename : str
        The path and the name of the file.
    
    sample_regex : str
        A regex string for isolating sample columns.
    
    ptm_type : str, None
        The type of PTM (i.e. 'P' for phosphoprotein).
    
    use_ruler : bool
        Whether to use the proteomic ruler method.
    
    Attributes
    ----------
    df : pandas.DataFrame
        A DataFrame containing the tabular experimental data.
"""


class MSTable(object):
    """Base class for all MS table types. This is an abstract class so it
    cannot be implemented. Please see CDAPiTraqTable, CDAPPrecursorAreaTable,
    MaxQuantProteomeTable, and/or MaxQuantPTMTable.
    """ + BASE_ATTR
    
    __metaclass__  = abc.ABCMeta
    
    def __init__(self, filename, sample_regex, ptm_type=None, use_ruler=False):
        # these are private attributes
        self.use_ruler = use_ruler
        self.self_path = __file__.split('/')[0] + '/'
        
        # these are public attributes
        self.df = pd.DataFrame()
        
        # read in data from csv file
        self.load_data(filename)
        
        # perform tests for data columns and data length
        self.check_data()
        
        # clean and format the data
        self.clean_rows()
        self.format_genes(ptm_type)
        self.clean_columns(sample_regex)
        self.clean_values()
        self.average_duplicates()
        
        # transform data
        if use_ruler:
            self.proteomic_ruler()
    
    def __repr__(self):
        return str(self.df)
    
    @property
    def columns(self):
        return self.df.columns
    
    @property
    def index(self):
        return self.df.index
    
    @property
    def shape(self):
        return self.df.shape
    
    def head(self, n=5):
        """A wrapper for Pandas' `head` method.
        
        Parameters
        ----------
        n : int
            Number of rows to return.
        
        Returns
        -------
        head : pandas.DataFrame
            The top `n` rows of the dataframe.
        """
        head = self.df.head(n)
        return head
    
    def load_data(self, filename):
        """Load in a tab-separated experimental results file. Attaches the
        DataFrame to the `df` attribute of the class.
        
        Parameters
        ----------
        filename : str
            The path and the name of the input file.
        
        Returns
        -------
        None
        """
        try:
            self.df = pd.read_csv(filename, sep='\t', header=0)
        except OSError:
            print('File {0} not found. Could not be converted into a Pandas dataframe.').format(filename)
            raise
        return
        
    def clean_columns(self, sample_regex):
        """Isolate sample columns from the larger tabular data using
        a regular expression. Modifies instance `df` attribute.
        
        Parameters
        ----------
        sample_regex : str
            A regex string for isolating sample columns.
        
        Returns
        -------
        None
        """
        sample_columns = [col for col in self.df.columns if re.match(sample_regex, col)]
        if not sample_columns:
            raise ValueError('`sample_regex` was not able to find any columns.')
        self.df = self.df[sample_columns]
        return
    
    @abc.abstractmethod
    def check_data(self):
        """Method containing checks that the data and the columns used
        for the processing exist. Requires method overloading."""
        return
    
    @abc.abstractmethod
    def clean_rows(self):
        """Method that removes extraneous rows depending on the experiment
        type. Requires method overloading."""
        return
    
    @abc.abstractmethod
    def format_genes(self, ptm_type):
        """Method that rewrites the gene column and sets it as the index.
        Requires method overloading."""
        return
    
    def clean_values(self):
        """Clean the data of `np.inf`.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        self.df.replace((np.inf, -np.inf), np.nan, inplace=True)
        return
    
    def log_transform(self):
        """Log2 transform data.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        self.df = np.log2(self.df + 1)
        return
    
    def proteomic_ruler(self):
        """Apply the proteomic ruler to the data to obtain copy number
        per cell estimates based on Wiśniewski et al. 2014. Uses the formula:
        
        Protein Copy Number = Protein MS Signal x (Avogadro's Number/Molar Mass)
                                                x (DNA Mass/Histone MS Signal)
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        
        # check for negative values
        min_val = np.nanmin(self.df.values)
        if min_val < 0:
            raise ValueError("""Negative values detected. Proteomic ruler only supports MS intensity values. 
                                If this is a mistake, set use_ruler=False.""")
        
        # load some constants
        molar_mass = pickle.load(open(self.self_path + 'molar_mass.pkl', 'rb'))
        avogadro = 6.022 * 10**23
        dna_mass = 6.5 * 10**-12 # 6.5 pg
        
        # get col of avogadro/molar mass
        molar_col = []
        for ind in self.df.index:
            try:
                molar_col.append(molar_mass[ind.split('|')[0]])
            except KeyError:
                molar_col.append(np.nan)
        
        molar_col = pd.Series(molar_col)
        molar_col.index = self.df.index
        avo_mol = molar_col.rtruediv(avogadro)
        
        # get row of sum of histones
        hist_row = self.df[self.df.index.to_series().str.contains('HIST', na=False)].sum()
        
        # apply equation
        self.df = self.df.div(hist_row).multiply(avo_mol, axis=0).multiply(dna_mass)
        return
    
    def average_duplicates(self):
        """Average any duplicated sample columns by taking the mean across
        the duplicate samples.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        nondup_df = self.df[self.df.columns[~self.df.columns.duplicated(keep=False)]]
        uniq_dup_cols = list(set(self.df.columns[self.df.columns.duplicated(keep=False)]))
        dedup_df = pd.DataFrame(index=nondup_df.index, columns=uniq_dup_cols)
        for col in uniq_dup_cols:
            dedup_df[col] = self.df[col].nanmean(axis=1)
        self.df = pd.concat((nondup_df, dedup_df), axis=1)
        return
    
    def rename_columns(self, renamer):
        """Rename all sample columns using a custom function.
        
        Parameters
        ----------
        renamer : function
            A function that applies a transformation on each column name.
        
        Returns
        -------
        None
        """
        self.df.columns = [renamer(col) for col in self.df.columns]
        return
    
    def vertical_concat(self, table):
        """Concatenate data from two experiments together by vertical concatenation.
        This is for combining proteins and PTMs from the same tissue type. It will
        not average duplicate genes.
        
        Parameters
        ----------
        table : ms2cbioportal.MSTable
            An instance derived from `MSTable`.
        
        Returns
        -------
        None
        """
        self.df = pd.concat((self.df, table.df), axis=0)
        return
    
    def horizontal_concat(self, table):
        """Concatenate data from two experiments together by horizontal concatenation.
        This is for combining the same experiment type but with two (or more, since
        you can chain the concatenation) sample sets.
        
        Parameters
        ----------
        table : ms2cbioportal.MSTable
            An instance derived from `MSTable`.
        
        Returns
        -------
        None
        """
        self.df = pd.concat((self.df, table.df), axis=1)
        return
    
    def write_csv(self, filename):
        """Write finished DataFrame to file.
        
        Parameters
        ----------
        filename : str
            The path and the name of the output file.
        
        Returns
        -------
        None
        """
        self.df.to_csv(filename, sep='\t', header=True, index=True)
        return


class CDAPTable(MSTable):
    """Base class for all CDAP table types. Please see CDAPiTraqTable, CDAPPrecursorAreaTable.
    """ + BASE_ATTR
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def check_data(self):
        """Method containing checks that the data and the columns used
        for the processing exist. Checks for 'Gene' columns and the presence
        of multiple rows.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        assert len(self.df.index) > 0
        assert 'Gene' in self.df.columns
        return
    
    def clean_rows(self):
        """Method that removes extraneous rows. For CDAP experiments, rows with 'Mean',
        'Median', and 'StdDev' are mixed in with genes in the 'Gene' column.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        self.df = self.df[~self.df['Gene'].str.contains('Mean|Median|StdDev', na=True)]
        return
    
    def format_genes(self, ptm_type):
        """Method that rewrites the gene column and sets it as the index. First checks
        to see if the table is a PTM by checking for the associated name (only phosphosites
        supported at this time for CDAP). PTMs are built using the format
        
        GENE|GENE_<ptm_type><amino acid><amino acid site>
        
        Multi-site PTMs are chained with underscores. Otherwise, the gene is
        formatted as
        
        GENE|GENE
        
        Parameters
        ----------
        ptm_type : str
            A string indicating the PTM type, i.e. 'P'
        
        Returns
        -------
        None
        """
        ptm_name = {
            'P' : 'Phosphosite',
        }
        
        if ptm_type and ptm_name[ptm_type] in self.df.columns:
            anno = []
            for g, ph in zip(self.df['Gene'], self.df[ptm_name[ptm_type]]):
                ptm_site = ph.split(':')[1]
                ptm_build = re.sub('[A-Z]+', lambda m: '_'+m.group(0), ptm_site.upper())[1:]
                anno.append('{0}|{0}_{1}{2}'.format(g, ptm_type, ptm_build))
        else:
            anno = ['{0}|{0}'.format(g) for g in self.df['Gene']]
        self.df.index = anno
        self.df.index.names = ['Hugo_Symbol']
        return


class CDAPiTraqTable(CDAPTable):
    """Class for CDAP experiments for proteome and PTM quantification using iTRAQ.
    """ + BASE_ATTR
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class CDAPPrecursorAreaTable(CDAPTable):
    """Class for CDAP experiments for proteome and PTM quantification using intensities.
    """ + BASE_ATTR
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not self.use_ruler:
            self.log_transform()


class MaxQuantProteomeTable(MSTable):
    """Class for MaxQuant experiments for proteome MS quantification.
    """ + BASE_ATTR
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def check_data(self):
        """Method containing checks that the data and the columns used
        for the processing exist. Checks for the presence of multiple rows.
        Also checks for 'Protein IDs', 'Gene names', and 'Q-value' columns.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        assert len(self.df.index) > 0
        assert all(col in self.df.columns for col in ('Protein IDs', 
                                                      'Gene names', 
                                                      'Q-value'))
        return
    
    def clean_rows(self):
        """Method that removes extraneous rows. For MaxQuant proteome quantification, 
        it takes out rows with 'REV' or 'CON' and rows with a 'Q-value' over 0.05.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        self.df = self.df[~self.df['Protein IDs'].str.contains('REV|CON', na=True)]
        self.df = self.df[pd.notnull(self.df['Gene names'])]
        self.df = self.df[self.df['Q-value'] < 0.05]
        return
    
    def format_genes(self, ptm_type):
        """Method that rewrites the gene column and sets it as the index. Format
        used is 'GENE|GENE'.
        
        Parameters
        ----------
        ptm_type : str
            A string indicating the PTM type, i.e. 'P'
        
        Returns
        -------
        None
        """
        anno = ['{0}|{0}'.format(g.split(';')[0]) for g in self.df['Gene names']]
        self.df.index = anno
        self.df.index.names = ['Hugo_Symbol']
        return


class MaxQuantPTMTable(MSTable):
    """Class for MaxQuant experiments for PTM MS quantification.
    """ + BASE_ATTR
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def check_data(self):
        """Method containing checks that the data and the columns used
        for the processing exist. Checks for the presence of multiple rows.
        Also checks for 'Protein', 'Amino acid', 'Positions within proteins',
        'Gene names', and 'Localization prob' columns.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        assert len(self.df.index) > 0
        assert all(col in self.df.columns for col in ('Protein', 
                                                      'Amino acid', 
                                                      'Positions within proteins', 
                                                      'Gene names', 
                                                      'Localization prob'))
        return
    
    def clean_rows(self):
        """Method that removes extraneous rows. For MaxQuant PTM quantification, 
        it takes out rows with 'REV' or 'CON', rows with a null 'Gene names' entry,
        and rows where 'Localization prob' is less than 0.75.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        self.df = self.df[~self.df['Protein'].str.contains('REV|CON', na=True)]
        self.df = self.df[pd.notnull(self.df['Gene names'])]
        self.df = self.df[self.df['Localization prob'] > 0.75]
        return
    
    def format_genes(self, ptm_type):
        """Method that rewrites the gene column and sets it as the index. Phosphosites 
        are built using the format 'GENE|GENE_P<ptm_type><amino acid site>'.
        
        Parameters
        ----------
        ptm_type : str
            A string indicating the PTM type, i.e. 'P'
        
        Returns
        -------
        None
        """
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
    """A class for building metadata for MSTable derived classes.
    
    Parameters
    ----------
    cancer_id : str
        Cancer ID string supplied by cBioPortal (i.e. 'brca_tcga').
    
    prot_file : str
        Path and filename of the file generated by MSTable derived classes.
    """
    
    def __init__(self, cancer_id, prot_file):
        self.text = META_TEXT.format(cancer_id, prot_file)
    
    def __repr__(self):
        return self.text
    
    def write(self, filename):
        """Write finished metadata to file.
        
        Parameters
        ----------
        filename : str
            The path and the name of the output file.
        
        Returns
        -------
        None
        """
        with open(filename, 'w') as f:
            f.write(self.text)
        return








