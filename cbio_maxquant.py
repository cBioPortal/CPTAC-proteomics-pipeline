import argparse
import re
import numpy as np
import pandas as pd


def load_transform_proteome(fname, pipeline, 
    sample_regex='[A-Za-z0-9\s]+'):
    
    df = pd.read_csv(fname, sep='\t', header=0)
    df = df[~df['Protein IDs'].str.contains('REV|CON', na=True)]
    df = df[pd.notnull(df['Gene names'])]
    df = df[df['Q-value'] < 0.05]
    
    if pipeline == 'unlabeled':
        samp_cols = [col for col in df.columns if re.match('Intensity ' + sample_regex, col)]
    elif pipeline == 'labeled':
        samp_cols = [col for col in df.columns if re.match('Reporter intensity ' + sample_regex, col)]
    
    genes = [g.split(';')[0] for g in df['Gene names']]
    df = df[samp_cols]
    df.columns = [re.search(sample_regex, col).group() for col in df.columns]
    df.index = [g+'|'+g for g in genes]
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df


def load_transform_ptm(fname, pipeline, ptm_prefix,
    sample_regex='[A-Za-z0-9\s]+'):
    
    df = pd.read_csv(fname, sep='\t', header=0)
    df = df[~df['Protein'].str.contains('REV|CON', na=True)]
    df = df[pd.notnull(df['Gene names'])]
    df = df[df['Localization prob'] > 0.75]
    
    if pipeline == 'unlabeled':
        samp_cols = [col for col in df.columns if re.match('Intensity ' + sample_regex, col)]
    elif pipeline == 'labeled':
        samp_cols = [col for col in df.columns if re.match('Reporter intensity ' + sample_regex, col)]
    
    prot_pos = [p.split(';')[0] for p in df['Positions within proteins']]
    aa_col = df['Amino acid'].tolist()
    genes = [g.split(';')[0] for g in df['Gene names']]
    
    df = df[samp_cols]
    df.columns = [re.search(sample_regex, col).group() for col in df.columns]
    
    new_ids = []
    for g, aa, pos in zip(genes, aa_col, prot_pos):
        new_id = '{0}|{0}_{1}{2}{3}'.format(g, ptm_prefix, aa, pos)
        new_ids.append(new_id)
    
    df.index = new_ids
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df


def average_duplicates(df):
    nondup_df = df[df.columns[~df.columns.duplicated(keep=False)]]
    uniq_dup_cols = list(set(df.columns[df.columns.duplicated(keep=False)]))
    dedup_df = pd.DataFrame(index=nondup_df.index, columns=uniq_dup_cols)
    for col in uniq_dup_cols:
        dedup_df[col] = df[col].mean(axis=1)
    df = pd.concat((nondup_df, dedup_df), axis=1)
    return df


def combine_data(df_list):
    df = pd.concat(df_list, axis=1)
    df = average_duplicates(df)
    df.index.names = ['Hugo_Symbol']
    return df


def write_meta(cancer_id, output_file):
    blurb = """cancer_study_identifier: {0}
genetic_alteration_type: PROTEIN_LEVEL
datatype: Z-SCORE
stable_id: protein_quantification
profile_description: Protein Quantification (Mass Spec)
show_profile_in_analysis_tab: true
profile_name: Protein levels (mass spectrometry by CPTAC)
data_filename: {1}""".format(cancer_id, output_file.split('/')[-1])
    return blurb


def main(args):
    prot_data = []
    for prot_fname in args.proteome_files.split(','):
        prot_df = load_transform_proteome(prot_fname, args.pipeline, sample_regex=args.sample_regex)
        prot_data.append(prot_df)
    ptm_params = (args.ptm_files, args.pipeline, args.ptm_prefixes)
    if all(ptm_params):
        ptm_data = []
        for ptm_fname, ptm_prefix in zip(args.ptm_files.split(','), args.ptm_prefixes.split(',')):
            ptm_df = load_transform_ptm(ptm_fname, args.pipeline, ptm_prefix, sample_regex=args.sample_regex)
            ptm_data.append(ptm_df)
        total_prot_df = combine_data(prot_data)
        total_ptm_df = combine_data(ptm_data)
        comb_df = pd.concat((total_prot_df, total_ptm_df), axis=0)
    elif any(ptm_params) and not all(ptm_params):
        raise RuntimeError('Attempting to process PTMs without all arguments set (--ptm-files, --pipeline, --ptm-prefixes).')
    else:
        comb_df = combine_data(prot_data)
    comb_df.to_csv(args.output_file, sep='\t', header=True, index=True)
    with open(args.meta_file, 'w') as f:
        f.write(write_meta(args.cancer_id, args.output_file))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This script converts MaxQuant quantitation data and PTM data into importable text files for cBioPortal. It will also output a meta file for use with cbioportalImporter.py. Use for one tissue type (i.e. breast cancer) at a time.')
    
    parser.add_argument('--proteome-files', 
        metavar='<filename_1,filename_2,filename_n>', 
        help='comma-separated list of paths to proteinGroup.txt files', 
        required=True)
    parser.add_argument('--pipeline', 
        choices=('labeled', 'unlabeled'), 
        help='processing pipeline (i.e. use "unlabeled" for LFQ/iBAQ and "labeled" for TMT)', 
        required=True)
    parser.add_argument('--sample-regex', 
        help='regular expression for sample names (i.e. "[A-Za-z0-9\s]+"', 
        default='[A-Za-z0-9\s]+',
        required=True)
    parser.add_argument('--ptm-files', 
        metavar='<filename_1,filename_2,filename_n>', 
        help='comma-separated list of PTM filenames', 
        required=False)
    parser.add_argument('--ptm-prefixes', 
        metavar='<prefix_1,prefix_2,prefix_n>',
        help='a comma-separated list of letters indicating the class of PTM of each PTM file, where p=phosphosite, etc. (i.e. p,p,p)', 
        required=False)
    parser.add_argument('--output-file', 
        metavar='<output-file>', 
        help='path and name given to output', 
        required=True)
    parser.add_argument('--meta-file', 
        metavar='<meta-file>', 
        help='path and name given to meta file', 
        required=True)
    parser.add_argument('--cancer-id', 
        metavar='<cancer-id>', 
        help='an official cancer study ID, i.e. "brca_tcga_pub"', 
        required=True)
    
    args = parser.parse_args()
    main(args)

