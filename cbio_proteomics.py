import argparse
import re
import numpy as np
import pandas as pd


def load_transform_data(fname, pipeline, 
    sample_regex='([A-Z0-9]{2}\-[A-Z0-9]{4}\-[A-Z0-9]{2})[A-Z0-9\-]+'):
    
    df = pd.read_csv(fname, sep='\t', header=0, index_col=0)
    for row in ['Mean', 'Median', 'StdDev']:
        if row in df.index:
            df.drop([row], inplace=True)
    
    if pipeline == 'itraq':
        new_cols = [col for col in df.columns if re.match(sample_regex + ' Log Ratio', col)]
        df = df[new_cols]
    elif pipeline == 'precursor_area':
        new_cols = [col for col in df.columns if re.match(sample_regex + ' Area', col)]
        df = np.log2(df[new_cols])
    
    df.columns = ['TCGA-'+re.match(sample_regex, col).group(1) for col in df.columns]
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df


def annotate_gene(df):
    df.index = [g+'|'+g for g in df.index]
    return df


def annotate_ptm(df, anno_file, ptm_prefix):
    anno = pd.read_csv(anno_file, sep='\t', header=0, dtype='str')
    anno = anno.drop_duplicates('Protein')
    anno_dict = dict(zip(anno['Protein'], anno['Gene']))
    
    new_index = []
    for item in df.index:
        if_anno, if_ptm = False, False
        try:
            lookup_id = item.split('.')[0]
            anno_id = anno_dict[lookup_id]
            if_anno = True
        except (IndexError, KeyError):
            pass
        
        try:
            ptm_site = item.split(':')[1]
            ptm_build = ptm_prefix + re.sub('[A-Z]+', 
                lambda m: '_'+m.group(0), ptm_site.upper())[1:]
            if_ptm = True
        except (IndexError, KeyError):
            pass
        
        if if_anno and if_ptm:
            new_index.append('{0}|{0}_{1}'.format(anno_id, ptm_build))
        else:
            new_index.append('NA|NA')
        
    df.index = new_index
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
stable_id: ms_abundance
profile_description: Protein levels (mass spectrometry)
show_profile_in_analysis_tab: true
profile_name: Protein levels (mass spectrometry)
data_filename: {1}""".format(cancer_id, output_file.split('/')[-1])
    return blurb


def main(args):
    prot_data = []
    for prot_fname in args.proteome_files.split(','):
        prot_df = load_transform_data(prot_fname, args.proteome_pipeline)
        prot_df = annotate_gene(prot_df)
        prot_data.append(prot_df)
    ptm_params = (args.ptm_files, args.ptm_pipeline, args.ptm_prefixes)
    if all(ptm_params):
        ptm_data = []
        for ptm_fname, ptm_prefix in zip(args.ptm_files.split(','), args.ptm_prefixes.split(',')):
            ptm_df = load_transform_data(ptm_fname, args.ptm_pipeline)
            ptm_df = annotate_ptm(ptm_df, args.annotation, ptm_prefix)
            ptm_data.append(ptm_df)
        total_prot_df = combine_data(prot_data)
        total_ptm_df = combine_data(ptm_data)
        comb_df = pd.concat((total_prot_df, total_ptm_df), axis=0)
    elif any(ptm_params) and not all(ptm_params):
        raise RuntimeError('Attempting to process PTMs without all arguments set (--ptm-files, --ptm_pipeline, --ptm-prefixes).')
    else:
        comb_df = combine_data(prot_data)
    comb_df.to_csv(args.output_file, sep='\t', header=True, index=True)
    with open(args.meta_file, 'w') as f:
        f.write(write_meta(args.cancer_id, args.output_file))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This script converts CPTAC proteome abundance data and PTM data into importable text files for cBioPortal. It will also output a meta file for use with cbioportalImporter.py. Use for one tissue type (i.e. breast cancer) at a time.')
    
    parser.add_argument('--proteome-files', 
        metavar='<filename_1,filename_2,filename_n>', 
        help='comma-separated list of protein abundance filenames', 
        required=True)
    parser.add_argument('--proteome-pipeline', 
        choices=('itraq', 'precursor_area'), 
        help='processing pipeline (i.e. "itraq" for Breast and Ovarian, "precursor_area for Colon")', 
        required=True)
    parser.add_argument('--ptm-files', 
        metavar='<filename_1,filename_2,filename_n>', 
        help='comma-separated list of PTM filenames', 
        required=False)
    parser.add_argument('--ptm-pipeline', 
        choices=('itraq', 'precursor_area'), 
        help='processing pipeline (i.e. "itraq" for Breast and Ovarian, "precursor_area for Colon")', 
        required=False)
    parser.add_argument('--ptm-prefixes', 
        metavar='<prefix_1,prefix_2,prefix_n>',
        help='a comma-separated list of letters indicating the class of PTM of each PTM file, where p=phosphosite, etc. (i.e. p,p,p)', 
        required=False)
    parser.add_argument('--annotation', 
        metavar='<annotation-file>', 
        help='annotation file generated by process_refseq.py', 
        required=True)
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
