import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

from typing import List
from functools import cache

@cache
def import_gene_mapping_table():
    mapping = pd.read_table('mart_export_grch38.txt')
    # Ensembl gene identifier mappings to gene names
    ensembl_gene_id_map = mapping.set_index('Gene stable ID').to_dict()['Gene name']

    # NCBI column has NaNs which turns the integers into floats, so drop NaNs, change floats --> int --> str, and create dict
    ncbi_gene_id_map = mapping.dropna(subset='NCBI gene (formerly Entrezgene) ID')
    ncbi_gene_id_map['NCBI gene (formerly Entrezgene) ID'] = ncbi_gene_id_map['NCBI gene (formerly Entrezgene) ID'].astype(int).astype(str)
    ncbi_gene_id_map = ncbi_gene_id_map.set_index('NCBI gene (formerly Entrezgene) ID').to_dict()['Gene name']

    # Combine the transcript IDs, gene IDs, and NCBI IDs dictionaries for all possible mappings
    mapping_dict = ensembl_gene_id_map | ncbi_gene_id_map 

    return mapping_dict


def harmonize_genes(matrix:pd.DataFrame):
    """
    matrix (pd.DataFrame): RNA-seq or Epigenetic feature matrix where the index is assigned as the gene identifiers for each respective dataset.
    """
    mapping_dict = import_gene_mapping_table()
    matrix['gene_name'] = matrix.index.astype(str).map(mapping_dict)
    matrix.dropna(subset='gene_name', inplace=True)
    matrix.set_index('gene_name', inplace=True)
    return matrix

def create_target_vector(rna_seq_df:pd.DataFrame, target_value:str="TPM"):
    """
    The target vector in this case is the accompanying RNA-seq data for each gene in the feature matrix.

    rna_seq_df (pd.DataFrame): RNA-seq dataframe
    target_value (str): the value used as the target variable to predict (i.e. TPM or FPKM)
    """
    rna_seq_df['gene_id_clean'] = rna_seq_df['gene_id'].apply(lambda x: x.split(".")[0]) # strip versions from ENSG identifiers
    rna_seq_df.set_index('gene_id_clean', inplace = True)
    rna_seq_df = harmonize_genes(rna_seq_df)
    y = np.log10(rna_seq_df[target_value]+0.001)
    return y

def create_feature_matrix(epigenetics_df:List[pd.DataFrame]):
    """
    epigenetics_df (List[pd.DataFrame]): dataframe containing tissue or cell line histone modification assay features from Roadmap Epigenomics.
    """
    # Take this as an example for now
    df_aggregated = pd.concat(epigenetics_df)

    # Get the number of peaks per gene region
    peak_count_df = df_aggregated.groupby(["geneId","annotation"],as_index=False).agg(
        peak_counts = ('annotation','size')
    )

    feature_matrix = peak_count_df.pivot_table(columns='annotation',index=['geneId'])['peak_counts']

    feature_matrix = harmonize_genes(feature_matrix)

    return feature_matrix

def create_train_test_matrices(epigenetics_df:List[pd.DataFrame], rna_seq_df:pd.DataFrame):
    feature_matrix = create_feature_matrix(epigenetics_df)
    target_vector = create_target_vector(rna_seq_df)
    data = pd.merge(feature_matrix,target_vector,left_index=True,right_index=True)
    X, y = data.iloc[:, :-1], data.iloc[:, -1]
    X.fillna(0, inplace = True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42) # may want to parameterize these values?
    return X_train, X_test, y_train, y_test

