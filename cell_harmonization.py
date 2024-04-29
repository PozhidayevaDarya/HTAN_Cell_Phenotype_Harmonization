import os
import synapseclient
from google.cloud import bigquery

import scanpy as sc
import pandas as pd
import re
import scipy
from scipy.sparse import csr_matrix
import celltypist

from celltypist import models
import glob
import subprocess

def get_ensembl_synonyms(ensembl_list):
    '''
    Get symbols for ensembl Ids 
    '''
    my_gene = mygene.MyGeneInfo()
    gene_list = list()
    #clean ensemble if needed
    for ensembl in ensembl_list:
        #remove .num
        ensembl_clean = ensembl.split('.', 1)[0]
        gene_list.append(ensembl_clean)

    saved_ensembl_ids = pd.DataFrame()
    saved_ensembl_ids['clean_ensembl_ids'] = gene_list
    saved_ensembl_ids['original_ensembl_ids'] = ensembl_list
    
    # query full list and save
    gene_synonym = my_gene.querymany(gene_list, scopes='ensembl.gene', fields='symbol', species='human')
    gene_symbols = list()
    query_list = list()
    for item in gene_synonym:
        if "symbol" in item and "query" in item:
            query_list.append(item['query'])
            gene_symbols.append(item['symbol'])
    query_return = pd.DataFrame()
    query_return['queried_ensembl'] = query_list
    query_return['returned_symbol'] = gene_symbols
    
    translated_gene_data = pd.merge(saved_ensembl_ids, query_return, left_on='clean_ensembl_ids', right_on='queried_ensembl')
    return translated_gene_data


def run_adata_processing(adata_object):
    '''
    Run Scanpy data processing on Ann Data object and then predict cell type based on low resolution immune model
    '''
    #Control cutoffs
    min_genes = 500
    min_cells = 3
    min_counts = 20
    cell_typist_model = 'Immune_All_Low.pkl'
    cell_typist_probability_threshold = 0.5
    # 
    sc.pp.calculate_qc_metrics(adata_object)
    sc.pp.filter_cells(adata_object, min_genes = min_genes)
    sc.pp.filter_genes(adata_object, min_cells = min_cells)
    sc.pp.filter_genes(adata_object, min_counts = min_counts)
    sc.pp.normalize_total(adata_object, target_sum=1e4)
    sc.pp.log1p(adata_object)
    sc.tl.pca(adata_object, svd_solver='arpack')
    sc.pp.neighbors(adata_object, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_object)
    # Impltemenet cell typist and save the outputs to ann data object
    predictions=celltypist.annotate(adata_object, model = cell_typist_model, mode = 'best match', p_thres = cell_typist_probability_threshold, majority_voting = True)
    adata_object = predictions.to_adata()
    
    return adata_object

if __name__ == '__main__':   
    #Directory Managment
    current_directory = os.getcwd()
    dir_10x = 'Temp_Data/10xformat'
    dir_csv = 'Temp_Data/csv'
    dir_final = 'Harmonized_h5ad'
    dirs = [dir_10x,dir_csv, dir_final]
    for dir in dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)
        
    #syn = synapseclient.Synapse()
    #try:
        #syn.login()
    #except synapseclient.core.exceptions.SynapseNoCredentialsError:
    #    print("Please ensure 'username' and 'password'/'api_key' \
    #        are filled out in .synapseConfig.")
    #except synapseclient.core.exceptions.SynapseAuthenticationError:
    #    print("Please make sure the credentials in the .synapseConfig file are correct.")
        
    # instantiate google bigquery client
    client = bigquery.Client.from_service_account_json("htan-dcc-ccc1ca71c3a3.json")
    
    data_query = client.query("""
    SELECT DISTINCT T1.entityId, T1.HTAN_Center, T1.Filename, T1.File_Format, T1.scRNAseq_Workflow_Parameters_Description, T1.Cell_Total, T1.Matrix_Type, T1.Component, T1.HTAN_Data_File_ID, T1.Data_Release, T2.HTAN_Participant_ID, T2.HTAN_Assayed_Biospecimen_ID, 
    T3.Tissue_or_Organ_of_Origin, T4.Tumor_Tissue_Type  FROM `htan-dcc.combined_assays.ScRNA-seqLevel3` AS T1
    #
    INNER JOIN `htan-dcc.ISB_CGC_r4.upstream_ids` AS T2
    ON T1.HTAN_Data_File_ID = T2.HTAN_Data_File_ID
    #
    LEFT JOIN `htan-dcc.ISB_CGC_r4.Diagnosis` AS T3
    ON T2.HTAN_Participant_ID = T3.HTAN_Participant_ID
    #
    LEFT JOIN `htan-dcc.ISB_CGC_r4.Biospecimen` AS T4
    ON HTAN_Assayed_Biospecimen_ID = T4.HTAN_Biospecimen_ID
    #
    WHERE T1.Data_Release IS NOT NULL 
    AND (T1.Matrix_Type = "Raw Counts")
    AND T1.Data_Category = "Gene Expression"
    AND (T3.Tissue_or_Organ_of_Origin = "Colon NOS" OR T3.Tissue_or_Organ_of_Origin = "Lung NOS")
    """).result().to_dataframe()

    # Collapse files by biospecimen id
    subset_query = data_query[['entityId', 'HTAN_Assayed_Biospecimen_ID', 'Filename', 'File_Format']]
    files_by_biospecimen = subset_query.groupby('HTAN_Assayed_Biospecimen_ID', as_index=False).agg(list)

    # Process all the files by type and biospecimen id
    for ids, biospecimen, filename, fileformat in zip(files_by_biospecimen.entityId, files_by_biospecimen.HTAN_Assayed_Biospecimen_ID, files_by_biospecimen.Filename, files_by_biospecimen.File_Format):
        
        if 'gzip' in fileformat and all(x == fileformat[0] for x in fileformat) and len(fileformat) >= 3:
            if [v for v in filename if 'barcodes' in v] and [v for v in filename if 'features' in v]:
                for id in ids:
                    syn.get(id, downloadLocation="Temp_Data/10xformat")
            #ASSUMPTION
            file_listing = os.listdir("Temp_Data/10xformat")
            prefix = re.match("(.*?)_", file_listing[0]).group()
            current_adata = sc.read_10x_mtx("Temp_Data/10xformat", prefix = prefix)
            current_adata.layers["raw"] = current_adata.X.copy()
            run_adata_processing(current_adata)
            current_adata.write_h5ad(prefix + "_harmonized_cell.h5ad")
            for f in file_listing:
                os.remove("Temp_Data/10xformat/"+f)
        
        if 'csv' in fileformat and all(x == fileformat[0] for x in fileformat):
            for id in ids:
                syn.get(id, downloadLocation="Temp_Data/csvformat")
            file_listing = os.listdir("Temp_Data/csvformat")
            for file in file_listing:
                current_adata = pd.read_csv("Temp_Data/csvformat/" + file, header=0, index_col=0)
                ens_dict = get_ensembl_synonyms(current_adata.index)
                ens_dict = ens_dict.set_index('original_ensembl_ids')
                mat = pd.merge(current_adata, ens_dict, left_index=True, right_index=True)
                mat = mat.set_index('returned_symbol')
                mat = mat.drop(columns=['clean_ensembl_ids', 'queried_ensembl'])
                adata = sc.AnnData(mat.transpose())
                # convert the pandas table to a sparse matrix
                adata.X = sparse.csr_matrix(adata.X)
                adata.layers["raw"] = adata.X.copy()
                run_adata_processing(adata)
                adata.write_h5ad(file + "harmonized_cell.h5ad")
                os.remove("Temp_Data/csvformat/" + file )