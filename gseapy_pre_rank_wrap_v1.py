# gseapy_wrap_v1.py
'''
gene_sets=['GO_Biological_Process_2025','GO_Cellular_Component_2025','GO_Molecular_Function_2025',
           'KEGG_2021_Human','Reactome_2022', 'Reactome_Pathways_2024','WikiPathways_2024_Human',
            'MSigDB_Computational','MSigDB_Hallmark_2020','MSigDB_Oncogenic_Signatures',
            ]

run_gseapy_prerank_multiple_term_collections(
    data_csv_path=deseq2_results_file,
    set_of_term_collections=gene_sets,
    rank_metric_calc_flavor='signed_neg_log10_pvalue',
    l2fc_column_name= 'log2FoldChange',
    pvalue_column_name='padj',
    permutation_num=1000,
    threads=32,
    seed=6,
    file_base_name='HOM_over_WT',
    results_dir='./',
    GSEA_run_filename_tag='GSEA_prerank_run_multi',
    )
'''
 
# config

####################################


# module imports
import os
import numpy as np
import pandas as pd
import gseapy as gp
from typing import Mapping, Sequence, Union, Optional

import logging
logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)


#####################################




def process_input_table(
        input_df:pd.DataFrame = None,
        feature_name_column_name: str ='feature_name',
        make_feature_names_unique: bool = True,
        keep_first_occurrence_for_rows_with_multiple_feature_names: bool = False,
        explode_rows_with_multiple_feature_names: bool = False,
        rows_with_multiple_feature_names_delimiter: Optional[str] = None,
        rank_metric_calc_flavor='signed_neg_log10_pvalue', # 'use_column_key' or 'signed_neg_log10_pvalue'
        existing_rank_metric_column_name: Optional[str] | None = None,
        new_rank_metric_column_name: Optional[str]  = 'rank',
        l2fc_column_name: Optional[str] | None = 'log2FoldChange',
        pvalue_column_name: Optional[str] | None ='pvalue',
        save_table_output: bool = True,
        results_dir: str ='./',
        file_base_name: str ='datasetID',
        GSEA_run_filename_tag: str ='GSEA_prerank_run',
        )->pd.DataFrame:
    '''
    '''
    _df = input_df.copy()
    # handle rows with multiple feature names
    if keep_first_occurrence_for_rows_with_multiple_feature_names and rows_with_multiple_feature_names_delimiter is not None:
        _df[feature_name_column_name] = _df[feature_name_column_name].str.split(rows_with_multiple_feature_names_delimiter).str[0]
    # explode multiple feature names into multiple rows if needed
    if explode_rows_with_multiple_feature_names and rows_with_multiple_feature_names_delimiter is not None:
        _df[feature_name_column_name] = _df[feature_name_column_name].str.split(rows_with_multiple_feature_names_delimiter)
        _df = _df.explode(feature_name_column_name, ignore_index=True)
    # calculate the ranking metric
    if rank_metric_calc_flavor=='signed_neg_log10_pvalue': 
        _df = _df.dropna(subset=[feature_name_column_name, l2fc_column_name, pvalue_column_name])
        _df['neg_log10_pvalue'] = -np.log10(_df[pvalue_column_name].replace(0, np.nextafter(0, 1)))
        _df[new_rank_metric_column_name]= _df[l2fc_column_name].apply(np.sign) * _df['neg_log10_pvalue']
        ranking_series = _df.set_index(feature_name_column_name)[new_rank_metric_column_name].sort_values(ascending=False)
    elif rank_metric_calc_flavor=='use_column_key':
        if existing_rank_metric_column_name is None:
            raise ValueError("existing_rank_metric_column_name must be provided when using 'use_column_key' flavor.")
        ranking_series = _df.set_index(feature_name_column_name)[existing_rank_metric_column_name].sort_values(ascending=False)
    else:
        raise ValueError(f"Unknown rank_metric_calc_flavor: {rank_metric_calc_flavor}")
    # make the rnk_df from the ranking series
    rnk_df = ranking_series.reset_index().rename(columns={ranking_series.index.name:'feature_name', ranking_series.name:new_rank_metric_column_name})
    # make feature names unique by keeping the highest absolute values rank value
    if make_feature_names_unique:
        rnk_df['absolute_value_rnk'] = rnk_df[new_rank_metric_column_name].abs()
        rnk_df = rnk_df.sort_values(by=['absolute_value_rnk'], ascending=False).drop_duplicates(subset=['feature_name'], keep='first').drop(columns=['absolute_value_rnk'])
        # re-sort by rank metric
        rnk_df = rnk_df.sort_values(by=[new_rank_metric_column_name], ascending=False)
    # save the rnk_df if needed
    if save_table_output:
        out_file = os.path.join(results_dir, f"{file_base_name}.{GSEA_run_filename_tag}.rnk_df.csv")
        rnk_df.to_csv(out_file, index=False)
        LOGGER.info(f"Ranked DataFrame saved at {out_file}")
    return rnk_df


def run_single_gsea_prerank(
        rnk_df,
        single_gene_set_term_collection,
        collection_name: Optional[str] | None = None,
        min_size=10,
        max_size=1000,
        permutation_num=100,
        threads=32,
        seed=6,
        verbose=True,
        **kwargs,
        ):
    single_pre_res = gp.prerank(
        rnk=rnk_df,
        gene_sets=single_gene_set_term_collection,
        min_size=min_size,
        max_size=max_size,
        threads=threads,
        permutation_num=permutation_num,
        #outdir='./',
        seed=seed,
        verbose=True,
    )
    # Attach term collection information to each result in the dictionary.
    for key, value in single_pre_res.results.items():
        if isinstance(value, dict):
            single_pre_res.results[key]['TermCollection'] = collection_name if collection_name is not None else single_gene_set_term_collection
    # **Update the res2d DataFrame to include 'TermCollection'.**
    #single_pre_res.res2d['TermCollection'] = collection_name if collection_name is not None else single_gene_set_term_collection
    return single_pre_res



def make_results_dataframe(
        single_pre_res,
        collection_name: Optional[str] | None = None,
        results_table_min_genes_matched: int = 10,
        save_table_output: bool = True,
        results_dir: str ='./',
        file_base_name: str ='datasetID',
        GSEA_run_filename_tag: str ='GSEA_prerank_run',
        ):
    """Convert the pre_res.results into a tidy DataFrame and save if needed.
    The 'collection_name' parameter is used to annotate the results.
    returns: pd.DataFrame with columns 'TermCollection', 'Term', 'name', 'es', 'nes', 'pval', 'fdr', 'fwerp'

    """
    if collection_name is None:
        collection_name = getattr(single_pre_res, 'gene_sets', 'Unknown_TermCollection')
    try:
        df = pd.DataFrame(single_pre_res.results).T
        # Filter for gene sets with enough matched genes.
        df = df[df['matched_genes'].str.split(';').str.len() >= results_table_min_genes_matched].copy()
        # Add custom columns based on Term and collection.
        df['Term'] = df.index
        df['TermCollection'] = collection_name
        df['TermCollection_Term'] = collection_name + '__' + df['Term']
        if 'tag %' in df.columns:
            df['leadhits_total'] = df['tag %'].apply(lambda x: x.split('/')[0] + '_' + x.split('/')[1] if isinstance(x, str) else x)
        # Rearrange columns.
        front_cols = ['TermCollection', 'Term', 'name', 'es', 'nes', 'pval', 'fdr', 'fwerp']
        if 'leadhits_total' in df.columns:
            front_cols.append('leadhits_total')
        remaining = [col for col in df.columns if col not in front_cols]
        df = df[front_cols + remaining]
        # Ensure numeric columns are floats.
        for col in ['es', 'nes', 'pval', 'fdr', 'fwerp']:
            df[col] = df[col].astype(float)
        df = df.sort_values(by=['TermCollection', 'fdr'], ascending=[True, True])
        if save_table_output:
            out_file = os.path.join(results_dir, f"{file_base_name}.{GSEA_run_filename_tag}.{collection_name}.csv")
            df.drop(columns=['RES'], inplace=True, errors='ignore')
            df.to_csv(out_file, index=False, )# quoting=csv.QUOTE_MINIMAL
            LOGGER.info(f"Results saved for {collection_name} at {out_file}")
    except Exception as e:
        LOGGER.info(f"Error making results DataFrame for {collection_name}: {e}")
    return df


def run_gseapy_prerank_multiple_term_collections(
        ## for input data
        data_df: Optional[pd.DataFrame] = None,
        data_csv_path: Optional[Union[str, os.PathLike]] = None,
        ##   process_input_table() # for preparing the ranking metric
        feature_name_column_name: str ='feature_name',
        make_feature_names_unique: bool = True,
        keep_first_occurrence_for_rows_with_multiple_feature_names: bool = False,
        explode_rows_with_multiple_feature_names: bool = False,
        rows_with_multiple_feature_names_delimiter: Optional[str] = None,
        rank_metric_calc_flavor: str ='signed_neg_log10_pvalue', # 'use_column_key' or 'signed_neg_log10_pvalue'
        existing_rank_metric_column_name: Optional[str] | None = None,
        new_rank_metric_column_name: Optional[str] | None = 'rank',
        l2fc_column_name: Optional[str] | None = 'log2FoldChange',
        pvalue_column_name: Optional[str] | None ='padj',
        ## for GSEApy prerank run_single_gsea_prerank() and make_results_dataframe() and main loop through term families
        set_of_term_collections: Union[str, Sequence[str],None] = None,
        min_size: int =10,
        max_size: int =1000,
        permutation_num: int = 100,
        results_table_min_genes_matched: int = 0,
        threads: int = 4,
        seed: int = 6,
        verbose: bool = True,
        save_table_output: bool = True,
        results_dir: str ='./',
        file_base_name: str ='datasetID',
        GSEA_run_filename_tag: str ='GSEA_prerank_run',
        )-> None:
        '''Run GSEApy prerank analysis on the given data.
        Either data_df or data_csv_path must be provided.
        '''   
        ## input verification and loading  
        if data_df is None and data_csv_path is None:
            raise ValueError("Either data_df or data_csv_path must be provided.")
        if data_df is None:
            data_df = pd.read_csv(data_csv_path)
        if set_of_term_collections is None:
            raise ValueError("set_of_term_collections must be provided.")
        if isinstance(set_of_term_collections, str):
            set_of_term_collections = [set_of_term_collections]
        # log all parameters for process input table
        LOGGER.info(f"Running GSEApy prerank with parameters:\n"
                    f"feature_name_column_name: {feature_name_column_name}\n"
                    f"make_feature_names_unique: {make_feature_names_unique}\n"
                    f"keep_first_occurrence_for_rows_with_multiple_feature_names: {keep_first_occurrence_for_rows_with_multiple_feature_names}\n"
                    f"explode_rows_with_multiple_feature_names: {explode_rows_with_multiple_feature_names}\n"
                    f"rows_with_multiple_feature_names_delimiter: {rows_with_multiple_feature_names_delimiter}\n"
                    f"rank_metric_calc_flavor: {rank_metric_calc_flavor}\n"
                    f"existing_rank_metric_column_name: {existing_rank_metric_column_name}\n"
                    f"new_rank_metric_column_name: {new_rank_metric_column_name}\n"
                    f"l2fc_column_name: {l2fc_column_name}\n"
                    f"pvalue_column_name: {pvalue_column_name}\n"
                    f"save_table_output: {save_table_output}\n"
                    f"results_dir: {results_dir}\n"
                    f"file_base_name: {file_base_name}\n"
                    f"GSEA_run_filename_tag: {GSEA_run_filename_tag}\n"
                    )
        
        rnk_df = process_input_table(
            input_df=data_df,
            feature_name_column_name=feature_name_column_name,
            make_feature_names_unique=make_feature_names_unique,
            keep_first_occurrence_for_rows_with_multiple_feature_names=keep_first_occurrence_for_rows_with_multiple_feature_names,
            explode_rows_with_multiple_feature_names=explode_rows_with_multiple_feature_names,
            rows_with_multiple_feature_names_delimiter=rows_with_multiple_feature_names_delimiter,
            rank_metric_calc_flavor=rank_metric_calc_flavor,
            existing_rank_metric_column_name=existing_rank_metric_column_name,
            new_rank_metric_column_name=new_rank_metric_column_name,
            l2fc_column_name=l2fc_column_name,
            pvalue_column_name=pvalue_column_name,
            save_table_output=save_table_output,
            results_dir=results_dir,
            file_base_name=file_base_name,
            GSEA_run_filename_tag=GSEA_run_filename_tag,
        )

        for single_term_collection in set_of_term_collections:
            collection_name = None
            LOGGER.info(f"\nRunning GSEApy prerank for term collection: {single_term_collection}")
            # log all parameters for run_single_gsea_prerank
            LOGGER.info(f"Parameters for run_single_gsea_prerank:\n"
                        f"rnk_df shape: {rnk_df.shape}\n"
                        f"single_term_collection: {single_term_collection}\n"
                        f"collection_name: {collection_name}\n"
                        f"permutation_num: {permutation_num}\n"
                        f"min_size: {min_size}\n"
                        f"max_size: {max_size}\n"
                        f"threads: {threads}\n"
                        f"seed: {seed}\n"
                        f"verbose: {verbose}\n"
                        )
            single_pre_res = run_single_gsea_prerank(
                rnk_df,
                single_term_collection,
                collection_name,
                permutation_num=permutation_num,
                min_size=min_size,
                max_size=max_size,
                threads=threads,
                seed=seed,
                verbose=verbose,
            )
            # log all parameters for make_results_dataframe
            LOGGER.info(f"Parameters for make_results_dataframe:\n"
                        f"collection_name: {collection_name}\n"
                        f"results_table_min_genes_matched: {results_table_min_genes_matched}\n"
                        f"save_table_output: {save_table_output}\n"
                        f"results_dir: {results_dir}\n"
                        f"file_base_name: {file_base_name}\n"
                        f"GSEA_run_filename_tag: {GSEA_run_filename_tag}\n"
                        )
            results_df = make_results_dataframe(
                single_pre_res, collection_name,
                results_table_min_genes_matched=results_table_min_genes_matched,
                save_table_output=save_table_output,
                results_dir=results_dir,
                file_base_name=file_base_name,
                GSEA_run_filename_tag=GSEA_run_filename_tag,)
            LOGGER.info(f"Processed term collection: {single_term_collection} with {results_df.shape[0]} significant terms with at least {results_table_min_genes_matched} matched genes.")
        # end of main loop through term collections
        LOGGER.info("GSEApy prerank analysis completed for all term collections.")
        return 


 