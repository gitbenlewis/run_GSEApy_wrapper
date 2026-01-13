#!/usr/bin/env python3
"""script doc string."""
# /home/ubuntu/projects/gitbenlewis/run_GSEApy_wrapper/examples/GSE68719/scripts/make_gsea_tables.py


# Adds the project root (two levels up) to the search path
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))


####################################
import sys
import os
from pathlib import Path
from typing import Mapping, Sequence, Union, Optional
import pandas as pd
import gseapy as gp
import numpy as np
import glob
import logging
from datetime import datetime
import gseapy.utils as gutils
import gseapy.base as gbase
import gseapy.biomart as gbiomart
import gseapy.enrichr as genrich

import gseapy_pre_rank_wrap
from gseapy_pre_rank_wrap import (
    run_gseapy_prerank_multiple_term_collections,
)
from gseapy_pre_rank_wrap import (
    make_gmts_from_enrichr_libraries,
    make_gmts_from_Msigdb,
    basic_set_of_enrichr_term_collections,
    expansion_set_of_enrichr_term_collections,
    all_enrichr_libraries_2025_sin_archive,
)

####################################



# config
####################################

from pathlib import Path
import yaml

REPO_CONFIG_YAML_PATH = Path(__file__).resolve().parent.parent / "config" / "config.yaml"
with REPO_CONFIG_YAML_PATH.open() as f:
    CFG = yaml.safe_load(f)


REPO_ROOT = Path(__file__).resolve().parent.parent
# GSEApy parameters
# general params
GSEA_RESULTS_DIR = REPO_ROOT / 'results' / 'GSEApy'
#DEFAULT_GSEA_RESULTS_DIR = REPO_ROOT / 'results' / 'GSEApy'
os.makedirs(GSEA_RESULTS_DIR, exist_ok=True)

GSEA_FEATURE_NAME_COLUMN_NAME= CFG['GSEApy_parameters_defaults']['feature_name_column_name']#'gene_name'
GSEA_PERMUTATION_NUM=  CFG['GSEApy_parameters_defaults']['gsea_permutation_num'] # 10000
GSEA_MIN_SIZE = CFG['GSEApy_parameters_defaults']['min_size'] # 10
GSEA_MAX_SIZE = CFG['GSEApy_parameters_defaults']['max_size'] # 3000
GSEA_RESULTS_TABLE_MIN_GENES_MATCHED = CFG['GSEApy_parameters_defaults']['results_table_min_genes_matched']# 0
GSEA_THREADS = CFG['GSEApy_parameters_defaults']['threads']# 30
GSEA_SEED = CFG['GSEApy_parameters_defaults']['seed']# 6  
GSEA_RUN_FILENAME_TAG= CFG['GSEApy_parameters_defaults']['GSEA_run_filename_tag']#'GSEA'
GSEA_RANK_METRIC_CALC_FLAVOR= CFG['GSEApy_parameters_defaults']['rank_metric_calc_flavor']#'signed_neg_log10_pvalue'
GSEA_L2FC_COLUMN_NAME= CFG['GSEApy_parameters_defaults']['l2fc_column_name']#'log2FoldChange'
GSEA_PVALUE_COLUMN_NAME = CFG['GSEApy_parameters_defaults']['pvalue_column_name']#'pvalue'
GENE_TERM_COLLECTIONS_2_USE= CFG['GSEApy_parameters_defaults']['gene_term_collections_2_use']

GENE_TERM_COLLECTIONS_2_USE=[
    "GO_Biological_Process_2025",
    "GO_Cellular_Component_2025",
    "GO_Molecular_Function_2025",
    "KEGG_2021_Human",
    "Reactome_2022",
    "Reactome_Pathways_2024",
    "WikiPathways_2024_Human",
    "WikiPathways_2024_Mouse",
    "MSigDB_Computational",
    "MSigDB_Hallmark_2020",
    "CellMarker_2024",
    "ARCHS4_TFs_Coexp",
    "ARCHS4_Tissues",
    "Disease_Perturbations_from_GEO_down",
    "Disease_Perturbations_from_GEO_up",
    'c1.all',
    'c2.all',
    'c2.cgp',
    'c2.cp.biocarta',
    'c2.cp.kegg_legacy',
    'c2.cp.kegg_medicus', 
    'c2.cp.pid',
    'c2.cp.reactome', 
    'c2.cp',
    'c2.cp.wikipathways',
    'c3.all', 
    'c3.mir.mir_legacy', 
    'c3.mir.mirdb', 
    'c3.mir', 
    'c3.tft.gtrd', 
    'c3.tft.tft_legacy',
    'c3.tft',
    'c4.3ca',
    'c4.all', 
    'c4.cgn', 
    'c4.cm',
    'c5.all', 
    'c5.go.bp',
    'c5.go.cc',
    'c5.go.mf', 
    'c5.go', 
    'c5.hpo', 
    'c6.all', 
    'c7.all', 
    'c7.immunesigdb',
    'c7.vax', 
    'c8.all',
    'h.all',
    'msigdb'
]

GENE_TERM_COLLECTIONS_2_USE=[
'GO_Biological_Process_2025','GO_Cellular_Component_2025',"GO_Molecular_Function_2025","Reactome_Pathways_2024",
 'c5.go.bp','c5.go.cc','c5.go.mf','c2.cp.reactome', 
 ]




#### start #### log file setup
# ---------- logging setup ----------
LOG_DIR = GSEA_RESULTS_DIR / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILENAME=f"make_gsea_tables_{datetime.now():%Y%m%d_%H%M%S}.log"
RESULTS_LOG_FILE  = LOG_DIR / LOG_FILENAME
SCRIPT_LOG_FILE = Path(__file__).resolve().parent / LOG_FILENAME
root = logging.getLogger()
root.setLevel(logging.INFO)
if not any(isinstance(h, logging.FileHandler) for h in root.handlers): # add file handler once
    for log_path in (RESULTS_LOG_FILE, SCRIPT_LOG_FILE):
        fh = logging.FileHandler(log_path)
        fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
        fh.setLevel(logging.INFO)
        root.addHandler(fh)
LOGGER = logging.getLogger(__name__)
LOGGER.info("Logging to %s", RESULTS_LOG_FILE)
LOGGER.info("Logging to %s", SCRIPT_LOG_FILE)
# #) make gseapy loggers use our handlers
# file_handler = next(h for h in root.handlers if isinstance(h, logging.FileHandler))
file_handlers = []
_seen_log_paths = set()
for handler in root.handlers:
    if isinstance(handler, logging.FileHandler):
        log_path = getattr(handler, "baseFilename", None)
        if log_path in _seen_log_paths:
            continue
        _seen_log_paths.add(log_path)
        file_handlers.append(handler)
def log_init_with_file(name, log_level=logging.INFO, filename=None):
    logger = logging.getLogger(name)
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)
    for handler in file_handlers:
        logger.addHandler(handler)
    logger.addHandler(logging.StreamHandler())
    logger.propagate = True
    return logger
gutils.log_init = log_init_with_file
gbase.log_init = log_init_with_file
gbiomart.log_init = log_init_with_file
genrich.log_init = log_init_with_file



# -----END----- logging setup ---------- 

####### GMT prep #########################
SET_OF_GMT_FILE_PATHS = []
LIST_OF_TERM_COLLECTION_NAMES= []
## loop over gmt file dirs
for gmt_file_dir in CFG['gene_term_collections_gmt_file_directory_list']:
    gmt_dir_path = Path(CFG['gene_term_collections_gmt_file_directory_list'][gmt_file_dir])
    LOGGER.info(f"Checking GMT file dir path: {gmt_dir_path}, exists: {gmt_dir_path.exists()}")
    gmt_file_paths= glob.glob(str(gmt_dir_path / '*.gmt'))
    LOGGER.info(f"Found {len(gmt_file_paths)} GMT files in {gmt_dir_path}: {gmt_file_paths}")
    kept_gmt_file_paths= [gmt_file_path for gmt_file_path in gmt_file_paths
                                if os.path.basename(gmt_file_path).split('.gmt')[0] in GENE_TERM_COLLECTIONS_2_USE]
    LOGGER.info(f"Keeping {len(kept_gmt_file_paths)} GMT files for {gmt_file_dir} term collections: {kept_gmt_file_paths}")
    SET_OF_GMT_FILE_PATHS += kept_gmt_file_paths
    ## #) now make LIST_OF_TERM_COLLECTION_NAMES
    kept_term_collection_names= [
        os.path.basename(gmt_file_path).split('.gmt')[0] for gmt_file_path in kept_gmt_file_paths ]
    kept_term_collection_names=[name.replace('.symbols','') for name in kept_term_collection_names] 
    LOGGER.info(f"Keeping {len(kept_term_collection_names)} term collection names for {gmt_file_dir} term collections: {kept_term_collection_names}")
    LIST_OF_TERM_COLLECTION_NAMES += kept_term_collection_names
####### end GMT prep #########################



####### end GMT prep #########################


#SET_OF_ENRICHR_TERM_COLLECTIONS=basic_set_of_enrichr_term_collections
gseapy_runs=CFG['gseapy_runs']
for run in gseapy_runs:
    # 
    run_name=run['gseapy_run_name']
    da_run_name=run['DA_run_name']
    da_table_path=run['DA_table']
    data_csv_path_sepstr=getattr(run,'DA_table_sepstr',',')

    #
    file_base_name = os.path.basename(da_table_path).replace('.deseq2.results.csv','')
    full_output_dir = os.path.join(GSEA_RESULTS_DIR,da_run_name, run_name)
    # 
    os.makedirs(full_output_dir, exist_ok=True)
    LOGGER.info(f"Processing file: {da_table_path} with base name: {file_base_name}, and output dir: {full_output_dir}")
    run_gseapy_prerank_multiple_term_collections(
        data_csv_path=da_table_path,
        data_csv_path_sepstr=data_csv_path_sepstr,
        feature_name_column_name=GSEA_FEATURE_NAME_COLUMN_NAME,
        list_of_term_collections=SET_OF_GMT_FILE_PATHS,
        list_of_term_collection_names=LIST_OF_TERM_COLLECTION_NAMES,
        rank_metric_calc_flavor=GSEA_RANK_METRIC_CALC_FLAVOR,
        l2fc_column_name= GSEA_L2FC_COLUMN_NAME,
        pvalue_column_name=GSEA_PVALUE_COLUMN_NAME,
        permutation_num=GSEA_PERMUTATION_NUM,
        min_size=GSEA_MIN_SIZE,
        max_size=GSEA_MAX_SIZE,
        results_table_min_genes_matched=GSEA_RESULTS_TABLE_MIN_GENES_MATCHED,
        threads=GSEA_THREADS,
        seed=GSEA_SEED,
        file_base_name=file_base_name,
        results_dir=full_output_dir,
        GSEA_run_filename_tag=GSEA_RUN_FILENAME_TAG,
        )
    


LOGGER.info(f"make_gsea_tables.py All done!")
