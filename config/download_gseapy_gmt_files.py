#!/usr/bin/env python3
"""Download Enrichr and MSigDB GMT files into data/ref."""
# config/download_gseapy_gmt_files.py


import sys
import os

# Adds the project root (two levels up) to the search path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))

from pathlib import Path
from datetime import datetime
import logging


from gseapy_pre_rank_wrap import (
    make_gmts_from_enrichr_libraries,
    make_gmts_from_Msigdb,
    convert_gmt_homologs_directory,
    basic_set_of_enrichr_term_collections,
    expansion_set_of_enrichr_term_collections,
    all_enrichr_libraries_2025_sin_archive,
    
)


REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_REF_DIR = REPO_ROOT / "data" / "ref"
ENRICHR_HU_DIR = DATA_REF_DIR / "gmt_files" / "enrichr_Hu"
MSIGDB_HU_DIR = DATA_REF_DIR / "gmt_files" / "msigdb_Hu"
MSIGDB_MM_DIR = DATA_REF_DIR / "gmt_files" / "msigdb_Mm"
# ortholog converted 
H2M_MAPPING_FILE = REPO_ROOT / "data" / "ref" / "h2m_agg.csv"
M2H_MAPPING_FILE = REPO_ROOT / "data" / "ref" / "m2h_agg.csv"
ENRICHR_HU_H2M_DIR = DATA_REF_DIR /"gmt_files" / "enrichr_Hu_h2m"
MSIGDB_HU_H2M_DIR = DATA_REF_DIR / "gmt_files" /"msigdb_Hu_h2m"
MSIGDB_MM_M2H_DIR = DATA_REF_DIR / "gmt_files" /"msigdb_Mm_m2h"


BASIC_SET_OF_TERM_COLLECTIONS = basic_set_of_enrichr_term_collections
EXPANSION_SET_OF_TERM_COLLECTIONS = expansion_set_of_enrichr_term_collections
ALL_SET_OF_TERM_COLLECTIONS = all_enrichr_libraries_2025_sin_archive    


LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"


def _setup_logging() -> None:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"download_gseapy_gmt_files_{timestamp}.log"
    gmt_log_path  = DATA_REF_DIR / "gmt_files" / log_filename
    script_log_path = Path(__file__).resolve().parent / log_filename
    log_paths = []
    for log_path in (gmt_log_path, script_log_path):
        if log_path not in log_paths:
            log_paths.append(log_path)

    for log_path in log_paths:
        log_path.parent.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO)
    root_logger = logging.getLogger()
    formatter = logging.Formatter(LOG_FORMAT)

    for log_path in log_paths:
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    for log_path in log_paths:
        root_logger.info("Logging to %s", log_path)


def main() -> None:
    _setup_logging()
    DATA_REF_DIR.mkdir(parents=True, exist_ok=True)
    make_gmts_from_enrichr_libraries(    list_of_term_collections=ALL_SET_OF_TERM_COLLECTIONS,
        outdir=str(ENRICHR_HU_DIR), organism="Human",)

    #make_gmts_from_enrichr_libraries(
    #    list_of_term_collections=BASIC_SET_OF_TERM_COLLECTIONS,
    #    outdir=str(ENRICHR_HU_DIR), organism="Human",
    #)
    #make_gmts_from_enrichr_libraries(
    #    list_of_term_collections=EXPANSION_SET_OF_TERM_COLLECTIONS,
    #    outdir=str(ENRICHR_HU_DIR), organism="Human",
    #)

    make_gmts_from_Msigdb(    dbver="2025.1.Hs",list_of_Msigdb_categories=None,
       outdir=str(MSIGDB_HU_DIR),)
    make_gmts_from_Msigdb(    dbver="2025.1.Mm", list_of_Msigdb_categories=None,
        outdir=str(MSIGDB_MM_DIR),)

    # ortholog converted H2M
    convert_gmt_homologs_directory(
        gmt_input_directory=str(ENRICHR_HU_DIR),
        gmt_output_directory=str(ENRICHR_HU_H2M_DIR),
        covert2organism="Mouse",
        homolog_mapping_h2m_file=str(H2M_MAPPING_FILE),
        #homolog_mapping_m2h_file=str(M2H_MAPPING_FILE),
        source_index_col="gene_name",
        #mouse_gene_name_col="mouse_gene_name",
        mouse_gene_name_col="mouse_gene_name",
    )
    # ortholog converted H2M
    convert_gmt_homologs_directory(
        gmt_input_directory=str(MSIGDB_HU_DIR),
        gmt_output_directory=str(MSIGDB_HU_H2M_DIR),
        covert2organism="Mouse",
        homolog_mapping_h2m_file=str(H2M_MAPPING_FILE),
        #homolog_mapping_m2h_file=str(M2H_MAPPING_FILE),
        source_index_col="gene_name",
        #mouse_gene_name_col="mouse_gene_name",
        human_gene_name_col="human_gene_name",
    )
    # ortholog converted M2H
    convert_gmt_homologs_directory(
        gmt_input_directory=str(MSIGDB_MM_DIR),
        gmt_output_directory=str(MSIGDB_MM_M2H_DIR),
        covert2organism="Human",
        homolog_mapping_m2h_file=str(M2H_MAPPING_FILE),
        source_index_col="gene_name",
        human_gene_name_col="human_gene_name",
    )


if __name__ == "__main__":
    main()
