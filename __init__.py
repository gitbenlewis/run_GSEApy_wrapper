"""  Bens wrapper to run GSEApy with reproducible FDR and save tables etc ... """

from . import gseapy_pre_rank_wrap as gpw

from .gseapy_pre_rank_wrap import (
    run_gseapy_prerank_multiple_term_collections,
    process_input_table,
    run_single_gsea_prerank,
    geneset_table_with_multiple_feature_label_split_to_gseapy_gmt,
    metabolon_chemical_annotation_to_pathway_table,
    make_gmts_from_enrichr_libraries,
    make_gmts_from_Msigdb,
    convert_gmt_homologs_directory,
    basic_set_of_enrichr_term_collections,
    expansion_set_of_enrichr_term_collections,
    all_enrichr_libraries_2025_sin_archive,
)