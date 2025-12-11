"""  Bens wrapper to run GSEApy with reproducible FDR and save tables etc ... """

from . import gseapy_pre_rank_wrap_v1 as rgw

from .gseapy_pre_rank_wrap_v1 import (
    run_gseapy_prerank_multiple_term_collections,
    process_input_table,
    run_single_gsea_prerank,
)