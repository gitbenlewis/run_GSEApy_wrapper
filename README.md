# run_GSEApy_wrapper
wrapper to run GSEApy with reproducible FDR
# set up
## make the conda enviroment
    conda env create -f ./config/env_GSEApy.yaml
## download the gmt files and convert orthologs
    # if not previously done
    conda activate env_GSEApy
    python config/download_gseapy_gmt_files.py

# directory structure
    repo_root/
    config/
        config.yaml
        env_GSEApy.yaml
        download_gseapy_gmt_files.py
    __init__.py
    gseapy_pre_rank_wrap.py
    data/
        ref/
            h2m_agg.csv
            m2h_agg.csv
            gmt_files/
                enrichr_Hu
                msigdb_Hu
                msigdb_Mm
                enrichr_Hu_h2m
                msigdb_Hu_h2m
                msigdb_Mm_m2h
    examples
        GSE68719
            config
            data
            scripts
            results
            notebooks
