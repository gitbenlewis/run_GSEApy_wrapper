# run_GSEApy_wrapper

Wrapper to run GSEApy with reproducible FDR.

## Set up
### Make the conda environment
```bash
conda env create -f ./config/env_GSEApy.yaml
```

### Download the GMT files and convert orthologs
```bash
# if not previously done
conda activate env_GSEApy
python config/download_gseapy_gmt_files.py
```

## Directory structure
```text
repo_root/
  LICENSE
  README.md
  __init__.py
  gseapy_pre_rank_wrap.py
  config/
    config.yaml
    env_GSEApy.yaml
    download_gseapy_gmt_files.py
  data/
    ref/
      h2m_agg.csv
      m2h_agg.csv
      gmt_files/
        download_gseapy_gmt_files_*.log
        enrichr_Hu/
          *.gmt
        enrichr_Hu_h2m/
          *.h2m.gmt
        msigdb_Hu/
          *.gmt
        msigdb_Hu_h2m/
          *.h2m.gmt
        msigdb_Mm/
          *.gmt
        msigdb_Mm_m2h/
          *.m2h.gmt
  examples/
    GSE68719/
      config/
      data/
      scripts/
        make_gsea_tables.py
        make_gsea_tables_*.log
      results/
        GSEApy/
          DESeq2_diffexp_pmi_age_rin/
            GSEApy/
              *.GSEA.*.csv
              *.rnk_df.csv
          logs/
      notebooks/
```
