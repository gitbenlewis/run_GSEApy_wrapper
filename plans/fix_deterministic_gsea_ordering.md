# Deterministic GSEA CSV Ordering

## Summary
Update `make_results_dataframe()` in `run_GSEApy_wrapper` so result-table serialization is fully deterministic for tied `fdr` rows by sorting on `TermCollection` ascending, `fdr` ascending, `nes` descending, then `Term` ascending. Do not change ranking inputs, filtering thresholds, seed, permutations, or any numeric GSEA calculations.

## Implementation Changes
- Change only `_gseapy_pre_rank_wrap.py` in the `make_results_dataframe()` path.
- Replace the current two-key sort `['TermCollection', 'fdr']` with `['TermCollection', 'fdr', 'nes', 'Term']`.
- Use ascending order for `TermCollection`, `fdr`, and `Term`, and descending order for `nes`.
- Keep the implementation aligned to the emitted lowercase column name `nes`.
- Do not modify `process_input_table()`, `run_single_gsea_prerank()`, config defaults, archived code, or any downstream plotting helpers.

## APIs / Output Contract
- No public API or function signature changes.
- The only behavior change is row ordering in emitted CSVs when multiple rows share the same `fdr`.
- Within a single collection, equal-`fdr` rows must now serialize in descending `nes` order, with `Term` as the final deterministic tiebreaker.

## Test Plan
- Use the existing GSE68719 example harness at `examples/GSE68719/scripts/make_gseapy_tables.py` with `examples/GSE68719/config/config.yaml`.
- Regenerate a tie-heavy table such as `GSE68719_mlpd_DESeq2_diffexp_pmi_age_rin_default.GSEA.ARCHS4_Cell-lines.csv` and confirm the change is ordering-only.
- Verify aligned rows retain identical `es`, `nes`, `pval`, `fdr`, and `fwerp`.
- Run the same GSEA generation twice on identical inputs and confirm the emitted CSVs are byte-identical.

## Assumptions
- The canonical fix remains in `/home/ubuntu/projects/gitbenlewis/run_GSEApy_wrapper`.
- The external `diff_test_INPUTnc_MAD20_Post_over_all_controls/...CHEM_ID.csv` example mentioned earlier is not available in this workspace, so repo-local verification should use the shipped GSE68719 example instead.
- `Term` is sufficient as the final tiebreaker because it is already materialized in the result table and makes the per-CSV ordering total.
