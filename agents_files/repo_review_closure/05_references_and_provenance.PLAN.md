# Species-aware reference acquisition and provenance

1. **Status:** Planned, not implemented. Created 2026-10-28. Findings #9, #10, #18, #30 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Oversight: current primary agent; implementation owner unassigned until dispatch.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** _gseapy_pre_rank_wrap.py reference functions; config/download_gseapy_gmt_files.py; reference manifest schema; DATA_LICENSES.md; focused reference tests. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Dependencies and decision gates:** Requires ortholog policy, species provenance and reference release decisions. No live replacement or data publication authorized by planning.

5. **Execution step 1:** Inventory source species, releases, mapping sources and available license/citation evidence. Mark unknown provenance explicitly; do not fabricate retrieval dates or original checksums.

6. **Execution step 2:** Approve species routing and ortholog mapping policy. Native mouse resources must remain native; human-to-mouse conversion applies only to declared human inputs.

7. **Execution step 3:** Detect conflicting source mappings and emit per-library mapped/unmapped/collision/empty-set QC without imposing snapshot-specific count guards.

8. **Execution step 4:** Implement a versioned manifest including source URL, species, release, file hashes, mapping hashes and tool versions; distinguish acquisition metadata from scientific contents.

9. **Execution step 5:** Stage downloads in a separate directory; validate schema and expected requested keys; atomically publish complete snapshots with interruption tests. Preserve existing checked-in snapshots during development.

10. **Execution step 6:** Review official collection terms and record source-specific citations/notices. Escalate unresolved redistribution questions before publication; software MIT is not a data license.

11. **Execution step 7:** Test mocked identical snapshots, failed downloads, conversion sentinels and many-to-one/one-to-many fixtures.

12. **Acceptance:** Species sentinel correctness; explicit conflict/QC accounting; reproducible manifests; interruption preserves prior snapshot; traceable license evidence or recorded unresolved disposition. Oversight reviews actual diffs and reruns focused checks; implementer delivery alone does not close findings.

13. **Verification:** Use the existing unittest runner with Agg and writable MPLCONFIGDIR; run focused new tests first, then the full suite, py_compile and git diff --check. Use synthetic fixtures, temporary outputs and mocked remote calls. Record runtime versions, commands not run, limitations and before/after behavior. Do not rerun private/full enrichment or refresh references as a substitute for targeted tests.

14. **Execution boundary:** This file is a plan, not a claim of completion or approval for policy choices, external publication, installs or destructive work. Dispatch implementation only after resolving its decision gates; retain separate ownership for independent review. No source-count snapshot guards for small files.

# Implementation scratchpad

15. **2026-10-28 — Created:** Mapped original findings to ordered work, files, policy gates and verification. Existing lollipop work remains uncommitted and must be preserved. No implementation or fresh reproduction performed for this phase. Next: revalidate findings and resolve the stated decisions before dispatch.

