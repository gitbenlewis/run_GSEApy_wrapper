# Efficient postprocessing and defensible networks

1. **Status:** Planned, not implemented. Created 2026-10-28. Findings #13, #14, #25, #26 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Oversight: current primary agent; implementation owner unassigned until dispatch.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** _gseapy_post_process_helpers.py; tests/test_postprocessing.py; docs/networks.md. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Dependencies and decision gates:** Requires network support decision, kNN and distance semantics; new NetworkX installation/declaration needs approval. Do not change scientific network metrics implicitly.

5. **Execution step 1:** First resolve whether network helpers are supported API or experiments; preserve exports until that decision.

6. **Execution step 2:** Fix leading-edge mapping to parse the source once and split the already-loaded data. Verify exact mapping parity and one read_csv call.

7. **Execution step 3:** Obtain union-kNN versus mutual-kNN policy and stable neighbor tie rule. Implement permutation-invariant named edge extraction with explicit symmetry.

8. **Execution step 4:** Obtain similarity-to-distance transform policy; maintain similarity weight for strength and a separate positive distance for shortest paths. Test the triangle counterexample.

9. **Execution step 5:** Keep incidence sparse, compute bounded candidate neighbors before full pairwise materialization, and characterize dense Pearson and exact centrality memory/runtime.

10. **Execution step 6:** Define resource limits from working-set estimates and available memory, not fixed source file record counts; document approximate versus exact methods separately.

11. **Execution step 7:** Run small dense/sparse parity, zero/one-node, row-shuffle and realistic synthetic scale tests with measured memory/runtime.

12. **Acceptance:** One-pass I/O parity; row-order invariant edge sets; approved distance metric fixture; sparse scale evidence and bounded unsupported operations. Oversight reviews actual diffs and reruns focused checks; implementer delivery alone does not close findings.

13. **Verification:** Use the existing unittest runner with Agg and writable MPLCONFIGDIR; run focused new tests first, then the full suite, py_compile and git diff --check. Use synthetic fixtures, temporary outputs and mocked remote calls. Record runtime versions, commands not run, limitations and before/after behavior. Do not rerun private/full enrichment or refresh references as a substitute for targeted tests.

14. **Execution boundary:** This file is a plan, not a claim of completion or approval for policy choices, external publication, installs or destructive work. Dispatch implementation only after resolving its decision gates; retain separate ownership for independent review. No source-count snapshot guards for small files.

# Implementation scratchpad

15. **2026-10-28 — Created:** Mapped original findings to ordered work, files, policy gates and verification. Existing lollipop work remains uncommitted and must be preserved. No implementation or fresh reproduction performed for this phase. Next: revalidate findings and resolve the stated decisions before dispatch.

