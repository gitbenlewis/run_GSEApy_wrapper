# Portable config-driven workflows

1. **Status:** Planned, not implemented. Created 2026-10-28. Findings #15, #16, #17, #27 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Oversight: current primary agent; implementation owner unassigned until dispatch.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** config/download_gseapy_gmt_files.py; examples/GSE68719/config/config.yaml; examples/GSE68719/scripts/*.py; focused subprocess/config tests. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Dependencies and decision gates:** Reference release selection requires approval; do not silently replace 2021 with 2024. Path changes must not rewrite existing outputs.

5. **Execution step 1:** Read config-driven-script skill and inspect all config consumers before editing. Map effective defaults and run modes.

6. **Execution step 2:** Repair downloader import topology with mocked remote calls; coordinate the final installed entrypoint with packaging plan 07.

7. **Execution step 3:** Resolve bundled assets relative to config/repository and accept explicit config location; remove the unrelated external output path from active example configuration without deleting its files.

8. **Execution step 4:** Align permutation configuration keys while retaining effective 1000; implement documented selected/all semantics or fail clearly on unsupported modes.

9. **Execution step 5:** Resolve every requested collection to exactly one file before any run. Report missing or ambiguous names rather than silently filtering.

10. **Execution step 6:** Keep 2021 versus 2024 reference choice pending explicit approval; run subprocess tests from unrelated directories and copied checkouts.

11. **Acceptance:** Offline mocked CLI smoke from arbitrary cwd; config overrides reach wrapper; requested collection inventory resolves exactly; no unrelated directories created. Oversight reviews actual diffs and reruns focused checks; implementer delivery alone does not close findings.

12. **Verification:** Use the existing unittest runner with Agg and writable MPLCONFIGDIR; run focused new tests first, then the full suite, py_compile and git diff --check. Use synthetic fixtures, temporary outputs and mocked remote calls. Record runtime versions, commands not run, limitations and before/after behavior. Do not rerun private/full enrichment or refresh references as a substitute for targeted tests.

13. **Execution boundary:** This file is a plan, not a claim of completion or approval for policy choices, external publication, installs or destructive work. Dispatch implementation only after resolving its decision gates; retain separate ownership for independent review. No source-count snapshot guards for small files.

# Implementation scratchpad

14. **2026-10-28 — Created:** Mapped original findings to ordered work, files, policy gates and verification. Existing lollipop work remains uncommitted and must be preserved. No implementation or fresh reproduction performed for this phase. Next: revalidate findings and resolve the stated decisions before dispatch.

