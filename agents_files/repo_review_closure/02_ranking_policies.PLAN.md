# Deterministic scientific ranking policies

1. **Status:** Planned, not implemented. Created 2026-10-28. Findings #6, #8, #21 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Oversight: current primary agent; implementation owner unassigned until dispatch.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** _gseapy_pre_rank_wrap.py; tests/test_ranking_policies.py; docs/ranking_policies.md. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Dependencies and decision gates:** Requires explicit tie, duplicate and zero-input policy approval. Defaults remain unchanged until approved.

5. **Execution step 1:** Reproduce shuffle-dependent tied-rank results and characterize duplicate symbols including opposite signs and equal absolute ranks.

6. **Execution step 2:** Prepare policy examples showing the consequences of lexical stable tie ordering, duplicate max-absolute selection, and the current subnormal zero-p replacement. Report effects without selecting a new scientific policy.

7. **Execution step 3:** After approval, implement the selected tie order and configure GSEApy to preserve it. Do not perturb rank values merely to break ties.

8. **Execution step 4:** Implement or retain the approved duplicate aggregation and zero-input policy; report conflicts and serialize policy metadata.

9. **Execution step 5:** Run deterministic small enrichments across shuffled rows, comparing ES, NES, nominal P, FDR and FWER plus rank/output ordering.

10. **Acceptance:** Repeated/shuffled fixture parity under approved policy; documented expected numeric differences from baseline; duplicate-direction and zero-neighbor tests. Oversight reviews actual diffs and reruns focused checks; implementer delivery alone does not close findings.

11. **Verification:** Use the existing unittest runner with Agg and writable MPLCONFIGDIR; run focused new tests first, then the full suite, py_compile and git diff --check. Use synthetic fixtures, temporary outputs and mocked remote calls. Record runtime versions, commands not run, limitations and before/after behavior. Do not rerun private/full enrichment or refresh references as a substitute for targeted tests.

12. **Execution boundary:** This file is a plan, not a claim of completion or approval for policy choices, external publication, installs or destructive work. Dispatch implementation only after resolving its decision gates; retain separate ownership for independent review. No source-count snapshot guards for small files.

# Implementation scratchpad

13. **2026-10-28 — Created:** Mapped original findings to ordered work, files, policy gates and verification. Existing lollipop work remains uncommitted and must be preserved. No implementation or fresh reproduction performed for this phase. Next: revalidate findings and resolve the stated decisions before dispatch.

