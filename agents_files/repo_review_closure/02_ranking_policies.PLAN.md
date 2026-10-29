# Deterministic scientific ranking policies

1. **Status:** Complete and independently accepted locally on 2026-10-28; publication authorized on 2026-10-29. Findings #6, #8, #21 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Implementation: `/root/plan01_implementation`; independent acceptance: `/root`; file ownership returned to oversight.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** _gseapy_pre_rank_wrap.py; tests/test_ranking_policies.py; docs/ranking_policies.md. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Approved policies:** User approved descending ranks with alphabetical feature-ID ties and preserving prepared order in GSEApy; retain largest-absolute duplicate selection while reporting opposing signs and rejecting equal-magnitude opposing-sign ambiguity; retain the existing subnormal replacement for zero input probabilities and document/report its use. Plan01 is pushed as 52d1cf5. No other scientific changes authorized.

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

14. **2026-10-28 — Approved and dispatched:** User annotations approved tied-rank and duplicate policies and approve_4 retained zero handling. Assigned implementation to the existing plan01 agent with wrapper/tests/ranking docs ownership; primary retains independent review and roadmap. Requested real shuffle-invariance experiments and explicit treatment of caller ascending overrides. No commit/push authorization inferred for this new phase from the annotation targeting plan01.

15. **2026-10-28 — Delivered and independently accepted:** Reviewed source, eight new tests and final documentation. Full unittest suite: 47 passed; compilation and tracked diff whitespace checks passed. Real three-shuffle signed-rank tests match all five enrichment metrics and every saved CSV/JSON byte. Independent all-equal-rank experiment confirmed old shuffle sensitivity and new exact invariance. Canonical lexical ordering can materially change enrichment direction when membership correlates with IDs; documented as a reproducibility convention, not biological evidence. Opposing-sign maximum ties fail only when collapsing duplicates; lower opposite-sign ties with a unique larger maximum remain valid. Zero probability replacement is unchanged and reported before name expansion. Direct calls preserve supplied order by default; explicit ascending overrides remain supported. Runtime: Python 3.10.20, GSEApy 1.1.13, pandas 2.3.3, NumPy 2.2.6, Matplotlib 3.10.8. Pinned GSEApy 1.1.11 and private/full analyses were not run. No new dependencies, commits or pushes; unrelated lollipop changes preserved. Next: Plan03 display-policy approval.

16. **2026-10-29 — Publication authorized:** User approved separate commits and push for Plans02–04. Final suite passed 62 tests; diff and blob-size checks passed. This phase is included in its own implementation commit; see Git history for commit identity. No further scientific policy changes or live analyses.
