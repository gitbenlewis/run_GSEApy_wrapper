# Existing dotplot correctness

1. **Status:** Complete and independently accepted locally 2026-10-28; publication authorized on 2026-10-29. Findings #11, #12, #19 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Implementation: `/root/plan01_implementation`; independent acceptance: `/root`; ownership returned to oversight.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** _gseapy_dot_plots.py; tests/test_gseapy_dotplots.py; docs/dotplots.md. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Approved policy and precedence:** User approved original probabilities for selection/order and display-only zero replacement by min(1/N, smallest positive probability / 2), retaining 1/N for all-zero data. With N=None, retain the existing all-zero error. Explicit column takes precedence over pvalue_col; preserve the chosen nominal or adjusted measure. Shorten tick labels only, retaining full term identities. New lollipop behavior does not close these findings.

5. **Execution step 1:** Freeze current signatures and ordinary valid calls; add regressions for zero FDR ordering, opposing nominal-P/FDR rankings, raw aliases, numeric strings and colliding truncated labels.

6. **Execution step 2:** Resolve the requested significance measure consistently; preserve an explicit column argument and document precedence when it conflicts with pvalue_col. Do not silently substitute nominal and adjusted measures.

7. **Execution step 3:** Separate selection/sorting values from render-only zero handling. Define and document an explicit floor behavior that cannot reverse significance order, including all-zero data; obtain agreement for changes to existing floor defaults.

8. **Execution step 4:** Keep stable categorical IDs through selection and plotting; forward custom y correctly and shorten rendered labels only.

9. **Execution step 5:** Exercise both single and multiple-table entrypoints, including mixed schemas and missing cells; preserve caller data and keep lollipop regressions green.

10. **Execution step 6:** Render focused before/after examples and review actual selected terms, labels, circle positions and legend meaning.

11. **Acceptance:** Zeros rank ahead of positive values; explicit nominal and adjusted selections differ correctly; aliases work; distinct labels remain distinct rows; visually accepted examples. Oversight reviews actual diffs and reruns focused checks; implementer delivery alone does not close findings.

12. **Verification:** Use the existing unittest runner with Agg and writable MPLCONFIGDIR; run focused new tests first, then the full suite, py_compile and git diff --check. Use synthetic fixtures, temporary outputs and mocked remote calls. Record runtime versions, commands not run, limitations and before/after behavior. Do not rerun private/full enrichment or refresh references as a substitute for targeted tests.

13. **Execution boundary:** This file is a plan, not a claim of completion or approval for policy choices, external publication, installs or destructive work. Dispatch implementation only after resolving its decision gates; retain separate ownership for independent review. No source-count snapshot guards for small files.

# Implementation scratchpad

14. **2026-10-28 — Created:** Mapped original findings to ordered work, files, policy gates and verification. Existing lollipop work remains uncommitted and must be preserved. No implementation or fresh reproduction performed for this phase. Next: revalidate findings and resolve the stated decisions before dispatch.

15. **2026-10-28 — Approved and dispatched:** User approve_3 resolves the zero-display policy gate. Read current source and confirmed early zero replacement, FDR substitution and truncated categorical identifiers. Reused implementation agent with ownership of plotting module, focused tests and dotplot documentation; preserve lollipop and all unrelated edits. Oversight owns plan/roadmap and independent numerical/visual acceptance. No new dependencies or publication authorized.

16. **2026-10-28 — Independent acceptance:** Reviewed full source diff, seven focused tests and docs/dotplots.md. Independently ran `MPLCONFIGDIR=/tmp/gseapy-review-mpl MPLBACKEND=Agg /Users/ben/miniconda3/envs/not_base/bin/python -m unittest discover -s tests`: all 54 tests passed. Compilation passed at implementation handoff; independent git diff --check passed. Compared ordinary positive-FDR scatter coordinates, colors and sizes to HEAD: equal; raw GSEApy default aliases work. Visually reviewed /tmp/plan03-dotplots/before.png, after.png and allzero-multiple.png: corrected nominal selection, separated colliding shortened labels, correct zero display. One obsolete mock expectation in the lollipop test now expects the prepared display column; lollipop implementation unchanged. Selected raw probabilities drive cutoff and top terms; direct -log10 display transformation avoids GSEApy subnormal clipping. Floor scope is the complete selected rendering table. Unrepresentable positive floors raise explicitly. Function signatures, caller inputs and save/return contracts preserved. Runtime remains GSEApy 1.1.13; pinned 1.1.11, private datasets and full biological reruns not tested. No commit/push or dependencies. Next: Plan04 reference-release decision.

17. **2026-10-29 — Publication authorized:** User approved separate commits and push for Plans02–04. Final suite passed 62 tests; diff and blob-size checks passed. This phase is included in its own implementation commit; see Git history for commit identity. No further scientific policy changes or live analyses.
