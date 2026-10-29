# Rank input and result correctness

1. **Status:** Complete and independently accepted locally, 2026-10-28; uncommitted. Findings #5, #7, #20, #22, #23, #24 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Oversight: `/root`; implementation delivered by `/root/plan01_implementation`. User-approved decisions 2–4 governed invalid inputs, permutation validation and opt-in all-results returns.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** _gseapy_pre_rank_wrap.py; tests/test_gseapy_pre_rank_wrap.py; docs/rank_inputs.md. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Dependencies and decision gates:** Invalid-row reject/drop policy, ES-only callers, and all-results return API require decisions. No new statistical defaults.

5. **Execution step 1:** Characterize current missing-value behavior in both rank modes before choosing reject-versus-drop; reproduce missing IDs, invalid probabilities, infinities, empty split tokens and conflicting multi-name flags.

6. **Execution step 2:** Validate required columns and numerical domains before converting identifiers. Preserve valid ranks and zero handling. Obtain the invalid-row policy decision before changing existing drop behavior.

7. **Execution step 3:** Trace kwargs into GSEApy: preserve the public signature and default call; forward supported caller options such as weight=0 and reject unknown options explicitly. Test exact forwarded arguments.

8. **Execution step 4:** Correct the reported/tested-terms log and document per-collection FDR. Record family and effective parameters in a sidecar manifest rather than silently changing the CSV schema.

9. **Execution step 5:** Reject nonpositive permutation counts for inferential runs after confirming callers do not depend on ES-only behavior; retain default 100.

10. **Execution step 6:** Characterize empty/single/multiple collection returns. Decide an opt-in all-results return contract before adding API surface; preserve the existing default.

11. **Execution step 7:** Remove genuinely missing metabolite IDs before conversion; preserve legitimate identifiers and deterministic order. Test missing values without treating every literal string nan as missing.

12. **Acceptance:** Golden valid-rank parity; invalid-input diagnostics; exact gp.prerank calls; zero/one/multi-return tests; missing-metabolite regression; existing naming/collision tests green. Oversight reviews actual diffs and reruns focused checks; implementer delivery alone does not close findings.

13. **Verification:** Use the existing unittest runner with Agg and writable MPLCONFIGDIR; run focused new tests first, then the full suite, py_compile and git diff --check. Use synthetic fixtures, temporary outputs and mocked remote calls. Record runtime versions, commands not run, limitations and before/after behavior. Do not rerun private/full enrichment or refresh references as a substitute for targeted tests.

14. **Execution boundary:** This file is a plan, not a claim of completion or approval for policy choices, external publication, installs or destructive work. Dispatch implementation only after resolving its decision gates; retain separate ownership for independent review. No source-count snapshot guards for small files.

# Implementation scratchpad

15. **2026-10-28 — Created:** Mapped original findings to ordered work, files, policy gates and verification. Existing lollipop work remains uncommitted and must be preserved. No implementation or fresh reproduction performed for this phase. Next: revalidate findings and resolve the stated decisions before dispatch.

16. **2026-10-28 — Execution started:** User authorized roadmap execution starting with this plan and required consultation before scientific-policy changes. Inspected current code and worktree. Synthetic read-only reproduction confirmed missing feature None becomes the string 'None'; missing signed-mode P is silently dropped; P=2 reverses a positive effect's rank to -0.3010299957; negative P produces NaN; precomputed NaN/infinite ranks pass through. Searched active examples/config/tests: configured permutations are 1000 and no explicit ES-only usage was found; this does not establish external caller behavior. No source edits yet. Decisions requested: missing-row policy, nonpositive permutation policy, and multi-result return contract. Preserve zero replacement, duplicate aggregation and tie ordering for plan 02.

17. **2026-10-28 — Decisions approved and implementation dispatched:** User replied `approve_2_through_4`: reject missing/blank IDs and invalid/missing/nonfinite metrics, reject nonpositive permutation counts without changing default 100, and add optional return_all_pre_res=False while preserving existing default/single return behavior. Use an ordered result list to preserve duplicate collection names; reject simultaneous return modes and empty batches before outputs. Existing multi-name option precedence remains unchanged pending any separate policy decision. Assigned wrapper, focused rank/result tests and docs/rank_inputs.md to `/root/plan01_implementation`; oversight owns this plan. Lollipop work and all unrelated files must remain untouched. Supported-kwargs forwarding, truthful logging, sidecar metadata and genuine missing-metabolite cleanup are included in authorized plan scope. No scientific changes beyond approved decisions, dependency installs, commits or publication.

18. **2026-10-28 — Independent review correction:** New sidecar initially failed on valid NumPy scalar settings after writing CSVs. Oversight reproduced results_table_min_genes_matched=np.int64(1) causing JSON TypeError and requested correction. Implementer added NumPy scalar serialization, removed inappropriate truncating metadata casts and added a regression. This issue is resolved; actual scientific parameter values are preserved.

19. **2026-10-28 — Acceptance evidence:** Oversight read actual diff, new tests and docs/rank_inputs.md. All 39 tests passed, changed Python files compiled, and git diff --check passed. Compared valid ranks with committed HEAD in both modes including exact zero, p=1 and duplicate IDs: exact DataFrame parity. A real two-collection GSEApy run returned ordered results, produced the mapped CSVs/metadata and matched direct GSEApy ES, NES, nominal P, FDR and FWER exactly. Installed environment is Python 3.10/GSEApy 1.1.13; pinned 1.1.11 and full dataset runs were not executed. No dependencies or publication changes. New tests live in tests/test_rank_inputs.py; existing output-file assertions now include the intentional JSON sidecar. No-save outputs remain absent.

20. **Closure and limits:** Findings #5/#7/#20/#22/#23/#24 accepted locally. Preserved legacy last-result return unless the new opt-in list mode is selected. CSV schema unchanged; sidecar is an intentional additional output and is not a complete provenance manifest. Prior batch files can remain after a failure; transactional output is not claimed. Plan02's tie/duplicate/zero scientific policies remain unresolved and unchanged. No commit or push performed.

