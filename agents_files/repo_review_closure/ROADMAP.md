# Repository review closure roadmap

1. **Status:** Closure planning complete on 2026-10-28; remediation remains open. Tracks all 30 findings (#4–#33) from [the original review](../sep2026_repo_review.md). Original observations remain historical, not rewritten as current evidence. Context, positive observations and recommendations (#1–#3, #34–#41) are not separate defects.

2. **Authorization and ownership:** The request is treated as authorization to create the closure roadmap and execution plans. Scientific decisions, dependencies, publication, data deletion and history rewriting remain gated. Oversight owns this roadmap; phase implementation owners are unassigned. Use separate implementation and independent acceptance roles on dispatch, as requested earlier. Do not overlap file ownership between plans sharing the wrapper or plotting module.

3. **Completed finding:** #4 was fixed and pushed in commit 27bcdef, with dictionary/GMT naming, propagated write failures, deterministic collision suffixes and regression coverage. Preserve that fix. Re-run its tests before accepting related changes.

4. **Partial finding:** #32 now has automated tests (25 passed at the preceding lollipop handoff), but CI, run manifests and tied-rank scientific reproducibility remain open. The lollipop feature is delivered/accepted locally in [its plan](../issue_1_lollipop.PLAN.md); it does not close original dotplot findings #11/#12/#19. Existing source/doc/example changes are uncommitted and belong to that feature.

5. **Rank input and result correctness:** [Execution plan](01_rank_and_results.PLAN.md). Findings #5, #7, #20, #22, #23, #24 accepted locally on 2026-10-28; uncommitted. Delivered by `/root/plan01_implementation`, independently reviewed by `/root`. Evidence: 39 tests, compilation/diff checks, baseline rank parity and real two-collection direct-GSEApy numerical parity. User approved strict invalid-input rejection, positive permutation validation and opt-in all-results list. Next: scientific decisions for plan02.

6. **Deterministic scientific ranking policies:** [Execution plan](02_ranking_policies.PLAN.md). Findings #6, #8, #21. Implementation: planned; independent review: pending; owner: unassigned. Next: inspect current source and satisfy the plan's gates.

7. **Existing dotplot correctness:** [Execution plan](03_dotplot_correctness.PLAN.md). Findings #11, #12, #19. Implementation: planned; independent review: pending; owner: unassigned. Next: inspect current source and satisfy the plan's gates.

8. **Portable config-driven workflows:** [Execution plan](04_portable_workflows.PLAN.md). Findings #15, #16, #17, #27. Implementation: planned; independent review: pending; owner: unassigned. Next: inspect current source and satisfy the plan's gates.

9. **Species-aware reference acquisition and provenance:** [Execution plan](05_references_and_provenance.PLAN.md). Findings #9, #10, #18, #30. Implementation: planned; independent review: pending; owner: unassigned. Next: inspect current source and satisfy the plan's gates.

10. **Efficient postprocessing and defensible networks:** [Execution plan](06_postprocessing_and_networks.PLAN.md). Findings #13, #14, #25, #26. Implementation: planned; independent review: pending; owner: unassigned. Next: inspect current source and satisfy the plan's gates.

11. **Packaging, operational failure visibility and reproducibility:** [Execution plan](07_packaging_and_reproducibility.PLAN.md). Findings #28, #29, #32, #33. Implementation: planned; independent review: pending; owner: unassigned. Next: inspect current source and satisfy the plan's gates.

12. **Recoverable data and repository hygiene:** [Execution plan](08_data_hygiene.PLAN.md). Findings #31. Implementation: planned; independent review: pending; owner: unassigned. Next: inspect current source and satisfy the plan's gates.

13. **Dependency order:** Begin characterization for plans 01 and 03; plan 06's one-read improvement can be independent. Resolve ranking policies before deterministic inference acceptance. Coordinate 04 with 07's package/import contract. Complete species/mapping decisions before reference refresh work in 05. Full CI/manifests/docs acceptance in 07 follows the affected correctness plans. Data migration 08 follows provenance and distribution decisions. Work sharing one file is serialized or explicitly handed off.

14. **Scientific decisions to resolve:** Choose invalid-row rejection versus documented dropping; tied-rank ordering; duplicate aggregation; zero-input ranking; species provenance and ortholog mapping policy; union versus mutual kNN; similarity-to-distance transform; and the intended 2021 versus 2024 reference. Preserve per-collection FDR and current scientific defaults unless explicitly changed. A planning recommendation is not approval.

15. **Engineering decisions to resolve:** Choose support status for private network helpers, supported platforms, optional dependencies, all-results return API, existing dotplot floor compatibility, and data distribution destination. No fixed source record-count guards; use schemas, keys, parseability, null/duplicate checks and domain invariants.

16. **Closure protocol:** planned → dispatched → in progress → delivered → independent review → accepted/closed. Findings needing corrections return to changes requested; blocked items record a specific decision and resume condition. Closure requires a reproduction or justified disposition, regression evidence, reviewed diff, documentation and numerical/visual validation appropriate to risk. Record implementation, acceptance, commit and publication separately. Never close a finding merely because a plan exists.

17. **Per-phase handoff:** Record changed paths/revision, tests and commands, actual environments, scientific impact, unresolved limitations and next action in that plan's append-only scratchpad. Update this index with accepted finding IDs and evidence links. Do not close or comment on GitHub issues without explicit instruction; the selected review's numbered findings are not GitHub issue numbers.

# Implementation scratchpad

18. **2026-10-28 — Roadmap created:** Read the entire original review and checked current Git state. Created eight plans and this index covering every finding #4–#33 exactly once (#4 closed separately; #32 partial). Preserved all existing lollipop code, tests, docs and artifacts. No remediation source edits, dependency installation, publication, or destructive cleanup performed in this turn. Next action: select/authorize implementation scope and resolve phase-specific policy choices; start with rank/result characterization and dotplot correctness.

19. **2026-10-28 — Execution authorization and first acceptance:** User authorized roadmap execution starting plan01 and consultation before scientific policy changes. Plan01 decisions subsequently approved and implemented; six findings accepted with independent evidence linked above. This supersedes the planning-only authorization description in item 2 for in-scope execution. Dependencies, external publication and destructive actions remain gated. Next unresolved scientific decisions are tie ordering, duplicate aggregation and zero-input rank handling; no changes to those are authorized yet.

