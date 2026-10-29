# Recoverable data and repository hygiene

1. **Status:** Planned, not implemented. Created 2026-10-28. Findings #31 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Oversight: current primary agent; implementation owner unassigned until dispatch.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** inventory/proposal documents; precise .gitignore entries if approved; data manifest; README acquisition instructions. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Dependencies and decision gates:** Depends on provenance/terms plan 05 and reproducibility plan 07. No deletion, LFS migration, release upload, force-push or history rewrite authorized.

5. **Execution step 1:** Recompute current tracked sizes and duplicate hashes; do not reuse historical file counts as integrity checks.

6. **Execution step 2:** Identify compact deterministic fixtures required by tests and immutable snapshots required for reproducibility.

7. **Execution step 3:** Propose versioned release/data-registry destinations with checksums, access requirements, retention and restore commands. Keep local originals intact.

8. **Execution step 4:** Obtain explicit approval for file removal, publication and any history rewrite separately; create and verify recoverable backups before destructive migration.

9. **Execution step 5:** After approved migration, test acquisition and full restore, update exact ignore rules, and validate outgoing Git blob sizes.

10. **Execution step 6:** Close only after retrieval/restore evidence and storage disposition are accepted; a written migration proposal alone is not closure.

11. **Acceptance:** Verified inventory, approved destination, checksum-preserving restore, small test fixtures, and explicit recovery/publication evidence. Oversight reviews actual diffs and reruns focused checks; implementer delivery alone does not close findings.

12. **Verification:** Use the existing unittest runner with Agg and writable MPLCONFIGDIR; run focused new tests first, then the full suite, py_compile and git diff --check. Use synthetic fixtures, temporary outputs and mocked remote calls. Record runtime versions, commands not run, limitations and before/after behavior. Do not rerun private/full enrichment or refresh references as a substitute for targeted tests.

13. **Execution boundary:** This file is a plan, not a claim of completion or approval for policy choices, external publication, installs or destructive work. Dispatch implementation only after resolving its decision gates; retain separate ownership for independent review. No source-count snapshot guards for small files.

# Implementation scratchpad

14. **2026-10-28 — Created:** Mapped original findings to ordered work, files, policy gates and verification. Existing lollipop work remains uncommitted and must be preserved. No implementation or fresh reproduction performed for this phase. Next: revalidate findings and resolve the stated decisions before dispatch.

