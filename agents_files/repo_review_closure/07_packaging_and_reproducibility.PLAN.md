# Packaging, operational failure visibility and reproducibility

1. **Status:** Planned, not implemented. Created 2026-10-28. Findings #28, #29, #32, #33 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Oversight: current primary agent; implementation owner unassigned until dispatch.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** pyproject.toml if approved; config/env_GSEApy.yaml; download/batch scripts; README.md; CI configuration; run-manifest code; tests. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Dependencies and decision gates:** Depends on plans 01/02/04/05 for full reproducibility claims. Supported platforms, dependencies and CI installation changes need approval; target environment validation cannot be inferred from local tests.

5. **Execution step 1:** Inventory direct imports and select supported platforms with the user. Preserve public package imports; design minimal packaging without unrelated layout changes.

6. **Execution step 2:** Declare only approved dependencies/extras. Keep existing environments until a canonical replacement and target-platform checks are available.

7. **Execution step 3:** Add bounded download timeouts, content/gzip/hash validation where trusted expected hashes exist, and preserve temporary-file replacement behavior.

8. **Execution step 4:** Make batch plotting failures visible with explicit accounting and nonzero failure exit; define expected-empty behavior without inventing success.

9. **Execution step 5:** Define run manifest sidecars with input/reference/mapping hashes, effective parameters, package versions, commit and scientific policies. Do not turn manifest timestamps into numerical determinism claims.

10. **Execution step 6:** Add offline CI/import/CLI smoke checks and tiny repeated-run enrichment checks after ranking policies are settled; keep numerical and serialized-order assertions distinct.

11. **Execution step 7:** Replace stale documentation with a tested end-to-end workflow and supported/private API statement; carry current lollipop documentation forward.

12. **Acceptance:** Installed import/CLI success outside checkout; corruption/failure exits tested; complete sidecars; repeated-run evidence; executable docs; platform-specific results recorded. Oversight reviews actual diffs and reruns focused checks; implementer delivery alone does not close findings.

13. **Verification:** Use the existing unittest runner with Agg and writable MPLCONFIGDIR; run focused new tests first, then the full suite, py_compile and git diff --check. Use synthetic fixtures, temporary outputs and mocked remote calls. Record runtime versions, commands not run, limitations and before/after behavior. Do not rerun private/full enrichment or refresh references as a substitute for targeted tests.

14. **Execution boundary:** This file is a plan, not a claim of completion or approval for policy choices, external publication, installs or destructive work. Dispatch implementation only after resolving its decision gates; retain separate ownership for independent review. No source-count snapshot guards for small files.

# Implementation scratchpad

15. **2026-10-28 — Created:** Mapped original findings to ordered work, files, policy gates and verification. Existing lollipop work remains uncommitted and must be preserved. No implementation or fresh reproduction performed for this phase. Next: revalidate findings and resolve the stated decisions before dispatch.

