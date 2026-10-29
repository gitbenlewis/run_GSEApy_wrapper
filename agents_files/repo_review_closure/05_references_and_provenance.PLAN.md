# Species-aware reference acquisition and provenance

1. **Status:** Complete: offline implementation independently accepted, 2026-10-29; publication authorized. Created 2026-10-28. Findings #9, #10, #18, #30 in [original review](../sep2026_repo_review.md). Repository: run_GSEApy_wrapper. Oversight: current primary agent; implementation delivered by `/root/plan01_implementation`; independently accepted by `/root`; ownership returned.

2. **Objective and scope:** Resolve the linked findings through evidence-backed corrections. Preserve public signatures, valid numerical behavior and unrelated working-tree changes. No new dependencies or scientific assumptions without explicit approval.

3. **Proposed diff:** _gseapy_pre_rank_wrap.py reference functions; config/download_gseapy_gmt_files.py; reference manifest schema; DATA_LICENSES.md; focused reference tests. Add only focused regressions and documentation needed for these findings. Reinspect current files before implementation because review line numbers are historical.

4. **Dependencies and decision gates:** User approved the four species/mapping decisions below. Unknown provenance remains explicitly unresolved; no species evidence is inferred. No live replacement or data publication authorized by planning.

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

# Approved decisions (2026-10-29)

1. **Species routing:** Declare both source-study organism (when known) and the organism of the actual gene identifiers. Route by verified identifier organism, not library name, folder, capitalization or the Enrichr endpoint alone. Native target-species sets stay unchanged; convert only explicitly declared opposite-species sets; retain native and converted collections separately. Unknown or mixed identifier species blocks automatic routing for that collection with an actionable report. Do not silently omit a requested library. Human analyses must not directly use mouse identifiers, and vice versa.

2. **Ambiguous orthologs:** Approved conservative default: retain only symbol pairs with one distinct target per source and one distinct source per target in the supplied mapping's complete relationship columns, after collapsing identical records. Exclude one-to-many and many-to-one relationships and report them. This is a uniqueness rule within these mapping files, not proof of evolutionary one-to-one orthology. Example: AOC1 lists Aoc1l1, Aoc1, Aoc1l3 and Aoc1l2; it would be excluded instead of selecting the existing scalar Aoc1l1. Alternative requiring explicit selection: expand to every listed target and deduplicate each term; this retains more genes but changes term size and composition. Do not infer a best ortholog from row order or scalar first choice.

3. **Conflicting mapping records:** Reject distinct contradictory records for the same source identifier before publishing any converted snapshot; identical duplicates may collapse. Report affected IDs for upstream resolution. The current scalar dict silently selects the last record. Proposal does not automatically merge conflicting records as if they were verified alternatives. Explicit multi-target relationships are handled by decision2, separately from contradictory records.

4. **Unmapped genes and empty sets:** Drop unmapped or policy-excluded genes with per-term and per-library counts; deduplicate retained targets within terms; omit terms that become empty and report them. No arbitrary minimum coverage cutoff. Native sets remain unchanged. Report source, mapped, unmapped, ambiguous, collided/deduplicated and empty-term counts so users can judge losses.

5. **Engineering boundary:** Preserve checked-in snapshots and current releases, including approved MGI2024 and configured MSigDB2025.1. Record hashes and available provenance; mark unavailable mapping release/source metadata unknown. Build and test routing, QC, staged snapshot publication and manifests offline. No live refresh, snapshot replacement, data redistribution or Git publication is authorized by these scientific approvals. Source-specific licenses and citations still require primary-source verification before closing Plan05.

# Implementation scratchpad

15. **2026-10-28 — Created:** Mapped original findings to ordered work, files, policy gates and verification. Existing lollipop work remains uncommitted and must be preserved. No implementation or fresh reproduction performed for this phase. Next: revalidate findings and resolve the stated decisions before dispatch.

16. **2026-10-29 — Decision preparation:** User requested starting Plan05 and presenting approval decisions. Clean Git state at start. Read downloader, conversion functions and mapping headers; current downloader fetches every Enrichr library into a human directory and converts them all. Current loader uses dict(zip(...)), ignores complete target lists and overwrites conflicting source symbols. Measured 14/2 conflicting scalar-source symbols and 1110/1198 multi-target rows in h2m/m2h respectively; these are descriptive counts, not integrity guards. Identifier-species labels are not verified by filename hints. Saved proposed routing, conservative unique-pair mapping, conflict rejection and loss/empty-set handling for user approval. No scientific source edits, downloads or converted snapshots. Mapping provenance remains unresolved. Historical item15's uncommitted lollipop note is superseded by published commit98b6812.

17. **2026-10-29 — Approved and dispatched:** User approve_2_through_5 authorizes the four proposed scientific policies. Reused implementation agent for reference code, downloader, manifest contract, QC, tests and reference documentation. Oversight owns DATA_LICENSES.md, this plan and independent acceptance. Existing snapshots remain unchanged; live downloads and publication are outside scope. Unknown/mixed species declarations must fail explicitly; contradictory checked-in mappings must not be silently repaired. Primary-source license/provenance review in progress.

18. **2026-10-29 — Independent acceptance:** Reviewed reference code, downloader, declaration example, ten reference tests, updated CLI smoke, docs/references.md and primary-source DATA_LICENSES.md. Independently ran full unittest discovery with not_base Python and Agg: 72 tests passed. Public wrapper signature AST comparison unchanged; compile and diff checks passed. Corrected review gaps: contradictory scalar records with same full list reject; positive MSigDB routes and release mismatch covered; endpoint retained in provenance. Native bytes preserved, both conversion directions tested, strict exclusion/QC and empty-term omission verified, identical source acquisitions yield deterministic manifests, interruptions leave no final partial snapshot. Final source-hash positioning correction independently verified for direct in-place conversion. No changes in data/ref.

19. **2026-10-29 — Limits and disposition:** Existing maps fail approved complete-relationship validation (1425 h2m and 41 m2h contradictory source IDs plus blank sources); these descriptive counts differ from prior scalar-only14/2 and are not integrity guards. Do not run live conversion until upstream records are resolved. Unknown/mixed species blocks snapshot requests, and the example declaration deliberately has no fabricated evidence. Mapping provider/release and collection-specific redistribution gaps remain explicitly unresolved in DATA_LICENSES.md; acceptance is of offline engineering and documented disposition, not provenance verification or clearance to publish data. Current primary MSigDB, Enrichr, Reactome and KEGG sources reviewed; production WikiPathways terms inaccessible (403). New snapshot publication is atomic to a previously absent directory; hard kills may leave hidden staging directories, concurrent publishers unsupported, direct low-level GMT/QC writes are not a transaction. Pinned GSEApy1.1.11 environment not recreated; no live downloads, mapping repairs, biological reruns, commits or pushes. Next: Plan06 graph policies.

20. **2026-10-29 — Publication authorized:** User requested Plan05 commit/push, followed by Plan06 decision preparation. Rechecked full suite (72 passing), diff and file sizes; publish this offline implementation only. Existing data snapshots and live-use prerequisites remain unchanged.
