# Paper edits, 18 Nov 2026 (Jeff's 17–18 Nov comments)

Edited copies, with tracked changes by "Claude":

- `paper/incidence-peer_review 4.docx` (from `~/claude_visible/incidence-peer_review 3.docx`)
- `paper/incidence-supplementary-peer_review 5.docx` (from `~/claude_visible/incidence-supplementary-peer_review 4.docx`)

The originals in `~/claude_visible/` are untouched. No models were refitted. Every number comes from the current
`outputs/paper_numbers.csv` and `outputs/tables/`.

## Done: main text

- Methods: "OVID" changed to "OVID Global Health" (comment 391, with a reply).
- "25+" changed to "25-49" throughout:
  - Abstract, Methods (cohort and PHIA age models, share of new infections), Results ¶3/¶4, Discussion, Figure 5 caption.
  - Siaya footnote: a real Word footnote on "older adults (25-49 years)" in the Methods, plus the same sentence at
    the end of the Figure 5 caption.
- Figure 5 caption, panel B: "PHIA surveys (national incidence rate)".
- Results numbers:
  - Rakai 73% (64, 79%)
  - Kyamulibwa 69% (51, 80%)
  - Siaya 71% (62, 79%)
  - 2022 share of new infections under 25: IQR 26-34%
  - 2022 share of new infections among women: IQR 67-70%
  - The PHIA 10% (−4, 24%) was already correct.
- Figures 2, 3, 4 and 5 replaced with the regenerated PNGs. Image swaps are not tracked by Word.
- Figure 4: Jeff had already deleted the "Namibia not shown" sentence, and Namibia (NAM) is already in the caption's
  country list, so no further change.
- Replies added to comments 391, 2475 ("Bayesian": no change, INLA), 2456, 2468, 2490, 2489 and 2531.

## Done: supplement

- Text S1: "OVID global health" changed to "OVID Global Health", in the prose and in the search-term heading.
- Text S5:
  - "25+" changed to "25-49" in the prose and in the equation (a ∈ {15-24, 25-49}).
  - New sentence: the age-specific cohort intercepts are constrained to sum to zero across cohorts within each age
    group.
- Table S9: the 25-49 row is now 1.06 (0.88, 1.30); row labels are 25-49.
- Tables S10 and S11: labels are 25-49. That is 248 cells in S11.
- Figure S4 caption: 25-49.
- Figures replaced: Text S4 KZN figure (larger text), Text S7 (with Namibia), Text S9 risk-of-bias figure (colour
  legend, no caption), and Figures S2, S3 and S4. Figure S1 is unchanged; the PNG is identical.
- Replies added to comments 25 (Table S9 constraint), 14, 4 and 11.
- No supplement text says Namibia was left out, so nothing needed removing.

## For Oliver to decide

1. **Abstract CrI.** Resolved: changed back to "(95%CrI 1, 6%)" to match `paper_numbers.csv`, as a tracked change
   with a comment to Jeff. "15-24 declined faster" is the inverse of the 25-49 interaction (1 − 1/1.067 = 6.3%),
   not aIRR − 1 (6.7%).
2. **Hayes 1807 ("randomised HIV prevention trials").** Draft options:
   - (a) "four population-representative community-randomised trials of universal HIV testing and treatment (UTT): …"
   - (b) Keep "prevention" but make the UTT focus explicit: "four community-randomised HIV prevention trials of
     universal testing and treatment (UTT): …"
3. **Grabowski 1821 (baseline).** An optional Discussion sentence after "…among young adults aged 15-24 years.":
   "These estimates provide an important baseline against which to assess the impact of new prevention
   technologies, including long-acting PrEP, as they are scaled up."
4. **Supplement comment 10 (ChatGPT risk-of-bias extraction).** Not answered. Suggested approach: hand-check a random
   ~10% of studies against the PDFs, report the item-level agreement in Text S9, and reply to Jeff with the result.
5. **Flowchart comments (2060, 2061, 2062, 2065, 2066).** Not touched; they belong to the separate flowchart handoff.

## Not done (needs refits or publishing)

- **Reproducible INLA sampling** (`num.threads = "1:1"`). This needs `R/04_fit_cohort_phia.R` re-run, which was out
  of scope ("don't refit"). If it is done later, re-check every number listed above.
- **Zenodo.** Needs a new version with Namibia plus the code archive, and the new DOI in the data-availability text.
  Nothing was published.
- Accepting tracked changes and resolving comments is left for the final pre-submission pass.

## Refit edits (after the Imbokodo refit, 18 Nov)

New copies, tracked changes by "Claude" (2026-11-18T12:00), built from the files above:

- `paper/incidence-peer_review 5.docx` (from `incidence-peer_review 4.docx`)
- `paper/incidence-supplementary-peer_review 6.docx` (from `incidence-supplementary-peer_review 5.docx`)

- Gray et al. 2024 (Lancet Infect Dis, doi 10.1016/S1473-3099(24)00358-X) was already reference 54; no CROI 2022 or
  natap citation remains, so no citation edit was needed.
- Main text: every number moved by the refit (see `paper/refit_numbers_2026-11-18.md`): Abstract and Results aIRRs
  (calendar time 1.07; trials 2.71 (1.70, 4.32); pregnant 3.68; other 3.37; FSW 6.30 (4.21, 9.41); MSM 28.86;
  retrospective 1.62), cohort declines (uMkhanyakude 47-61, Rakai 64-80, Kyamulibwa 68%, Siaya 61-79), IQR 27-34%
  and 66-69%, PHIA women 25-49 (5, 28%). Table 1: 8 cells, including RW2 0.03 (0.00, 0.10).
- Figure 1 replaced with the revised flowchart (`outputs/figures/figure_1.png`), with a comment on the caption
  (image swaps are not tracked).
- Supplement: Tables S6 (46 cells), S7 (7), S8 (68; some were already stale before the refit), S12 (5) now match
  `outputs/tables/`; the Table S6 footnote's negative binomial size parameter is now 11.15.
- Checks: accepted text differs from the input only in these numbers; reject-all text is identical to the input.
