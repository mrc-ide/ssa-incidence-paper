# Pre-submission review: main text + supplement (17 Nov 2026)

Files reviewed: `paper/HIV incidence declines ... sub_refs-fixed.docx` (latest main text),
`paper/incidence-supplementary-peer_review.docx`, against `outputs/paper_numbers.csv`,
`outputs/fits/*`, `outputs/tables/*` and `outputs/figures/*`.

Severity: **[A]** wrong number or claim, fix before submission; **[B]** inconsistency a reviewer would
notice; **[C]** copy edit.

Status (fourth pass, 17 Nov 2026, docx saved 16:38): ✅ fixed · ◐ partly fixed · ❌ still open.
The main text was read with all tracked changes accepted.

---

## 0. Status summary (fourth pass; fifth-pass fixes marked in §1–§3)

- **Figures:** every embedded figure is byte-identical to the current `outputs/figures/*.png`, including the
  regenerated Fig 5, S2, S3, S4, Text S7 and Text S9 figures. **Figure 1 (flowchart)** is the exception: it is
  still the old version and now reads "N=179, n observations = 3546" (item 1).
- **Tracked changes [B]:** the main docx still has about 1,080 insertions and 1,090 deletions, plus 26 open
  comments (Hayes, Grabowski, Imai-Eaton, Stevens). Accept the changes and resolve the comments before submission.
  The affiliation fix depends on accepting them (see §2).
- **Numbers changed since the third pass:** Wahome 2020 participation is now set to NA, so S9 reads "78 of 175"
  and the key-population correlations become n=10 / n=11. The Kyamulibwa sex ratio is now 3.0:1, an incidence rate
  ratio. The supplement text matches both.
- **Still open, [A]:** 1 (flowchart), 3 (2022 shares extrapolated), 4 ("South African" one-in-three),
  16 / refs 1–2 (2024 decline cited to 2021/2022 reports).
- **New issues** are listed in §1 (items 33–38).

---

## 1. Numeric consistency

### Main text

| # | Status | Where | Issue | Fix |
|---|---|---|---|---|
| 1 | ✅ | Figure 1 | **[A]** Flowchart is from the old analysis. It now reads "Analysis dataset (N=179, n observations = 3546)", but the text, Table S5 and outputs give 175 studies / 3,509 obs. The PHIA box says 2016–2023 (data 2016–2022). The AHRI box says 2005–2021 (data midpoints 2004–2022) and uses "AHRI" and "KEMRI, Kenya" rather than uMkhanyakude / Siaya. Pooled box N=368, but 286+78 = 364 and `flowchart_counts.csv` gives 353. FSW+MSM box N=52, but 32+24 = 56 (Table 1: 28+22 = 50; `flowchart_counts.csv`: 27+22 = 49). Exclusion counts (99/16/16/39/31/22) don't match `flowchart_counts.csv`. | Regenerate from `outputs/data/flowchart_counts.csv` **Fifth pass:** flowchart.pptx regenerated from `flowchart_counts.csv` (pooled 353; 8 exclusion steps 137/15/19/12/17/24/33/4; FSW+MSM 49; uMkhanyakude 2004–2022, Siaya; PHIA 2016–2022; 175 / 3,509). Exported to `outputs/figures/figure_1.png` and replaced in the paper. |
| 2 | ✅ | Results ¶3 | "(Figure 5B)" for shares → now "Figure 5C". | |
| 3 | ✅ | Results ¶3/¶4 | **[A]** The 2022 medians (32% under-25; 68% female) are taken across all 7 cohorts in 2022. Only uMkhanyakude and Kyamulibwa have data to 2022. The other five are extrapolated (Karonga ends 2010, Manicaland 2013, Kisesa 2015, Rakai 2019, Siaya 2021), and Fig 5C/S4 only show the observed period. | Restrict to each cohort's latest observed year or say "projected to 2022" **Fifth pass:** "Projected to 2022" added to both sentences. |
| 4 | ✅ | Discussion ¶3 | **[A]** "as few as one in three infections occurred among young adults within South African, Ugandan, and Kenyan cohorts". uMkhanyakude is at 45% under-25 in 2022. | Drop "South African" or rephrase **Fifth pass:** "South African" dropped. |
| 5 | ➖ | Discussion ¶3 | **[B]** "up to half of all new HIV infections within cohort studies now occur in women over age 25". The maximum in any cohort's latest year is 43%. | "over 40%" **Fifth pass:** kept as "up to half" (Oliver). Cohort maxima: Kyamulibwa 43%, Siaya 42%, Rakai 40%. |
| 6 | ◐ | Results ¶4 | Kyamulibwa now "about 3:1 (2020-2022)", which matches the regenerated `paper_numbers.csv` (3.0:1 incidence rate ratio) ✅. **Still open:** "latest data from uMkhanyakude, Kyamulibwa, Siaya, and Rakai ... indicated a rapid increase in the sex ratio", but uMkhanyakude's last point in Fig S3B falls. | Drop uMkhanyakude from the list |
| 7 | ✅ | Table 1 | Convenience-sample and "Matched modelled estimates" rows are no longer bold. | |
| 8 | ✅ | Table 1 | ‡ is now on the calendar-time and interaction rows ✅. **Still open:** the Convenience, FSW and MSM rows are indented one level less than the other study types, and there is no "Key populations" sub-header. | Fix indentation; add sub-header **Fifth pass:** Convenience/FSW/MSM rows indented; "Key populations" sub-header row added. |
| 9 | ❌ | Results ¶1 | **[C]** Per-country Ns sum to 127 (>109) because multi-country studies count in each. | Add "(studies may contribute to more than one country)" |
| 10 | ❌ | Results ¶1 | **[C]** "20,400 infections"; output is 20,351. | Consider "20,351" |
| 11 | ❌ | Abstract / Discussion | **[C]** "measurements 1990-2023" and "up to 2023": observation midpoint years run 1990–2022. | Check |
| 12 | ✅ | Results ¶2 | **[C]** "empirical estimates of incidence among men **from population-representative studies**": the male coefficient is pooled over all study types. | Drop the qualifier **Fifth pass:** Qualifier dropped. |
| 13 | ✅ | Results ¶3 | Now "women aged 15-24" / "women aged 25+". | |
| 14 | ❌ | Results ¶4 | **[C]** "and faster in PHIA surveys (interaction aIRR per year 0.92, 95%CrI 0.82, 1.04)": the CrI includes 1, but the clause follows "significantly". | Add "though not significantly" |
| 15 | ✅ | Methods | Ref 40 is now UNAIDS AIDSinfo (accessed Nov 2026), which supports "Global HIV Epidemic Estimates published in July 2026". | |
| 16 | ✅ | Introduction ¶1 | **[B]** "declined by nearly 70% between 2010-2024 ^(1,2)". Refs 1 (2021 strategy) and 2 (2022 update) can't support a 2024 figure. | Cite the 2025/2026 Global AIDS Update, or ref 40 **Fifth pass:** Refs 1–2 replaced by AIDSinfo + Spectrum files (2026 estimates); pooled SSA 15–49 incidence fell 70% 2010–2024 (now computed in `paper_numbers.R`). References renumbered 1–67. |
| 17 | ◐ | Figure 4 caption | Now gives the reason Namibia isn't plotted ✅. **Still open:** the caption doesn't say that NAMPHIA 2017 *is* included in R = 0.91 and 19/23. | Add "but is included in the correlation and the 19/23 count" |
| 18 | ✅ | Figure 3 caption | **[C]** "cohort data (orange)": the points are yellow, and the Fig 4C / S3 captions say "yellow". | "yellow" **Fifth pass:** "yellow". |
| 33 | new | Abstract | **Verified:** the new "declined 4% faster per year (95%CrI 1, 6%)" matches `paper_numbers.csv` (4%, 1–6%). | |
| 34 | new | Results ¶3 | **Verified:** the new cohort-decline CrIs (55% [47, 62], 73% [64, 80], 68% [52, 80], 71% [61, 79]) match `paper_numbers.csv`. | |
| 35 | ❌ | Methods, Regression | **[C]** Empty Word equation objects remain inline after tracked edits. "sex interacted with , year" (a comma inside an equation object), and invisible ones before "relative to modelled estimates" and after "age group". | Delete the equation objects; → "sex interacted with year" |

**Verified (no change needed):** 3,509 obs / 175 studies / 25 countries; 2.5 m py; n=2,556 (73%) / n=634 (18%);
42 cross-sectional studies (23 PHIA); 29 WCA (12 total population); 99% (3,466/3,509) sex-stratified;
median year 2014 (IQR 2006–2016); R = 0.61 (0.57, 0.64); every Table 1 IRR, N and random-effect SD; "7–32% higher";
men −12% and −2.3%/year; cohort-only 1.10 / 0.90 / 0.88; measurement-method 0.71 / 0.89; median 70%; cohort age
model 0.94, 1.04, 0.96; PHIA 0.93, 0.92, 10%, 17%; PHIA R = 0.91 / 0.94 [14] / 0.46 [9]; 19/23; 58% / 8% per year.

### Supplement

| # | Status | Where | Issue | Fix |
|---|---|---|---|---|
| 19 | ✅ | Table S6 fn 7 | Size parameter now 11.13 (9.39, 13.23). | |
| 20 | ✅ | Title | Now matches the main title. | |
| 21 | ✅ | Text S1 | November 2026; Joshi (2021); search terms → S8. | |
| 22 | ✅ | Text S5 (primary) | **[B]** The equation still omits the **sex × year** term (β₄S·t), and ι_o has no distribution line. Unbalanced bracket in "(log σ ~ N(−1.97, 0.32²)." "weighted by the observed person-years" (see §2). | Add β₄S·t and ι_o ∼ N(0, σ_ι); close the bracket **Fifth pass:** Equations rewritten in Word linear format (select each, Equation → Professional). β₄S·t, ι_o line, t_s definition, closed bracket, offset wording. |
| 23 | ✅ | Text S5 (age model) | **[B]** The equation omits the cohort-specific random slope (now described in the text and in Table S9). It still has "ϵ_c ∼ N(0, σ_ϵ)" for ϵ_(ca), and η has no sex index, but Table S9 reports it by sex. Bracket "a∈{15-24, 25-49)"; the label should also now be 25+. | Fix the equation **Fifth pass:** Linear-format equations: ϵ_(ca), η_(tcS) with sex-specific AR2, γ_(ta), cohort slope δ_c x_a t with x = ±1/2, a∈{15-24, 25+}. |
| 24 | ✅ | Text S5 (PHIA) + main Methods | **[C]** "random slopes for country over year **and binary age group**". Table S10 reports only country intercepts and a country time trend. | Align **Fifth pass:** "and binary age group" dropped in main Methods and S5. |
| 25 | ➖ | Text S6 | **[C]** "a smoothing random effect" in the text vs "Year (AR1)" in the table. | Name AR1 in the text **Fifth pass:** No change (Oliver). |
| 26 | ◐ | Table S3 | **[B]** uMkhanyakude "2005-2023" vs 2004–2022 elsewhere; Siaya "2008-2021" vs midpoints 2011–2021. | Align **Fifth pass:** uMkhanyakude → 2004-2022. Siaya "2008-2021" kept as is (Oliver). |
| 27 | ✅ | Table S5 footnote | **[C]** Still says "study method" sums to >100% (it sums exactly). | Drop "study method" **Fifth pass:** "study method" dropped. |
| 28 | ◐ | Text S9 | Wahome 2020 participation is now NA (179/172 was not a participation rate) ✅. The text "78 of 175" and the figure are regenerated ✅. **Still open:** the figure has no caption or number, colours are explained only in text, and "Supplementary File S2" must exist. | Add caption/number; confirm File S2 |
| 29 | ◐ | PRISMA tables | Items 7, 11, 18, 20a, 13f, 20d and 23d now point to the right places ✅. Line numbers have been refreshed, but they can only be checked in Word with line numbering on. **Still open:** 24a–c say "N/A", but PRISMA asks you to *state* that the review was not registered or no protocol was prepared, and **Text S1 says "As defined in the protocol"**, which contradicts 24b. Abstract checklist item 9 (limitations) is "N/A"; it should be "No". 16b now cites "Supplementary File S3"; confirm it exists. | "Not registered" / point 24b to the protocol or drop "As defined in the protocol" **Fifth pass:** 24a–c → "Not registered"; "As defined in the protocol" dropped from S1. Abstract item 9 "N/A" and the line numbers not changed. |
| 30 | ❌ | TOC | **[C]** Supplementary Figure S1 is missing from the TOC. Text S4 still reads "1990-2024", in the TOC and in the heading. | |
| 31 | ❌ | Text S8 | **[C]** "equitorial guinea" (Embase) and "Democratic Repubic of Congo" (Scopus) are verbatim typos in the executed strings; keep them but be ready to explain. The Scopus/OVID import date is "feb 23" vs 22 Feb in the main text. The PubMed range reads "12/31/3000 (feb 22)". | |
| 32 | ❌ | Fig S2 | **[C]** Panel labels ("ANC", "Intervention Trials", "Non-representative - Other", "No inclusion criteria", "UTT/Cohort") still don't match Table 1 names. The caption also says "ANC regression coefficient". | Relabel in `paper_figures.R` (`loo_panel` labels) |
| 36 | new | Table S7 | **[C]** FSW and MSM rows are absent while the header matches Table 1. Check that this is intended (measurement method fitted to total-population studies only) and say so in a footnote. | |

---

## 2. Copy edit

### Main text
- ✅ **Affiliations:** the blank entries #19 and #21 are tracked deletions. Once accepted, the list renumbers to
  19 SACEMA, 20 JHU Pathology (Grabowski), 21 Harvard CCDD (Imai-Eaton), 22 Deakin. **Confirm** that Bock "4,22"
  = Deakin University is intended.
- ✅ "include"; "MeSH"; no U+2010 hyphens left in the main text.
- ✅ Methods now says "resolved through consensus", matching S1.
- ✅ "Population-based HIV Impact Assessment" used throughout.
- ✅ Kassanjee 2012 estimator paper now cited (ref 38).
- ◐ Methods (Regression): offset wording fixed in main text and S5 (fifth pass) ✅. Still open: "weighted by observed person-years" → "with log person-years as an offset" (same in
  Text S5). "This captured both structured and random heterogeneity between studies." still follows the
  slope-prior sentence.
- ◐ Sensitivity analyses: "used to derive" ✅ and "2015-present; (4)" ✅. Still "(4) ... cohort study data only, (5)":
  comma → semicolon.
- ✅ Results ¶2: parentheses balanced, aIRR defined once, "FSW aIRR", "other inclusion criteria",
  "leave-one-out analysis".
- ✅ Limitations: "leave-one-out analysis indicates".
- ✅ Rate notation now "aIRR per year" / "aIRR per five calendar years" throughout; "95%CrI" everywhere and
  "95%CI 0.57, 0.64".
- ✅ "person-years" and "mixed-effects" consistent.
- ◐ "Côte d'Ivoire" fixed in Fig 4 caption ✅. "eSwatini" → "Eswatini" still open in the Fig 4 caption, the
  Fig 4C panel (`country_name()` in `R/functions.R`), Text S4 tables and Fig S7. Ref 65's title is a quotation; leave it.
- ◐ Cohort naming: text consistent; Figure 1 still says "AHRI" and "KEMRI, Kenya" (item 1).
- ✅ "Reference Group on Estimates, Modelling, and Projections" consistent; ref 23's title quotes the original.
- ✅ "Figure 1:" colon. ✅ "postpartum". ✅ "3,509". ✅ "ageing".
- ❌ Conclusion: "the control arms of **most** HIV prevention/vaccine efficacy trials" → drop "most".
- ❌ **Figure 4 alt text** is still auto-generated ("A graph of different sizes and colors…").

### Supplement
- ✅ S1: "standardised", "CI" for cumulative incidence, "K.P.", "MeSH". ("synthesize" remains only in the
  verbatim PRISMA checklist wording, which is fine.)
- ✅ S2: "were available", "impute them", "five-year".
- ✅ S3: "more than one location".
- ◐ S4: "Guinea-Bissau" consistent ✅. "Multiply **Naomi** 2025 estimates": Naomi is still undefined ❌.
- ✅ S9: "follow-up"; model now dated ("GPT-5.6 Terra, October 2026").
- ✅ Table S6: title now above the table.
- ✅ Table S8: two decimals throughout; thousands separators; "population-representative" in the title.
- ◐ Tables S9/S10: "Sex" header ✅; "25+" now matches the text ✅. "Male: 25+" vs "Male : year" spacing ❌.
- ✅ Fig S4 caption: "sex". Fig S3B caption now "Female:male HIV incidence rate ratio", matching the axis.
- ❌ Tables S7/S12: header "(N studies)" but no Ns.

---

## 3. References

Fifth pass: refs 1–2 (2021 strategy, 2022 update) removed; AIDSinfo and the Spectrum files are now refs 1–2 and everything from old ref 42 moves down by two (67 refs). Citation order checked by script: first appearances run 1–67. Reference metadata edits and main-text edits are tracked changes by "Claude"; the renumbering is not tracked. Mendeley shows the bibliography as manually edited, so a Mendeley refresh would overwrite these fixes.

### Citation–content mismatches
- ❌ **[A] Refs 1, 2:** see item 16.
- ◐ **[A] Ref 40 (now ref 1):** now AIDSinfo ✅. **Fifth pass, new [A]:** in the 2026 Spectrum files the per-country median 15–49 incidence decline is 65% for 2010–2023, 69% for 2010–2024 and 75% for 2010–2025, so the Discussion's "median 75% ... between 2010 and 2023" should read 2010–2025 (or 65%). ✅ Changed to 2010–2025 (tracked); `paper_numbers.R` now calculates it. The Discussion's "median 75% decline across countries between 2010 and 2023"
  (refs 40, 41) still reads like your own calculation from 2026 Spectrum files; say so, and check whether it should be 2010–2025.
- ✅ **Ref 38:** Kassanjee 2012.
- ➖ **Ref 42:** left as is at Oliver's request.
- ◐ **Ref 46:** metadata fixed (UNAIDS, access date; still no year). It still sits at the end of the
  geocode-and-assign sentence, after "(Supplementary Text S3)".
- ✅ **Ref 56:** Gray 2024 *Lancet Infect Dis*.
- ✅ (fifth pass: text now "MSM ... and FSW") **[C] "FSW ... and MSM ... ^(36,37)":** 36 = Stannah (MSM), 37 = Jones (FSW), the reverse of the text order.
- ✅ **Refs 41, 45:** access dates Nov 2026.
- Borgdorff 2018 (now ref 68) could be co-cited with ref 29 for Siaya (optional).

### Metadata fixes
| Ref | Status | Fix |
|---|---|---|
| 1, 2, 4, 6, 8, 18, 19, 21, 22, 23, 25, 39, 44, 48, 51, 54, 59, 63 | ✅ | Done |
| 9, 16, 17 | ✅ | Add issue: 33(Suppl 3) |
| 12 | ✅ | e0133255 |
| 13 | ✅ (PubMed) | 66(Suppl 1):S13-S26† |
| 15 | ✅ | e34410 |
| 27 | ✅ | Abstract 150, CROI 2020, Boston (per Healio report of the presentation) |
| 28 | ✅ | dyae058 |
| 31, 33, 35, 55, 58 | ✅ | N Engl J Med; Havlir DV; DOI suffixes stripped |
| 36 | ✅ | DOI, Lam JKS, volume/pages; journal italicised |
| 46 | ◐ | Add a year: not done, the ADR dataset is updated annually and the version year couldn't be confirmed |
| 49, 67 (was 65), 6 | ✅ | Italics now on the journal name only |
| 53, 64 | ✅ | Verified on PubMed |
| 56 | ✅ | Journal italicised |
| 62 | ✅ | "Beamon T" correct (given name Ta'Adhmeeka) |

Duplicate records in the Mendeley library that could confuse the citation plugin: ten Brink 2025 (×2),
DREAMS/Saul 2018 (×2), subnational viewer (×2), Havlir JIAS 2020 (×2), Facente (preprint plus published).

---

## 4. Oliver's review list (17 Nov 2026, second pass): status

- ❌ Methods, Regression: the placement of "This captured both structured and random heterogeneity…".
- ❌ Tables S7 and S12: "(N studies)" header without Ns.
- Supplement TOC: ❌ Figure S1 missing; ❌ Text S4 "1990-2024"; refresh page numbers at the end.
- ❌ Figure 4 alt text.
- Other open items:
  - ◐ PRISMA 16b/24a–c: see item 29 (the protocol contradiction is new).
  - ❌ Confirm Supplementary Files S1, S2 **and S3** exist.
  - ❌ Text S9 figure caption/number.
  - ❌ Table S5 footnote "study method".
- ✅ Reference issues (third pass): refs 2, 21, 41, 45, 46, 63 fixed; ref 62 correct; ref 42 left.
