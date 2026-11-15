# Selected target journal — Analytica Chimica Acta (Elsevier)

Decision date: 2026-11-15. Decided **before** the research topic was finalised (Phase 0), then re-checked
against the chosen topic (retention-time uncertainty for LC–HRMS non-target screening; see `../../PROJECT_STATE.md`).

| Field | Value | Evidence / status |
|---|---|---|
| TARGET JOURNAL | **Analytica Chimica Acta** (ACA) | — |
| PUBLISHER | Elsevier B.V. (journal founded 1947) | Official journal page (shop.elsevier.com/journals/analytica-chimica-acta/0003-2670) |
| ISSN | 0003-2670 (print); 1873-4324 (online) | JCR export; official page |
| JCR STATUS | Indexed, SCIE | JCR export "Selected Categories: CHEMISTRY, ANALYTICAL; Selected JCR Year: 2025" (`data/JCR2025_Chemistry_Analytical_export_UCM_08-2026.pdf`, sha256 d31d577a…) |
| JCR QUARTILE | **Q1** — Chemistry, Analytical; rank 23 of 111 | same JCR export (latest release, June 2026, JCR year 2025) |
| JOURNAL IMPACT FACTOR | **6.1** (2025 JIF); 5-year JIF 5.8; JCI 1.43 | JCR export + official Elsevier page (both report 6.1) |
| JCR SUBJECT CATEGORY | Chemistry, Analytical | JCR export |
| VERIFIED REVIEW/DECISION TIMELINE | First-decision time **not independently verified** (Elsevier "journal insights" page blocked, HTTP 403). | — |
| VERIFIED ACCEPTANCE TIMELINE | Submission(received)→acceptance, computed from PubMed article-history dates: **all research articles published Nov 2025–Oct 2026: median 60 d (IQR 43–81), 90.5 % ≤100 d, 50.2 % ≤60 d (n = 400)**. **Topic-relevant articles (ML/chemometrics/retention/uncertainty) in the same window: median 72 d (IQR 53.5–96), 77.5 % ≤100 d (n = 173).** | `scripts/pubmed_timelines.py`, `scripts/topic_timelines_window.py`; `data/*.csv`. These are **accepted** papers only (survivorship), so they say nothing about rejection rates and cannot guarantee any individual manuscript's timeline. |
| APC | Hybrid journal: **no charge** for subscription publication; optional gold OA ≈ USD 4,030 (third-party directory listing updated 2024-01-17; **not verified on the official site**) | Cirad "Où publier" listing |
| SCOPE FIT | Official scope: "fundamental research … including chemometric techniques"; "data treatment that truly demonstrate new principles and/or important applications" is especially encouraged; requires "proper statistical treatment of data"; new methods must show "obvious and quantifiable advantages over existing approaches". | Official Elsevier journal page (fetched 2026-11-15) |
| RESEARCH FIT | ACA published 11 retention-time-prediction papers (2023–2026), 39 non-target-screening papers and 63 chemometrics papers (2024–2026) (PubMed counts, `data/topic_fit_counts_2024-2026.json`, `data/rt_prediction_papers_by_journal_2023-2026.json`). No conformal/coverage-guaranteed retention-time paper found in any journal (PubMed query, 2026-11-15). | — |
| KEY AUTHOR GUIDELINES (as verifiable) | Structured abstract Background / Results / Significance (used by ~48 % of 150 most recent ACA articles; median 225 words, max ≈300); ~5 keywords (median 5, max 8); Highlights and graphical abstract (standard ACA article furniture); Elsevier declarations: competing interests, CRediT author contributions, data availability statement, **Declaration of generative AI and AI-assisted technologies** (mandatory under Elsevier policy). Official Guide for Authors could not be fetched (ScienceDirect 403) → re-check before submission. | PubMed XML of 150 recent ACA records; Elsevier policy |
| KEY REJECTION RISKS | (1) Desk rejection for insufficient novelty/impact ("particular scrutiny"); (2) perceived as "application of an existing ML tool"; (3) lack of experimental (wet-lab) validation; (4) need to demonstrate quantifiable advantage over existing RT-window practice; (5) scope perception that pure data-treatment papers belong in chemometrics journals (mitigated: chemometrics journals are Q2/Q3 and ACA scope explicitly invites data treatment). | — |

## Why ACA and not the alternatives (short form; full table in `journal_comparison.md`)
* **Analytical Chemistry (ACS, JIF 7.3, Q1, rank 12/111)** — best topical home for retention-time work, but its verified
  topic-specific median is **95 d** with only **56.6 %** of topic papers accepted within 100 d (overall: 97 d; 55.0 %).
  Under the user's strict 100-day ceiling this is a coin flip, so it was not selected (kept as a fallback).
* **Talanta (Elsevier, JIF 6.7, Q1, rank 15/111)** — same speed as ACA (topic median 71 d, 77.3 % ≤100 d) and higher JIF,
  but its official scope emphasises sensors/instrumentation and "demonstrated analytical applicability" and does not
  mention chemometrics/data treatment, whereas ACA's scope names both explicitly. Fit is therefore **not** comparable for
  a computational data-treatment paper, so the higher JIF does not dominate. Talanta is the first fallback.
* **All Oxford University Press candidates** (Bioinformatics, Briefings in Bioinformatics, JAMIA, GigaScience, NAR Genomics
  & Bioinformatics, Bioinformatics Advances, EHJ Digital Health, PNAS Nexus) have verified medians of **118–198 d** — all
  above the 100-day ceiling.

## Requirements that could NOT be fully verified or satisfied
* No journal with a verified ~30-day median submission-to-acceptance time and adequate scope/JIF was found
  (fastest Q1 candidate with fit: ACA, 60 d overall / 72 d topic-specific).
* First-decision time, acceptance rate and current Guide-for-Authors text for ACA were not accessible (HTTP 403) — **unverified**.
* APC amount is from a third-party directory (2024) — **unverified**.
* Acceptance-time statistics describe accepted papers; **they do not guarantee acceptance or a timeline for this manuscript.**
