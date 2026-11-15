# Journal verification log (Phase 0)

## 1. Clarivate JCR status, quartile and JIF of the selected journal
* Source: Clarivate Journal Citation Reports export "Journal Data Filtered By: Selected Categories: CHEMISTRY, ANALYTICAL;
  Selected Editions: SCIE;SSCI;AHCI;ESCI; Selected JCR Year: 2025" hosted by Universidad Complutense de Madrid
  (https://www.ucm.es/quimica_analitica/file/anal-chem-jcr-journalresults-08-2026), archived here as
  `data/JCR2025_Chemistry_Analytical_export_UCM_08-2026.pdf` (sha256 d31d577aa78d4443005fbc45a14d9134cb0c946880f830b3eb4c848dd03363aa).
* Extracted rows: Analytica Chimica Acta — ELSEVIER, ISSN 0003-2670 / eISSN 1873-4324, SCIE, total cites 63,897, **2025 JIF 6.1, Q1**, JCI 1.43, 15.50 % OA.
  Category size 111 journals → ACA rank 23 (Q1 = ranks 1-28).
* Cross-check: the official Elsevier journal page (shop.elsevier.com/journals/analytica-chimica-acta/0003-2670, fetched
  2026-11-15) reports "Impact factor: 6.1; 5-year impact factor: 5.8". Consistent.
* Clarivate's own JCR portal requires a subscription and could not be accessed directly; the export above is a
  university-hosted copy of the official JCR output (latest release, June 2026).

## 2. Acceptance timeline (submission → acceptance)
* Method: PubMed E-utilities, `esearch` for `"Anal Chim Acta"[ta] AND journal article[pt]` excluding reviews, editorials,
  comments, letters and errata; publication date window 2025/11/01–2026/10/31; most recent 400 records; `efetch` XML;
  interval = `accepted` − `received` history dates (0–2000 d kept).
* Result (all research articles): n = 400 with both dates; **median 60 d; IQR 43–81 d; mean reported in CSV; 50.2 % ≤60 d; 90.5 % ≤100 d**.
* Result (topic-relevant articles, same window; title/abstract containing machine learning / deep learning / neural network /
  chemometric* / partial least squares / retention time / uncertainty / prediction interval*): n = 173;
  **median 72 d; IQR 53.5–96 d; 35.3 % ≤60 d; 77.5 % ≤100 d; 90th percentile 137 d**.
* 2024-01 … 2026-10 topic window: n = 299; median 82 d; 67.6 % ≤100 d (slower in earlier years).
* These are **medians of accepted papers**, not an acceptance guarantee; first-decision time and acceptance rate are **unverified**.

## 3. Scope (official text, quoted fragments, Elsevier page fetched 2026-11-15)
* "An International Journal Devoted to All Branches of Analytical Chemistry."
* Fundamental research "including chemometric techniques" is encouraged; "Data treatment that truly demonstrate new
  principles and/or important applications" is especially encouraged; "proper statistical treatment of data" required.
* "particular scrutiny" of "the degree of novelty and impact"; routine methods rejected unless they offer "obvious and
  quantifiable advantages over existing approaches"; "does not focus on application papers".

## 4. Publisher legitimacy / integrity
* Elsevier B.V.; journal published since 1947; indexed in SCIE; not on any predatory list; no Clarivate suppression or
  "on hold" notice found for this title. (No integrity concern identified.)

## 5. Recent-paper fit (PubMed counts, 2024–2026)
* ACA: 3,304 records total; ML/DL/NN 177; chemometric* 63; non-target screening 39; retention-time prediction 8 (11 in 2023–2026);
  calibration transfer/domain adaptation 8; conformal/prediction interval/UQ 0.
* Full lists: `data/topic_fit_counts_2024-2026.json`, `data/rt_prediction_papers_by_journal_2023-2026.json`.

## 6. Reproducibility of this verification
```
python scripts/pubmed_timelines.py data/pubmed_timelines_2025-11_2026-10.csv
python scripts/pubmed_timelines.py data/pubmed_timelines_batch2_2025-11_2026-10.csv 2025/11/01 2026/10/31 batch2
python scripts/crossref_timelines.py data/crossref_timelines_2025-11_2026-10.csv
python scripts/topic_timelines.py data/topic_specific_timelines_2024-2026.csv
python scripts/topic_timelines_window.py data/topic_specific_timelines_2025-11_2026-10.csv 2025/11/01 2026/10/31
python scripts/topic_fit_counts.py data/topic_fit_counts_2024-2026.json
```
(Network queries to NCBI/Crossref; results drift as new papers are indexed.)
