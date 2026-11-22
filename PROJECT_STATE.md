# PROJECT_STATE.md

Last updated: 2026-11-22.

## Status
Submission-ready for **Analytica Chimica Acta** (Elsevier). Final package:
`Target_Journal_Submission/submission_package/` (and the zip next to it).

Title: *Laboratory-level guarantees for retention-time filtering in non-target screening by cross-system conformal calibration*

## Author-side items still to fill in before upload
* Second author name, affiliations of both authors, corresponding e-mail (manuscript title page, cover letter, declarations).
* Date in the cover letter.
* E-mail addresses of the suggested reviewers (Editorial Manager requires a working address).

## Rebuilding the documents
```
cd Target_Journal_Submission/code
pip install -r requirements.txt
python3 scripts/build_manuscript.py
python3 scripts/build_supplementary.py
python3 scripts/build_cover_letter.py
python3 scripts/make_submission_package.py
```
All numbers in the documents are resolved from the result files in `results/` and `tables/`; unresolved tokens abort the build.

## Data
RepoRT v1.0 (Zenodo 10.5281/zenodo.16267768), NORMAN interlaboratory RT data (Malm & Kruve), PubChemLite for Exposomics
2026-09-25 (Zenodo 10.5281/zenodo.22953038). Download: `python3 code/scripts/download_data.py`.
