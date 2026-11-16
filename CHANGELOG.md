# CHANGELOG

Created 2026-11-15.

## 2026-11-16
* Pipeline completed; analysis (tables 1–12, S1–S17; figures 1–7, S1–S2; metrics/statistics JSON).
* Manuscript v1 assembled from token templates (no hand-typed results); ACA structure (structured abstract, highlights,
  declarations incl. generative-AI statement); DOCX via pandoc + reference style; PDF via LibreOffice.
* Corrections found during assembly: share of seen test compounds (73 %, not 81 %), HILIC training-system counts, GP miss factor,
  α-dependence of parametric miscalibration, decoy counts (per system, not per compound), unsupported literature figures removed.
* Rendering: portable math (no \lceil, \Big, bare ± superscripts), typographic minus, numbered sections, page numbers,
  content-proportional table widths, deterministic figure PDFs.
* NORMAN per-stratum coverage now reported per calibration scenario (previously pooled).
* Supplementary Material builder; supplementary tables renumbered S1–S26 in document order; manuscript references updated.
* README: pandoc/LibreOffice requirements for the document stage; run_all builds the Supplementary Material.
* Audit-driven revision (13:30): see REVIEW/claim_evidence_audit.md and REVIEW/citation_audit.md (resolution sections).
* Comparator B7 (meta-learned GP prior) added post hoc (Addendum E) and reported in all main tables and figures.
