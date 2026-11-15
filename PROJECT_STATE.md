# PROJECT_STATE.md — persistent research checkpoint

Last updated: 2026-11-15 (session 1)

## Environment
- 4 CPU cores, 15 GB RAM, no GPU (no CUDA, no MPS). ~31 GB disk.
- Python 3.13; pip installs work; pandoc + LibreOffice available.
- Reachable: Crossref API, PubMed E-utilities, OpenML, UCI, PhysioNet, Zenodo, arXiv, UCR archive, HuggingFace.
- ScienceDirect HTML pages return HTTP 403 (cannot scrape Elsevier "journal insights" directly).

## Phase status
- [ ] Phase 0 — journal selection (IN PROGRESS)
- [ ] Phases 1–4 — topic discovery / novelty
- [ ] Phases 5–10 — design / implementation
- [ ] Phases 11–20 — experiments
- [ ] Phases 21–23 — manuscript / audits
- [ ] Review cycles (10 agents x 10 cycles)
- [ ] Final packaging

## Exact next step
Compute verified received->accepted intervals for candidate journals from PubMed/Crossref history metadata.
