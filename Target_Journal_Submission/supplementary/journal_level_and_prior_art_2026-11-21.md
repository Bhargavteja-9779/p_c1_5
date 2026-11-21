# Journal-level calibration and prior-art checks (2026-11-21)

All records come from PubMed E-utilities, Crossref, Europe PMC and arXiv, queried on 2026-11-21. The coding below is based on abstracts only. Scripts, query files and raw JSON/XML are in `/tmp/claude-0/-home-user-p-c1-5/106b1790-cbf3-5c4e-815f-e1b9280f47cc/scratchpad/litscan/` (`scripts/`, `raw/`).

Main Task A searches:
- PubMed: journal filter (Anal Chim Acta, Anal Chem, Talanta, J Chromatogr A, J Cheminform) AND date 2025/01/01–2026/11/30. This was combined with three topic blocks: RT prediction (89 hits), NTS or annotation with ML (130 hits), and uncertainty or conformal (31 hits).
- Crossref: the same ISSNs, searched by relevance.

## A. Journal-level calibration (2025–2026)

Codes used in the table:
- **Data:** C = computational on public data only; O = the authors' own measurements; R = a real-sample case study.
- **Novelty:** Alg = new algorithm; Mod = new model or descriptor; WF/Tool = new workflow or software; Bench = benchmark; Data = new dataset.
- **ID gain:** a practical identification gain shown in the abstract.

| # | Paper (short) | J, yr | DOI | Data | Scale | Novelty | ID gain shown |
|---|---|---|---|---|---|---|---|
| 1 | Graphormer-RT, method-independent RT | AC 2025 | 10.1021/acs.analchem.4c05859 | C | RepoRT: 191 RP + 49 HILIC methods | Mod | No (MAE only) |
| 2 | MDL-TL: multi-dataset transfer learning | AC 2025 | 10.1021/acs.analchem.5c01703 | C | 28 datasets | Mod | No |
| 3 | Uni-RT multitask RT | AC 2026 | 10.1021/acs.analchem.5c07973 | C | 28 datasets | Mod | No |
| 4 | MC-Retention, uncertainty-aware multi-condition RT | AC 2026 | 10.1021/acs.analchem.5c06729 | O | 480 analytes × 8 conditions | Mod+Data | Method development, not ID |
| 5 | ROASMI retention-order model with applicability domain | JCh 2025 | 10.1186/s13321-025-00968-8 | C+R | 71 RPLC sets + 4 real datasets | Mod | Yes: isomers separated, peaks annotated |
| 6 | TC-TopoRT | Tal 2027 | 10.1016/j.talanta.2026.130639 | C | SMRT + 6 PredRet | Mod + calibrated filter | Yes: −41 %/−30 % candidates, 93 %/99 % of true candidates kept, top-k kept |
| 7 | ML for identification probability (RTI + MS2) | AC 2025 | 10.1021/acs.analchem.5c01873 | O+R | 4,713 calibrants; 1.7 M matches; spiked tea | WF/Mod | Yes: IP +47–55 % in matrices |
| 8 | LC RI with two-phase selectivity filtration | AC 2025 | 10.1021/acs.analchem.4c01784 | O+R | 116 compounds, 2 columns | WF | Yes: >70 % of impossible structures removed; drinking water |
| 9 | Retip RT + HRMS, Schisandra–Evodia | AC 2025 | 10.1021/acs.analchem.5c05620 | O+R | 165 compounds, in vivo | WF | Yes: fewer false positives at ±1 min |
| 10 | AI-SNPS2 (NPS screening, RT module) | AC 2025 | 10.1021/acs.analchem.5c02830 | O+R | 164 substances, spiked matrices | Tool | Yes: analogues identified |
| 11 | ML-enhanced pesticide suspect/NTS workflow | AC 2026 | 10.1021/acs.analchem.5c06726 | O+R | Fortified + real samples | WF (uncertainty budget) | Yes: >6 % more identifications |
| 12 | MetaboAnnotate (LLM multi-agent) | AC 2026 | 10.1021/acs.analchem.6c01717 | C+R | CASMI 2016/22 + faecal set | Tool | Yes: lower FDR by consensus |
| 13 | ML RT for lipid identification | JCA 2025 | 10.1016/j.chroma.2024.465650 | O | Lipidomics, external sets | Mod | Claimed annotation gain |
| 14 | RT-Pred web server | JCA 2025 | 10.1016/j.chroma.2025.465816 | C | >40 CMs | Tool | No quantitative gain |
| 15 | Pharma RT across 20 RP columns (HSM) | JCA 2025 | 10.1016/j.chroma.2024.465628 | O | 51 compounds × 20 SPs | Mod+Data | Method transfer |
| 16 | Prior knowledge / Bayesian D-optimal design | JCA 2025 | 10.1016/j.chroma.2025.465787 | O/C | Multi-analyte | Alg (design) | Precision |
| 17 | Bayesian multilevel GP (Wiczling) | JCA 2026 | 10.1016/j.chroma.2025.466503 | C | 1,026 analytes | Alg | No |
| 18 | GA-MLR QSRR, plant bioactives | JCA 2025 | 10.1016/j.chroma.2025.466194 | C/O | 3 LC systems | Mod | No |
| 19 | Chiral cellulose RT (PLS/NN/GNN) | ACA 2025 | 10.1016/j.aca.2025.344733 | O | 535 mol., 1,414 RTs | Mod+Tool | Method development |
| 20 | Multi-condition SFC RT | ACA 2026 | 10.1016/j.aca.2025.345026 | O | 1,217 cpds × 51 conditions | Mod+mechanism | Indirect |
| 21 | Project-specific QSRR from miniature datasets | ACA 2026 | 10.1016/j.aca.2026.346198 | O | 3 × 17–21 cpds | WF | Impurity RTs |
| 22 | IAM QSRR: descriptors vs GNN | ACA 2026 | 10.1016/j.aca.2026.346201 | O | 2,003 analytes measured | Bench+Data | No |
| 23 | XGBoost QSRR with stereo descriptors | JCA 2026 | 10.1016/j.chroma.2026.467239 | O | 43,329 RTs, 13 SPs | Mod (descriptors) | Stereoisomers 2.5× |
| 24 | "Identification probability" (perspective) | AC 2025 | 10.1021/acs.analchem.4c04060 | – | – | Concept | – |
| 25 | Measurement precision vs ID probability | AC 2025 | 10.1021/acs.analchem.5c01067 | C | Metabolite space | Framework | Quantifies ID probability |
| 26 | TIMS improves exposomics annotation | AC 2025 | 10.1021/acs.analchem.5c04665 | O+R | 769 cpds, 5 matrices | Experimental eval. | Yes: fewer inconclusive assignments |
| 27 | INTERPRET NTA candidate prioritisation | AC 2025 | 10.1021/acs.analchem.5c02223 | O+R | 77 known chemicals | Tool | Ranking of known chemicals |
| 28 | Measurable-feature prediction (chemical space) | AC 2026 | 10.1021/acs.analchem.5c07705 | O+C | IS + CompTox | Framework | Coverage estimate |
| 29 | Active learning for IE prediction | AC 2025 | 10.1021/acs.analchem.5c00816 | C/O | IE sets | WF | Quantification |
| 30 | IE from cumulative neutral losses | JCh 2025 | 10.1186/s13321-025-01129-7 | C | – | Mod | No |
| 31 | Hierarchical Bayesian CCS cross-platform correction | AC 2026 | 10.1021/acs.analchem.5c06667 | O (multi-lab) | 840 meas., 347 cpds | Alg | 94.6 % coverage of 95 % PI; ≥3 reference cpds |
| 32 | HyperCCS (LLM CCS) | AC 2025 | 10.1021/acs.analchem.5c03492 | C | METLIN-CCS, AllCCS2 | Mod | No |
| 33 | CCSBase2 | AC 2026 | 10.1021/acs.analchem.6c03711 | C | 66,153 CCS | Data+Mod | No |
| 34 | Group-conformal NIR intervals under shift | Tal 2027 | 10.1016/j.talanta.2026.130941 | C | 146,873 spectra, 12 shifts | Alg | Yes: fewer false specification acceptances |
| 35 | Calibrated conformal UQ, Raman SOC | AC 2026 | 10.1021/acs.analchem.5c04616 | O | SERDS field data | Framework | Coverage |
| 36 | Cross-platform LC-MS classifier transfer (honey) | ACA 2026 | 10.1016/j.aca.2026.345986 | O | 262 samples, 4 platforms | Alg | Yes: cross-instrument MCC |

Journal abbreviations: AC = Anal Chem; ACA = Anal Chim Acta; Tal = Talanta; JCA = J Chromatogr A; JCh = J Cheminform.

**What sets accepted papers apart from a "benchmark of an existing ML method on public data":**
1. **Own data.** About 23 of the 36 papers measure their own data or run a real-sample case study. Every ACA paper in the set has its own measurements.
2. **Public-data-only papers that were accepted bring something new.** In Anal Chem these are new architectures at RepoRT scale (#1–3). Talanta and J Cheminform also accept public-data papers, but only with a measured identification endpoint (#5, #6, #34).
3. **An identification endpoint.** RT/NTS papers in Anal Chem report one of: fewer candidates at a stated retention of true structures (#6: −30–41 % at 93–99 %), more identifications or higher identification probability (#7, #8, #11: +6 % to +55 %), or fewer false positives in real matrices (#9, #10, #12, #26).
4. **Pure benchmarks are rare.** The only one (#22) adds a new measured dataset of 2,003 analytes.
5. **The uncertainty papers that were accepted propose an algorithm adapted to the analytical shift.** Examples are group-conformal (#34) and hierarchical Bayesian with ≥3 reference compounds (#31). Neither applies an off-the-shelf method. Each also ties coverage to a decision outcome, such as false acceptances or probability-based matching.

## B. Prior-art checks

| Item | Queries (PubMed / arXiv / Europe PMC / Crossref; files `scripts/qB*.tsv`) | Closest hits | Verdict for analytical chemistry / RT |
|---|---|---|---|
| **B1** Weighted / covariate-shift conformal; density-ratio weights from structures | conformal AND (covariate shift OR weighted conformal OR density ratio) AND (chromatograph* OR retention OR mass spectrometry OR spectroscop* OR molecul*); arXiv abs:conformal AND "covariate shift" AND molecul*/drug; abs:jackknife AND conformal AND "covariate shift" (0) | Tibshirani et al. 2019 (arXiv 1904.06019); Prinster et al. 2022 JAWS (arXiv 2207.10716); **Laghuvarapu et al. 2023 CoDrug** (arXiv 2310.12033): KDE / energy-based density weights on molecules for drug-property CP; Fannjiang et al. 2022 PNAS (10.1073/pnas.2204569119, PMID 36256807); Jin & Candès weighted conformal p-values (arXiv 2307.09291); Lee et al. 2026 label-shift molecular CP (arXiv 2608.17678); scaffold-shift CP benchmark (10.1016/j.compbiolchem.2026.109371) | **None found** for RT, LC or analytical calibration. The method exists in cheminformatics (CoDrug), so the new part would be the application: calibrant set to suspect list with JAW/weighted jackknife+ |
| **B2** Borrowing calibration across systems/tasks (few-shot, meta, hierarchical) | conformal AND (hierarchical OR meta-learning OR few-shot OR transfer OR auxiliary) AND (chromatograph* OR retention OR spectroscop* OR chemometric*); arXiv abs:conformal AND "retention time" (0); ti:conformal AND (few-shot OR meta OR hierarchical OR transfer OR multi-source) | Fisch et al. 2021 (arXiv 2102.08898); Park et al. 2022 meta-XB (arXiv 2210.03067); Dunn, Wasserman & Ramdas 2023 JASA (10.1080/01621459.2022.2060112); Lee, Barber & Willett (arXiv 2306.06342); **Mallick, Tchetgen Tchetgen, Dobriban & Lee 2026 "Generalized hierarchical CP"** (arXiv 2608.15500): exactly "few observations from the new group", but illustrated only on census data; multi-source CP (arXiv 2405.09331, 2609.14531). Closest analytical analogues: group-conformal NIR (10.1016/j.talanta.2026.130941); hierarchical Bayesian CCS (10.1021/acs.analchem.5c06667) | **None found** in chemistry or LC. Directly targets the manuscript's 1/α−1 limit |
| **B3** Conformal p/e-values or conformal selection for candidate ranking, evidence fusion, annotation FDR | conformal AND (metabolom* OR mass spectr* OR non-target* OR suspect screening OR proteom*) AND (annotation OR candidate* OR false discovery OR retrieval); conformal selection / conformal p-value / e-value AND molecul*/metabol*; decoy AND (conformal OR e-value OR knockoff) | Conformal MS/MS retrieval sets (Rakhshaninejad 2026 JCIM, 10.1021/acs.jcim.6c00727, PMID 42113637); conformal CCS prioritisation (Ferrando Plo & Kruve 2026 JCIM, 10.1021/acs.jcim.6c02492, PMID 42789453: 16 % of wrong candidates removed at 90 % retention); CP for NTS feature prioritisation (10.1007/s00216-025-06303-2); Jin & Candès conformal selection (arXiv 2210.01408) plus later variants (multivariate 2505.00917, multi-condition 2510.08075); non-conformal target-decoy FDR (Group-walk, 10.1093/bioinformatics/btac471; 10.1111/biom.13830) | **None found** for conformal FDR or selection in metabolomics/NTS, or for a guaranteed combination of RT + CCS + MS2 |
| **B4** Calibrant (anchor) set design for RT projection | (retention time/index) AND (calibrant* OR anchor/reference compound*) AND (select* OR optim* OR design); "calibrant selection/set/mix"; active learning / D-optimal AND retention; Crossref "optimal selection of calibrant compounds…" | Aalizadeh et al. 2021 Anal Chem (10.1021/acs.analchem.1c02348): NORMAN RTI calibrants (2 × 18) chosen by RT overlap and chemical similarity out of 2,123; FastRet 2026 JCIM (10.1021/acs.jcim.6c01344): "selective measuring" chooses metabolites to remeasure for model adjustment; Wiczling 2025 Bayesian D-optimal experiments (10.1016/j.chroma.2025.465787); CiRT peptide standards (10.1074/mcp.O114.042267); AL for IE (10.1021/acs.analchem.5c00816) | **Partly open.** Heuristic selection rules and an adjustment-oriented choice exist. No optimisation of calibrant sets against transfer error, window width or coverage was found |
| **B5** RT prediction intervals with coverage on many systems (2025–26) | RT/RI/order AND (prediction interval* OR uncertainty OR coverage OR credible OR conformal OR tolerance window) 2025–26; author search Kruve/Kretschmer/Wiczling/Böcker; Europe PMC RepoRT + interval terms | Kretschmer et al. 2026 Nat Methods (10.1038/s41592-026-03243-2; point predictions only); MC-Retention (#4, uncertainty on 8 own conditions); Wiczling (#17, one dataset); hybrid GNN with 23 forensic labs (10.1021/acs.jpcb.6c03944; no coverage); Malm & Kruve 2025 Analyst (10.1039/d5an00323g; accuracy only); TC-TopoRT (#6; calibrated thresholds on 2 sets) | **None found** with interval coverage across many systems. The manuscript's niche is still unoccupied |

## Implications

- **Open in analytical chemistry:**
  - **B2** (cross-system or hierarchical conformal calibration to get below 1/α−1 calibrants): no application to chemistry or LC was found. The statistics are ready (GHCP 2026, Dunn 2023, meta-XB), and this addresses the manuscript's own stated limitation. It is the strongest candidate.
  - **B3** (conformal selection/FDR, or guaranteed fusion of RT + CCS + MS2): nothing in metabolomics or NTS.
  - **B1** (JAW / weighted jackknife+ with structure-based density ratios between the calibrant set and the suspect list): open for RT, but CoDrug already exists in drug discovery. It maps directly onto the manuscript's documented failure mode: interlaboratory calibrant mixes gave 85–87 % coverage.
  - **B4**: only partly open. Heuristic and adjustment-oriented selection exists (Aalizadeh 2021, FastRet 2026). Optimal design against window width or coverage has not been done.
- **Evidence level the target journals expect:** accepted Anal Chem/ACA papers in this niche show an identification endpoint, usually on the authors' own measurements or real samples:
  - candidate-space reduction at a stated retention of true structures (TC-TopoRT: −30–41 % at 93–99 %);
  - an increase in identification probability or in the number of identifications (+47–55 %; +6 %);
  - removal of >70 % of impossible structures, with a drinking-water case;
  - fewer false positives in real matrices.

  Even the conformal CCS paper reports the fraction of wrong candidates removed (16 %, 8 % above random). A resubmission should therefore add:
  1. Real candidate lists (e.g. PubChemLite isomers per formula), reporting the fraction of wrong candidates removed and the rank of the true structure at guaranteed retention, against ±k·RMSE.
  2. At least one real-sample case study with the laboratory's own calibrants, or a NORMAN dataset with real suspect lists.
  3. A method that answers a failure documented in the manuscript, B1 or B2, rather than a plain jackknife+.
