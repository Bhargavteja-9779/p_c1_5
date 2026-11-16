export const meta = {
  name: 'aca-ten-reviewer-cycle',
  description: 'One review cycle: 10 independent specialist reviewers (ACA submission) + consolidation into a prioritised action list',
  phases: [
    { title: 'Review', detail: '10 reviewer agents with distinct perspectives' },
    { title: 'Consolidate', detail: 'deduplicate, prioritise, map to concrete fixes and experiments' },
    { title: 'Verify', detail: 'one adversarial skeptic per fatal/critical/major action tries to refute it against the files' },
  ],
}

const A = args || {}
const CYCLE = A.cycle || 1
const ROOT = '/home/user/p_c1_5/Target_Journal_Submission'
const PACKET = `REVIEW PACKET (read these first; open other repository files only when needed for your lens):
- Manuscript (Markdown, built from templates; every number comes from result files): ${ROOT}/manuscript/Final_Manuscript.md
- Supplementary Material (Sections S1-S6, Tables S1-S26; built from result files): ${ROOT}/supplementary/Supplementary_Material.md ; table index: ${ROOT}/supplementary/table_index.csv
- Supplementary proof of the coverage guarantee (also Section S2 of the Supplementary Material): ${ROOT}/supplementary/S2_coverage_proof.md
- Formatted outputs (for layout checks only): ${ROOT}/manuscript/Final_Manuscript.docx and .pdf ; ${ROOT}/supplementary/Supplementary_Material.pdf
- Manuscript section templates with {{token}} placeholders (what authors edit): ${ROOT}/manuscript/src/*.md ; token values: ${ROOT}/code/scripts/manuscript_facts.py
- Pre-registered plan and addenda: ${ROOT}/code/EXPERIMENT_PLAN.md
- Machine-generated tables: ${ROOT}/tables/*.md ; statistics: ${ROOT}/results/statistics.json ; metrics: ${ROOT}/results/metrics.json
- Code: ${ROOT}/code/src/conformrt/ (evaluation/system_eval.py, evaluation/intervals.py, evaluation/protocol.py, evaluation/driver.py, models/*.py), experiments/*.py, tests/test_core.py
- Journal selection and ACA benchmark: ${ROOT}/journal_selection/selected_journal.md, aca_journal_quality_benchmark.md
- Prior-art map and verified literature: ${ROOT}/supplementary/prior_art_novelty_report.md, literature_matrix.csv, citation_verification.csv, claim_source_verification.md
- Claim-evidence and citation audits, previous review cycles and responses (do NOT repeat issues that were already raised and fixed; DO check whether fixes are adequate): ${ROOT}/REVIEW/
Target journal: Analytica Chimica Acta (Elsevier; Q1, JIF 6.1). Review cycle ${CYCLE} of 10.`

const RULES = `Rules: be hostile but scientifically fair; find NEW, concrete, verifiable problems (cite file/section/table and quote the text or code). Do not invent problems to reach a number; if a concern is speculative, say so and lower its severity. Do not praise. Never put any person's name from the author team into web queries. Keep tool use focused (read the files you need; at most ~15 tool calls; web searches only if your lens requires them). Severity: fatal (invalidates the central claim), critical, major, moderate, minor. For each issue give a precise corrective action and, if needed, the exact additional experiment.`

const SCHEMA = {
  type: 'object',
  properties: {
    reviewer: { type: 'string' },
    recommendation: { type: 'string', enum: ['acceptable', 'minor revision', 'major revision', 'reject'] },
    summary: { type: 'string' },
    issues: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          id: { type: 'string' }, severity: { type: 'string', enum: ['fatal', 'critical', 'major', 'moderate', 'minor'] },
          title: { type: 'string' }, evidence: { type: 'string' }, location: { type: 'string' },
          corrective_action: { type: 'string' }, required_experiment: { type: 'string' }, confidence: { type: 'string' },
        },
        required: ['id', 'severity', 'title', 'evidence', 'location', 'corrective_action', 'required_experiment', 'confidence'],
      },
    },
    verified_previous_fixes: { type: 'array', items: { type: 'string' } },
  },
  required: ['reviewer', 'recommendation', 'summary', 'issues', 'verified_previous_fixes'],
}

const REVIEWERS = [
  { key: 'R1_novelty', lens: 'AGENT 1 - NOVELTY AND PRIOR-ART REVIEWER. Try to invalidate the novelty claims: search for similar publications, methods, patents (PubMed E-utilities and Crossref via curl; WebSearch via ToolSearch "select:WebSearch") for conformal / distribution-free / coverage-guaranteed retention-time or retention-index windows, few-calibrant RT transfer with intervals, meta-learned normalisation in chemistry. Check that the manuscript states honestly what is not new.' },
  { key: 'R2_methodology', lens: 'AGENT 2 - METHODOLOGY REVIEWER. Check mathematical correctness (jackknife+ quantile indices, normalised jackknife+ proof in S2, monotone map, finite-window condition n >= 1/alpha - 1, interval score), assumptions (exchangeability, symmetry of the LOO procedure), and whether the implementation matches the equations (read evaluation/intervals.py, evaluation/system_eval.py, models/difficulty.py, models/mapping.py).' },
  { key: 'R3_experimental_design', lens: 'AGENT 3 - EXPERIMENTAL DESIGN REVIEWER. Check fairness of baselines (same base score/calibrants/maps?), hyper-parameter selection (pre-registration vs post hoc), evaluation protocol (leave-source-out folds, draws, n grid), whether the fixed-window and NORMAN analyses are properly labelled as post hoc, and reproducibility.' },
  { key: 'R4_statistics', lens: 'AGENT 4 - STATISTICAL REVIEWER. Check confidence intervals, unit of analysis (systems within sources are correlated), Wilcoxon/Holm usage, effect sizes, multiple comparisons across n and alpha, claims of significance vs equivalence, coverage validity claims (mean coverage vs per-system), and whether statistical language in the text matches statistics.json.' },
  { key: 'R5_dataset_leakage', lens: 'AGENT 5 - DATASET AND LEAKAGE REVIEWER. Check RepoRT cleaning (InChIKey connectivity, stereoisomer handling), duplicate or near-duplicate systems across folds, compound overlap between training and test systems (seen/unseen), the inner cross-fitting of the difficulty model, NORMAN overlap with RepoRT, licensing (CC BY-SA), and whether any target information leaks into models.' },
  { key: 'R6_robustness_generalisation', lens: 'AGENT 6 - ROBUSTNESS AND GENERALISATION REVIEWER. Try to show the method fails: non-random calibrant mixes, calibrant errors, very small n, HILIC, different base model, novel compounds, ionisable compounds, extrapolation outside the calibrant range, the external NORMAN data. Judge whether the claims are bounded correctly by the evidence.' },
  { key: 'R7_computation', lens: 'AGENT 7 - COMPUTATIONAL EFFICIENCY AND DEPLOYMENT REVIEWER. Check reported cost (training, calibration, per-candidate), memory, scalability to large candidate lists, practicality for laboratories (software, inputs needed), and whether the improvement is worth the added complexity compared with plain jackknife+.' },
  { key: 'R8_literature_citations', lens: 'AGENT 8 - LITERATURE AND CITATION REVIEWER. Verify a sample of references and DOIs (Crossref API via curl), check that statements attributed to cited papers are supported (use citation_verification.csv and claim_source_verification.md), identify missing important competing studies, and check the bibliography format for Elsevier.' },
  { key: 'R9_associate_editor', lens: 'AGENT 9 - ANALYTICA CHIMICA ACTA ASSOCIATE EDITOR. Evaluate scope fit, novelty/impact expectations of ACA (see the ACA benchmark), structured abstract (Background/Results/Significance), highlights, keywords, length, figure/table quality, declarations (including generative-AI declaration), data availability, and desk-rejection risks. Decide whether the paper would be sent to review.' },
  { key: 'R10_hostile_competitor', lens: 'AGENT 10 - HOSTILE COMPETING RESEARCHER. Argue that the work is incremental or unnecessary: is there a simpler existing method (plain jackknife+, ±2·RMSE with LOO, larger fixed windows, GP projection) that gives equivalent practical results? Is any comparison unfair? Does the dataset favour the proposed method? Use the numbers in the tables to make your case.' },
]

phase('Review')
const reviews = await parallel(REVIEWERS.map(r => () => agent(`${r.lens}\n\n${PACKET}\n\n${RULES}\nReturn your review via the structured output (reviewer = "${r.key}").`, { label: `c${CYCLE}:${r.key}`, phase: 'Review', schema: SCHEMA })))
const ok = reviews.filter(Boolean)
log(`cycle ${CYCLE}: ${ok.length}/10 reviews returned`)

phase('Consolidate')
const CONS = {
  type: 'object',
  properties: {
    actions: { type: 'array', items: { type: 'object', properties: {
      action_id: { type: 'string' }, priority: { type: 'string', enum: ['fatal', 'critical', 'major', 'moderate', 'minor'] },
      problem: { type: 'string' }, raised_by: { type: 'array', items: { type: 'string' } }, fix: { type: 'string' },
      experiment_needed: { type: 'string' }, scientific_hypothesis: { type: 'string' }, feasible_now: { type: 'boolean' } },
      required: ['action_id', 'priority', 'problem', 'raised_by', 'fix', 'experiment_needed', 'scientific_hypothesis', 'feasible_now'] } },
    recommendations: { type: 'object', additionalProperties: { type: 'string' } },
    overall: { type: 'string' },
  },
  required: ['actions', 'recommendations', 'overall'],
}
const consolidated = await agent(`You consolidate review cycle ${CYCLE} for an Analytica Chimica Acta submission. Below are ${ok.length} independent reviews (JSON). Merge duplicate issues (keep all reviewers in raised_by), drop issues that are clearly wrong given the manuscript files (check ${ROOT}/manuscript/Final_Manuscript.md when unsure), prioritise (fatal > critical > major > moderate > minor), and for every remaining issue write a concrete fix, whether an experiment is needed, and the scientific hypothesis behind a methodological change. Mark feasible_now=false only if the fix needs resources unavailable here (wet-lab data, private data).\n\nREVIEWS:\n${JSON.stringify(ok)}`, { label: `c${CYCLE}:consolidate`, phase: 'Consolidate', schema: CONS })
phase('Verify')
const VERDICT = {
  type: 'object',
  properties: {
    action_id: { type: 'string' }, verdict: { type: 'string', enum: ['real', 'partly', 'refuted'] },
    evidence: { type: 'string' }, minimal_correct_fix: { type: 'string' },
  },
  required: ['action_id', 'verdict', 'evidence', 'minimal_correct_fix'],
}
const heavy = ((consolidated && consolidated.actions) || []).filter(a => ['fatal', 'critical', 'major'].includes(a.priority))
log(`cycle ${CYCLE}: verifying ${heavy.length} fatal/critical/major actions adversarially`)
const verdicts = (await parallel(heavy.map(a => () => agent(`You are an adversarial VERIFIER for review cycle ${CYCLE} of an Analytica Chimica Acta submission. A reviewer panel raised the issue below. Try to REFUTE it: open the cited files (manuscript ${ROOT}/manuscript/Final_Manuscript.md, templates ${ROOT}/manuscript/src/*.md, Supplementary ${ROOT}/supplementary/Supplementary_Material.md, results ${ROOT}/results/*.json and ${ROOT}/tables/*.csv, code ${ROOT}/code/) and check whether the problem really exists as described and whether the proposed fix is scientifically right. Verdict 'real' = the problem exists and matters; 'partly' = exists but is overstated, or the proposed fix is wrong/excessive; 'refuted' = factually wrong given the files. Quote file evidence. Give the minimal correct fix (wording, analysis or experiment).\n\nISSUE (JSON):\n${JSON.stringify(a)}\n\nKeep tool use focused (at most ~12 tool calls).`, { label: `c${CYCLE}:verify:${a.action_id}`, phase: 'Verify', schema: VERDICT })))).filter(Boolean)

// full reviews are kept in the workflow journal (journal.jsonl) and rendered by journal_to_md.py; return a compact summary
return { cycle: CYCLE, n_reviews: ok.length, recommendations: ok.map(r => ({ reviewer: r.reviewer, recommendation: r.recommendation, n_issues: (r.issues || []).length })), consolidated, verdicts }
