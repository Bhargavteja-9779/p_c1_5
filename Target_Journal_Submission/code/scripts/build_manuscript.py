"""Assemble the manuscript from section templates, resolving every {{token}} from machine-generated result files.
Unresolved tokens abort the build (no hand-typed numbers). Produces manuscript/Final_Manuscript.{md,docx,pdf}.
Usage: python scripts/build_manuscript.py [--sections-dir ../manuscript/src] [--no-pdf]"""
import argparse, os, sys, glob

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
from docbuild import resolve, pandoc_docx, postprocess_docx, docx_to_pdf  # noqa: E402
from manuscript_facts import build_facts  # noqa: E402  project-specific: derives every reported number from result files

ap = argparse.ArgumentParser(); ap.add_argument("--sections-dir", default=os.path.join(ROOT, "manuscript", "src"))
ap.add_argument("--out-dir", default=os.path.join(ROOT, "manuscript")); ap.add_argument("--no-pdf", action="store_true")
ap.add_argument("--name", default="Final_Manuscript"); a = ap.parse_args()
RES, TAB, FIG = (os.path.join(ROOT, d) for d in ("results", "tables", "figures"))

# ---------------------------------------------------------------- assemble
facts, tables = build_facts(RES, TAB)
parts = sorted(glob.glob(os.path.join(a.sections_dir, "*.md")))
body = "\n\n".join(open(p).read() for p in parts)
for tname, tmd in tables.items():
    body = body.replace("{{table:" + tname + "}}", tmd)
body, missing = resolve(body, facts)
if missing:
    sys.exit("UNRESOLVED TOKENS: " + ", ".join(sorted(set(missing))))
body = body.replace("](figures/", "](" + FIG + "/")
md_path = os.path.join(a.out_dir, a.name + ".md"); open(md_path, "w").write(body)
print("markdown written:", md_path)

# ---------------------------------------------------------------- pandoc -> docx -> pdf
build = os.path.join(a.out_dir, "build")
docx = os.path.join(a.out_dir, a.name + ".docx")
pandoc_docx(md_path, docx, os.path.join(build, "reference_aca.docx"), ROOT, bib=os.path.join(a.out_dir, "references.bib"),
            csl=os.path.join(build, "elsevier-with-titles.csl"), number_sections=True)
postprocess_docx(docx, line_numbers=True, table_pt=9)
print("docx written:", docx)
if not a.no_pdf:
    ok, err = docx_to_pdf(docx, a.out_dir)
    print("pdf:", "ok" if ok else err)
