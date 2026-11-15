"""Assemble the manuscript from section templates, resolving every {{token}} from machine-generated result files.
Unresolved tokens abort the build (no hand-typed numbers). Produces manuscript/Final_Manuscript.{md,docx,pdf}.
Usage: python scripts/build_manuscript.py [--sections-dir ../manuscript/src] [--no-pdf]"""
import argparse, os, re, json, subprocess, sys, glob
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
ap = argparse.ArgumentParser(); ap.add_argument("--sections-dir", default=os.path.join(ROOT, "manuscript", "src"))
ap.add_argument("--out-dir", default=os.path.join(ROOT, "manuscript")); ap.add_argument("--no-pdf", action="store_true")
ap.add_argument("--name", default="Final_Manuscript"); a = ap.parse_args()
RES, TAB, FIG = (os.path.join(ROOT, d) for d in ("results", "tables", "figures"))
sys.path.insert(0, os.path.join(HERE, "facts")) if os.path.isdir(os.path.join(HERE, "facts")) else None

# ---------------------------------------------------------------- facts
def load_facts():
    sys.path.insert(0, HERE)
    from manuscript_facts import build_facts  # project-specific: derives every reported number from result files
    return build_facts(RES, TAB)

def fmt(v, spec):
    if spec is None or spec == "":
        return str(v)
    if spec == "pct":
        return f"{100 * float(v):.1f} %"
    if spec == "pct0":
        return f"{100 * float(v):.0f} %"
    if spec == "int":
        return f"{int(round(float(v))):,}"
    if spec.startswith("p"):  # p-value formatting
        x = float(v); return "< 0.001" if x < 0.001 else f"{x:.{int(spec[1:] or 3)}f}"
    return format(float(v), spec)

def resolve(text, facts):
    missing = []
    def rep(m):
        key, _, spec = m.group(1).partition("|")
        key = key.strip()
        if key not in facts:
            missing.append(key); return m.group(0)
        return fmt(facts[key], spec.strip() or None)
    out = re.sub(r"\{\{([^{}]+)\}\}", rep, text)
    return out, missing

# ---------------------------------------------------------------- assemble
facts, tables = load_facts()
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

# ---------------------------------------------------------------- pandoc -> docx
build = os.path.join(a.out_dir, "build")
docx = os.path.join(a.out_dir, a.name + ".docx")
cmd = ["pandoc", md_path, "-o", docx, "--from", "markdown+tex_math_dollars+pipe_tables+implicit_figures+table_captions",
       "--citeproc", "--bibliography", os.path.join(a.out_dir, "references.bib"), "--csl", os.path.join(build, "elsevier-with-titles.csl"),
       "--reference-doc", os.path.join(build, "reference_aca.docx"), "--resource-path", ROOT, "--metadata", "link-citations=true"]
r = subprocess.run(cmd, capture_output=True, text=True)
if r.returncode != 0:
    sys.exit("pandoc failed: " + r.stderr)
if r.stderr.strip():
    print("pandoc warnings:\n" + r.stderr)

# ---------------------------------------------------------------- docx post-processing
from docx import Document
from docx.shared import Pt, Cm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
d = Document(docx)
for sec in d.sections:  # continuous line numbering for review, A4
    sp = sec._sectPr; ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); ln.set(qn("w:distance"), "283")
    sp.append(ln)
    pg = OxmlElement("w:pgNumType"); pg.set(qn("w:start"), "1"); sp.append(pg)
for t in d.tables:  # three-line style tables, 9 pt, single spaced
    t.autofit = True
    tbl = t._tbl; tblPr = tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge, val in [("top", "single"), ("bottom", "single"), ("insideH", "nil"), ("insideV", "nil"), ("left", "nil"), ("right", "nil")]:
        e = OxmlElement(f"w:{edge}"); e.set(qn("w:val"), val); e.set(qn("w:sz"), "6"); e.set(qn("w:color"), "000000"); borders.append(e)
    tblPr.append(borders)
    for ri, row in enumerate(t.rows):
        for c in row.cells:
            for p in c.paragraphs:
                p.paragraph_format.line_spacing = 1.0; p.paragraph_format.space_after = Pt(0)
                for run in p.runs:
                    run.font.size = Pt(9); run.font.name = "Times New Roman"
                    if ri == 0: run.font.bold = True
        if ri == 0:
            for c in row.cells:
                tcPr = c._tc.get_or_add_tcPr(); b = OxmlElement("w:tcBorders"); bt = OxmlElement("w:bottom")
                bt.set(qn("w:val"), "single"); bt.set(qn("w:sz"), "4"); bt.set(qn("w:color"), "000000"); b.append(bt); tcPr.append(b)
for shp in d.inline_shapes:  # cap figure width to text width (16 cm)
    if shp.width > Cm(16):
        ratio = Cm(16) / shp.width; shp.width = Cm(16); shp.height = int(shp.height * ratio)
d.save(docx); print("docx written:", docx)

if not a.no_pdf:
    r = subprocess.run(["soffice", "--headless", "--convert-to", "pdf", "--outdir", a.out_dir, docx], capture_output=True, text=True, timeout=600)
    print("pdf:", "ok" if r.returncode == 0 else r.stderr[-500:])
