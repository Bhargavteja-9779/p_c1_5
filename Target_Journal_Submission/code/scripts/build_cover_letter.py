"""Build the cover letter from manuscript/cover_letter/cover_letter_template.md; every number is a token resolved from the result
files (unresolved tokens abort). Outputs manuscript/Cover_Letter.{md,docx,pdf}. Usage: python scripts/build_cover_letter.py [--no-pdf]"""
import argparse, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
from docbuild import resolve, pandoc_docx, postprocess_docx, docx_to_pdf  # noqa: E402
from manuscript_facts import build_facts  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--no-pdf", action="store_true"); a = ap.parse_args()
M = os.path.join(ROOT, "manuscript")
facts, _ = build_facts(os.path.join(ROOT, "results"), os.path.join(ROOT, "tables"))
text, missing = resolve(open(os.path.join(M, "cover_letter", "cover_letter_template.md")).read(), facts)
if missing:
    sys.exit("UNRESOLVED TOKENS: " + ", ".join(sorted(set(missing))))
md = os.path.join(M, "Cover_Letter.md"); open(md, "w").write(text)
docx = os.path.join(M, "Cover_Letter.docx")
pandoc_docx(md, docx, os.path.join(M, "build", "reference_aca.docx"), ROOT)
postprocess_docx(docx, line_numbers=False)
print("cover letter written:", docx)
if not a.no_pdf:
    ok, err = docx_to_pdf(docx, M); print("pdf:", "ok" if ok else err)
