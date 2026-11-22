"""Assemble the Analytica Chimica Acta upload package from the built documents.
Output: submission_package/ and ACA_submission_package.zip (both in Target_Journal_Submission/).
Usage: python scripts/make_submission_package.py"""
import os, re, shutil, zipfile
from PIL import Image
from docx import Document
from docx.shared import Pt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MS, SUP, FIG = (os.path.join(ROOT, d) for d in ("manuscript", "supplementary", "figures"))
OUT = os.path.join(ROOT, "submission_package")
TITLE = ("Laboratory-level guarantees for retention-time filtering in non-target screening "
         "by cross-system conformal calibration")

shutil.rmtree(OUT, ignore_errors=True); os.makedirs(OUT)


def cp(src, name):
    shutil.copy2(src, os.path.join(OUT, name))


def doc(name, heading, paras):
    d = Document(); st = d.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(12)
    d.add_heading(heading, level=1)
    for p in paras:
        if p.startswith("- "):
            d.add_paragraph(p[2:], style="List Bullet")
        else:
            d.add_paragraph(p)
    d.core_properties.author = ""; d.core_properties.title = heading; d.core_properties.comments = ""
    d.save(os.path.join(OUT, name))


# 1-2 cover letter and manuscript
cp(os.path.join(MS, "Cover_Letter.docx"), "01_Cover_Letter.docx")
cp(os.path.join(MS, "Cover_Letter.pdf"), "01_Cover_Letter.pdf")
cp(os.path.join(MS, "Final_Manuscript.docx"), "02_Manuscript.docx")
cp(os.path.join(MS, "Final_Manuscript.pdf"), "02_Manuscript.pdf")

# 3 highlights (taken from the built manuscript, so they always match)
md = open(os.path.join(MS, "Final_Manuscript.md")).read()
hl = re.search(r"# Highlights.*?\n(.*?)\n# Abstract", md, re.S).group(1)
hl = [l[2:].strip() for l in hl.splitlines() if l.startswith("- ")]
assert 3 <= len(hl) <= 5 and all(len(h) <= 85 for h in hl), hl
doc("03_Highlights.docx", "Highlights", ["- " + h for h in hl])

# 4 graphical abstract
cp(os.path.join(FIG, "graphical_abstract.tiff"), "04_Graphical_Abstract.tiff")
cp(os.path.join(FIG, "graphical_abstract.pdf"), "04_Graphical_Abstract.pdf")

# 5 figures (vector PDF + 300 dpi TIFF, RGB)
for i, stem in enumerate(["fig1_workflow", "figX1_lab_level", "figX2_selection_rules"], 1):
    cp(os.path.join(FIG, stem + ".pdf"), f"05_Fig{i}.pdf")
    im = Image.open(os.path.join(FIG, stem + ".png")).convert("RGB")
    im.save(os.path.join(OUT, f"05_Fig{i}.tiff"), dpi=(300, 300), compression="tiff_lzw")

# 6 supplementary material
cp(os.path.join(SUP, "Supplementary_Material.pdf"), "06_Supplementary_Material.pdf")
cp(os.path.join(SUP, "Supplementary_Material.docx"), "06_Supplementary_Material.docx")

# 7 declaration of interests (Elsevier standard wording)
doc("07_Declaration_of_Interests.docx", "Declaration of interests", [
    "Manuscript title: " + TITLE,
    "☒ The authors declare that they have no known competing financial interests or personal relationships that could "
    "have appeared to influence the work reported in this paper.",
    "☐ The authors declare the following financial interests/personal relationships which may be considered as "
    "potential competing interests: (none)",
    "Funding: This research did not receive any specific grant from funding agencies in the public, commercial, or "
    "not-for-profit sectors.",
    "Corresponding author: P N Bhargav Teja, [e-mail address]"])

# 8 suggested reviewers (independent experts in the field; e-mail addresses must be added by the authors)
doc("08_Suggested_Reviewers.docx", "Suggested reviewers", [
    "The following researchers are experts in retention-time prediction, candidate filtering and non-target screening. "
    "None of them has collaborated with the authors. E-mail addresses must be taken from their most recent publications "
    "or institutional pages before submission, and current affiliations re-checked.",
    "- Prof. Anneli Kruve (Stockholm University, Sweden) — retention-time projection across laboratories, NORMAN "
    "interlaboratory comparison, uncertainty-aware candidate prioritisation in NTS. E-mail: [to be added]",
    "- Dr. Emma L. Schymanski (Luxembourg Centre for Systems Biomedicine, University of Luxembourg) — identification "
    "confidence in NTS, PubChemLite. E-mail: [to be added]",
    "- Prof. Michael Witting (Helmholtz Munich, Germany) — RepoRT, retention-time data for metabolomics. "
    "E-mail: [to be added]",
    "- Dr. Jan Stanstrup (University of Copenhagen, Denmark) — PredRet, retention-time projection between systems. "
    "E-mail: [to be added]",
    "- Dr. Eric Bach (Aalto University, Finland) — retention-order prediction for metabolite identification. "
    "E-mail: [to be added]",
    "Note: Prof. Kruve and Dr. Malm authored the NORMAN data set used here; the editor may prefer to treat them as "
    "data providers rather than reviewers. Do not suggest anyone with whom an author has a personal or professional "
    "relationship."])

# 9 checklist (not uploaded)
open(os.path.join(OUT, "00_SUBMISSION_CHECKLIST.txt"), "w").write(f"""ANALYTICA CHIMICA ACTA - SUBMISSION CHECKLIST (Elsevier Editorial Manager)
Article type: Research Paper
Title: {TITLE}

Upload as (Editorial Manager item type -> file):
  Cover Letter                     -> 01_Cover_Letter.docx
  Manuscript (editable, with tables, figures embedded) -> 02_Manuscript.docx   (02_Manuscript.pdf for checking only)
  Highlights                       -> 03_Highlights.docx
  Graphical Abstract               -> 04_Graphical_Abstract.tiff
  Figure                           -> 05_Fig1.tiff, 05_Fig2.tiff, 05_Fig3.tiff (vector PDFs also provided)
  Supplementary Material           -> 06_Supplementary_Material.pdf
  Declaration of Interest          -> 07_Declaration_of_Interests.docx (or regenerate with Elsevier's declarations tool)
  Suggested reviewers              -> enter in the submission form (08_Suggested_Reviewers.docx)

BEFORE UPLOAD - fill in the placeholders (search for "[" in each file):
  [ ] Second author name, affiliations a/b, corresponding e-mail - 02_Manuscript.docx title page + CRediT
  [ ] [Date] and signature line - 01_Cover_Letter.docx
  [ ] Corresponding e-mail - 07_Declaration_of_Interests.docx
  [ ] Reviewer e-mail addresses - 08_Suggested_Reviewers.docx
  [ ] Every author approves the final version and the CRediT roles
  [ ] Keep the "Declaration of generative AI and AI-assisted technologies" section unchanged (mandatory Elsevier policy)
  [ ] Online form: data availability = public data (RepoRT, NORMAN ESI, PubChemLite); code = on request
""")

zp = os.path.join(ROOT, "ACA_submission_package.zip")
with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
    for f in sorted(os.listdir(OUT)):
        z.write(os.path.join(OUT, f), os.path.join("ACA_submission_package", f))
print("package:", OUT); print("zip:", zp, os.path.getsize(zp) // 1024, "kB")
