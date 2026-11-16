"""Shared helpers for building the manuscript and the Supplementary Material:
token resolution ({{key|fmt}} -> value from machine-generated result files), pandoc conversion to DOCX,
DOCX post-processing (line numbers, three-line tables, figure width) and LibreOffice PDF export."""
import os, re, subprocess

TYPO_MINUS = re.compile(r"(?<![\w.])-(?=\d)")  # hyphen-minus before a digit that is not part of a word/range


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


def fmt_text(v, spec):
    """fmt() with a typographic minus sign (U+2212) for negative values."""
    out = fmt(v, spec)
    return "\u2212" + out[1:] if out.startswith("-") and len(out) > 1 and (out[1].isdigit() or out[1] == ".") else out


def resolve(text, facts):
    """Replace every {{key|fmt}} by the formatted fact; return (text, list of unresolved keys)."""
    missing = []

    def rep(m):
        key, _, spec = m.group(1).partition("|")
        key = key.strip()
        if key not in facts:
            missing.append(key); return m.group(0)
        return fmt_text(facts[key], spec.strip() or None)
    return re.sub(r"\{\{([^{}]+)\}\}", rep, text), missing


def mdtab(header, rows, max_w=36):
    """Pipe table with a typographic minus sign for negative numbers. The dash counts of the separator line are made proportional
    to the content length of each column, which pandoc uses as relative column widths when the table wraps."""
    cell = lambda v: TYPO_MINUS.sub("\u2212", str(v))
    H = [cell(h) for h in header]; R = [[cell(x) for x in r] for r in rows]
    w = []
    for j, h in enumerate(H):
        content = max([len(r[j]) for r in R] + [1])
        longest_word = max([len(t) for t in h.split()] + [1])
        w.append(max(5, min(max_w, max(content, longest_word + 2, len(h) // 2))))
    return "\n".join(["| " + " | ".join(H) + " |", "|" + "|".join("-" * k for k in w) + "|"] + ["| " + " | ".join(r) + " |" for r in R])


def pandoc_docx(md_path, docx, reference_doc, resource_path, bib=None, csl=None, number_sections=False):
    cmd = ["pandoc", md_path, "-o", docx, "--from", "markdown+tex_math_dollars+pipe_tables+implicit_figures+table_captions",
           "--reference-doc", reference_doc, "--resource-path", resource_path]
    if bib:
        cmd += ["--citeproc", "--bibliography", bib, "--csl", csl, "--metadata", "link-citations=true"]
    if number_sections:
        cmd += ["--number-sections"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("pandoc failed: " + r.stderr)
    if r.stderr.strip():
        print("pandoc warnings:\n" + r.stderr)


def postprocess_docx(docx, line_numbers=True, table_pt=9, max_fig_cm=16.0):
    from docx import Document
    from docx.shared import Pt, Cm
    from docx.oxml.ns import qn
    from docx.oxml import OxmlElement
    d = Document(docx)
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    for sec in d.sections:  # continuous line numbering for review, centred page number in the footer
        sp = sec._sectPr
        footer = sec.footer; footer.is_linked_to_previous = False
        fp = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()
        if "PAGE" not in footer._element.xml:
            fp.alignment = WD_ALIGN_PARAGRAPH.CENTER; run = fp.add_run()
            fp._p.get_or_add_pPr().append(OxmlElement("w:suppressLineNumbers"))  # no line number on the footer line
            for kind, text in [("begin", None), (None, "PAGE"), ("end", None)]:
                if kind:
                    el = OxmlElement("w:fldChar"); el.set(qn("w:fldCharType"), kind)
                else:
                    el = OxmlElement("w:instrText"); el.set(qn("xml:space"), "preserve"); el.text = text
                run._r.append(el)
            run.font.size = Pt(10); run.font.name = "Times New Roman"
        if line_numbers:
            ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); ln.set(qn("w:distance"), "283")
            sp.append(ln)
        pg = OxmlElement("w:pgNumType"); pg.set(qn("w:start"), "1"); sp.append(pg)
    for t in d.tables:  # three-line tables, single spaced
        t.autofit = True
        tblPr = t._tbl.tblPr
        borders = OxmlElement("w:tblBorders")
        for edge, val in [("top", "single"), ("bottom", "single"), ("insideH", "nil"), ("insideV", "nil"), ("left", "nil"), ("right", "nil")]:
            e = OxmlElement(f"w:{edge}"); e.set(qn("w:val"), val); e.set(qn("w:sz"), "6"); e.set(qn("w:color"), "000000"); borders.append(e)
        tblPr.append(borders)
        mar = OxmlElement("w:tblCellMar")  # narrow left/right cell margins (0.1 cm) so that narrow numeric columns do not wrap
        for side in ("left", "right"):
            e = OxmlElement(f"w:{side}"); e.set(qn("w:w"), "57"); e.set(qn("w:type"), "dxa"); mar.append(e)
        tblPr.append(mar)
        for ri, row in enumerate(t.rows):
            for c in row.cells:
                for p in c.paragraphs:
                    p.paragraph_format.line_spacing = 1.0; p.paragraph_format.space_after = Pt(0)
                    for run in p.runs:
                        run.font.size = Pt(table_pt); run.font.name = "Times New Roman"
                        if ri == 0:
                            run.font.bold = True
            if ri == 0:
                for c in row.cells:
                    tcPr = c._tc.get_or_add_tcPr(); b = OxmlElement("w:tcBorders"); bt = OxmlElement("w:bottom")
                    bt.set(qn("w:val"), "single"); bt.set(qn("w:sz"), "4"); bt.set(qn("w:color"), "000000"); b.append(bt); tcPr.append(b)
    for shp in d.inline_shapes:  # cap figure width to the text width
        if shp.width > Cm(max_fig_cm):
            ratio = Cm(max_fig_cm) / shp.width; shp.width = Cm(max_fig_cm); shp.height = int(shp.height * ratio)
    d.save(docx)


def docx_to_pdf(docx, out_dir):
    r = subprocess.run(["soffice", "--headless", "--convert-to", "pdf", "--outdir", out_dir, docx], capture_output=True, text=True, timeout=900)
    return r.returncode == 0, r.stderr[-500:]
