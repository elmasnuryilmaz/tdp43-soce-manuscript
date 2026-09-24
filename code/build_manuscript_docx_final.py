#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the submission .docx, with Tables 1-5, from the manuscript Markdown.

1. Writes Tables 1-5, formatted for publication, from tables/*.csv into
   manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL.md between the markers <!-- table:N --> and
   <!-- /table:N -->. The captions above the markers are part of the manuscript and are left
   alone; the CSV files remain the machine-readable versions of the tables.
2. Converts the Markdown with pandoc (default reference document).
3. Post-processes word/document.xml: A4 pages with 2.54 cm margins; the Tables section starts
   on a new page and every table sits in its own section, Tables 2-5 in landscape; table text
   and table notes are 9 pt; tables get top and bottom rules; column widths are written in
   DXA on the grid and on every cell.
4. Submission layout: text double-spaced (tables, their captions and notes single-spaced),
   continuous line numbers in the main-text section, page numbers in the footer.

Usage:  /usr/bin/python3 code/build_manuscript_docx.py [--check]
        --check only reports whether the tables in the Markdown are up to date.
"""
import io, math, os, re, subprocess, sys, zipfile
import pandas as pd

P = "/Users/elmas/Desktop/MAKALE/11_NEUROCHEMISTRY_INTERNATIONAL_FIGURE_REVISION"
MD = f"{P}/manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.md"
DOCX = f"{P}/manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx"
TAB = f"{P}/tables"

MINUS = "−"
SUPER = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")

# A4 in twentieths of a point, 1-inch margins
PAGE_W, PAGE_H, MARGIN = 11906, 16838, 1440
TEXT_W = {"portrait": PAGE_W - 2 * MARGIN, "landscape": PAGE_H - 2 * MARGIN}
ORIENT = {1: "portrait", 2: "landscape", 3: "landscape", 4: "landscape", 5: "landscape"}


# ------------------------------------------------------------------ formatting helpers
def signed(x, d=3):
    return "—" if pd.isna(x) else f"{float(x):+.{d}f}".replace("-", MINUS)


def pval(x):
    """p and q values: two significant digits, scientific notation below 0.01."""
    if pd.isna(x):
        return "—"
    x = float(x)
    if x == 0:
        return "< 1 × 10⁻¹⁶"
    if x >= 0.01:
        return f"{x:.{max(1, 1 - math.floor(math.log10(x)))}f}"
    m, e = f"{x:.1e}".split("e")
    return f"{m} × 10{str(int(e)).translate(SUPER)}"


def integer(x):
    return "—" if pd.isna(x) else f"{int(round(float(x))):,}"


def gene(g):
    return f"*{g}*"


def md_table(header, rows, align, widths):
    """Pipe table; widths are dash counts, which pandoc turns into relative column widths
    because the separator line is longer than its 72-column limit."""
    esc = lambda c: str(c).replace("|", "\\|")
    line = lambda cells: "| " + " | ".join(esc(c) for c in cells) + " |"
    sep = "|" + "|".join((":" if a == "l" else "") + "-" * w + (":" if a == "r" else "")
                         for a, w in zip(align, widths)) + "|"
    assert len(sep) > 72
    return "\n".join([line(header), sep] + [line(r) for r in rows])


# ------------------------------------------------------------------ the five tables
def table1():
    t = pd.read_csv(f"{TAB}/Table2_coverage_prefilter.csv").set_index("dataset")
    order = ["SH-SY5Y (GSE296712)", "iPSC colonies (GSE230647)",
             "iPSC-derived motor neurons (GSE77702)", "Mouse striatum (GSE27394)",
             "C2C12 (GSE171714)", "NSC34 (GSE171714)"]
    rows = [[d, integer(r.events_tested), integer(r.events_after_filter),
             f"{r.events_removed_pct:.1f}", integer(r.significant_before_filter),
             integer(r.significant_after_filter), f"{r.significant_lost_pct:.1f}"]
            for d, r in t.loc[order].iterrows()]
    head = ["Dataset", "Events tested", "Events after filter", "Removed (%)",
            "Significant before filter", "Significant after filter", "Significant calls lost (%)"]
    note = ("*Note.* Significant: FDR < 0.05 and |ΔPSI| ≥ 0.10. After the filter, "
            "Benjamini–Hochberg q values were recomputed within the retained events (Methods 2.3). "
            "FDR, false discovery rate; PSI, percent spliced in.")
    return md_table(head, rows, "lrrrrrr", [26, 9, 9, 8, 9, 9, 10]), note


def table2():
    t = pd.read_csv(f"{TAB}/Table3_robust_SOCE_splicing_events.csv").set_index("gene")
    rows = []
    for g in ["STIMATE", "ORAI3", "STIM2", "STIM1"]:
        r = t.loc[g]
        assert r.event_class == "SE", "Table 2 legend states that all events are skipped exons"
        rows.append([gene(g),
                     f"{r.chrom}:{int(r.start_1based):,}–{int(r.end):,} ({str(r.strand).replace('-', MINUS)})",
                     integer(r.exon_bp), str(r.reading_frame).capitalize(),
                     signed(r.delta_PSI), pval(r.FDR),
                     f"{signed(r.CI95_low)} to {signed(r.CI95_high)}",
                     str(r.PSI_knockdown).replace(";", "; "), str(r.PSI_control).replace(";", "; "),
                     f"{r.mean_reads_per_sample:.1f}", integer(r.min_informative_reads)])
    head = ["Gene", "Exon, GRCh38 (strand)", "Length (bp)", "Reading frame", "ΔPSI",
            "FDR", "Bootstrap 95% CI", "PSI, knockdown replicates", "PSI, control replicates",
            "Mean reads per sample", "Minimum reads in a sample"]
    note = ("*Note.* SH-SY5Y, 0 versus 75 ng/mL doxycycline, three libraries per group. ΔPSI is "
            "knockdown minus control from rMATS junction counts; the confidence interval is a "
            "replicate-level bootstrap (10,000 resamples). Coordinates are 1-based and inclusive. "
            "CI, confidence interval; FDR, false discovery rate; PSI, percent spliced in.")
    return md_table(head, rows, "lllrrrrlrrr", [14, 43, 12, 17, 11, 17, 17, 19, 19, 12, 13]), note


def table3():
    t = pd.read_csv(f"{TAB}/Table4_cryptic_events_eleven_comparisons.csv").set_index("comparison")
    order = ["SH-SY5Y 75 ng/mL", "SH-SY5Y 25 ng/mL", "iPSC colonies", "iPSC-MN, TDP-43 KD",
             "iPSC-MN, FUS KD", "iPSC-MN, TAF15 KD", "K562 total RNA", "K562 poly(A)+ mRNA",
             "C2C12", "NSC34", "Mouse striatum"]

    def genes(v):
        v = str(v)
        if v in ("—", "nan"):
            return "—"
        if v.startswith("n/a"):
            return "n/a"
        return ", ".join(gene(x.strip()) for x in v.split(","))

    def count(v):
        return "n/a" if str(v).startswith("n/a") else integer(v)

    rows = []
    for c in order:
        r = t.loc[c]
        null = "—" if pd.isna(r.null_calls_high_confidence) else \
            f"{integer(r.null_calls_high_confidence)} ({r.null_to_real_ratio:.2f})"
        rows.append([c, f"{integer(r.permissive_events)} ({integer(r.permissive_genes)})",
                     count(r.positive_controls_permissive),
                     f"{integer(r.high_confidence_events)} ({integer(r.genes)})",
                     count(r.positive_controls_high_confidence),
                     genes(r.positive_control_genes_high_confidence),
                     genes(r.Tier1_genes), genes(r.SOCE_machinery_genes), null])
    head = ["Comparison", "Permissive calls (genes)", "Positive controls, permissive",
            "High-confidence calls (genes)", "Positive controls, high-confidence",
            "Positive-control genes, high-confidence", "Tier 1 genes", "SOCE-machinery genes",
            "Split-control null: calls (ratio)"]
    note = ("*Note.* Calls are unannotated splicing changes in the regtools junction set, with the "
            "number of genes carrying them in parentheses (Methods 2.5). Positive controls are "
            "counted among the sixteen literature cryptic genes (Supplementary Table S6); n/a, not "
            "assessed in mouse. The null column gives the high-confidence calls of the "
            "control-versus-control split and their ratio to the real calls; —, fewer than four "
            "control replicates. iPSC-MN, iPSC-derived motor neurons; KD, knockdown.")
    return md_table(head, rows, "lrrrrllll", [18, 10, 9, 10, 9, 32, 10, 9, 11]), note


def table4():
    t = pd.read_csv(f"{TAB}/Table1_transcript_family_abundance.csv")
    rows, last = [], None
    for _, r in t.iterrows():
        fam = str(r.Family).replace("Ca2+", "Ca²⁺")
        total = r.Gene == "FAMILY TOTAL"
        name = "Family total" if total else (f"***{r.Gene}***" if r.dominant == "dominant"
                                             else gene(r.Gene))
        rows.append([fam if fam != last else "", name,
                     f"{r.TPM_control:.2f}", f"{r.TPM_KD:.2f}",
                     f"{r.TPM_control_adjusted:.2f}", f"{r.TPM_KD_adjusted:.2f}",
                     f"{r.share_of_family_control_pct:.1f}", signed(r.change_pct_adjusted, 1),
                     signed(r.log2FC, 3), pval(r.padj)])
        last = fam
    head = ["Family", "Gene", "TPM, control", "TPM, knockdown", "Adjusted TPM, control",
            "Adjusted TPM, knockdown", "Share of family, control (%)", "Adjusted change (%)",
            "log2FC", "p_adj"]
    # two lines at most: a third pushes the note of this full-page table onto its own page
    note = ("*Note.* TPM (transcripts per million), mean of three libraries per group; adjusted "
            "TPM, after per-library median-of-ratios scaling for library composition (Methods "
            "2.10). Share, percentage of the family total in control cells; bold, dominant member "
            "(> 50%). log2FC and p_adj from DESeq2 on the gene-level Salmon counts; —, not "
            "computed; CRAC, Ca²⁺ release-activated Ca²⁺.")
    return md_table(head, rows, "llrrrrrrrr", [38, 16, 17, 18, 17, 18, 20, 15, 11, 16]), note


def table5():
    t = pd.read_csv(f"{TAB}/Table5_cross_disease_comparison.csv")
    COH = {"ALS (NYGC GSE153960)": "ALS (NYGC, GSE153960)",
           "GSE123496": "Multiple sclerosis (GSE123496)",
           "GSE138614": "Multiple sclerosis (GSE138614)"}
    REG = {"All five regions pooled (centred within region)": "Five regions pooled",
           "Normal-appearing white matter vs control white matter": "Normal-appearing white matter",
           "MS lesions vs control white matter": "Lesions",
           "Donor level (10 MS vs 5 control donors)": "All samples, averaged per donor"}
    rows, last = [], None
    keys = t[["cohort", "region"]].drop_duplicates().itertuples(index=False)
    for coh, reg in keys:
        sub = t[(t.cohort == coh) & (t.region == reg)].set_index("gene")
        first = sub.iloc[0]
        cells = []
        for g in ["TRPC1", "SARAF", "CBARP"]:
            if g not in sub.index:
                cells.append("—")
                continue
            r = sub.loc[g]
            c = f"{signed(r.cliffs_delta)} ({pval(r.q_value)})"
            cells.append(f"**{c}**" if float(r.q_value) < 0.05 else c)
        name = COH.get(coh, coh)
        rows.append([name if name != last else "", REG.get(reg, reg),
                     f"{integer(first.n_case)}/{integer(first.n_control)}"] + cells)
        last = name
    head = ["Cohort", "Region", "n, case/control", "*TRPC1* δ (q)", "*SARAF* δ (q)",
            "*CBARP* δ (q)"]
    note = ("*Note.* δ, Cliff's delta, case minus control, with the Benjamini–Hochberg q value in "
            "parentheses; bold, q < 0.05; —, not tested. n counts samples, or donors in the "
            "donor-level row. The NYGC comparisons of each region share its non-neurological "
            "controls. BA9, Brodmann area 9; NYGC, New York Genome Center.")
    return md_table(head, rows, "lllrrr", [22, 24, 10, 16, 16, 16]), note


TABLES = {1: table1, 2: table2, 3: table3, 4: table4, 5: table5}


def blocks():
    out = {}
    for n, fn in TABLES.items():
        tab, note = fn()
        out[n] = f"<!-- table:{n} -->\n\n{tab}\n\n{note}\n\n<!-- /table:{n} -->"
    return out


def write_markdown(text):
    """Insert or refresh each table block after its caption paragraph."""
    for n, blk in blocks().items():
        pat = re.compile(rf"<!-- table:{n} -->.*?<!-- /table:{n} -->", re.S)
        if pat.search(text):
            text = pat.sub(lambda m: blk, text)
        else:
            cap = re.search(rf"^\*\*Table {n}\.\*\*[^\n]*\n", text, re.M)
            assert cap, f"caption of Table {n} not found"
            text = text[:cap.end()] + "\n" + blk + "\n" + text[cap.end():]
    return text


# ------------------------------------------------------------------ Word post-processing
FOOTER_ID = "rIdPageNumberFooter"
FOOTER = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
          '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
          '<w:p><w:pPr><w:suppressLineNumbers />'
          '<w:spacing w:before="0" w:after="0" w:line="240" w:lineRule="auto" />'
          '<w:jc w:val="center" /></w:pPr>'
          '<w:r><w:fldChar w:fldCharType="begin" /></w:r>'
          '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
          '<w:r><w:fldChar w:fldCharType="separate" /></w:r>'
          '<w:r><w:t>1</w:t></w:r>'
          '<w:r><w:fldChar w:fldCharType="end" /></w:r></w:p></w:ftr>')
SINGLE = 'w:line="240" w:lineRule="auto"'
# Word numbers lines per section, LibreOffice for the whole document; paragraphs outside the
# main text say so explicitly, so both number the main text only
NOLN = "<w:suppressLineNumbers />"


def sect_pr(orient, line_numbers=False):
    w, h = (PAGE_W, PAGE_H) if orient == "portrait" else (PAGE_H, PAGE_W)
    o = ' w:orient="landscape"' if orient == "landscape" else ""
    ln = ('<w:lnNumType w:countBy="1" w:distance="500" w:restart="continuous" />'
          if line_numbers else "")
    return (f'<w:sectPr><w:footerReference w:type="default" r:id="{FOOTER_ID}" />'
            '<w:footnotePr><w:numRestart w:val="eachSect" /></w:footnotePr>'
            '<w:type w:val="nextPage" />'
            f'<w:pgSz w:w="{w}" w:h="{h}"{o} />'
            f'<w:pgMar w:top="{MARGIN}" w:right="{MARGIN}" w:bottom="{MARGIN}" w:left="{MARGIN}" '
            f'w:header="708" w:footer="708" w:gutter="0" />{ln}</w:sectPr>')


def add_to_ppr(par, xml, last=True):
    """Put xml into the paragraph properties: last child (sectPr) or right after pStyle."""
    if "<w:pPr>" in par:
        if last:
            return par.replace("</w:pPr>", xml + "</w:pPr>", 1)
        m = re.search(r"<w:pStyle [^>]*/>", par)
        at = m.end() if m else par.index("<w:pPr>") + len("<w:pPr>")
        return par[:at] + xml + par[at:]
    at = re.match(r"<w:p(?: [^>]*)?>", par).end()
    return par[:at] + "<w:pPr>" + xml + "</w:pPr>" + par[at:]


def size_runs(xml, half_points=18):
    sz = f'<w:sz w:val="{half_points}" /><w:szCs w:val="{half_points}" />'
    xml = re.sub(r"(<w:r>)(?!<w:rPr>)", r"\1<w:rPr>" + sz + "</w:rPr>", xml)
    return re.sub(r"(<w:r><w:rPr>(?:(?!</w:rPr>).)*?)(</w:rPr>)",
                  lambda m: m.group(1) + (sz if "<w:sz " not in m.group(1) else "") + m.group(2),
                  xml)


def format_table(tbl, width):
    grid = [int(v) for v in re.findall(r'<w:gridCol w:w="(\d+)" />', tbl)]
    cols = [round(width * g / sum(grid)) for g in grid]
    cols[-1] += width - sum(cols)
    tbl = re.sub(r"<w:tblW [^>]*/>",
                 f'<w:tblW w:type="dxa" w:w="{width}" /><w:tblBorders>'
                 '<w:top w:val="single" w:sz="8" w:space="0" w:color="000000" />'
                 '<w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000" />'
                 '</w:tblBorders>', tbl, count=1)
    tbl = re.sub(r"<w:tblGrid>.*?</w:tblGrid>",
                 "<w:tblGrid>" + "".join(f'<w:gridCol w:w="{c}" />' for c in cols) + "</w:tblGrid>",
                 tbl, count=1, flags=re.S)

    def row(m):
        i = iter(cols)
        return re.sub(r"<w:tcPr />", lambda _: f'<w:tcPr><w:tcW w:w="{next(i)}" w:type="dxa" /></w:tcPr>',
                      m.group(0))
    tbl = re.sub(r"<w:tr>.*?</w:tr>", row, tbl, flags=re.S)
    tbl = tbl.replace('<w:pStyle w:val="Compact" />',
                      f'<w:pStyle w:val="Compact" />{NOLN}<w:spacing w:before="0" w:after="0" {SINGLE} />')
    # the paragraph mark sets the height of an empty cell, so it is 9 pt as well
    tbl = tbl.replace("</w:pPr>", '<w:rPr><w:sz w:val="18" /><w:szCs w:val="18" /></w:rPr></w:pPr>')
    return size_runs(tbl)


def postprocess(docx):
    with zipfile.ZipFile(docx) as z:
        files = {n: z.read(n) for n in z.namelist()}
    x = files["word/document.xml"].decode("utf-8")

    # the main text ends with the paragraph before the "Tables" heading
    h = re.search(r"<w:p><w:pPr><w:pStyle w:val=\"Heading2\" /></w:pPr>"
                  r"(?:(?!</w:p>).)*?<w:t xml:space=\"preserve\">Tables</w:t>", x, re.S)
    assert h, "Tables heading not found"
    ps = x.rfind("<w:p>", 0, h.start())
    pe = x.index("</w:p>", ps) + len("</w:p>")
    x = x[:ps] + add_to_ppr(x[ps:pe], sect_pr("portrait", line_numbers=True)) + x[pe:]
    hs = x.index("<w:p>", ps + 1)                  # the "Tables" heading itself
    he = x.index("</w:p>", hs) + len("</w:p>")
    assert ">Tables</w:t>" in x[hs:he]
    x = x[:hs] + add_to_ppr(x[hs:he], NOLN, last=False) + x[he:]

    tables = list(re.finditer(r"<w:tbl>.*?</w:tbl>", x, re.S))
    assert len(tables) == 5, f"expected 5 tables, found {len(tables)}"
    for n in range(5, 0, -1):                      # right to left keeps offsets valid
        m = list(re.finditer(r"<w:tbl>.*?</w:tbl>", x, re.S))[n - 1]
        orient = ORIENT[n]
        # note paragraph after the table: 9 pt, and it closes the table's section
        ns = x.index("<w:p>", m.end())
        ne = x.index("</w:p>", ns) + len("</w:p>")
        note = add_to_ppr(size_runs(x[ns:ne]), f'{NOLN}<w:spacing w:before="120" w:after="0" {SINGLE} />',
                          last=False)
        if n < 5:
            note = add_to_ppr(note, sect_pr(orient))
        x = x[:ns] + note + x[ne:]
        x = x[:m.start()] + format_table(m.group(0), TEXT_W[orient]) + x[m.end():]
        # the caption keeps with its table
        cs = x.rfind("<w:p>", 0, m.start())
        ce = x.index("</w:p>", cs) + len("</w:p>")
        cap = x[cs:ce]
        assert f">Table {n}.</w:t>" in cap, f"caption of Table {n} not found before it"
        cap = add_to_ppr(cap, f'<w:keepNext />{NOLN}<w:spacing w:before="0" w:after="120" {SINGLE} />',
                         last=False)
        x = x[:cs] + cap + x[ce:]

    # the final section (Table 5)
    x = re.sub(r"<w:sectPr>(?:(?!<w:sectPr>).)*?</w:sectPr>\s*</w:body>",
               sect_pr(ORIENT[5]) + "</w:body>", x, count=1, flags=re.S)
    files["word/document.xml"] = x.encode("utf-8")

    # double spacing is the document default; tables, captions and notes set single spacing
    st = files["word/styles.xml"].decode("utf-8")
    st, n = re.subn(r"(<w:pPrDefault>\s*<w:pPr>\s*<w:spacing )([^>]*?)\s*/>",
                    lambda m: m.group(1) + m.group(2) + ' w:line="480" w:lineRule="auto" />',
                    st, count=1)
    assert n == 1, "default paragraph spacing not found"
    files["word/styles.xml"] = st.encode("utf-8")

    # page numbers: one footer part, referenced from every section
    assert "word/footer1.xml" not in files
    files["word/footer1.xml"] = FOOTER.encode("utf-8")
    rel = files["word/_rels/document.xml.rels"].decode("utf-8")
    files["word/_rels/document.xml.rels"] = rel.replace(
        "</Relationships>",
        '<Relationship Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/'
        f'footer" Id="{FOOTER_ID}" Target="footer1.xml" /></Relationships>').encode("utf-8")
    ct = files["[Content_Types].xml"].decode("utf-8")
    files["[Content_Types].xml"] = ct.replace(
        "</Types>",
        '<Override PartName="/word/footer1.xml" ContentType="application/'
        'vnd.openxmlformats-officedocument.wordprocessingml.footer+xml" /></Types>').encode("utf-8")
    tmp = docx + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for n, b in files.items():
            z.writestr(n, b)
    os.replace(tmp, docx)


# ------------------------------------------------------------------ main
if __name__ == "__main__":
    text = io.open(MD, encoding="utf-8").read()
    new = write_markdown(text)
    if "--check" in sys.argv:
        print("tables up to date" if new == text else "tables OUT OF DATE")
        sys.exit(0 if new == text else 1)
    if new != text:
        io.open(MD, "w", encoding="utf-8").write(new)
        print("Markdown tables refreshed")
    subprocess.run(["pandoc", MD, "-o", DOCX], check=True)
    postprocess(DOCX)
    # Keep all manuscript headings black for a neutral submission layout.
    from docx import Document
    from docx.shared import RGBColor
    doc = Document(DOCX)
    for style_name in ("Heading 1", "Heading 2", "Heading 3", "Heading 4", "Title"):
        if style_name in doc.styles:
            doc.styles[style_name].font.color.rgb = RGBColor(0, 0, 0)
    doc.save(DOCX)
    print("written:", DOCX)
