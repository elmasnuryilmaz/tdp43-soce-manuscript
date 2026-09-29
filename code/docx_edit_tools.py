"""Helpers for editing the manuscript DOCX without disturbing Zotero citation fields.

A paragraph is read as a sequence of *atoms*: ordinary text runs and citation groups. A citation
group is the five-run Zotero structure (begin, instrText, separate, visible result, end); it is
copied verbatim when a paragraph is rebuilt, so the field code and its visible text can never drift
apart. Text is given as light markup: ``**bold**``, ``*italic*`` and ``[[cite:LABEL]]`` where LABEL is a
substring of the visible citation text (it must identify one citation of the original paragraph).
Gene symbols in ``GENES`` are italicised automatically.
"""
from __future__ import annotations

import copy
import re
from pathlib import Path

from docx import Document
from docx.oxml.ns import qn
from docx.text.paragraph import Paragraph

GENES = [
    "TARDBP", "STMN2", "UNC13A", "ACTL6B", "PFKP", "HDGFL2", "AGRN", "ARHGAP32", "ATG4B",
    "ELAVL3", "SETD5", "RSF1", "GPSM2", "KALRN", "CAMK2B", "POLDIP3", "SYNJ2", "TRPM3",
    "STIM1", "STIM2", "STIMATE", "SARAF", "CRACR2A", "CRACR2B", "ORAI1", "ORAI2", "ORAI3",
    "TRPC1", "TRPC3", "TRPC4", "TRPC5", "TRPC6", "ATP2A1", "ATP2A2", "ATP2A3", "CBARP", "SELENOK",
    "MCU", "MCUB", "MICU1", "MICU2", "MICU3", "MCUR1", "ITPR1", "ITPR2", "ITPR3", "RYR2", "RYR3",
    "SNAP25", "GFAP", "MBP", "PLP1", "MOG", "MAG", "XRN1", "UPF1", "SMG6", "ATP2B1", "ATP2B2",
    "ATP2B3", "ATP2B4", "CHGA", "RBFOX3", "TBP", "B2M", "GAPDH", "AIF1",
    "Stmn2", "Saraf", "Atp2a2", "Trpc1", "Cbarp",
]
_GENE_RE = re.compile(r"(?<![A-Za-z0-9])(" + "|".join(sorted(map(re.escape, GENES), key=len, reverse=True)) + r")(?![A-Za-z0-9])")
_TOKEN_RE = re.compile(r"(__|\*\*|\*|\[\[cite:.+?\]\])")


# ------------------------------------------------------------------ reading
def is_citation_begin(run) -> bool:
    return any(c.get(qn("w:fldCharType")) == "begin" for c in run.findall(qn("w:fldChar")))


def _fld_type(run):
    for c in run.findall(qn("w:fldChar")):
        return c.get(qn("w:fldCharType"))
    return None


def atoms(paragraph):
    """Split a paragraph into text-run atoms and complete field groups (any complex field)."""
    out, group, depth = [], None, 0
    for child in paragraph._p:
        if child.tag == qn("w:pPr"):
            continue
        if child.tag != qn("w:r"):
            out.append({"kind": "other", "el": child})
            continue
        ft = _fld_type(child)
        if ft == "begin" and group is None:
            group, depth = [child], 1
        elif group is not None:
            group.append(child)
            if ft == "begin":
                depth += 1
            elif ft == "end":
                depth -= 1
                if depth == 0:
                    text = "".join(t.text or "" for r in group for t in r.findall(qn("w:t")))
                    instr = "".join(t.text or "" for r in group for t in r.findall(qn("w:instrText")))
                    out.append({"kind": "field", "els": group, "label": text, "instr": instr})
                    group = None
        else:
            out.append({"kind": "run", "el": child,
                        "text": "".join(t.text or "" for t in child.findall(qn("w:t")))})
    if group is not None:  # a field that continues into the next paragraph (bibliography)
        out.append({"kind": "other", "el": group[0]})
        out.extend({"kind": "other", "el": g} for g in group[1:])
    return out


def describe(paragraph, width=100000):
    parts = []
    for a in atoms(paragraph):
        if a["kind"] == "field":
            parts.append("{{" + a["label"] + "}}")
        elif a["kind"] == "run":
            parts.append(a["text"])
        else:
            parts.append("<" + a["el"].tag.split("}")[1] + ">")
    return "".join(parts)[:width]


# ------------------------------------------------------------------ rebuilding
def _clean_rpr(run):
    rpr = run.find(qn("w:rPr"))
    if rpr is None:
        return None
    rpr = copy.deepcopy(rpr)
    for tag in ("w:b", "w:bCs", "w:i", "w:iCs"):
        for el in rpr.findall(qn(tag)):
            rpr.remove(el)
    return rpr


def _make_run(template_rpr, text, bold=False, italic=False):
    from docx.oxml import OxmlElement
    r = OxmlElement("w:r")
    if template_rpr is not None or bold or italic:
        rpr = copy.deepcopy(template_rpr) if template_rpr is not None else OxmlElement("w:rPr")
        # OOXML order: rFonts, b, i ... insert b/i before other properties that must follow them
        if italic:
            rpr.insert(0, OxmlElement("w:i"))
        if bold:
            rpr.insert(0, OxmlElement("w:b"))
        r.append(rpr)
    t = OxmlElement("w:t")
    t.text = text
    t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    r.append(t)
    return r


def _pieces(markup):
    """Yield (kind, payload, bold, italic, upright) from the light markup.

    ``__text__`` marks upright text (for example protein names) that is never auto-italicised.
    """
    bold = italic = upright = False
    for tok in _TOKEN_RE.split(markup):
        if tok == "":
            continue
        if tok == "**":
            bold = not bold
        elif tok == "*":
            italic = not italic
        elif tok == "__":
            upright = not upright
        elif tok.startswith("[[cite:") and tok.endswith("]]"):
            yield ("cite", tok[7:-2], bold, italic, upright)
        else:
            yield ("text", tok, bold, italic, upright)


def italicise_genes(text, bold, italic):
    """Split text on gene symbols; yield (chunk, italic_flag). Explicit *italic* wins."""
    if italic:
        yield text, True
        return
    pos = 0
    for m in _GENE_RE.finditer(text):
        if m.start() > pos:
            yield text[pos:m.start()], False
        yield m.group(0), True
        pos = m.end()
    if pos < len(text):
        yield text[pos:], False


def rebuild(paragraph, markup, template_rpr=None, auto_genes=True, extra=()):
    """Replace a paragraph's content by ``markup`` keeping its citation groups verbatim.

    ``extra`` lists further paragraphs whose citation groups may also be used (for merged paragraphs).
    """
    at = atoms(paragraph)
    assert not any(a["kind"] == "other" for a in at), "paragraph has non-run content"
    for a in at:
        if a["kind"] == "run":
            assert not a["el"].findall(qn("w:tab")) and not a["el"].findall(qn("w:br")) \
                and not a["el"].findall(qn("w:drawing")), "paragraph has tabs, breaks or drawings"
    fields = [a for a in at if a["kind"] == "field"]
    extra_fields = [a for q in extra for a in atoms(q) if a["kind"] == "field"]
    if template_rpr is None:
        base = next((a["el"] for a in at if a["kind"] == "run" and a["text"].strip()), None)
        template_rpr = _clean_rpr(base) if base is not None else None
    used = set()
    new_nodes = []
    for kind, payload, bold, italic, upright in _pieces(markup):
        if kind == "cite":
            hits = [f for f in fields + extra_fields if payload in f["label"]]
            assert len(hits) == 1, f"citation label {payload!r} matches {len(hits)} groups: {[f['label'] for f in fields + extra_fields]}"
            f = hits[0]
            used.add(id(f))
            new_nodes.extend(copy.deepcopy(r) for r in f["els"])
        else:
            chunks = italicise_genes(payload, bold, italic) if (auto_genes and not upright) else [(payload, italic)]
            for chunk, it in chunks:
                new_nodes.append(_make_run(template_rpr, chunk, bold=bold, italic=it))
    dropped = [f["label"] for f in fields if id(f) not in used]
    assert not dropped, f"rebuild would drop citations: {dropped}"
    p = paragraph._p
    for a in at:
        for el in (a["els"] if a["kind"] == "field" else [a["el"]]):
            p.remove(el)
    for n in new_nodes:
        p.append(n)


def edit_text(paragraph, old, new):
    """Replace ``old`` by ``new`` (plain text) across ordinary text nodes of one paragraph.

    Works on the individual ``w:t`` nodes, so tabs and line breaks inside runs stay in place.
    Text inside citation fields, hyperlinks and other non-run content is never touched, and a
    match may not span such an element; use ``rebuild`` for paragraphs that need larger changes.
    """
    segments, cur = [], []
    for a in atoms(paragraph):
        if a["kind"] == "run":
            cur.extend(a["el"].findall(qn("w:t")))
        else:
            segments.append(cur)
            cur = []
    segments.append(cur)
    hits = []
    for nodes in segments:
        whole = "".join(t.text or "" for t in nodes)
        start = whole.find(old)
        while start != -1:
            hits.append((nodes, start))
            start = whole.find(old, start + 1)
    assert len(hits) == 1, f"span {old[:70]!r} found {len(hits)} times outside citations"
    nodes, start = hits[0]
    end = start + len(old)
    pos, done = 0, False
    for t in nodes:
        txt = t.text or ""
        a0, a1 = pos, pos + len(txt)
        pos = a1
        if a1 <= start or a0 >= end:
            continue
        lo, hi = max(start, a0) - a0, min(end, a1) - a0
        t.text = txt[:lo] + (new if not done else "") + txt[hi:]
        t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        done = True


# ------------------------------------------------------------------ structure helpers
def italicise_in_place(paragraph):
    """Italicise gene symbols that sit in ordinary, non-italic runs by splitting the run around each symbol.

    Only runs made of one text node are split; code-styled runs, citation fields and runs that are
    already italic are left alone, and the paragraph text is unchanged. Returns the number of symbols
    made italic. ``STIM2.1`` (the isoform) is not treated as the gene ``STIM2``.
    """
    from docx.oxml import OxmlElement
    changed = 0
    for a in atoms(paragraph):
        if a["kind"] != "run":
            continue
        run = a["el"]
        rpr = run.find(qn("w:rPr"))
        if _on(rpr, "w:i") or (rpr is not None and rpr.find(qn("w:rStyle")) is not None):
            continue
        kids = [c for c in run if c.tag != qn("w:rPr")]
        if len(kids) != 1 or kids[0].tag != qn("w:t"):
            continue
        text = kids[0].text or ""
        hits = [m for m in _GENE_RE.finditer(text) if not re.match(r"\.\d", text[m.end():m.end() + 2])]
        if not hits:
            continue
        pieces, pos = [], 0
        for m in hits:
            if m.start() > pos:
                pieces.append((text[pos:m.start()], False))
            pieces.append((m.group(0), True))
            pos = m.end()
        if pos < len(text):
            pieces.append((text[pos:], False))
        parent = run.getparent()
        at = parent.index(run)
        for j, (chunk, ital) in enumerate(pieces):
            new = copy.deepcopy(run)
            node = new.find(qn("w:t"))
            node.text = chunk
            node.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
            if ital:
                nrpr = new.find(qn("w:rPr"))
                if nrpr is None:
                    nrpr = OxmlElement("w:rPr")
                    new.insert(0, nrpr)
                pos_i = 0
                for k, child in enumerate(nrpr):  # after rStyle, rFonts, b, bCs (OOXML order)
                    if child.tag in (qn("w:rStyle"), qn("w:rFonts"), qn("w:b"), qn("w:bCs")):
                        pos_i = k + 1
                nrpr.insert(pos_i, OxmlElement("w:i"))
            parent.insert(at + j, new)
        parent.remove(run)
        changed += len(hits)
    return changed


def find(doc, startswith=None, contains=None, style=None, nth=0):
    hits = [p for p in doc.paragraphs
            if (startswith is None or p.text.startswith(startswith))
            and (contains is None or contains in p.text)
            and (style is None or p.style.name == style)]
    assert len(hits) > nth, f"paragraph not found: {startswith or contains!r}"
    if nth == 0:
        assert len(hits) == 1, f"paragraph not unique ({len(hits)}): {startswith or contains!r}"
    return hits[nth]


def insert_after(paragraph, markup, style=None, like=None):
    """Insert a new paragraph after ``paragraph``, formatted like ``like`` (default: paragraph)."""
    from docx.oxml import OxmlElement
    src = like if like is not None else paragraph
    node = OxmlElement("w:p")
    ppr = src._p.find(qn("w:pPr"))
    if ppr is not None:
        node.append(copy.deepcopy(ppr))
    paragraph._p.addnext(node)
    new = Paragraph(node, paragraph._parent)
    if style:
        new.style = style
    base = next((a["el"] for a in atoms(src) if a["kind"] == "run" and a["text"].strip()), None)
    rpr = _clean_rpr(base) if base is not None else None
    for kind, payload, bold, italic, upright in _pieces(markup):
        assert kind == "text", "no citations in inserted paragraphs"
        for chunk, it in ([(payload, italic)] if upright else italicise_genes(payload, bold, italic)):
            node.append(_make_run(rpr, chunk, bold=bold, italic=it))
    return new


def delete(paragraph, citations_reused=False):
    """Remove a paragraph; one with citations is deleted only when they were copied into a merged paragraph."""
    assert citations_reused or not paragraph._p.xpath(".//w:instrText"), "refusing to delete a paragraph with citations"
    paragraph._p.getparent().remove(paragraph._p)


def move_after(paragraph, anchor):
    anchor._p.addnext(paragraph._p)


def retarget_hyperlinks(doc, paragraph, old, new):
    """Replace ``old`` by ``new`` in the visible text and in the URL of hyperlinks inside ``paragraph``."""
    for t in paragraph._p.iter(qn("w:t")):
        if old in (t.text or ""):
            t.text = t.text.replace(old, new)
    for h in paragraph._p.iter(qn("w:hyperlink")):
        rid = h.get(qn("r:id"))
        if rid:
            rel = doc.part.rels[rid]
            if old in rel.target_ref:
                rel._target = rel.target_ref.replace(old, new)


def field_count(doc):
    return len(doc._element.xpath(".//w:instrText"))


def replace_image_before(doc, caption_prefix, png, width_in=None):
    """Swap the picture in the paragraph that precedes the caption starting with ``caption_prefix``."""
    from docx.shared import Inches
    from PIL import Image
    idx = next(i for i, p in enumerate(doc.paragraphs) if p.text.startswith(caption_prefix))
    cands = [p for p in doc.paragraphs[max(0, idx - 3):idx] if p._p.xpath(".//a:blip")]
    assert len(cands) == 1, (caption_prefix, len(cands))
    pp = cands[0]
    blip = pp._p.xpath(".//a:blip")[0]
    rid = blip.get(qn("r:embed"))
    pp.part.related_parts[rid]._blob = Path(png).read_bytes()
    exts = pp._p.xpath(".//wp:extent") + pp._p.xpath(".//a:ext")
    cx = int(Inches(width_in)) if width_in else int(exts[0].get("cx"))
    w, h = Image.open(png).size
    cy = int(cx * h / w)
    for e in exts:
        e.set("cx", str(cx))
        e.set("cy", str(cy))
    return pp


def set_alt_text(paragraph, text):
    for dp in paragraph._p.xpath(".//wp:docPr"):
        dp.set("descr", text)


def _on(rpr, tag):
    """True when the run property ``tag`` (w:b, w:i) is present and not switched off with w:val="0"."""
    if rpr is None:
        return False
    el = rpr.find(qn(tag))
    return el is not None and el.get(qn("w:val")) not in ("0", "false", "off")


def to_markup(paragraph):
    """Return the paragraph as light markup (inverse of ``rebuild``); gene italics are left implicit."""
    out = []
    for a in atoms(paragraph):
        if a["kind"] == "field":
            out.append("[[cite:" + a["label"] + "]]")
        elif a["kind"] == "run":
            el = a["el"]
            rpr = el.find(qn("w:rPr"))
            b = _on(rpr, "w:b")
            i = _on(rpr, "w:i")
            t = a["text"]
            if not t:
                continue
            if b:
                t = "**" + t + "**"
            if i and t.strip() and t.strip() not in GENES:
                t = ("*" + t + "*") if not b else t
            out.append(t)
    return "".join(out)


def load(path):
    return Document(str(path))
