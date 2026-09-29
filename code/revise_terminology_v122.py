"""Rename analysis labels in the data files and in the code that reads them (release v1.0.8).

'High-confidence' implied a validated specificity that the control-versus-control null does not
support (null-to-real call ratios of 0.64, 2.17 and 0.83), so the calls are named for what they
are: unannotated splicing changes that pass a stringent filter. 'SOCE machinery' is replaced by
the fixed gene lists it stood for. The numbers in every file are untouched; only labels change.

Run once from the repository root:  /usr/bin/python3 code/revise_terminology_v122.py
"""
import subprocess
from pathlib import Path

import pandas as pd

R = Path(__file__).resolve().parents[1]

COLUMNS = {
    "high_confidence_events": "stringent_filter_events",
    "positive_controls_high_confidence": "positive_controls_stringent_filter",
    "positive_control_genes_high_confidence": "positive_control_genes_stringent_filter",
    "SOCE_machinery_genes": "SOCE_associated_set_genes",
    "null_calls_high_confidence": "null_calls_stringent_filter",
}

for rel in ("tables/Table4_cryptic_events_eleven_comparisons.csv",
            "supplementary/S12_cryptic_counts_by_dataset.csv"):
    path = R / rel
    text = path.read_text(encoding="utf-8")
    header, _, body = text.partition("\n")
    cols = header.split(",")
    assert all(c in cols for c in COLUMNS), rel
    header = ",".join(COLUMNS.get(c, c) for c in cols)
    path.write_text(header + "\n" + body, encoding="utf-8")

s14 = R / "supplementary/S14_control_vs_control_null_test.csv"
t = s14.read_text(encoding="utf-8")
assert t.count("tier 2 (high-confidence)") == 3
s14.write_text(t.replace("tier 2 (high-confidence)", "tier 2 (stringent filter)"), encoding="utf-8")

old = R / "supplementary/S5_high_confidence_cryptic_events.csv"
new = R / "supplementary/S5_stringent_filter_unannotated_splicing_candidates.csv"
if old.exists():
    subprocess.run(["git", "mv", str(old), str(new)], check=True, cwd=R)


# generators and readers of the renamed files (bare identifiers cover quoted and attribute use)
path = R / "code/build_tables_v3.py"
text = path.read_text(encoding="utf-8")
for a, b in COLUMNS.items():
    assert a in text, a
    text = text.replace(a, b)
text = text.replace('"S5_high_confidence_cryptic_events"', '"S5_stringent_filter_unannotated_splicing_candidates"')
text = text.replace('"tier 2 (high-confidence)"', '"tier 2 (stringent filter)"')
assert "high_confidence" not in text.replace("S13_yuksek_guven", "")
path.write_text(text, encoding="utf-8")

path = R / "code/fig_nmd_descriptive_v118.py"
text = path.read_text(encoding="utf-8")
assert "S5_high_confidence_cryptic_events.csv" in text
path.write_text(text.replace("S5_high_confidence_cryptic_events.csv",
                             "S5_stringent_filter_unannotated_splicing_candidates.csv"), encoding="utf-8")

# the historical baseline builders still regenerate the v4 Markdown tables, which keep the old
# wording; they read the renamed CSV columns through this alias
ALIAS = '''
    t = t.rename(columns={"stringent_filter_events": "high_confidence_events",
                          "positive_controls_stringent_filter": "positive_controls_high_confidence",
                          "positive_control_genes_stringent_filter": "positive_control_genes_high_confidence",
                          "SOCE_associated_set_genes": "SOCE_machinery_genes",
                          "null_calls_stringent_filter": "null_calls_high_confidence"})
'''
for rel in ("code/build_manuscript_docx.py", "code/build_manuscript_docx_final.py"):
    path = R / rel
    text = path.read_text(encoding="utf-8")
    anchor = 't = pd.read_csv(f"{TAB}/Table4_cryptic_events_eleven_comparisons.csv").set_index("comparison")'
    assert text.count(anchor) == 1
    text = text.replace(anchor, 't = pd.read_csv(f"{TAB}/Table4_cryptic_events_eleven_comparisons.csv")' + ALIAS.rstrip("\n")
                        + '\n    t = t.set_index("comparison")')
    path.write_text(text, encoding="utf-8")
print("terminology renamed in data files and readers")
