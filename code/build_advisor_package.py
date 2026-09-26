#!/usr/bin/env python3
"""Build the package sent to the supervisor: the manuscript to read, everything else zipped.

Layout produced in ~/Desktop/MAKALE/HOCAYA_GONDERIM_<date>/
    MAKALE_HOCAYA_REVIZYON.docx   the manuscript, loose, figures embedded in the text
    EKLER_<date>.zip              supplementary material, highlights, figures, tables and data

Run: /usr/bin/python3 code/build_advisor_package.py 2026-09-26
"""
import hashlib
import shutil
import sys
import zipfile
from pathlib import Path

DATE = sys.argv[1] if len(sys.argv) > 1 else "2026-09-26"
ROOT = Path("/Users/elmas/Desktop/MAKALE")
PKG = ROOT / "09_YAYIN_PAKETI"
OUT = ROOT / f"HOCAYA_GONDERIM_{DATE}"
STAGE = OUT / f"EKLER_{DATE}"

if OUT.exists():
    shutil.rmtree(OUT)
STAGE.mkdir(parents=True)

# the manuscript stays outside the archive so it can be opened directly
shutil.copy2(PKG / "manuscript/MANUSCRIPT_NEUROCHEMISTRY_INTERNATIONAL_FINAL.docx",
             OUT / "MAKALE_HOCAYA_REVIZYON.docx")

copies = [
    (PKG / "supplementary/SUPPLEMENTARY_MATERIAL.docx", "01_EK_MATERYAL/SUPPLEMENTARY_MATERIAL.docx"),
    (PKG / "highlights_Neurochemistry_International.docx", "01_EK_MATERYAL/HIGHLIGHTS.docx"),
    (PKG / "figures/graphical_abstract.png", "02_GORSELLER/Grafik_Ozet/graphical_abstract.png"),
    (PKG / "figures/graphical_abstract.svg", "02_GORSELLER/Grafik_Ozet/graphical_abstract.svg"),
]
# main figures keep the numbering used in the manuscript
for n, stem in enumerate(sorted(p.stem for p in (PKG / "figures/main").glob("Figure*.png")), 1):
    for ext in ("png", "svg"):
        copies.append((PKG / f"figures/main/{stem}.{ext}", f"02_GORSELLER/Ana_Figurler/Figure{n}.{ext}"))
for p in sorted((PKG / "figures/supplementary").glob("Supplementary_Figure_S*")):
    if p.suffix in (".png", ".svg"):
        copies.append((p, f"02_GORSELLER/Ek_Figurler/{p.name}"))
for p in sorted((PKG / "tables").glob("Table*.csv")):
    copies.append((p, f"03_TABLOLAR/{p.name}"))
for p in sorted((PKG / "supplementary").glob("S*")):
    if p.suffix in (".csv", ".xlsx", ".gz") and p.name[1].isdigit():
        copies.append((p, f"04_EK_VERI/{p.name}"))

for src, rel in copies:
    dst = STAGE / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

readme = STAGE / "00_OKU_BENI.txt"
readme.write_text(f"""MAKALE EKLERI — {DATE[8:]}.{DATE[5:7]}.{DATE[:4]}

Makale bu arsivin disinda, ayri bir dosyadir: MAKALE_HOCAYA_REVIZYON.docx
Alti ana figur metnin icine gomulu oldugu icin makale tek basina okunabilir; bu arsiv
figurlerin ayri dosyalarini, ek materyali ve veri dosyalarini icerir.

01_EK_MATERYAL
  SUPPLEMENTARY_MATERIAL.docx  Ek Tablo S1-S17 aciklamalari ve Ek Figur S1-S8 (gorsellerle)
  HIGHLIGHTS.docx              Dort maddelik one cikanlar listesi

02_GORSELLER
  Ana_Figurler   Figure1-6, PNG (baski cozunurlugu) ve SVG (duzenlenebilir)
  Ek_Figurler    Supplementary Figure S1-S8, PNG ve SVG
  Grafik_Ozet    Grafiksel ozet

03_TABLOLAR      Makaledeki Tablo 1-5'in CSV surumleri
04_EK_VERI       Ek Tablo S1-S17 veri dosyalari (S1 laboratuvar ham verisi XLSX)

DOSYA_LISTESI_SHA256.txt her dosyanin saglama toplamini verir.
Kod ve veri deposu: github.com/elmasnuryilmaz/tdp43-soce-manuscript (surum v1.0.5)
""", encoding="utf-8")

files = sorted(p for p in STAGE.rglob("*") if p.is_file() and p.name != "DOSYA_LISTESI_SHA256.txt")
lines = []
for p in files:
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    lines.append(f"{h}  {p.relative_to(STAGE)}")
(STAGE / "DOSYA_LISTESI_SHA256.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

archive = OUT / f"EKLER_{DATE}.zip"
with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
    for p in sorted(STAGE.rglob("*")):
        if p.is_file():
            z.write(p, Path(f"EKLER_{DATE}") / p.relative_to(STAGE))
shutil.rmtree(STAGE)

with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    print(f"{archive.name}: {len(z.namelist())} dosya, {archive.stat().st_size / 1e6:.1f} MB")
print(f"{OUT.name}/ hazir")
