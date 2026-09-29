#!/usr/bin/env python3
"""Build the package sent to the supervisor: the manuscript to read, everything else zipped.

Layout produced in ~/Desktop/MAKALE/HOCAYA_GONDERIM_<date>/
    MAKALE_HOCAYA_REVIZYON.docx   the manuscript, loose, figures embedded in the text
    EKLER_<date>.zip              supplementary material, highlights, figures, tables and data

Run: python code/build_advisor_package.py 2026-09-29
"""
import hashlib
import shutil
import sys
import zipfile
from pathlib import Path

DATE = sys.argv[1] if len(sys.argv) > 1 else "2026-09-29"
ROOT = Path("/Users/elmas/Desktop/MAKALE")
PKG = ROOT / "09_YAYIN_PAKETI"
SUFFIX = sys.argv[2] if len(sys.argv) > 2 else ""
OUT = ROOT / f"HOCAYA_GONDERIM_{DATE}{SUFFIX}"
STAGE = OUT / f"EKLER_{DATE}"

if OUT.exists():
    raise FileExistsError(f"Refusing to overwrite an existing package: {OUT}")
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

for name in ['STIM1_transcript_test_eligibility.csv', 'STIM1_isoformswitch_DEXSeq.csv', 'nmd_descriptive_all_genes.csv.gz',
             'svaseq_sensitivity_SHSY5Y.csv', 'qpcr_biological_replicate_tests.csv', 'power_simulation_S1.csv',
             'DESeq2_featureCounts_vs_Salmon_calcium_genes.csv']:
    copies.append((PKG / 'source_data' / name, f'05_KAYNAK_KONTROLLER/{name}'))
for name in ['QPCR_CORRECTION_2026-09-29.md', 'CORRECTIONS_2026-09-29_v1.0.8.md']:
    copies.append((PKG / name, f'05_KAYNAK_KONTROLLER/{name}'))

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
  SUPPLEMENTARY_MATERIAL.docx  Ek Tablo S1-S18d aciklamalari, Ek Sonuc 1-6 ve Ek Figur S1-S9 (gorsellerle)
  HIGHLIGHTS.docx              Dort maddelik one cikanlar listesi

02_GORSELLER
  Ana_Figurler   Figure1-6, PNG (baski cozunurlugu) ve SVG (duzenlenebilir)
  Ek_Figurler    Supplementary Figure S1-S9, PNG ve SVG
  Grafik_Ozet    Istege bagli dergi grafik ozeti; ayni zamanda depo ozetidir. Makalenin numarali bir figuru
                 degildir ve ana metinde ya da ekte gosterilmeyen bir bulgu tasimaz.

03_TABLOLAR      Makaledeki Tablo 1-5'in CSV surumleri
04_EK_VERI       Ek Tablo S1-S18d veri dosyalari (S1 laboratuvar ham verisi XLSX)
05_KAYNAK_KONTROLLER  STIM1, NMD, svaseq, guc simulasyonu ve featureCounts/Salmon karsilastirma tablolari; duzeltme ozetleri

Bu surum (v1.0.8) bir onceki surumun hakem bicimli degerlendirmesine yanit verir: ana iddia kanit duzeyine
indirildi (Fura-2 tek plaka, betimsel ve RNA verisiyle eslestirilmemis gozlem; RNA analizleri aday uretir),
"high-confidence" adi "stringent-filter" olarak degistirildi, guc simulasyonunun varsayimlari yazildi,
Sekil 2, 4, 5, 6 ve grafik ozet yeniden cizildi, Ek Figur S9 eklendi, DESeq2 sayim kaynaklari (Salmon/featureCounts)
Yontem 2.2'de acikca ayrildi. Ayrintilar: 05_KAYNAK_KONTROLLER/CORRECTIONS_2026-09-29_v1.0.8.md.
Ek Figur aciklamalarinin bagimsiz kopyasi depoda: supplementary/Supplementary_Figure_Legends.md.

RT-qPCR grup basina dort biyolojik tekrardir. Fura-2 grup basina
tek kultur plakasindaki uc kuyudan, mevcut WST-1 verisi bir deneydeki dort kuyudan gelir.
RT-qPCR icin Delta Ct uzerinde Welch testi ve Holm duzeltmesi kullanilir.
Fura-2 ve WST-1 karsilastirmalari betimseldir; bu iki analizde cikarimsal p degeri verilmez.

DOSYA_LISTESI_SHA256.txt her dosyanin saglama toplamini verir.
Kod ve veri deposu: github.com/elmasnuryilmaz/tdp43-soce-manuscript
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
