# MANUSCRIPT_v4 — değişiklik kaydı (gönderime hazır sürüm)

**Tarih:** 18 Eylül 2026
**Kaynak:** `MANUSCRIPT_v3.md` → `MANUSCRIPT_v4_SUBMISSION.md` (+ `.docx`)
**Paket:** `09_YAYIN_PAKETI/` (şekiller, tablolar, ekler, kod, kaynak veri, QA raporu)
**Otomatik denetim:** `09_YAYIN_PAKETI/logs/qa_manuscript_vs_data.txt` — **94 kontrol, hepsi geçti**

v3'te "kapatılamadı" diye işaretlenen 13 maddenin **11'i kapatıldı**. Kalan ikisi harici disk
(`/Volumes/10TBElmas`) veya yazar kararı gerektiriyor; aşağıda §4'te.

---

## 1. Yeniden çalıştırılan analizler

### 1.1 Salmon toplamada %26'lık transkript kaybı bulundu ve düzeltildi

Salmon indeksi **kapsamlı** GENCODE v47 transkript kümesinden kurulmuş (384.354 transkript),
ama gen düzeyine toplama **basic** anotasyonun haritasıyla yapılmıştı; bu harita
transkriptlerin yalnızca 157.588'ini kapsıyor ve **TPM'in %26'sı** toplama sırasında düşüyordu.
Etkisi gen bazında büyük: ATP2A2 sayımı 2,39×, ORAI1 2,71×, CBARP 1,96×, STIM1 1,45× artıyor.

Yapılan: tam haritayla yeniden toplama (`code/deseq_full.R` + `gene_TPM_full.csv`,
`gene_counts_full.csv`) ve DESeq2'nin yeniden çalıştırılması. Bütün yönler korundu:

| Gen | v3 (basic harita) | v4 (tam harita) |
|---|---|---|
| STIM1 | +0,803 (5,0 × 10⁻¹⁹) | **+0,929 (4,3 × 10⁻⁴⁷)** |
| TRPC1 | +0,967 (5,7 × 10⁻⁹) | **+0,958 (1,3 × 10⁻¹²)** |
| ORAI1 | +0,364 (anlamsız) | **+0,433 (2,4 × 10⁻⁴, anlamlı)** |
| ATP2A3 | +1,202 (2,7 × 10⁻¹⁸) | **+1,306 (1,4 × 10⁻²⁹)** |
| CBARP | −1,746 (2,6 × 10⁻¹⁵) | **−1,254 (3,1 × 10⁻²⁵)** |

ORAI1'in anlamlı hâle gelmesi RT-qPCR ile uyumu **güçlendiriyor** (1,7 kat artış, p < 0,01).

### 1.2 Aile bolluğu TPM ile yeniden hesaplandı (rapor §13.1/1)

`Table4_transcript_family_abundance.csv` ve Figure 7 artık uzunluk düzeltilmiş TPM'e dayanıyor.
Nitel sonuçların hepsi ayakta: ORAI2 (%77) ve ATP2A2 (%98) baskın ve azalıyor; ORAI1 (%15) ve
ATP2A3 (%1,5) minör ve artıyor; SARAF SOCE düzenleyici havuzunun %82'si ve artıyor.
Aile net değişimleri: SERCA −%32, mitokondriyal −%43, STIM +%13, ORAI −%3, SOCE düzenleyicileri −%1.

### 1.3 Panel düzeyi NMD testleri dört koşullu tabloyla yeniden hesaplandı (rapor §13.1/3)

Sonuç neredeyse aynı çıktı: Tier1 p = 0,435; Tier2 0,905; Tier3 0,964; Tier4 0,911
(makalede 0,44–0,96). Ek olarak on altı pozitif kontrol geni **grup olarak da NMD-duyarlı değil**
(p = 0,62); bu cümle metne eklendi. v3'teki "eski tablodan geliyor" uyarısı kaldırıldı.
Çıktı: `S10b_NMD_panel_level_tests.csv`.

### 1.4 MS'te donör bağımlılığı modellendi (rapor §13.1/8)

`code/ms_donor_level.py`: her testin örnek düzeyi ve **donör düzeyi** sürümü.

| Karşılaştırma | örnek düzeyi | donör düzeyi |
|---|---|---|
| NAWM vs kontrol WM, TRPC1 | δ = −0,463 (p = 0,0076) | **δ = −0,771 (p = 0,030)**, 7 vs 5 donör |
| Miyelin+gliya düzeltilmiş | δ = −0,395 (p = 0,0034) | **δ = −0,640 (p = 0,055)** |
| Remiyelinize / inaktif lezyon | — | δ = −1,000 (p = 0,016) |

Metin artık donör düzeyi sayıları veriyor; düzeltilmiş analizin sınırda kaldığı açıkça yazılı.
Çıktı: `S16b_multiple_sclerosis_donor_level.csv`.

### 1.5 SUPPA2 sonucu metne eklendi

v3'te "yöntem uyumu için çalıştırıldı" deyip sonuç verilmiyordu. Tezden (Şekil E.13):
31.474 olay / 5.345 gen test edildi, **ortak eşiklerden geçen olay kalmadı**; eşleşen 9.369
ekzon atlama olayında Spearman ρ = 0,319, yön uyumu %62,4. §3.1'e eklendi.

## 2. Üretilen dosyalar

**Şekiller (9/9, PNG 300 dpi + PDF).** Numaralandırma ilk anılma sırasına göre; numara artık
görselin içine gömülü değil, dosya adında.
- **Figure 6 sıfırdan üretildi** (laboratuvar verisi): TARDBP susturma, dört hedef mRNA,
  Fura-2 ER salıverilişi/SOCE, WST-1 24 ve 48 saat. Bütün istatistikler ham veriden yeniden
  hesaplandı ve tezdeki değerlerle **birebir** çıktı (SOCE p = 0,0115; canlılık 61,5 ± 0,8).
- **Figure 7** TPM ile yeniden çizildi; başlıktaki "Absolute stoichiometry" ifadesi kaldırıldı.
- **Figure 5C**'deki hata düzeltildi: ölçülmüş sıfır artık `NaN` yapılıp gizlenmiyor, açık
  daire ve "measured zero" etiketiyle gösteriliyor.
- **Figure 8**'e hipokampus, oksipital korteks ve üç omurilik seviyesi eklendi; MS satırları
  donör düzeyine çevrildi.
- **Figure 5A** başlığı "Cryptic calls are TDP-43-specific" → "Positive-control recovery by
  comparison (permissive definition)" (gevşek tanımdan geldiği için).

**Tablolar (5/5, İngilizce başlıklarla).** Table 3 (on bir karşılaştırma) ve Table 5
(hastalıklar arası) **ilk kez üretildi**; Table 2'ye 1-tabanlı koordinat sütunu eklendi.

**Ek dosyalar (S1–S17).**
- **S1 sıfırdan üretildi**: sekiz sayfalık Excel — ham Ct'ler, replika düzeyi değerler,
  Fura-2 ve WST-1 ölçümleri, **primer dizileri**, termal profil, her panelin istatistiği.
- **S3**: 23 GB'lık rMATS çıktısından eşikleri geçen bütün olaylar, ham birleşim sayılarıyla
  (6,1 MB gzip) — böylece kapsam süzgeci dışarıdan yeniden üretilebilir.
- **S6b**'nin bozuk sütun başlıkları (`562_mRNA`, `H_SY5Y`…) yeniden üretildi.
- **S17** yeni: bütün erişim numaraları, tasarım ve örnek→grup eşlemesi.
- Geçersiz eski S11 (kriptik + NMD-duyarlı genler) **pakete alınmadı**.

**Kaynakça.** Beş kayıttan **48 kayda** çıkarıldı; hepsi Europe PMC üzerinden veya tezin
kaynakçasından doğrulandı. Giriş, yöntem ve tartışmadaki atıflar metne yerleştirildi
(TDP-43 patolojisi, kriptik hedefler, ASO çalışmaları, SOCE, SARAF, VAPB–PTPIP51 ve
kullanılan bütün yazılımlar).

**Declarations bölümü eklendi:** etik/biyogüvenlik beyanı (yeni insan/hayvan materyali yok,
yerleşik hücre hattı, kurumsal biyogüvenlik kuralları), finansman, yazar katkısı, çıkar
çatışması, teşekkür. Finansman ve katkı satırları yazar tarafından doldurulacak yer tutucu.

## 3. Metindeki diğer düzeltmeler

- §2.4: Salmon indeksinin kapsamlı küme olduğu ve %26'lık kaybın düzeltildiği yazıldı.
- §2.10: yöntem TPM'e göre yeniden yazıldı.
- §2.13: iTaq katalog numarası, 10 µL reaksiyon hacmi, erime eğrisi kontrolü, iki ayrı RNA
  setinin tarihleri ve primerlerin S1'de olduğu eklendi.
- §2.14: CPA katalog numarası eklendi.
- §3.7: WST-1'in **24. saatte %116,8'e çıktığı** (p = 2,4 × 10⁻⁴) eklendi — v3'te yalnızca
  48. saat vardı; iki zaman noktasının ters yönde olması Limitations'ta da belirtildi.
- §5: n = 3'ün tek plakadan geldiği ve ikinci susturma aracı + kurtarma deneyi gerektiği;
  APA düzeltmesinin yalnızca SH-SY5Y'ye uygulandığı.
- §2.17: erişilebilirlik metni paketin gerçek içeriğine göre yeniden yazıldı (S3, S17, 23 GB
  rMATS arşivi).
- Özet 465 → **404 kelime**.

## 4. Kapatılamayanlar

1. **iPSC-MN ve fare APA düzeltmesi** — BAM'ler harici diskte, disk bağlı değil.
   `code/run_bedcov4.sh` hazır; disk takılınca 10. adım çalıştırılmalı.
2. **Yazar kararı gerekenler:** finansman (BAP proje numarası), yazar katkısı, teşekkür,
   GitHub/Zenodo bağlantıları, hedef dergi ve özet sınırı.

## 5. Sonradan kapatılanlar (18 Eylül, yazar teyidiyle)

**CPA konsantrasyonu 10 µM olarak teyit edildi.** Prism'den çıkan temsili trase görselindeki
"10⁻³ M" etiketi hatalıdır. Görsel `code/relabel_traces.py` ile yeniden etiketlendi (İngilizce
eksen ve başlıklar, "CPA (10 µM)"); **yalnızca metin ek açıklamaları değiştirildi, traseler ve
eksenler olduğu gibi bırakıldı.** Orijinal dosya kanıt olarak
`source_data/representative_Fura2_traces_from_thesis.png` içinde duruyor.

**Figure 6 beş panele çıkarıldı:** A TARDBP susturma · B hedef mRNA'lar · C WST-1 (24 ve 48 saat)
· **D temsili Fura-2 traseleri** · E ER salıverilişi ve SOCE grup verisi. Makaledeki efsane
buna göre güncellendi; `.docx` yeniden üretildi; QA denetimi yeniden çalıştırıldı (94/94).
