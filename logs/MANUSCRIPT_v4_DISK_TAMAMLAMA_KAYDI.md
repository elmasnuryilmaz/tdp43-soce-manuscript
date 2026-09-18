# Harici disk bağlandıktan sonra tamamlananlar

**Tarih:** 18 Eylül 2026 · **Disk:** `/Volumes/10TBElmas` bağlı
**Etkilenen dosyalar:** `MANUSCRIPT_v4_SUBMISSION.md/.docx`, `09_YAYIN_PAKETI/{supplementary,code}`
**Otomatik denetim:** 94 → **120 kontrol, hepsi geçti** (`09_YAYIN_PAKETI/logs/qa_manuscript_vs_data.txt`)

---

## 1. iPSC-türevi motor nöronlarda APA düzeltmesi çalıştırıldı (v4'ün 1 numaralı açık maddesi)

`samtools bedcov -Q 30 **-j**` ile 4 BAM üzerinde yeniden hesaplandı
(`code/run_bedcov_v3_remaining.sh` → `code/recalc_apa_remaining_R.R`).
Derinlik süzgecini **59 birim / 29 gen** geçti.

**En önemli sonuç — ölçümün pozitif kontrolü bu modelde çalışıyor:**

| Birim | kontrol | KD | Δ | aralık |
|---|---|---|---|---|
| **STMN2 intron 2** (kriptik poliadenilasyon) | 0,570 | 0,819 | **+0,249** | +0,208 … +0,290 |
| UNC13A intron 31 | 0,691 | 0,545 | −0,145 | −0,199 … −0,092 |

SH-SY5Y'de `-j` düzeltmesinden sonra STMN2 pozitif kontrolü zayıflamış ve derinlik süzgecinin
altında kalmıştı; bu yüzden oradaki negatif sonuç "duyarlılığı gösterilememiş bir ölçüm"
üzerine kuruluydu. Artık **hastalıkla ilgili insan modelinde ölçümün çalıştığı gösterildi** ve
o modelde de çekirdek SOCE genleri az hareket ediyor:

- en büyük: ATP2A2 intron 3 (−0,164), MICU3 intron 8 (−0,260)
- 0,10'un altında: SARAF intron 5 (+0,097), TRPC1 intron 1 (+0,083), STIM2 terminal ekzon (−0,074)
- SH-SY5Y'nin iki adayından STIM2 intron 13 tekrarlanmadı (−0,024, aralık sıfırı içeriyor);
  STIM1 intron 17 bu veri setinde derinlik süzgecini geçmedi.

**Sınır:** grup başına iki replika olduğu için tam sayım bootstrap'ı yalnızca 16 çekiliş
veriyor; aralıklar kabadır ve metinde bu açıkça yazılıdır. Bilgi taşıyan şey nokta tahminidir.

Metinde güncellenenler: §3.5 (yeni paragraf), §4 "RNA-processing analyses are informative but
not exclusionary", §5 Limitations. Ek dosya **S11 artık dört veri setini birden içeriyor** (290 satır).

## 2. Sürüm ve parametre doğrulamaları (rapor §13.1/12–13)

Disk açılınca hattın kendi kayıtlarından doğrulandı:

| İddia | Kaynak | Sonuç |
|---|---|---|
| Cutadapt v4.6 | `trimmed/*_trimming_report.txt` | **YANLIŞ → v5.2** (tezdeki 4.6 hatalı). Metin düzeltildi; adaptörün otomatik saptanan Illumina TruSeq olduğu da eklendi. |
| LeafCutter v0.2.9 | `tez_duzeltmeler/tools/leafcutter_repo/leafcutter/DESCRIPTION` | **DOĞRU** (Version: 0.2.9) |
| SUPPA2 v2.3 | `TEZ/output/extended_computational_layers_2026-05-16/suppa2/SUPPA-2.3/suppa.py` | **DOĞRU** (conda'da değil, kaynak kopyası) |
| HISAT2 2.2.2, `--dta`, bilinen kırpılma bölgesi dosyası | BAM `@PG` satırı | **DOĞRU** |
| SAMtools 1.23 (sıralama) | BAM `@PG` satırı | **DOĞRU** |
| featureCounts v2.1.1 | `analysis_results/featurecounts_countReadPairs_2026-06-27/.../*.log` | **DOĞRU**; parametreler de okundu ve §2.2'ye eklendi: `-s 0 -g gene_name`, çift uçlu için `-p --countReadPairs`, çoklu eşleşen ve çoklu örtüşen okumalar sayılmıyor, meta-özellik düzeyi |

## 3. §3.5'te kesinlik düzeltmesi

"No robust index change was found for ORAI1–3, TRPC1, SARAF, STIMATE, CBARP or the SERCA genes"
ifadesi, SARAF terminal ekzonu (−0,036) ve ORAI2 terminal ekzonu (+0,012) aralıkları sıfırı
dışladığı için **yanıltıcıydı**. Yeni ifade: "0,05 veya daha büyük hiçbir indeks değişimi
bulunmadı; sıfırı dışlayan diğer iki aralık ihmal edilebilir büyüklükte" + değerler verildi.

## 4. Fare veri setleri de tamamlandı

C2C12 (74 birim) ve NSC34 (131 birim) `-j` ile yeniden hesaplandı; böylece **hizalaması olan
dört karşılaştırmanın dördü de** düzeltilmiş APA analizinden geçmiş oldu.

- **Pozitif kontrol yalnızca motor nöron benzeri hatta bilgilendirici:** NSC34'te Stmn2 intron 2
  beklenen yönde (+0,148; aralık +0,005…+0,290). Myoblast hattı C2C12'de o intron derinlik
  süzgecini geçmiyor.
- Çekirdek SOCE'de |Δ| = 0,30'u aşan birim yok. En büyük ikisi C2C12'de (Atp2a2 intron 6 +0,285;
  Trpc1 intron 7 +0,260) ve ikisinin de aralığı sıfırı içeriyor, diğer hatta tekrarlanmıyor.
- **İki nöronal modelde birden pozitif olan tek birim: SARAF intron 5** (iPSC-MN +0,097;
  NSC34 +0,204; SH-SY5Y'de ölçülemiyor). Aday olarak kaydedildi — NSC34 aralığı geniş,
  myoblast hattında yok ve kapsam eğimi poli(A) bölgesi demek değil.

§2.7, §3.5 ve Limitations buna göre güncellendi; S11 dört veri setini birden içeriyor (290 satır).
Denetim **120 kontrole** çıktı, hepsi geçiyor.

## 5. Değişmeyen açık maddeler

BAP proje numarası, yazar katkısı/teşekkür, GitHub/Zenodo bağlantıları — yazar dolduracak.
