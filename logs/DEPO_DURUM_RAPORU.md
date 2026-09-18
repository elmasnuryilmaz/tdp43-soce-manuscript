# GitHub deposu — durum ve eksikler

**Tarih:** 18 Eylül 2026
**Depo:** `https://github.com/elmasnuryilmaz/tdp43-thesis-reproducibility`
**Yerel kopya:** `~/Desktop/TEZ/tdp43-thesis-reproducibility` (641 MB; `.git` 259 MB)

## 1. Depo güncel mi?

| Kontrol | Sonuç |
|---|---|
| Yerel `main` ↔ `origin/main` | **Eşit** — ileride 0, geride 0 commit |
| Son commit | `da0f170`, **9 Temmuz 2026**, "Add reusable workflow engine layer" |
| İçerik | 17 tez analiz klasörü, `workflows/` (Snakemake + Nextflow), `docs/`, `MANIFEST.sha256` |
| Çalışma ağacı | **Temiz değil:** 11 dosya yerelde silinmiş, commit edilmemiş |

**11 silinmiş dosya** (hepsi büyük, toplam ~135 MB; yer açmak için silinmiş görünüyor):
`salmon_isoform_tpm_raw.csv` (20,7 MB), `salmon_isoform_counts_raw.csv` (18,2 MB),
`apa_pas_switch_sites_all.csv` (16,7 MB), `TRPC1_splicing_loop_black.gif` (15,7 MB),
`drimseq_isoform_proportions.csv` (11,4 MB), `figure_4_4_compact_volcano_all_datasets.pdf`
(10,7 MB), `psi_coregulation_edges.csv` ×2 (9,9 MB), `Mouse_PE_C2C12_vs_NSC34_significant_DEG_counts…`
(9,5 MB), iki `DESeq2_all_genes_recomputed.csv` (6,1 MB).

Bunlar **GitHub'da duruyor**, kayıp yok. Ama dikkat: ileride `git commit -a` yapılırsa silinmeleri
uzaktaki depoya da işlenir. Geri almak için (indirme gerektirmez, nesneler yerelde):

```bash
cd ~/Desktop/TEZ/tdp43-thesis-reproducibility && git restore .
```

## 2. Makale çalışmasının hiçbiri depoda yok

Deponun son commit'i 9 Temmuz; makale çalışmasının tamamı Eylül'de `~/Desktop/MAKALE` altında
yapıldı. Depoda şu anahtar kelimelerin **hiçbiri geçmiyor**: `bedcov`, `recalc_nmd`, `lib_junc`,
`A13_kriptik`, `ms_donor`, `fig_v3`, `YUKSEK_GUVEN`, `recount3`.

Eksik olanlar:

| Eksik | Nerede |
|---|---|
| Sağlamlık analizleri (TRPC1 adli incelemesi, kapsam süzgeci, STIM2.1, güç, eşlenmiş permütasyon) | `04_KOD/` |
| Birleşim/kriptik analiz, FRASER, APA, NMD, NYGC recount3, boş test | `07_DISK_ANALIZLERI/kod/` |
| MS kohort analizi | `06_MS_ANALIZI/` |
| v4 düzeltmeleri: tam haritayla DESeq2, TPM aile bolluğu, donör düzeyi MS, dört veri setinde APA, QA denetimi | `04_KOD/v3/` |
| Makale, şekiller, tablolar, S1–S17, kaynak veri | `09_YAYIN_PAKETI/` |
| Tekrarlanabilirlik raporu ve denetim belgeleri | `00_TEKRARLANABILIRLIK_RAPORU.md`, `08_HAKEM_DEGERLENDIRMESI/` |

## 3. Depodaki bazı sonuçlar makale tarafından geçersiz kılındı

Bir hakem makaledeki bağlantıyı takip edip depoya girerse, makaleyle **çelişen** sayılarla
karşılaşır. Yayımlanan her şeyde bunun açıkça işaretlenmesi gerekir:

| Depodaki | Makaledeki |
|---|---|
| `-j` olmadan hesaplanmış APA kapsamı (STMN2 +0,489) | `-j` ile +0,145; eski değer geçersiz |
| Sekiz kontrastlı NMD testi (CBARP q = 0,0049) | Dört koşullu: q = 0,145, hiçbir gen FDR'yi geçmiyor |
| Basic haritayla toplanmış Salmon sayımları (CPM stokiyometrisi) | Tam harita + TPM; ORAI1 artık anlamlı (q = 2,4 × 10⁻⁴) |
| TRPC1 ekzon atlama olayı bulgu olarak | Makalede geri çekildi, yöntemsel örnek olarak sunuluyor |
| Örnek düzeyi MS testleri | Donör düzeyi (NAWM p = 0,030; miyelin düzeltmesi p = 0,055) |

## 4. Yayına hazır paket

`09_YAYIN_PAKETI/` artık doğrudan depo kökü olarak kullanılabilecek durumda (21 MB, `qa/` hariç):

```
README.md              paketin haritası, üretim komutları, kalan açık maddeler
DATA_AVAILABILITY.md   bütün erişim numaraları, çalışma anında indirilen kaynaklar, paylaşılmayanlar
environment.yml        doğrulanmış araç sürümleri (conda)
requirements.txt       Python paketleri
.gitignore             qa/ ve geçici dosyalar hariç tutuluyor
manuscript/ figures/ tables/ supplementary/ code/ source_data/ logs/
```

GitHub sınırları açısından sorun yok: 40 MB'ı aşan tek dosya yok, toplam 21 MB.

## 5. Karar ve sonraki adım

İki seçenek var:

**(a) Makale için ayrı depo** (önerilen) — `tdp43-soce-manuscript` gibi. Makalenin bağlantısı
yalnızca makaleye ait materyale gider; tez deposu olduğu gibi kalır ve README'lerden
karşılıklı bağlantı verilir. §3'teki çelişki riski en aza iner.

**(b) Tez deposuna alt klasör** — `manuscript_2026/` olarak eklenir. Tek yerde toplanır ama
okuyucu 17 tez klasörünün arasına düşer; geçersiz kılınan sonuçlar için depo kökündeki
README'ye büyük bir uyarı notu koymak şart olur.

Her iki durumda da **push işlemi dışa dönüktür ve sizin onayınızla yapılır.** Onay verirseniz:
depoyu hazırlar, commit mesajını yazar, `git push` öncesi ne gideceğini size gösteririm.
Makalenin 2.17 bölümündeki `[GitHub repository — to be created]` yer tutucusu da o adreste
güncellenecek.
