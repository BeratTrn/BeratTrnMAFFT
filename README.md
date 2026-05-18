# BeratTrnMAFFT 🧬

[![PyPI version](https://img.shields.io/pypi/v/BeratTrnMAFFT.svg)](https://pypi.org/project/BeratTrnMAFFT/)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**BeratTrnMAFFT**, biyoinformatik çalışmaları için geliştirilmiş, Hızlı Fourier Dönüşümü (FFT) tabanlı bir Çoklu Dizi Hizalaması (Multiple Sequence Alignment - MSA) Python kütüphanesidir. Geleneksel dinamik programlama algoritmalarının yüksek zaman karmaşıklığını aşmak için sinyal işleme tekniklerini ve aşamalı hizalama (progressive alignment) stratejilerini kullanır.

Bu kütüphane, İstanbul Rumeli Üniversitesi Biyoinformatik dersi dönem projesi kapsamında geliştirilmiştir.

---

## ✨ Özellikler

* **Sinyal Dönüşümü (One-Hot Encoding):** DNA dizilerini (A, C, G, T) matematiksel analiz için 4 kanallı dijital sinyallere dönüştürür.
* **FFT Tabanlı Çapraz Korelasyon:** Diziler arasındaki homolog (benzer) bölgeleri bulmak için frekans uzayında hızlı çapraz korelasyon (cross-correlation) hesaplaması yapar.
* **UPGMA Kılavuz Ağaç:** FFT mesafe matrisinden hiyerarşik kümeleme ile kılavuz ağaç (guide tree) oluşturur; hangi dizilerin önce hizalanacağını belirler.
* **Aşamalı Hizalama (Progressive Alignment):** Kılavuz ağaca göre Needleman-Wunsch ve profil-profil hizalaması ile dizileri gruplar ve birleştirir.
* **FASTA Desteği:** Biyolojik dizileri standart `.fasta` formatındaki dosyalardan kolayca okuyabilmeniz için dahili araçlar sunar.

---

## 📦 Kurulum

Kütüphaneyi PyPI üzerinden pip kullanarak kolayca sisteminize kurabilirsiniz:

```bash
pip install BeratTrnMAFFT
```

---

## 🚀 Kullanım

```python
from BeratTrnMAFFT import MAFFTHizalayici
from BeratTrnMAFFT.utils import fasta_oku, hizalamali_yazdir

# FASTA dosyasını oku
veriler = fasta_oku("ornek.fasta")
isimler = [v[0] for v in veriler]
diziler = [v[1] for v in veriler]

# Hizala
hizalayici = MAFFTHizalayici()
son_isimler, hizalanmis = hizalayici.coklu_hizala(diziler, isimler)

# Sonucu yazdır
hizalamali_yazdir(son_isimler, hizalanmis)
```

---

## 🔬 Algoritma

| Aşama | Yöntem | Açıklama |
|-------|--------|----------|
| 1 | FFT Korelasyon | One-hot encoding + çapraz korelasyon ile benzerlik skoru |
| 2 | UPGMA | Mesafe matrisinden hiyerarşik kılavuz ağaç |
| 3 | Progressive Alignment | Needleman-Wunsch + Profil-Profil hizalaması |

---

## 📁 Paket Yapısı

```
BeratTrnMAFFT/
├── __init__.py        # Paket girişi
├── core.py            # Ana MAFFT motoru (3 aşama)
├── alignment.py       # Needleman-Wunsch + Profil hizalaması
├── guide_tree.py      # UPGMA kılavuz ağaç
├── utils.py           # FASTA okuma, yazdırma
└── visualization.py   # Grafik ve görselleştirme
```

---

## 📋 Gereksinimler

- Python 3.7+
- numpy >= 1.20
- matplotlib
- seaborn
- scipy
