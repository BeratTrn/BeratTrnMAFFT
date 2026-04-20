# BeratTrnMAFFT 🧬

[![PyPI version](https://img.shields.io/pypi/v/BeratTrnMAFFT.svg)](https://pypi.org/project/BeratTrnMAFFT/)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**BeratTrnMAFFT**, biyoinformatik çalışmaları için geliştirilmiş, Hızlı Fourier Dönüşümü (FFT) tabanlı bir Çoklu Dizi Hizalaması (Multiple Sequence Alignment - MSA) Python kütüphanesidir. Geleneksel dinamik programlama algoritmalarının yüksek zaman karmaşıklığını ($O(nmk)$) aşmak için sinyal işleme tekniklerini ve aşamalı hizalama (progressive alignment) stratejilerini kullanır.

Bu kütüphane, İstanbul Rumeli Üniversitesi Biyoinformatik dersi dönem projesi kapsamında geliştirilmiştir.

## ✨ Özellikler

* **Sinyal Dönüşümü (One-Hot Encoding):** DNA dizilerini (A, C, G, T) matematiksel analiz için 4 kanallı dijital sinyallere dönüştürür.
* **FFT Tabanlı Çapraz Korelasyon:** Diziler arasındaki homolog (benzer) bölgeleri bulmak için frekans uzayında hızlı çapraz korelasyon (cross-correlation) hesaplaması yapar.
* **Aşamalı Hizalama (Progressive Alignment):** FFT skorlarından elde edilen benzerlik matrisini kullanarak, dizileri "Kılavuz Ağaç" (Guide Tree) mantığıyla aşamalı olarak gruplar ve hizalar.
* **FASTA Desteği:** Biyolojik dizileri standart `.fasta` formatındaki dosyalardan kolayca okuyabilmeniz için dahili araçlar sunar.

## 📦 Kurulum

Kütüphaneyi PyPI üzerinden `pip` kullanarak kolayca sisteminize kurabilirsiniz:

```bash
pip install BeratTrnMAFFT