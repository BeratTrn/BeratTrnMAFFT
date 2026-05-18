"""
PyPI Paket Kurulum Dosyası
pip install BeratTrnMAFFT   komutuyla kurulabilmesi için bu dosya gerekli.
"""

from setuptools import setup, find_packages
import os

here = os.path.abspath(os.path.dirname(__file__))
readme_yolu = os.path.join(here, "README.md")

if os.path.exists(readme_yolu):
    with open(readme_yolu, "r", encoding="utf-8") as f:
        uzun_aciklama = f.read()
else:
    uzun_aciklama = "MAFFT tabanlı çoklu dizi hizalama (MSA) paketi."

setup(
    # ── Paket Kimliği ────────────────────────────────────────────
    name="BeratTrnMAFFT",
    version="1.0.0",

    # ── Yazar Bilgileri ──────────────────────────────────────────
    author="Berat Turan",
    author_email="turanberatr@gmail.com",

    # ── Açıklama ─────────────────────────────────────────────────
    description=(
        "MAFFT (Multiple Alignment using Fast Fourier Transform) "
        "algoritmasının Python implementasyonu. "
        "Biyoinformatik dersi dönem projesi — Öğrenci No: 221201018"
    ),
    long_description=uzun_aciklama,
    long_description_content_type="text/markdown",

    # ── URL ──────────────────────────────────────────────────────
    url="https://github.com/BeratTuran/BeratTrnMAFFT", 

    # ── Paket Bulma ──────────────────────────────────────────────
    # find_packages() BeratTrnMAFFT/ dizinini otomatik bulur
    packages=find_packages(),

    # ── Bağımlılıklar ────────────────────────────────────────────
    # Sadece numpy gerekiyor; standart kütüphane dışında başka bağımlılık yok
    install_requires=[
        "numpy>=1.20",
    ],

    # ── Python Sürümü ────────────────────────────────────────────
    python_requires=">=3.7",

    # ── PyPI Kategorileri (Classifier) ───────────────────────────
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Education",
    ],

    # ── Anahtar Kelimeler ────────────────────────────────────────
    keywords="bioinformatics msa mafft sequence alignment fft",
)
