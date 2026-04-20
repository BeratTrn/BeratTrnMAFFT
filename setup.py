from setuptools import setup, find_packages
import os

# README dosyasının içeriğini okutuyoruz
with open("README.md", "r", encoding="utf-8") as fh:
    uzun_aciklama = fh.read()

setup(
    name="BeratTrnMAFFT", 
    version="0.1.1", # VERSİYONU 0.1.1 YAPTIK!
    author="Berat Turan",
    description="Biyoinformatik dersi MAFFT algoritması dönem projesi",
    long_description=uzun_aciklama, # README'deki yazıları buraya aktarıyor
    long_description_content_type="text/markdown", # Yazının formatını belirtiyor
    packages=find_packages(),
    install_requires=[
        "numpy", 
    ],
    python_requires='>=3.6',
)