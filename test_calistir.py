"""
BeratTrnMAFFT Test / Demo Scripti
==============================
Bu scripti çalıştırarak paketin doğru kurulup kurulmadığını
ve algoritmanın çalışıp çalışmadığını test edebilirsin.

Çalıştırma:
    python test_calistir.py
"""

import sys
import io

# Türkçe karakter sorunu yaşamamak için stdout'u UTF-8'e zorla
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from BeratTrnMAFFT import MAFFTHizalayici
from BeratTrnMAFFT.utils import fasta_oku, hizalamali_yazdir, kimlik_yuzdesi_hesapla


print("=" * 55)
print("  BeratTrnMAFFT  —  MAFFT Demo Testi")
print("  Öğrenci No: 221201018   (221201018 % 4 = 2 → MAFFT)")
print("=" * 55)


# ─── TEST 1: FASTA dosyasını okuma ─────────────────────────────
print("\n[TEST 1] FASTA dosyası okunuyor: ornek.fasta")

try:
    kayitlar = fasta_oku("ornek.fasta")
except FileNotFoundError:
    print("HATA: ornek.fasta bulunamadı. Scripti proje klasöründen çalıştır.")
    sys.exit(1)

isimler  = [k[0] for k in kayitlar]
diziler  = [k[1] for k in kayitlar]

print(f"  {len(diziler)} dizi okundu:")
for isim, dizi in zip(isimler, diziler):
    print(f"  {isim.ljust(20)} : {dizi}  (uzunluk: {len(dizi)})")


# ─── TEST 2: FFT çapraz korelasyon ─────────────────────────────
print("\n[TEST 2] İki dizi arasında FFT benzerlik skoru:")

hizalayici = MAFFTHizalayici()

for i in range(len(diziler)):
    for j in range(i + 1, len(diziler)):
        kaydirma, skor = hizalayici.fft_benzerlik_skoru(diziler[i], diziler[j])
        print(f"  {isimler[i]}  vs  {isimler[j]}")
        print(f"    → Kaydırma: {kaydirma:+d}  |  Skor: {skor:.2f}")


# ─── TEST 3: Mesafe matrisi ─────────────────────────────────────
print("\n[TEST 3] Mesafe matrisi:")

mesafe_mat = hizalayici.mesafe_matrisi_olustur(diziler)
print("  (küçük değer = diziler daha benzer)")

# Sütun başlıklarını yaz
baslik = "              " + "  ".join(f"{isim[:8]:>8}" for isim in isimler)
print(baslik)
for i in range(len(diziler)):
    satir = f"  {isimler[i][:12]:<12}  "
    satir += "  ".join(f"{mesafe_mat[i][j]:8.4f}" for j in range(len(diziler)))
    print(satir)


# ─── TEST 4: Basit iki dizi hizalaması (Needleman-Wunsch) ───────
print("\n[TEST 4] İki dizi arası Needleman-Wunsch hizalaması:")

from BeratTrnMAFFT.alignment import needleman_wunsch
h1, h2, skor = needleman_wunsch(diziler[0], diziler[1])
print(f"  {isimler[0]}: {h1}")
print(f"  {isimler[1]}: {h2}")
print(f"  NW Skoru : {skor:.2f}")
print(f"  Kimlik   : %{kimlik_yuzdesi_hesapla(h1, h2):.1f}")


# ─── TEST 5: Çoklu dizi hizalaması (MAFFT tam pipeline) ─────────
print("\n[TEST 5] MAFFT çoklu dizi hizalaması:")

son_isimler, hizalanmis = hizalayici.coklu_hizala(diziler, isimler)

# Sonuçları yazdır
hizalamali_yazdir(son_isimler, hizalanmis)

# Tüm çift kimlik yüzdelerini hesapla
print("  Çift Kimlik Yüzdeleri (Pairwise Identity):")
for i in range(len(hizalanmis)):
    for j in range(i + 1, len(hizalanmis)):
        yuzde = kimlik_yuzdesi_hesapla(hizalanmis[i], hizalanmis[j])
        print(f"    {son_isimler[i]}  vs  {son_isimler[j]} : %{yuzde:.1f}")


print("\n[TAMAMLANDI] Tüm testler başarıyla çalıştı.")
print("Paket kurulumu için: pip install BeratTrnMAFFT")


# ─── TEST 6: Grafik / Görselleştirme ────────────────────────────────
print("\n[TEST 6] Görselleştirme grafikleri oluşturuluyor...")

try:
    from BeratTrnMAFFT.visualization import tum_grafikleri_kaydet

    grafik_dosyasi = tum_grafikleri_kaydet(
        mesafe_mat=mesafe_mat,
        isimler=son_isimler,
        hizalanmis_diziler=hizalanmis,
        cikti_dosyasi="mafft_sonuclari.png"
    )
    print(f"  Grafik kaydedildi: {grafik_dosyasi}")
    print("  İçerik: Isı haritası | Dendrogram | Hizalama ızgarası | Kimlik bar grafiği")
except ImportError as e:
    print(f"  Görselleştirme için 'pip install matplotlib scipy seaborn' gerekli: {e}")
