import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

from BeratTrnMAFFT import MafftCore
from BeratTrnMAFFT.utils import fasta_oku
import numpy as np

# fasta dosyasini okuyorum
diziler = fasta_oku("ornek.fasta")

print("okunan diziler:")
for i in range(len(diziler)):
    print(i, "-->", diziler[i])

# mafft nesnesini olusturdum
mafft = MafftCore()

# iki dizi karsilastiriyorum
kaydirma, skor = mafft.fft_ile_karsilastir(diziler[0], diziler[1])
print("\ndizi0 ile dizi1 karsilastirmasi:")
print("kaydirma:", kaydirma, "  skor:", round(skor, 2))

# benzerlik matrisini hesapliyorum
matris = mafft.benzerlik_matrisi_olustur(diziler)
print("\nbenzerlik matrisi:")
print(np.round(matris, 2))

# hizalama sirasini buluyorum
print("\nassamali hizalama:")
siralama = mafft.asamali_hizalama(diziler)
print("hizalama sirasi:", siralama)
