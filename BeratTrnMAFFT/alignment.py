"""
BeratTrnMAFFT - Hizalama Modülü

Bu modül şunları içerir:
  1. needleman_wunsch()    → İki dizi arasında global hizalama (DP algoritması)
  2. profil_olustur()      → Hizalanmış dizi grubundan frekans profili çıkar
  3. profil_sutun_skoru()  → İki profil sütununu karşılaştır
  4. profil_profil_hizala()→ İki grubu profil üzerinden NW ile hizala

Tüm bunlar MAFFT'ın "progressive alignment" (aşamalı hizalama) aşamasında kullanılır.
"""

import numpy as np

# Sabit skor parametreleri

ESLESME = 1        # Aynı nükleotid eşleştiğinde verilen puan
ESLESMEME = -1     # Farklı nükleotidler eşleştiğinde verilen ceza
BOSLUK = -2        # Her gap karakteri için uygulanan lineer ceza
# Not: Daha gerçekçi bir model için "affine gap" (gap açma + uzatma ayrımı)
# kullanılabilir, ama basitlik açısından lineer gap cezasıyla devam ediyoruz.

NUKLEOTIDLER = ['A', 'C', 'G', 'T']


# Yardımcı fonksiyon

def nukleotid_skoru(k1, k2):
    """
    İki karakter (nükleotid veya gap) arasındaki eşleşme skorunu döndürür.

    Parametreler:
        k1, k2 (str): Tek karakter; 'A','C','G','T' veya '-' (gap)

    Döndürür:
        int: ESLESME, ESLESMEME veya BOSLUK
    """
    if k1 == '-' or k2 == '-':
        return BOSLUK
    if k1 == k2:
        return ESLESME
    return ESLESMEME


# 1. Needleman-Wunsch Algoritması

def needleman_wunsch(dizi1, dizi2):
    """
    İki dizi arasında global hizalama yapar (Needleman-Wunsch algoritması, 1970).

    Algoritma nasıl çalışır?
    Dinamik Programlama (DP) prensibini kullanır. İki aşamadan oluşur:

    Aşama 1 — DP Matrisini Doldurma:
        dp[i][j], dizi1'in ilk i karakteriyle dizi2'nin ilk j karakterini
        hizalamanın en iyi skorudur. Her hücreyi 3 seçenekten en iyisiyle doldururuz:

        Seçenek A (köşegen):  dp[i-1][j-1] + score(dizi1[i-1], dizi2[j-1])
            -> İki karakteri eşleştir (match ya da mismatch)

        Seçenek B (yukarı):   dp[i-1][j] + BOSLUK
            -> dizi1'den bir karakter al, dizi2'ye gap ekle

        Seçenek C (sol):      dp[i][j-1] + BOSLUK
            -> dizi2'den bir karakter al, dizi1'e gap ekle

    Aşama 2 — Traceback (Geriye İz Sürme):
        dp[n][m]'den başlayarak hangi hücreye nereden gelindiğini takip ederek
        hizalanmış dizileri oluşturuz.

    Neden ayrı bir 'yön' matrisi kullanıyoruz?
        Floating-point karşılaştırmasında sayısal sapma olabilir.
        Örneğin: 1.9999999 == 2.0 → False! Bunu önlemek için DP doldurulurken
        hangi yönden geldiğimizi ayrı bir matrise ('D','U','L') yazıyoruz.
        Bu sayede traceback'te tekrar hesaplama yapmaya gerek kalmıyor.

    Parametreler:
        dizi1 (str): Birinci DNA dizisi
        dizi2 (str): İkinci DNA dizisi

    Döndürür:
        tuple: (hizali_dizi1, hizali_dizi2, skor)
    """
    n = len(dizi1)
    m = len(dizi2)

    # dp[i][j] = dizi1[:i] ile dizi2[:j]'nin en iyi hizalama skoru
    dp = np.zeros((n + 1, m + 1), dtype=float)

    # yön[i][j]: bu hücreye hangi yönden gelindi?
    # 'D' = diagonal (köşegen), 'U' = up (yukarı), 'L' = left (sol)
    yon = np.full((n + 1, m + 1), '', dtype=object)

    # İlk satır ve sütunu gap cezasıyla başlat
    # (Birini sıfırdan hizalamak için hep gap eklenmesi gerekir)
    for i in range(1, n + 1):
        dp[i][0] = i * BOSLUK
        yon[i][0] = 'U'
    for j in range(1, m + 1):
        dp[0][j] = j * BOSLUK
        yon[0][j] = 'L'

    # DP matrisini satır satır doldur
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            # 3 seçeneğin skorlarını hesapla
            kosegen = dp[i-1][j-1] + nukleotid_skoru(dizi1[i-1], dizi2[j-1])
            yukari  = dp[i-1][j]   + BOSLUK
            sol     = dp[i][j-1]   + BOSLUK

            # En iyi skoru al
            en_iyi = max(kosegen, yukari, sol)
            dp[i][j] = en_iyi

            # Yönü kaydet (eşitlik durumunda öncelik sırası: köşegen > yukarı > sol)
            if en_iyi == kosegen:
                yon[i][j] = 'D'
            elif en_iyi == yukari:
                yon[i][j] = 'U'
            else:
                yon[i][j] = 'L'

    # Traceback: matrisin sağ alt köşesinden başla
    hizali1 = []
    hizali2 = []

    i, j = n, m

    while i > 0 or j > 0:
        adim = yon[i][j]

        if adim == 'D':
            # Köşegenden: her iki dizi de bir karakter veriyor
            hizali1.append(dizi1[i-1])
            hizali2.append(dizi2[j-1])
            i -= 1
            j -= 1

        elif adim == 'U':
            # Yukarıdan: dizi1'den karakter al, dizi2'ye gap ekle
            hizali1.append(dizi1[i-1])
            hizali2.append('-')
            i -= 1

        else:  # 'L'
            # Soldan: dizi2'den karakter al, dizi1'e gap ekle
            hizali1.append('-')
            hizali2.append(dizi2[j-1])
            j -= 1

    # Geriye doğru oluşturduğumuz için listeyi ters çevir
    hizali1.reverse()
    hizali2.reverse()

    return ''.join(hizali1), ''.join(hizali2), float(dp[n][m])


# 2. Profil Matrisi Oluşturma

def profil_olustur(hizalanmis_diziler):
    """
    Hizalanmış bir dizi grubundan "profil matrisi" çıkarır.

    Profil nedir?
    Her hizalama sütunu için, o sütunda hangi karakterin kaç kez göründüğünü
    oransal olarak tutan bir tablodur. Örneğin:

        Hizalı diziler:  ACG-
                         ACT-
                         AGT-

        Profil (her sütun için frekans sözlüğü):
          Sütun 0: {'A':1.0, 'C':0.0, 'G':0.0, 'T':0.0, '-':0.0}
          Sütun 1: {'A':0.0, 'C':0.67, 'G':0.33, 'T':0.0, '-':0.0}
          Sütun 2: {'A':0.0, 'C':0.0, 'G':0.33, 'T':0.67, '-':0.0}
          Sütun 3: {'A':0.0, 'C':0.0, 'G':0.0, 'T':0.0, '-':1.0}

    Bu yapı, ilerleyen adımda iki grubu birleştirirken "profil-profil hizalaması"
    yapmamızı sağlar.

    Parametreler:
        hizalanmis_diziler (list of str): Eşit uzunluklu, hizalanmış dizi listesi

    Döndürür:
        list of dict: Uzunluk L liste; her eleman bir sütunun frekans sözlüğü
    """
    if not hizalanmis_diziler:
        return []

    uzunluk = len(hizalanmis_diziler[0])
    adet = len(hizalanmis_diziler)
    profil = []

    for pos in range(uzunluk):
        # Bu pozisyondaki her karakteri say
        sayac = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0}

        for dizi in hizalanmis_diziler:
            k = dizi[pos] if pos < len(dizi) else '-'
            if k in sayac:
                sayac[k] += 1
            # Bilinmeyen karakter (N, R vb.) gelirse gap gibi davran

        # Sayıları frekansa (0-1 arası) normalize et
        frekans = {k: sayi / adet for k, sayi in sayac.items()}
        profil.append(frekans)

    return profil


# 3. İki Profil Sütunu Arasındaki Skor
def profil_sutun_skoru(sutun1, sutun2):
    """
    İki profil sütunu arasındaki beklenen hizalama skorunu hesaplar.

    Yöntem:
    Her nükleotid çifti (b1, b2) için:
        olasılık = sutun1[b1] * sutun2[b2]   (bağımsız olay çarpımı)
        katkı    = olasılık * score(b1, b2)

    Tüm çiftlerin katkısını toplarsak "beklenen skor" çıkar. Bu yaklaşım
    Thompson et al. (1994)'ün ClustalW makalesindeki ağırlıklı profil
    skorlamasıyla benzerdir.

    Gap frekansı da dikkate alınır: her iki taraftaki gap oranı arttıkça
    ek bir gap cezası uygularız.

    Parametreler:
        sutun1, sutun2 (dict): {'A': frek, 'C': frek, 'G': frek, 'T': frek, '-': frek}

    Döndürür:
        float: İki sütun arasındaki beklenen skor
    """
    skor = 0.0

    # Nükleotid × nükleotid katkılarını topla
    for b1 in NUKLEOTIDLER:
        for b2 in NUKLEOTIDLER:
            olasilik = sutun1.get(b1, 0.0) * sutun2.get(b2, 0.0)
            if b1 == b2:
                skor += olasilik * ESLESME
            else:
                skor += olasilik * ESLESMEME

    # Gap cezasını ekle: her iki sütundaki gap frekansından küçük bir ceza
    # (tam BOSLUK değildir, çünkü profil ortalama bir temsil)
    g1 = sutun1.get('-', 0.0)
    g2 = sutun2.get('-', 0.0)
    # "Herhangi birinde gap varsa ceza" formülü: 1 - (gap olmama olasılığı)
    gap_orani = 1.0 - (1.0 - g1) * (1.0 - g2)
    skor -= gap_orani * abs(BOSLUK) * 0.5

    return skor


# 4. Profil-Profil Hizalaması

def profil_profil_hizala(profil1, profil2, diziler1, diziler2):
    """
    İki hizalanmış dizi grubunu profil-profil yöntemiyle birleştirir.

    Bu fonksiyon esasen Needleman-Wunsch algoritmasının aynısıdır;
    fark şu: karakter-karakter skor yerine PROFIL_SUTUN skoru kullanırız.

    Neden profil kullanıyoruz?
    Progressive alignment'ın ortasında artık elimizde tek tek diziler değil,
    birkaç diziyi kapsayan hizalama grupları var. Örneğin Grup 1'de Seq1 ve Seq2
    zaten hizalanmış; Grup 2'de Seq3 var. Grup 1'i Seq3 ile hizalamak için
    Grup 1'i bir "profil"e (sütun bazlı frekans tablosuna) indirgeyip
    Seq3'ün profiline karşı NW çalıştırıyoruz. Bu sayede önceki hizalamayı
    bozmadan yeni diziyi ekliyoruz.

    Parametreler:
        profil1 (list of dict): Grup 1'in profil matrisi
        profil2 (list of dict): Grup 2'nin profil matrisi
        diziler1 (list of str): Grup 1'deki hizalanmış diziler
        diziler2 (list of str): Grup 2'deki hizalanmış diziler

    Döndürür:
        list of str: Tüm dizileri kapsayan yeni hizalanmış dizi listesi
    """
    n = len(profil1)
    m = len(profil2)

    # DP matrisi ve yön matrisi
    dp  = np.zeros((n + 1, m + 1), dtype=float)
    yon = np.full((n + 1, m + 1), '', dtype=object)

    # Başlangıç değerleri
    for i in range(1, n + 1):
        dp[i][0] = i * BOSLUK
        yon[i][0] = 'U'
    for j in range(1, m + 1):
        dp[0][j] = j * BOSLUK
        yon[0][j] = 'L'

    # DP matrisini doldur (profil sütun skoru kullanarak)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            kosegen = dp[i-1][j-1] + profil_sutun_skoru(profil1[i-1], profil2[j-1])
            yukari  = dp[i-1][j]   + BOSLUK
            sol     = dp[i][j-1]   + BOSLUK

            en_iyi = max(kosegen, yukari, sol)
            dp[i][j] = en_iyi

            if en_iyi == kosegen:
                yon[i][j] = 'D'
            elif en_iyi == yukari:
                yon[i][j] = 'U'
            else:
                yon[i][j] = 'L'

    # Traceback: her adımı kaydet
    # 'M' = her iki gruptan karakter al
    # 'G2' = sadece grup 1'den karakter al → grup 2'ye gap ekle
    # 'G1' = sadece grup 2'den karakter al → grup 1'e gap ekle
    yol = []
    i, j = n, m

    while i > 0 or j > 0:
        adim = yon[i][j]
        if adim == 'D':
            yol.append(('M', i - 1, j - 1))
            i -= 1
            j -= 1
        elif adim == 'U':
            yol.append(('G2', i - 1, -1))
            i -= 1
        else:  # 'L'
            yol.append(('G1', -1, j - 1))
            j -= 1

    yol.reverse()

    # Yola göre her diziye karakter veya gap ekle
    # Başlangıçta her dizi için boş karakter listesi oluştur
    yeni1 = [[] for _ in diziler1]
    yeni2 = [[] for _ in diziler2]

    for tip, pos1, pos2 in yol:
        if tip == 'M':
            # Her iki gruptan da o pozisyondaki karakteri al
            for k, d in enumerate(diziler1):
                yeni1[k].append(d[pos1])
            for k, d in enumerate(diziler2):
                yeni2[k].append(d[pos2])

        elif tip == 'G2':
            # Grup 1'den karakter al; grup 2'deki tüm dizilere gap ekle
            for k, d in enumerate(diziler1):
                yeni1[k].append(d[pos1])
            for k in range(len(diziler2)):
                yeni2[k].append('-')

        else:  # 'G1'
            # Grup 2'den karakter al; grup 1'deki tüm dizilere gap ekle
            for k in range(len(diziler1)):
                yeni1[k].append('-')
            for k, d in enumerate(diziler2):
                yeni2[k].append(d[pos2])

    # Karakter listelerini string'e çevir ve grupları birleştir
    sonuc1 = [''.join(chars) for chars in yeni1]
    sonuc2 = [''.join(chars) for chars in yeni2]

    return sonuc1 + sonuc2
