"""
BeratTrnMAFFT - UPGMA Kılavuz Ağaç Modülü

Bu modül UPGMA (Unweighted Pair Group Method with Arithmetic mean)
algoritmasını uygular. UPGMA bir hiyerarşik kümeleme yöntemidir ve
çoklu dizi hizalamasında hangi dizilerin hangi sırayla birleştirileceğini
(yani "kılavuz ağacı") belirlemek için kullanılır.
"""


def upgma(mesafe_matrisi, isimler):
    """
    UPGMA algoritması ile hiyerarşik kümeleme yapar ve
    birleştirme sırasını (kılavuz ağacı) döndürür.

    Algoritma nasıl çalışır?
    1. Başlangıçta her dizi kendi başına bir küme oluşturur.
       Örn: N=4 dizi → 4 küme: {0}, {1}, {2}, {3}

    2. Mesafe matrisine bakarak en yakın iki kümeyi bul.
       "En yakın" demek: iki kümedeki tüm dizi çiftleri arasındaki
       mesafelerin ORTALAMASI en küçük olan çift.

    3. Bu iki kümeyi birleştir → yeni bir küme oluşur.
       Birleştirme adımını kaydet (progressive alignment'ta lazım olacak).

    4. Toplam N-1 birleştirme işleminden sonra tek bir kök kümede bitecek.

    UPGMA vs. Neighbor-Joining:
        UPGMA daha basit ama "moleküler saat" varsayımı yapar (evrim hızı
        tüm kollarda eşit). NJ daha gerçekçidir ama daha karmaşık. Biz
        öğrenci projesi olduğu için UPGMA kullanıyoruz.

    Parametreler:
        mesafe_matrisi (np.ndarray): NxN simetrik mesafe matrisi
                                    (köşegen = 0, küçük değer = yakın)
        isimler (list): Dizi isimleri (sadece hata mesajlarında kullanılır)

    Döndürür:
        list of tuple: Her eleman (grup1_indeksler, grup2_indeksler) şeklinde.
                       indeksler = orijinal dizi pozisyonları (0-tabanlı)
                       Örn: [([0], [1]), ([0, 1], [2])]
                       -> Önce 0. ve 1. diziler birleşiyor,
                       -> sonra bu grup ile 2. dizi birleşiyor.
    """
    n = len(isimler)

    if n < 2:
        # Tek dizi varsa birleştirme gerekmez
        return []

    # Her kümenin içerdiği orijinal dizi indekslerini tut.
    # Başlangıçta her dizi kendi başına bir küme.
    kumeler = [[i] for i in range(n)]

    birlestirme_sirasi = []

    # N-1 adımda tek bir kümeye indirge
    for _ in range(n - 1):
        aktif_kume_sayisi = len(kumeler)

        if aktif_kume_sayisi < 2:
            break  # Sadece 1 küme kaldı, bitti

        # En yakın iki kümeyi bul (en küçük ortalama mesafe)
        en_kucuk_mesafe = float('inf')
        birinci_idx = 0
        ikinci_idx  = 1

        for i in range(aktif_kume_sayisi):
            for j in range(i + 1, aktif_kume_sayisi):

                # Küme i ile küme j arasındaki ortalama mesafeyi hesapla.
                # Bunu yapmak için küme i'deki her a ve küme j'deki her b için
                # mesafe_matrisi[a][b]'yi toplayıp dizi sayısına bölüyoruz.
                toplam = 0.0
                karsılasma = 0

                for a in kumeler[i]:
                    for b in kumeler[j]:
                        toplam += mesafe_matrisi[a][b]
                        karsılasma += 1

                if karsılasma == 0:
                    ort = float('inf')
                else:
                    ort = toplam / karsılasma

                if ort < en_kucuk_mesafe:
                    en_kucuk_mesafe = ort
                    birinci_idx = i
                    ikinci_idx  = j

        # Bu adımda birleşen grupları kaydet
        birlestirme_sirasi.append(
            (kumeler[birinci_idx].copy(), kumeler[ikinci_idx].copy())
        )

        # İki kümeyi birleştir: yeni küme her iki kümenin üyelerini içerir
        yeni_kume = kumeler[birinci_idx] + kumeler[ikinci_idx]

        # Önce büyük indeksi sil; yoksa küçük indeks kayar ve yanlış silinir!
        # (İlk sürümde bu sırayı karıştırmıştım, yanlış kümeler siliniyordu.)
        kumeler.pop(ikinci_idx)   # büyük indeksi önce sil
        kumeler.pop(birinci_idx)  # sonra küçük indeksi sil
        kumeler.append(yeni_kume) # yeni birleşik kümeyi ekle

    return birlestirme_sirasi
