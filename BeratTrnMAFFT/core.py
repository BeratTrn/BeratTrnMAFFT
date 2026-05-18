"""
BeratTrnMAFFT - Ana MAFFT Motoru
=============================
Bu modül MAFFT (Multiple Alignment using Fast Fourier Transform) algoritmasının
üç ana aşamasını bir arada yönetir:

  Aşama 1 — FFT ile Benzerlik Skoru
      Her dizi, nükleotid bazlı sinyal vektörlerine (one-hot encoding) çevrilir.
      İki dizi arasındaki çapraz korelasyon FFT ile hesaplanır. Bu, iki sinyalin
      birbirine ne kadar "benzediğini" hızlı bulmak için ses işlemeden ödünç alınan
      bir tekniktir.

  Aşama 2 — Mesafe Matrisi ve UPGMA Kılavuz Ağaç
      FFT skorları mesafeye çevrilir (yüksek benzerlik = düşük mesafe).
      UPGMA algoritması bu mesafelere bakarak hangi dizilerin önce birleşeceğini
      belirleyen bir kılavuz ağaç üretir.

  Aşama 3 — Progressive Alignment
      Kılavuz ağacın söylediği sırayla diziler ikişer ikişer hizalanır.
      İki tekil dizi → Needleman-Wunsch.
      İki grup → Profil-Profil hizalaması.
      Her adımda bir önceki hizalama bozulmaz; sadece gap eklenerek genişler.
"""

import numpy as np
from .alignment import needleman_wunsch, profil_olustur, profil_profil_hizala
from .guide_tree import upgma

# One-hot encoding için nükleotid listesi
NUKLEOTIDLER = ['A', 'C', 'G', 'T']


class MAFFTHizalayici:
    """
    MAFFT algoritmasını uygulayan ana sınıf.

    Kullanım:
        from BeratTrnMAFFT import MAFFTHizalayici
        hizalayici = MAFFTHizalayici()
        isimler, hizalanmis = hizalayici.coklu_hizala(diziler, isimler)
    """

    def __init__(self):
        """
        Şu an herhangi bir parametre almıyor; skor parametreleri
        alignment.py'deki sabitlerden gelir. İleride geliştirme için
        buraya gap cezası gibi parametreler eklenebilir.
        """
        pass


    # Aşama 1 Yardımcıları: One-Hot Encoding + FFT Korelasyon

    def dizi_sinyal_vektoru(self, dizi):
        """
        Bir DNA dizisini dört ayrı sinyal vektörüne (one-hot) çevirir.

        One-hot encoding nedir?
        ───────────────────────
        Her pozisyon için dört nükleotidin (A, C, G, T) birine karşılık gelen
        bir "0/1 dizisi" oluştururuz. Bunu yaparken her nükleotid için ayrı
        bir vektör tutarız:

            Dizi: "ACGT"
              A sinyal: [1, 0, 0, 0] 
              C sinyal: [0, 1, 0, 0] 
              G sinyal: [0, 0, 1, 0] 
              T sinyal: [0, 0, 0, 1] 

        Bu temsil sayesinde her nükleotid kendi "kanalında" ayrı ayrı ele alınır.
        Harflere doğrudan FFT uygulanamayacağı için bu adım zorunludur.

        Parametreler:
            dizi (str): DNA dizisi (büyük harf, sadece A/C/G/T)

        Döndürür:
            dict: {'A': np.array, 'C': np.array, 'G': np.array, 'T': np.array}
        """
        sinyaller = {}

        for nk in NUKLEOTIDLER:
            # Sıfırlarla dolu bir vektör başlat
            v = np.zeros(len(dizi))

            # Dizide bu nükleotid geçen her pozisyona 1 yaz
            for i, karakter in enumerate(dizi):
                if karakter == nk:
                    v[i] = 1.0

            sinyaller[nk] = v

        return sinyaller

    def fft_benzerlik_skoru(self, dizi1, dizi2):
        """
        İki dizi arasındaki benzerliği FFT çapraz korelasyonu ile hesaplar.

        Çapraz korelasyon (cross-correlation) nedir?
        İki sinyali üst üste koyup birini kaydırdığımızda en iyi örtüşmeyi bulmak
        istediğimizi düşün. Her kaydırma değeri için "ne kadar örtüşüyor" diye
        bir skor hesaplayabiliriz. Bu işleme korelasyon denir.

        Naif yöntemle O(N²) zaman alır. FFT bunu O(N log N)'e indirir çünkü
        frekans uzayında korelasyon sadece bir çarpım işlemine karşılık gelir:

            correlation(x, y)[k] = IFFT(FFT(x) · conj(FFT(y)))[k]

        Burada 'conj' kompleks eşlenik (conjugate) anlamındadır.

        Her nükleotid kanalı (A, C, G, T) için ayrı ayrı korelasyon hesaplayıp
        hepsini topluyoruz. Böylece tüm nükleotid benzerliklerini birden yakalıyoruz.

        Döndürülen max_skor: iki dizi ne kadar benziyorsa o kadar büyük.
        en_iyi_kaydirma: hangi offset'te en iyi örtüşme var (hizalama için fikir verir).

        Parametreler:
            dizi1 (str): Birinci DNA dizisi
            dizi2 (str): İkinci DNA dizisi

        Döndürür:
            tuple: (en_iyi_kaydirma: int, max_skor: float)
        """
        sig1 = self.dizi_sinyal_vektoru(dizi1)
        sig2 = self.dizi_sinyal_vektoru(dizi2)

        # Zero-padding boyutu: iki uzunluğun toplamı.
        # Bu boyutu kullanmazsak dairesel korelasyon (circular correlation) yaparız,
        # yani sonun başa sarması söz konusu olur. Toplamı kullanarak bunu önleriz.
        pad = len(dizi1) + len(dizi2)

        toplam_korelasyon = np.zeros(pad)

        for nk in NUKLEOTIDLER:
            # Zaman → frekans uzayı
            F1 = np.fft.fft(sig1[nk], n=pad)
            F2 = np.fft.fft(sig2[nk], n=pad)

            # Frekans uzayında korelasyon = F1 * conj(F2)
            korelasyon_freq = F1 * np.conj(F2)

            # Geri frekans → zaman uzayı; her kaydırma için bir skor aldık
            korelasyon = np.real(np.fft.ifft(korelasyon_freq))

            toplam_korelasyon += korelasyon

        # En yüksek skorun olduğu kaydırma
        en_iyi_kaydirma = int(np.argmax(toplam_korelasyon))
        max_skor = float(toplam_korelasyon[en_iyi_kaydirma])

        # Büyük indeks aslında negatif kaydırmayı temsil edebilir.
        # (FFT çıktısı 0..pad-1 arasında; yarısından büyükse negatife çevir)
        if en_iyi_kaydirma > len(dizi1):
            en_iyi_kaydirma -= pad

        return en_iyi_kaydirma, max_skor


    # Aşama 2: Mesafe Matrisi

    def mesafe_matrisi_olustur(self, diziler):
        """
        Tüm dizi çiftleri için FFT benzerlik skoru hesaplar ve
        bir NxN mesafe matrisi oluşturur.

        Benzerlik → Mesafe dönüşümü:
        FFT yüksek skor → diziler birbirine çok benziyor → DÜŞÜK mesafe.
        Dönüşüm formülü: mesafe = 1 / (skor + epsilon)

        epsilon küçük bir sabittir (1e-9); skor sıfır veya negatif olduğunda
        sıfıra bölme hatasını önler.

        Parametreler:
            diziler (list of str): DNA dizi listesi

        Döndürür:
            np.ndarray: NxN simetrik mesafe matrisi (köşegen = 0)
        """
        n = len(diziler)
        matris = np.zeros((n, n))
        epsilon = 1e-9

        for i in range(n):
            for j in range(i + 1, n):
                _, skor = self.fft_benzerlik_skoru(diziler[i], diziler[j])

                # Skor negatif gelebilir (korelasyon ters fazda olabilir).
                # Bu durumda sıfır benzerlik gibi davran.
                skor = max(skor, 0.0)

                mesafe = 1.0 / (skor + epsilon)

                # Matris simetrik: i→j ve j→i aynı mesafe
                matris[i][j] = mesafe
                matris[j][i] = mesafe

        return matris


    # Aşama 3: Ana Hizalama Fonksiyonu

    def coklu_hizala(self, dizi_listesi, isim_listesi=None):
        """
        MAFFT algoritmasının tüm üç aşamasını sırayla çalıştırır ve
        çoklu dizi hizalamasını (MSA) döndürür.

        Algoritma akışı:
        1. FFT ile mesafe matrisi → her dizi çiftinin uzaklığını bul
        2. UPGMA → mesafe matrisinden kılavuz ağaç çıkar (birleştirme sırası)
        3. Progressive alignment → kılavuz ağaca göre dizileri ikişer ikişer hizala:
             • İki tekil dizi → Needleman-Wunsch
             • İki grup veya dizi+grup → Profil-Profil hizalaması

        Parametreler:
            dizi_listesi (list of str): Hizalanacak DNA dizileri
            isim_listesi (list of str): Dizi isimleri (opsiyonel)

        Döndürür:
            tuple: (son_isimler, hizalanmis_diziler)
                   Her iki liste de aynı sırada, aynı uzunlukta.
        """
        n = len(dizi_listesi)

        if n == 0:
            raise ValueError("Hizalanacak hiç dizi yok! FASTA dosyasını kontrol et.")

        if n == 1:
            isimler = isim_listesi or ['Dizi_0']
            return isimler, dizi_listesi[:]

        # İsimler verilmediyse otomatik ata
        if isim_listesi is None:
            isim_listesi = [f'Dizi_{i}' for i in range(n)]

        print(f"\n{'='*55}")
        print(f"  BeratTrnMAFFT / MAFFT  |  {n} dizi hizalanıyor...")
        print(f"{'='*55}")

        # AŞAMA 1 
        print("\n[AŞAMA 1] FFT tabanlı mesafe matrisi oluşturuluyor...")
        mesafe_mat = self.mesafe_matrisi_olustur(dizi_listesi)

        print("  Mesafe matrisi (küçük değer = diziler daha benzer):")
        for i in range(n):
            satir = "  " + isim_listesi[i].ljust(12) + " | "
            satir += "  ".join(f"{mesafe_mat[i][j]:.4f}" for j in range(n))
            print(satir)

        # AŞAMA 2
        print("\n[AŞAMA 2] UPGMA ile kılavuz ağaç oluşturuluyor...")
        birlestirme_sirasi = upgma(mesafe_mat, isim_listesi)

        print("  Birleştirme adımları:")
        for adim_no, (g1, g2) in enumerate(birlestirme_sirasi):
            g1_isimleri = [isim_listesi[i] for i in g1]
            g2_isimleri = [isim_listesi[i] for i in g2]
            print(f"  Adım {adim_no + 1}: {g1_isimleri}  +  {g2_isimleri}")

        # AŞAMA 3
        print("\n[AŞAMA 3] Progressive alignment (aşamalı hizalama) başlıyor...")

        # Her dizi kendi başına bir "grup" olarak başlar.
        # Anahtar: o grubun ilk üyesinin orijinal indeksi (temsilci).
        # Değer: o gruba ait hizalanmış dizi listesi.
        gruplar       = {i: [dizi_listesi[i]] for i in range(n)}
        isim_gruplari = {i: [isim_listesi[i]] for i in range(n)}

        for adim_no, (indeksler1, indeksler2) in enumerate(birlestirme_sirasi):
            # Her grubun "temsilcisi" listedeki ilk orijinal dizi indeksidir.
            # (UPGMA adımları sırasında bu temsilci sabit kalır.)
            t1 = indeksler1[0]
            t2 = indeksler2[0]

            grup1 = gruplar[t1]
            grup2 = gruplar[t2]

            isim1_listesi = [isim_listesi[i] for i in indeksler1]
            isim2_listesi = [isim_listesi[i] for i in indeksler2]

            print(f"\n  Adım {adim_no + 1}: {isim1_listesi}  +  {isim2_listesi}")

            if len(grup1) == 1 and len(grup2) == 1:
                # ─ İki tekil dizi: standart Needleman-Wunsch ─
                h1, h2, skor = needleman_wunsch(grup1[0], grup2[0])
                yeni_hizalama = [h1, h2]
                print(f"    → NW hizalaması  |  skor: {skor:.2f}")
                print(f"    {isim1_listesi[0]}: {h1}")
                print(f"    {isim2_listesi[0]}: {h2}")

            else:
                # ─ En az biri grup: Profil-Profil hizalaması ─
                profil1 = profil_olustur(grup1)
                profil2 = profil_olustur(grup2)
                yeni_hizalama = profil_profil_hizala(profil1, profil2, grup1, grup2)
                print(f"    → Profil-Profil hizalaması  |  {len(grup1)} + {len(grup2)} dizi")

            # Yeni birleşik grubu t1 temsilcisi altında topla
            yeni_isimler = isim_gruplari[t1] + isim_gruplari[t2]

            # t2 grubunu sil (artık t1 altında birleşti)
            del gruplar[t2]
            del isim_gruplari[t2]

            gruplar[t1]       = yeni_hizalama
            isim_gruplari[t1] = yeni_isimler

        # Tüm birleştirmeler bitti; artık tek bir grup var
        son_anahtar     = list(gruplar.keys())[0]
        son_hizalama    = gruplar[son_anahtar]
        son_isimler     = isim_gruplari[son_anahtar]

        print(f"\n{'='*55}")
        print("  Hizalama tamamlandı!")
        print(f"{'='*55}\n")

        return son_isimler, son_hizalama
