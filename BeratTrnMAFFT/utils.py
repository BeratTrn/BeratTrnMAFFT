"""
BeratTrnMAFFT - Yardımcı Fonksiyonlar
Bu modül FASTA dosyası okuma ve hizalama sonuçlarını yazdırma gibi
tekrar kullanılan küçük yardımcı araçları içerir.
"""


def fasta_oku(dosya_adi):
    """
    FASTA formatındaki bir dosyayı okur ve içindeki dizi bilgilerini döndürür.

    FASTA formatı şu şekilde görünür:
        >DizininIsmi
        ACGTACGT...
        >IkinciDizininIsmi
        GCGTACGT...

    '>' ile başlayan satırlar isim/başlık satırıdır.
    Diğer satırlar nükleotid dizisidir (birden fazla satır da olabilir,
    hepsini birleştiriyoruz).

    Parametreler:
        dosya_adi (str): Okunacak FASTA dosyasının tam yolu

    Döndürür:
        list of tuple: [(isim1, dizi1), (isim2, dizi2), ...] şeklinde bir liste
                       isimler baştaki '>' işareti olmadan döner
                       diziler büyük harfe çevrilmiş olarak döner
    """
    sonuc = []
    gecici_isim = None
    gecici_satırlar = []

    with open(dosya_adi, 'r', encoding='utf-8') as f:
        for satir in f:
            satir = satir.strip()

            # Boş satırları atla
            if not satir:
                continue

            if satir.startswith('>'):
                # Yeni bir dizi başlıyor.
                # Eğer daha önce bir dizi toplanıyorsa onu listeye ekle.
                if gecici_isim is not None and gecici_satırlar:
                    dizi = ''.join(gecici_satırlar).upper()
                    sonuc.append((gecici_isim, dizi))

                # Başlık satırından '>' işaretini kaldır, geri kalanı isim olarak al
                gecici_isim = satir[1:].strip()
                gecici_satırlar = []
            else:
                # Dizi satırı, biriktiriyoruz
                gecici_satırlar.append(satir)

    # Döngü bitince dosyadaki son diziyi de eklemeyi unutma!
    # (İlk sürümde bunu unutmuştum, son dizi hep kayboluyordu.)
    if gecici_isim is not None and gecici_satırlar:
        dizi = ''.join(gecici_satırlar).upper()
        sonuc.append((gecici_isim, dizi))

    return sonuc


def hizalamali_yazdir(isimler, hizalanmis_diziler):
    """
    Hizalanmış dizileri düzgün sütunlar halinde terminale yazdırır.

    Örnek çıktı:
        ÇOKLU DİZİ HİZALAMASI SONUCU
        Seq1  : ACG-TACGT
        Seq2  : ACGTACG-T
        ...

    Ayrıca basit bir istatistik olarak tamamen korunan (conserved) pozisyon
    sayısını da gösterir.

    Parametreler:
        isimler (list): Dizi isimlerinin listesi
        hizalanmis_diziler (list): Hizalanmış dizi stringleri ('-' gap karakteridir)
    """
    if not isimler or not hizalanmis_diziler:
        print("Gösterilecek hizalama yok.")
        return

    # En uzun ismin kaç karakter olduğunu bul; hizalı yazdırmak için lazım
    max_uzunluk = max(len(isim) for isim in isimler)

    print()
    print("=" * 62)
    print("  ÇOKLU DİZİ HİZALAMASI SONUCU  (BeratTrnMAFFT / MAFFT)")
    print("=" * 62)

    for isim, dizi in zip(isimler, hizalanmis_diziler):
        # İsmi sola doldur ki tüm ':' işaretleri aynı sütunda dursun
        bosluk_doldurma = " " * (max_uzunluk - len(isim))
        print(f"  {isim}{bosluk_doldurma} : {dizi}")

    print("-" * 62)

    # Hizalama uzunluğu (tüm diziler aynı uzunlukta olmalı)
    hizalama_uzunlugu = len(hizalanmis_diziler[0])

    # Tamamen korunan pozisyonları say:
    # Bir pozisyon "korunmuş" demek → tüm dizilerde aynı karakter var (gap olmadan)
    korunan = 0
    for pos in range(hizalama_uzunlugu):
        karakterler = set()
        tamamı_gap_mi = True
        for dizi in hizalanmis_diziler:
            if pos < len(dizi):
                k = dizi[pos]
                if k != '-':
                    karakterler.add(k)
                    tamamı_gap_mi = False
        # Eğer o pozisyonda yalnızca tek tip nükleotid varsa → korunmuş
        if not tamamı_gap_mi and len(karakterler) == 1:
            korunan += 1

    print(f"  Hizalama uzunluğu  : {hizalama_uzunlugu} sütun")
    print(f"  Korunan pozisyon   : {korunan} / {hizalama_uzunlugu}")
    print("=" * 62)
    print()


def kimlik_yuzdesi_hesapla(dizi1, dizi2):
    """
    İki hizalanmış dizi arasındaki yüzde kimliği (percent identity) hesaplar.

    Formül:  (aynı olan pozisyon sayısı) / (gap olmayan toplam pozisyon) * 100

    Gap içeren pozisyonlar sayıma dahil edilmez; sadece iki tarafında da
    nükleotid olan sütunlar değerlendirilir.

    Parametreler:
        dizi1 (str): Hizalanmış birinci dizi
        dizi2 (str): Hizalanmış ikinci dizi

    Döndürür:
        float: 0.0 ile 100.0 arasında kimlik yüzdesi
    """
    eslesenler = 0
    toplam = 0

    for k1, k2 in zip(dizi1, dizi2):
        if k1 == '-' or k2 == '-':
            continue  # gap olan pozisyonu atla
        toplam += 1
        if k1 == k2:
            eslesenler += 1

    if toplam == 0:
        return 0.0

    return round((eslesenler / toplam) * 100, 2)
