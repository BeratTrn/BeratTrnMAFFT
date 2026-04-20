import numpy as np

class MafftCore:
    def __init__(self):
        # DNA nükleotidlerimiz
        self.harfler = ['A', 'C', 'G', 'T']
        
    def sinyale_cevir(self, dizi):
        # Harfleri 1 ve 0'lardan oluşan matrislere çeviriyoruz (one-hot)
        sinyaller = {h: np.zeros(len(dizi)) for h in self.harfler}
        
        for i, harf in enumerate(dizi):
            if harf in self.harfler:
                sinyaller[harf][i] = 1
                
        return sinyaller

    def fft_ile_karsilastir(self, dizi1, dizi2):
        # İki diziyi frekans uzayında karşılaştırıp en iyi skoru bulur
        sig1 = self.sinyale_cevir(dizi1)
        sig2 = self.sinyale_cevir(dizi2)
        
        # FFT'de boyut uyuşmazlığı olmaması için uzunlukları topluyoruz
        n = len(dizi1) + len(dizi2)
        toplam_skor = np.zeros(n)
        
        for harf in self.harfler:
            # Sinyalleri Fourier uzayına aktar
            fft1 = np.fft.fft(sig1[harf], n=n)
            fft2 = np.fft.fft(sig2[harf], n=n)
            
            # Çapraz korelasyon formülü (ikincinin eşleniği ile çarpım)
            korelasyon = np.fft.ifft(fft1 * np.conj(fft2))
            toplam_skor += np.real(korelasyon)
            
        # En yüksek benzerliğin olduğu kaydırma miktarını bul
        en_iyi_kaydirma = np.argmax(toplam_skor)
        
        # Negatif kaydırma kontrolü
        if en_iyi_kaydirma > len(dizi1):
            en_iyi_kaydirma -= n
            
        return en_iyi_kaydirma, toplam_skor[en_iyi_kaydirma]

    def benzerlik_matrisi_olustur(self, diziler):
        # Tüm dizileri ikili ikili karşılaştırıp skor matrisi çıkarıyoruz
        adet = len(diziler)
        matris = np.zeros((adet, adet))
        
        for i in range(adet):
            for j in range(i + 1, adet): # Kendisiyle karşılaştırmamak için i+1
                kaydirma, skor = self.fft_ile_karsilastir(diziler[i], diziler[j])
                matris[i][j] = skor
                matris[j][i] = skor # Matris simetrik olduğu için kopyalıyoruz
                
        return matris

    def asamali_hizalama(self, diziler):
        # Benzerlik matrisine göre kılavuz ağaç oluşturma
        matris = self.benzerlik_matrisi_olustur(diziler)
        adet = len(diziler)
        
        kalanlar = list(range(adet))
        siralamamiz = []
        
        # En yüksek skorlu iki diziyi bul (Merkez grup)
        max_skor = -1
        en_iyi_cift = (0, 0)
        
        for i in range(adet):
            for j in range(i + 1, adet):
                if matris[i][j] > max_skor:
                    max_skor = matris[i][j]
                    en_iyi_cift = (i, j)
                    
        # İlk iki diziyi sıraya ekle
        siralamamiz.append(en_iyi_cift[0])
        siralamamiz.append(en_iyi_cift[1])
        kalanlar.remove(en_iyi_cift[0])
        kalanlar.remove(en_iyi_cift[1])
        
        # Dışarıda kalan dizileri ortalama benzerliklerine göre gruba dahil et
        while kalanlar:
            siradaki_dizi = -1
            en_yuksek_ort = -1
            
            for k in kalanlar:
                ort_benzerlik = sum(matris[k][a] for a in siralamamiz) / len(siralamamiz)
                if ort_benzerlik > en_yuksek_ort:
                    en_yuksek_ort = ort_benzerlik
                    siradaki_dizi = k
            
            siralamamiz.append(siradaki_dizi)
            kalanlar.remove(siradaki_dizi)
            
        print("Bulunan Hizalama Sırası (İndekslere göre):", siralamamiz)
        return siralamamiz
    
    

"""
Ne Yaptık?
Bu kodu çalıştırdığında algoritma harfleri hiç klasik yollarla 
(+1/-1 diye) alt alta koymaya çalışmadan, doğrudan FFT ile matematikteki 
"ses eşleşmesi" mantığını kullanarak GCGT dizisinin, büyük dizinin içinde 
tam olarak 2 adım kaydırıldığında mükemmel uyum sağladığını şak diye bulacak.
Sinyallere dönüştürme (One-hot encoding) ve FFT ile iki diziyi üst üste bindirip 
tepe noktasını bulma mantığı kafana yattı mı?
"""