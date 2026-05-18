"""
BeratTrnMAFFT - Görselleştirme Modülü
MAFFT hizalama sonuçlarını grafik olarak çizen fonksiyonlar:

  1. mesafe_isi_haritasi()   → Mesafe matrisi ısı haritası (heatmap)
  2. upgma_dendrogrami()     → UPGMA kılavuz ağaç dendrogramı
  3. hizalama_gorseli()      → Renkli hizalama ızgarası
  4. kimlik_bar_grafigi()    → Çift kimlik yüzdesi bar grafiği
  5. tum_grafikleri_kaydet() → Hepsini tek PNG'ye kaydeder
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')   # ekran gerektirmeyen (headless) mod
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


# Nükleotid renk paleti
NÜKL_RENK = {
    'A': '#4ade80',   # yeşil
    'C': '#60a5fa',   # mavi
    'G': '#fbbf24',   # sarı
    'T': '#f87171',   # kırmızı
    '-': '#e2e8f0',   # açık gri (gap)
}


# 1. Mesafe Isı Haritası

def mesafe_isi_haritasi(mesafe_mat, isimler, ax=None):
    """
    NxN mesafe matrisini renkli bir ısı haritası (heatmap) olarak çizer.

    Koyu renk = diziler birbirine yakın (mesafe küçük)
    Açık renk = diziler uzak (mesafe büyük)

    Parametreler:
        mesafe_mat (np.ndarray): NxN mesafe matrisi
        isimler (list): Dizi isimleri
        ax: matplotlib Axes nesnesi (None ise yeni oluşturulur)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))

    # Köşegeni 0 yap (kendi kendine mesafe)
    mat = mesafe_mat.copy()
    np.fill_diagonal(mat, 0)

    # Kısa isimler (grafik için)
    kisa_isimler = [s[:12] for s in isimler]

    sns.heatmap(
        mat,
        annot=True,
        fmt=".3f",
        xticklabels=kisa_isimler,
        yticklabels=kisa_isimler,
        cmap="YlOrRd_r",       # ters: koyu = yakın
        linewidths=0.5,
        linecolor='white',
        ax=ax,
        annot_kws={"size": 8}
    )

    ax.set_title("FFT Mesafe Matrisi", fontsize=12, fontweight='bold', pad=12)
    ax.tick_params(axis='x', rotation=30, labelsize=8)
    ax.tick_params(axis='y', rotation=0,  labelsize=8)

    return ax


# 2. UPGMA Dendrogramı

def upgma_dendrogrami(mesafe_mat, isimler, ax=None):
    """
    UPGMA mesafe matrisinden scipy ile dendrogram çizer.

    Dendrogram, kılavuz ağacı görselleştirir: hangi diziler önce
    birleşiyor, hangiler daha uzak — bunu hiyerarşik ağaç olarak gösterir.

    Parametreler:
        mesafe_mat (np.ndarray): NxN mesafe matrisi
        isimler (list): Dizi isimleri
        ax: matplotlib Axes nesnesi
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))

    # squareform: kare matris → kondanse (üst üçgen) mesafe vektörü
    # Scipy linkage UPGMA için 'average' metodunu kullanır
    kondanse = squareform(mesafe_mat)
    baglanma  = linkage(kondanse, method='average')

    kisa_isimler = [s[:14] for s in isimler]

    dendrogram(
        baglanma,
        labels=kisa_isimler,
        orientation='left',
        ax=ax,
        leaf_font_size=9,
        color_threshold=0.7 * max(baglanma[:, 2]),
        above_threshold_color='#94a3b8'
    )

    ax.set_title("UPGMA Kılavuz Ağaç (Dendrogram)", fontsize=12,
                 fontweight='bold', pad=12)
    ax.set_xlabel("Mesafe", fontsize=9)
    ax.tick_params(axis='y', labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return ax


# 3. Renkli Hizalama Izgarası

def hizalama_gorseli(isimler, hizalanmis_diziler, ax=None):
    """
    Her nükleotidi farklı renkle boyayan hizalama ızgarası çizer.

      A -> Yeşil   C -> Mavi   G -> Sarı   T -> Kırmızı   - -> Gri

    Bu görselleştirme, korunan (conserved) pozisyonları ve
    gap bölgelerini bir bakışta gösterir.

    Parametreler:
        isimler (list): Dizi isimleri
        hizalanmis_diziler (list): Hizalanmış dizi stringleri
        ax: matplotlib Axes nesnesi
    """
    if ax is None:
        n_dizi = len(hizalanmis_diziler)
        uzunluk = len(hizalanmis_diziler[0]) if hizalanmis_diziler else 1
        genislik = max(10, uzunluk * 0.35)
        yukseklik = max(3, n_dizi * 0.7 + 1.5)
        fig, ax = plt.subplots(figsize=(genislik, yukseklik))

    n_dizi  = len(hizalanmis_diziler)
    uzunluk = len(hizalanmis_diziler[0]) if hizalanmis_diziler else 0

    for dizi_idx, dizi in enumerate(hizalanmis_diziler):
        for pos, karakter in enumerate(dizi):
            renk = NÜKL_RENK.get(karakter.upper(), '#e2e8f0')
            # Her hücre için renkli dikdörtgen çiz
            rect = mpatches.FancyBboxPatch(
                (pos, n_dizi - dizi_idx - 1),   # x, y
                0.92, 0.88,                       # genişlik, yükseklik
                boxstyle="round,pad=0.03",
                facecolor=renk,
                edgecolor='white',
                linewidth=0.4
            )
            ax.add_patch(rect)

            # Hücreye karakteri yaz (çok geniş hizalamada font küçülür)
            font_size = max(4, min(8, 120 // max(uzunluk, 1)))
            ax.text(
                pos + 0.46, n_dizi - dizi_idx - 0.56,
                karakter,
                ha='center', va='center',
                fontsize=font_size,
                fontweight='bold',
                color='#1e293b'
            )

    # Y eksenine dizi isimleri
    kisa_isimler = [s[:14] for s in isimler]
    ax.set_yticks([i + 0.44 for i in range(n_dizi)])
    ax.set_yticklabels(reversed(kisa_isimler), fontsize=8)

    # X eksenine pozisyon numaraları
    ax.set_xlim(0, uzunluk)
    ax.set_ylim(0, n_dizi)
    x_ticks = list(range(0, uzunluk, max(1, uzunluk // 10)))
    ax.set_xticks([x + 0.5 for x in x_ticks])
    ax.set_xticklabels([str(x + 1) for x in x_ticks], fontsize=7)

    ax.set_title("Çoklu Dizi Hizalaması (MSA)", fontsize=12,
                 fontweight='bold', pad=12)
    ax.set_xlabel("Pozisyon", fontsize=9)

    # Renk açıklaması (legend)
    legend_parcalari = [
        mpatches.Patch(facecolor=renk, edgecolor='#94a3b8', label=harf)
        for harf, renk in NÜKL_RENK.items()
    ]
    ax.legend(handles=legend_parcalari, loc='upper right',
              fontsize=7, ncol=5, framealpha=0.9,
              bbox_to_anchor=(1.0, -0.08))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return ax


# 4. Kimlik Yüzdesi Bar Grafiği

def kimlik_bar_grafigi(isimler, hizalanmis_diziler, ax=None):
    """
    Her dizi çifti için kimlik yüzdesini (percent identity) bar grafik ile gösterir.

    Parametreler:
        isimler (list): Dizi isimleri
        hizalanmis_diziler (list): Hizalanmış dizi stringleri
        ax: matplotlib Axes nesnesi
    """
    from .utils import kimlik_yuzdesi_hesapla

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))

    cift_isimleri = []
    yüzdeler      = []

    for i in range(len(hizalanmis_diziler)):
        for j in range(i + 1, len(hizalanmis_diziler)):
            yüzde = kimlik_yuzdesi_hesapla(hizalanmis_diziler[i],
                                           hizalanmis_diziler[j])
            isim1 = isimler[i][:10]
            isim2 = isimler[j][:10]
            cift_isimleri.append(f"{isim1}\nvs\n{isim2}")
            yüzdeler.append(yüzde)

    # Renk: yüksek kimlik → yeşil, düşük → turuncu
    bar_renkleri = [
        '#4ade80' if y >= 90 else '#60a5fa' if y >= 75 else '#fbbf24'
        for y in yüzdeler
    ]

    barlar = ax.bar(range(len(yüzdeler)), yüzdeler,
                    color=bar_renkleri, edgecolor='white',
                    linewidth=0.8, width=0.6)

    # Barların üstüne yüzde değerini yaz
    for bar, yuzde in zip(barlar, yüzdeler):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.8,
            f"%{yuzde:.1f}",
            ha='center', va='bottom', fontsize=8, fontweight='bold',
            color='#1e293b'
        )

    ax.set_xticks(range(len(cift_isimleri)))
    ax.set_xticklabels(cift_isimleri, fontsize=7)
    ax.set_ylim(0, 115)
    ax.set_ylabel("Kimlik Yüzdesi (%)", fontsize=9)
    ax.set_title("Çift Kimlik Yüzdeleri (Pairwise Identity)", fontsize=12,
                 fontweight='bold', pad=12)

    # Referans çizgileri
    for ref in [80, 90, 100]:
        ax.axhline(ref, color='#cbd5e1', linewidth=0.8, linestyle='--')

    # Renk açıklaması
    legend_parcalari = [
        mpatches.Patch(facecolor='#4ade80', label='≥ %90 (Çok Benzer)'),
        mpatches.Patch(facecolor='#60a5fa', label='%75–90'),
        mpatches.Patch(facecolor='#fbbf24', label='< %75'),
    ]
    ax.legend(handles=legend_parcalari, fontsize=7, loc='upper right')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return ax


# 5. Tüm Grafikleri Tek PNG'ye Kaydet

def tum_grafikleri_kaydet(mesafe_mat, isimler, hizalanmis_diziler,
                          cikti_dosyasi="mafft_sonuclari.png"):
    """
    Dört grafiği 2×2 düzeninde tek bir PNG dosyasına kaydeder.

    Parametreler:
        mesafe_mat (np.ndarray): NxN mesafe matrisi
        isimler (list): Dizi isimleri
        hizalanmis_diziler (list): Hizalanmış dizi stringleri
        cikti_dosyasi (str): Kaydedilecek PNG dosyasının yolu
    """
    plt.style.use('seaborn-v0_8-whitegrid')
    fig = plt.figure(figsize=(16, 13))
    fig.patch.set_facecolor('#f8fafc')

    # Başlık
    fig.suptitle(
        "BeratTrnMAFFT — MAFFT Algoritması Sonuçları\n"
        "Öğrenci No: 221201018  |  İstanbul Rumeli Üniversitesi",
        fontsize=14, fontweight='bold', color='#1a2744', y=0.98
    )

    # 2×2 grid: sol üst, sağ üst, sol alt, sağ alt
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    mesafe_isi_haritasi(mesafe_mat, isimler, ax=ax1)
    upgma_dendrogrami(mesafe_mat, isimler, ax=ax2)
    hizalama_gorseli(isimler, hizalanmis_diziler, ax=ax3)
    kimlik_bar_grafigi(isimler, hizalanmis_diziler, ax=ax4)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(cikti_dosyasi, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()

    print(f"Grafik kaydedildi: {cikti_dosyasi}")
    return cikti_dosyasi
