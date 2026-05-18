"""
Kullanim:
    from BeratTrnMAFFT import MAFFTHizalayici
    from BeratTrnMAFFT.utils import fasta_oku, hizalamali_yazdir

    kayitlar = fasta_oku("diziler.fasta")
    isimler  = [x[0] for x in kayitlar]
    diziler  = [x[1] for x in kayitlar]

    hizalayici = MAFFTHizalayici()
    son_isimler, hizalanmis = hizalayici.coklu_hizala(diziler, isimler)
    hizalamali_yazdir(son_isimler, hizalanmis)
"""

from .core import MAFFTHizalayici

__version__ = "1.0.0"
__author__  = "Berat Turan"
__email__   = "turanberatr@gmail.com"

__all__ = ["MAFFTHizalayici"]
