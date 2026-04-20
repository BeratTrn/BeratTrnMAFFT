def fasta_oku(dosya_adi):
    # Fasta formatındaki dosyayı okuyup dizileri liste yapar
    diziler = []
    gecici_dizi = []
    
    with open(dosya_adi, 'r') as dosya:
        for satir in dosya:
            satir = satir.strip()
            
            if satir.startswith(">"): # Başlık satırı geldiyse
                if gecici_dizi: # Eskisini listeye ekle
                    diziler.append("".join(gecici_dizi))
                    gecici_dizi = []
            else:
                gecici_dizi.append(satir)
                
        # Dosya bittiğinde son diziyi de ekle
        if gecici_dizi:
            diziler.append("".join(gecici_dizi))
            
    return diziler