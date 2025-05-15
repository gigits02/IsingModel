import numpy as np

#ESEGUE BINNING SUI DATI GREZZI (A DIMENSIONE FISSATA, GIRANDO SU BETA)
#E SALVA LE GRANDEZZE SU FILE DI TESTO "dim.txt"
# Carica dati da file di testo

for i in range(410, 472, 2):
    dati = np.loadtxt(f"L180/beta0{i}.txt")

    # Termalizzazione time  
    t = 100000

    # Nome colonne
    En = dati[t:, 0]
    m = dati[t:, 1]

    # Definisci binsize e binnumber
    dim = len(En)
    binsize = 2000
    binnumber = dim // binsize

    # Reshape i dati per calcolare la media e la deviazione standard
    En_reshaped = En[:binnumber * binsize].reshape(binnumber, binsize)
    m_reshaped = m[:binnumber * binsize].reshape(binnumber, binsize)

    # Calcola la media e la deviazione standard per energia e magnetizzazione
    mediaEn = np.mean(En_reshaped, axis=1)
    quadratoEn = np.mean(En_reshaped**2, axis=1)
    mediaM = np.mean(m_reshaped, axis=1)
    quadratoM = np.mean(m_reshaped**2, axis=1)
    moduloM = np.mean(np.abs(m_reshaped), axis=1)

    # Calcola le grandezze finali
    mediaEn_media = np.mean(mediaEn)
    mediaM_media = np.mean(mediaM)
    quadratoEn_media = np.mean(quadratoEn)
    quadratoM_media = np.mean(quadratoM)
    moduloM_media = np.mean(moduloM)

    # Calcola deviazione standard
    s_q_m_En = quadratoEn_media - mediaEn_media**2
    s_q_m_M = quadratoM_media - mediaM_media**2
    s_q_m_moduloM = quadratoM_media - moduloM_media**2

    # Leggi le vecchie soluzioni se presenti
    old_solutions = []
    try:
        with open("180.txt", "r") as file:
            old_solutions = file.readlines()
    except FileNotFoundError:
        pass

    # Scrivi le soluzioni su file
    with open("180.txt", "w") as file:
        # Scrivi le vecchie soluzioni
        for solution in old_solutions:
            file.write(solution)

        # Scrivi le nuove soluzioni
        file.write(f"{i/1000} {mediaEn_media} {s_q_m_En} {mediaM_media} {s_q_m_M} {moduloM_media} {s_q_m_moduloM}\n")