import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.loadtxt("fit.txt")
Ll = [10., 20., 30., 40., 60., 80., 100., 120., 140., 160., 180.]
gammaSUnu = [0., 0., 0.]
unoSUnu = [0., 0., 0.]
erGammaSUnu = [0., 0., 0.]
erUnoSUnu = [0., 0., 0.]
beta = [0., 0., 0.]
erbeta = [0., 0., 0.]

for i in range(3):
    
    L = Ll[i:]
    betaC = data[i:,0]
    chiM = data[i:,1]
    erBeta = data[i:,2]
    erChi = data[i:,3]

    #FIT LEGGE DI POTENZA PER CHIM E PER BETAC
    def Chipower_func(x, c1, c2, esponente):
        return c1 + c2*((x)**esponente)

    def Betapower_func(x, betacr, c2, esponente):
        return betacr + c2*((x)**(-esponente))

    p0 = (0., 1., 1.)
    poptChi, pcovChi = curve_fit(Chipower_func, L, chiM, sigma = erChi, p0 = p0)
    p0 = (0.4410, 1.0, 0.9638)
    poptBeta, pcovBeta = curve_fit(Betapower_func, L, betaC, sigma = erBeta, p0 = p0)

    print('L A FINE CICLO:', L)

    #PRINTING RESULTS
    print('GAMMA/NU:')
    print(poptChi[2])
    print('PCOV:')
    print(pcovChi)
    print('BETAC:')
    print(poptBeta[0])
    print('1/NU:')
    print(poptBeta[2])
    print('PCOV:')
    print(pcovBeta)
    
    #save vectors for systematics
    gammaSUnu[i] = poptChi[2]   
    unoSUnu[i] = poptBeta[2]
    beta[i] = poptBeta[0]
    erGammaSUnu[i] = np.sqrt(pcovChi[2][2])
    erUnoSUnu[i] = np.sqrt(pcovBeta[2][2])
    erbeta[i] = np.sqrt(pcovBeta[0][0])    

    #PLOTTINGS
    plt.errorbar(L, chiM, yerr=erChi, fmt='v', label = 'χ_max stimati', linestyle='None', capsize=3, markersize = 4)
    x = np.linspace(10, 180, 1000)
    plt.errorbar(x, Chipower_func(x, *poptChi), fmt='-', label = 'Power law fit' , linestyle='--', capsize=3, markersize = 4)
    plt.legend()
    plt.xlabel('L')
    plt.ylabel('χ_max')
    plt.title('L_min = 10')
    plt.grid(True)
    plt.show()

    plt.errorbar(L, betaC, yerr=erBeta, fmt='v', label = 'β_PC stimati', linestyle='None', capsize=3, markersize = 4)
    plt.errorbar(x, Betapower_func(x, *poptBeta), fmt='-', label = 'Power law fit', linestyle='--', capsize=3, markersize = 4)
    plt.legend()
    plt.xlabel('L')
    plt.ylabel('β_PC')
    plt.title('L_min = 10')
    plt.grid(True)
    plt.show()
    
#PLOTTING SYSTEMATICS ON LMIN
Lmin = [10.,20.,30.]
#GAMMA/NU
print('gammaSUnu =', gammaSUnu)
plt.errorbar(Lmin, gammaSUnu, yerr=erGammaSUnu, fmt='v', label = 'Stime', linestyle='None', capsize=3, markersize = 4)
plt.legend()
plt.xlabel('L_min')
plt.ylabel(' γ/ν')
plt.axhline(y=1.75, color='r', linestyle='--', label='Valore teorico')
plt.grid(True)
plt.show()
#1/NU
plt.errorbar(Lmin, unoSUnu, yerr=erUnoSUnu, fmt='v', label = 'Stime', linestyle='None', capsize=3, markersize = 4)
plt.legend()
plt.xlabel('Lmin')
plt.ylabel('1/ν')
plt.axhline(y=1.0, color='r', linestyle='--', label='Valore teorico')
plt.grid(True)
plt.show()
#BETAC
plt.errorbar(Lmin, beta, yerr=erbeta, fmt='v', label = 'Stime', linestyle='None', capsize=3, markersize = 4)
plt.legend()
plt.xlabel('Lmin')
plt.ylabel('βC')
plt.axhline(y=0.4406, color='r', linestyle='--', label='Valore teorico')
plt.grid(True)
plt.show()
