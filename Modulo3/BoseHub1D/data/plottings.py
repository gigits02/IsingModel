import numpy as np
import matplotlib.pyplot as plt

#PLOTTINGS DI CONVERGENZE DI DMRG, ENERGIE PER SITI, DIAGRAMMA DI FASE
#Load data from text file
dataMOTT = np.loadtxt("./infMI2.txt")

#Name columns
L_mott = dataMOTT[:,0]
Energy_mott = dataMOTT[:,1]
Ebond_mott = dataMOTT[:,2]
Edens_mott = dataMOTT[:,3]
TruncErr_mott = dataMOTT[:,4]
Ent_mott = dataMOTT[:,5]

#PLOTTINGS ENERGY
plt.errorbar(L_mott, Energy_mott, fmt='.', label = 'MI (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('L')
plt.ylabel('E')
plt.title('Energia')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTINGS ENERGY DENSITY
plt.errorbar(L_mott, Edens_mott, fmt='.', label = 'MI (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('L')
plt.ylabel('E/L')
plt.title('Densità di energia')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTINGS ENERGY PER BOND
plt.errorbar(L_mott, Ebond_mott,fmt='.', label = 'MI (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('L')
plt.ylabel('E_bond')
plt.title('Variazione di energia tra update')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTINGS TRUNC. ERROR
plt.errorbar(L_mott, TruncErr_mott, fmt='.', label = 'inf-DMRG (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('L')
plt.ylabel('trunc_err')
plt.title('Errore di troncamento')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTINGS BIPARTITE ENTANGLEMENT
plt.errorbar(L_mott/2, Ent_mott, fmt='.', label = 'inf-DMRG (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel(r'$L_A$')
plt.ylabel('Ent')
plt.title('Entanglement bipartito')
plt.grid(True)
plt.legend()
plt.tight_layout()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTINGS ENERGIES PER SITE FOR EACH MU (J = 0.0) 
#Load data from text file
data0 = np.loadtxt("./mu0/dataMu0.txt")
data1 = np.loadtxt("./mu1/dataMu1.txt")
data2 = np.loadtxt("./mu2/dataMu2.txt")

'''
#Sample data with given frequency
w = 1  # (Modifica questo valore a tuo piacimento)
sample = data[::w]
'''
#Name columns
n0 = data0[:,0]
E0 = data0[:,1]
n1 = data1[:,0]
E1 = data1[:,1]
n2 = data2[:,0]
E2 = data2[:,1]
n00 = 0.5 + 0
n01 = 0.5 + 1
n02 = 0.5 + 2

plt.errorbar(n0, E0/2, fmt='.', label = 'μ =0', linestyle='--', color = 'blue', capsize=3, markersize = 4)
plt.errorbar(n1, E1/2, fmt='.', label = 'μ =1', linestyle='--', color = 'orange', capsize=3, markersize = 4)
plt.errorbar(n2, E2/2, fmt='.', label = 'μ =2', linestyle='--', color = 'green', capsize=3, markersize = 4)
plt.axvline(x=n00, color='r', linestyle='--',  ymin = -1.00, ymax = 0.2, label = 'valori minimi previsti')
plt.axvline(x=n01, color='r', linestyle='--', ymin = -1.00, ymax = 0.15)
plt.axvline(x=n02, color='r', linestyle='--', ymin = -1.00, ymax = 0.05)
plt.xlabel('n')
plt.ylabel('E1')
plt.title('Energia per sito (J=0)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTINGS FINITE-SYSTEM CONVERGENCES

#Load data from text file
dataMOTT = np.loadtxt("./finMI.txt")
#Name columns
L_mott = dataMOTT[:,0]
Energy_mott = dataMOTT[:,1]
Ebond_mott = dataMOTT[:,2]
Edens_mott = dataMOTT[:,3]
TruncErr_mott = dataMOTT[:,4]
Ent_mott = dataMOTT[:,5]
x = np.linspace(1, 99, num=99)
step = x.astype(int)

plt.errorbar(step, Energy_mott, fmt='.', label = 'MI (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('#step')
plt.ylabel('E')
plt.title('Finite-system DMRG')
plt.grid(True)
plt.axvline(x=50, color='r', linestyle='--', label = 'inizio del finite-system')
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

plt.errorbar(step, Ebond_mott, fmt='.', label = 'MI (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('#step')
plt.ylabel('E_bond')
plt.title('Finite-system DMRG')
plt.grid(True)
plt.axvline(x=50, color='r', linestyle='--', label = 'inizio del finite-system')
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#COMPARING NRG AND DMRG

#Load data from text file
dataNRG = np.loadtxt("./nrg.txt")
dataDMRG = np.loadtxt("./inf.txt")
x = np.linspace(1, 30, num=30)
step = x.astype(int)
plt.errorbar(step[:5], dataNRG[:5,1], fmt='.', label = 'NRG (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(step[:5], dataDMRG[:5,1], fmt='.', label = 'DMRG (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('#step')
plt.ylabel('E')
plt.title('Confronto NRG/DMRG')
plt.grid(True)
plt.legend()
plt.show()

plt.errorbar(step[:], dataNRG[:,2], fmt='.', label = 'NRG (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(step[:], dataDMRG[:,2], fmt='.', label = 'DMRG (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('#step')
plt.ylabel('E_bond')
plt.title('Confronto NRG/DMRG')
plt.grid(True)
plt.legend()
plt.ylim(-1.010,-1.004)
plt.tight_layout()
plt.show()

plt.errorbar(step[:], dataNRG[:,3], fmt='.', label = 'NRG (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(step[:], dataDMRG[:,3], fmt='.', label = 'DMRG (μ,J) = (0.5, 0.05)', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('#step')
plt.ylabel('E/L')
plt.title('Confronto NRG/DMRG')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

