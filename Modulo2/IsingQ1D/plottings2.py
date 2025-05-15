import numpy as np
import matplotlib.pyplot as plt

#Load data from text file
data3 = np.loadtxt("./dinamica/suscLong3.txt")
data4 = np.loadtxt("./dinamica/suscLong4.txt")
data5 = np.loadtxt("./dinamica/suscLong5.txt")
data6 = np.loadtxt("./dinamica/suscLong6.txt")
data7 = np.loadtxt("./dinamica/suscLong7.txt")
data8 = np.loadtxt("./dinamica/suscLong8.txt")
data9 = np.loadtxt("./dinamica/suscLong9.txt")
data10 = np.loadtxt("./dinamica/suscLong10.txt")
'''
#Sample data with given frequency
w = 1  # (Modifica questo valore a tuo piacimento)
sample = data[::w]
'''


#Name columns
g = data3[:,0]
h = np.linspace(-0.02+(0.02/52), 0.02-(0.02/52), 50)

#PLOTTING CHI LONG PER UN SITO AL VARIARE DI g (h = piccolo) (usando Mz sensibile ad h)
plt.errorbar(g, data3[:,25],fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, data4[:,25],fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, data5[:,25],fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, data6[:,25],fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, data7[:,25],fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, data8[:,25],fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g[1:], data9[1:,25],fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g[1:], data10[1:,25],fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('g')
plt.ylabel('χz')
plt.title('χz (h ⟶ 0)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()


dataDEF3 = np.loadtxt("./dinamica/suscLongDEF3.txt")
dataDEF4 = np.loadtxt("./dinamica/suscLongDEF4.txt")
dataDEF5 = np.loadtxt("./dinamica/suscLongDEF5.txt")
dataDEF6 = np.loadtxt("./dinamica/suscLongDEF6.txt")
dataDEF7 = np.loadtxt("./dinamica/suscLongDEF7.txt")
dataDEF8 = np.loadtxt("./dinamica/suscLongDEF8.txt")
dataDEF9 = np.loadtxt("./dinamica/suscLongDEF9.txt")
dataDEF10 = np.loadtxt("./dinamica/suscLongDEF10.txt")

#FSS for Chi (Mz)
b = 1./8.
nu = 1.
gc = 1.0
gamma = 7./4.

x3 = (g-gc)*(3**(1/nu))
y3 = data3[:,25]/(3**(gamma/nu))
x4 = (g-gc)*(4**(1/nu))
y4 = data4[:,25]/(4**(gamma/nu))
x5 = (g-gc)*(5**(1/nu))
y5 = data5[:,25]/(5**(gamma/nu))
x6 = (g-gc)*(6**(1/nu))
y6 = data6[:,25]/(6**(gamma/nu))
x7 = (g-gc)*(7**(1/nu))
y7 = data7[:,25]/(7**(gamma/nu))
x8 = (g-gc)*(8**(1/nu))
y8 = data8[:,25]/(8**(gamma/nu))
x9 = (g-gc)*(9**(1/nu))
y9 = data9[:,25]/(9**(gamma/nu))
x10 = (g-gc)*(10**(1/nu))
y10 = data10[:,25]/(10**(gamma/nu))

plt.errorbar(x3, y3,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x4, y4,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x5, y5,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x6, y6,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x7, y7,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x8, y8,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x9[1:], y9[1:],fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x10[1:], y10[1:],fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('(g-gc)N^(1/ν)')
plt.ylabel('χz/(N^(γ/ν))')
plt.title('Finite size scaling per χz (h ⟶ 0)')
plt.grid(True)
plt.legend()
#plt.ylim(-0.25,0.2)
#plt.xlim(-0.0, 7.0)
plt.show()


#FSS for Chi (MzDEF)
x3 = (g-gc)*(3**(1/nu))
y3 = dataDEF3[:,25]/(3**(gamma/nu))
x4 = (g-gc)*(4**(1/nu))
y4 = dataDEF4[:,25]/(4**(gamma/nu))
x5 = (g-gc)*(5**(1/nu))
y5 = dataDEF5[:,25]/(5**(gamma/nu))
x6 = (g-gc)*(6**(1/nu))
y6 = dataDEF6[:,25]/(6**(gamma/nu))
x7 = (g-gc)*(7**(1/nu))
y7 = dataDEF7[:,25]/(7**(gamma/nu))
x8 = (g-gc)*(8**(1/nu))
y8 = dataDEF8[:,25]/(8**(gamma/nu))
x9 = (g-gc)*(9**(1/nu))
y9 = dataDEF9[:,25]/(9**(gamma/nu))
x10 = (g-gc)*(10**(1/nu))
y10 = dataDEF10[:,25]/(10**(gamma/nu))

Kdata3 = np.loadtxt("./N3/primaSp3.txt") 
Kdata4 = np.loadtxt("./N4/primaSp4.txt") 
Kdata5 = np.loadtxt("./N5/primaSp5.txt") 
Kdata6 = np.loadtxt("./N6/primaSp6.txt") 
Kdata7 = np.loadtxt("./N7/primaSp7.txt") 
Kdata8 = np.loadtxt("./N8/primaSp8.txt") 
Kdata9 = np.loadtxt("./N9/primaSp9.txt") 
Kdata10 = np.loadtxt("./N10/primaSp10.txt") 

'''
#Sample data with given frequency
w = 1  # (Modifica questo valore a tuo piacimento)
sample = data[::w]
'''

#PLOTTING FSS DI Mz ALLA TRANSIZIONE DI PRIMO ORDINE (g = 0.5)
plt.errorbar(Kdata3[:,0], Kdata3[:,1],fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(Kdata4[:,0], Kdata4[:,1],fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(Kdata5[:,0], Kdata5[:,1],fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(Kdata6[:,0], Kdata6[:,1],fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(Kdata7[:,0], Kdata7[:,1],fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(Kdata8[:,0], Kdata8[:,1],fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(Kdata9[:,0], Kdata9[:,1],fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(Kdata10[:,0], Kdata10[:,1],fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('K')
plt.ylabel('Mz')
plt.title('Finite size scaling per Mz (First Order Trans. g = 0.5)')
plt.grid(True)
plt.legend()
#plt.ylim(-0.25,0.2)
#plt.xlim(-0.0, 7.0)
plt.show()