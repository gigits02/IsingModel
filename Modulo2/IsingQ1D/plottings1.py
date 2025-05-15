import numpy as np
import matplotlib.pyplot as plt
'''
#Load data from text file
data3 = np.loadtxt("./N3/0solutions.txt")
data4 = np.loadtxt("./N4/0solutions.txt")
data5 = np.loadtxt("./N5/0solutions.txt")
data6 = np.loadtxt("./N6/0solutions.txt")
data7 = np.loadtxt("./N7/0solutions.txt")
data8 = np.loadtxt("./N8/0solutions.txt")
data9 = np.loadtxt("./N9/0solutions.txt")
data10 = np.loadtxt("./N10/0solutions.txt")

gdata3 = np.loadtxt("./N3/gsolutions.txt")
gdata4 = np.loadtxt("./N4/gsolutions.txt")
gdata5 = np.loadtxt("./N5/gsolutions.txt")
gdata6 = np.loadtxt("./N6/gsolutions.txt")
gdata7 = np.loadtxt("./N7/gsolutions.txt")
gdata8 = np.loadtxt("./N8/gsolutions.txt")
gdata9 = np.loadtxt("./N9/gsolutions.txt")
gdata10 = np.loadtxt("./N10/gsolutions.txt")


'''
'''
#Sample data with given frequency
w = 1  # (Modifica questo valore a tuo piacimento)
sample = data[::w]
'''

'''

#Name columns
g = data3[:,4]
h = gdata3[:,4]

mx3 = data3[:,0]
mz3 = data3[:,1]
gap3 = data3[:,2]
ggap3 = data3[:,3]
mx3g = gdata3[:,0]
mz3g = gdata3[:,1]
gap3g = gdata3[:,2]
ggap3g = gdata3[:,3]

mx4 = data4[:,0]
mz4 = data4[:,1]
gap4 = data4[:,2]
ggap4 = data4[:,3]
mx4g = gdata4[:,0]
mz4g = gdata4[:,1]
gap4g = gdata4[:,2]
ggap4g = gdata4[:,3]

mx5 = data5[:,0]
mz5 = data5[:,1]
gap5 = data5[:,2]
ggap5 = data5[:,3]
mx5g = gdata5[:,0]
mz5g = gdata5[:,1]
gap5g = gdata5[:,2]
ggap5g = gdata5[:,3]

mx6 = data6[:,0]
mz6 = data6[:,1]
gap6 = data6[:,2]
ggap6 = data6[:,3]
mx6g = gdata6[:,0]
mz6g = gdata6[:,1]
gap6g = gdata6[:,2]
ggap6g = gdata6[:,3]

mx7 = data7[:,0]
mz7 = data7[:,1]
gap7 = data7[:,2]
ggap7 = data7[:,3]
mx7g = gdata7[:,0]
mz7g = gdata7[:,1]
gap7g = gdata7[:,2]
ggap7g = gdata7[:,3]

mx8 = data8[:,0]
mz8 = data8[:,1]
gap8 = data8[:,2]
ggap8 = data8[:,3]
mx8g = gdata8[:,0]
mz8g = gdata8[:,1]
gap8g = gdata8[:,2]
ggap8g = gdata8[:,3]

mx9 = data9[:,0]
mz9 = data9[:,1]
gap9 = data9[:,2]
ggap9 = data9[:,3]
mx9g = gdata9[:,0]
mz9g = gdata9[:,1]
gap9g = gdata9[:,2]
ggap9g = gdata9[:,3]

mx10 = data10[:,0]
mz10 = data10[:,1]
gap10 = data10[:,2]
ggap10 = data10[:,3]
mx10g = gdata10[:,0]
mz10g = gdata10[:,1]
gap10g = gdata10[:,2]
ggap10g = gdata10[:,3]

#PLOTTING MAGNETIZZAZIONE TRASVERSA PER UN SITO AL VARIARE DI g (h = 0.00)
plt.errorbar(g, mx3,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mx4,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mx5,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mx6,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mx7,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mx8,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mx9,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mx10,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('g')
plt.ylabel('Mx')
plt.title('Mx (h=0.0)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTING MAGNETIZZAZIONE LONGITUDINALE PER UN SITO AL VARIARE DI g (h = 0.00)
plt.errorbar(g, mz3,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mz4,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mz5,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mz6,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mz7,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mz8,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mz9,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, mz10,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('g')
plt.ylabel('Mz')
plt.title('Mz (h=0.0)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTING GAP ENERGETICO (E1-E0) AL VARIARE DI g (h = 0.00)
plt.errorbar(g, gap3,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, gap4,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, gap5,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, gap6,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, gap7,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, gap8,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, gap9,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, gap10,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('g')
plt.ylabel('ΔE')
plt.title('ΔE (h=0.0)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTING FUNZIONE UNIVERSALE ALPHA(g) (E1-E0) AL VARIARE DI g (h = 0.00) per g < gc
plt.errorbar(g, -np.log(gap3)/3,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, -np.log(gap4)/4,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, -np.log(gap5)/5,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, -np.log(gap6)/6,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, -np.log(gap7)/7,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, -np.log(gap8)/8,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, -np.log(gap9)/9,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, -np.log(gap10)/10,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('g')
plt.ylabel('α(g)')
plt.title('Funzione universale α(g) (g < 1, h = 0)')
plt.grid(True)
plt.legend()
plt.xlim(0.2,0.9)
plt.ylim(0.0,2.0)
plt.show()

#PLOTTING GAP ENERGETICO (E2-E1) AL VARIARE DI g (h = 0.00)
plt.errorbar(g, ggap3,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, ggap4,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, ggap5,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, ggap6,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, ggap7,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, ggap8,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, ggap9,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, ggap10,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('g')
plt.ylabel('ΔE')
plt.title('E2-E0 (h=0.0)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTING MAGNETIZZAZIONE TRASVERSA PER UN SITO AL VARIARE DI h (g = 0.90)
plt.errorbar(h, mx3g,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mx4g,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mx5g,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mx6g,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mx7g,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mx8g,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mx9g,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mx10g,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('h')
plt.ylabel('Mx')
plt.title('Mx (g=0.9)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()


#PLOTTING MAGNETIZZAZIONE LONGITUDINALE PER UN SITO AL VARIARE DI h (g = 0.90)
plt.errorbar(h, mz3g,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mz4g,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mz5g,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mz6g,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mz7g,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mz8g,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mz9g,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, mz10g,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('h')
plt.ylabel('Mz')
plt.title('Mz (g=0.9)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTING GAP ENERGETICO (E1-E0) AL VARIARE DI h (g = 0.90)
plt.errorbar(h, gap3g,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, gap4g,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, gap5g,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, gap6g,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, gap7g,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, gap8g,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, gap9g,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, gap10g,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('h')
plt.ylabel('ΔE')
plt.title('ΔE (g=0.9)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)  
plt.show()

#PLOTTING GAP ENERGETICO (E2-E1) AL VARIARE DI h (g = 0.90)
plt.errorbar(h, ggap3g,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, ggap4g,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, ggap5g,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, ggap6g,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, ggap7g,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, ggap8g,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, ggap9g,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(h, ggap10g,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('h')
plt.ylabel('delta(E)')
plt.title('E2-E1 (g=0.9)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()


#FSS for Mz (h=0.0)
#(for pasting: βγνχΔ)
b = 1./8.
nu = 1.
gc = 1.0
gamma = 7./4.

x3 = (g-gc)*(3**(1/nu))
y3 = mz3*(3**(b/nu))
x4 = (g-gc)*(4**(1/nu))
y4 = mz4*(4**(b/nu))
x5 = (g-gc)*(5**(1/nu))
y5 = mz5*(5**(b/nu))
x6 = (g-gc)*(6**(1/nu))
y6 = mz6*(6**(b/nu))
x7 = (g-gc)*(7**(1/nu))
y7 = mz7*(7**(b/nu))
x8 = (g-gc)*(8**(1/nu))
y8 = mz8*(8**(b/nu))
x9 = (g-gc)*(9**(1/nu))
y9 = mz9*(9**(b/nu))
x10 = (g-gc)*(10**(1/nu))
y10 = mz10*(10**(b/nu))

#PLOTTINGS 
plt.errorbar(x3, y3, fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x4, y4, fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x5, y5, fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x6, y6, fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x7, y7, fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x8, y8, fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x9, y9, fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x10, y10, fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('(g-gc)N^(1/ν)')
plt.ylabel('Mz N^(β/ν)')
plt.title('Finite Size Scaling per Mz (h=0.0)')
plt.grid(True)
plt.legend()
#plt.xlim(-5000,5000)
plt.show()

#PLOTTING  FFS GAP ENERGETICO (E1-E0) AL VARIARE DI g (h = 0.00)
plt.errorbar(x3, gap3*3,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x4, gap4*4,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x5, gap5*5,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x6, gap6*6,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x7, gap7*7,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x8, gap8*8,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x9, gap9*9,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(x10, gap10*10,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('(g-gc)N^(1/ν)')
plt.ylabel('ΔE N^(-z)')
plt.title('Finite size scaling per ΔE(h=0.0)')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

#DINAMICA

#QUENCH (gQuench = 1.0, h = 0.5)
#Load data from text file
dinamica3 = np.loadtxt("./dinamica/dinamica3.txt")
dinamica4 = np.loadtxt("./dinamica/dinamica4.txt")
dinamica5 = np.loadtxt("./dinamica/dinamica5.txt")
dinamica6 = np.loadtxt("./dinamica/dinamica6.txt")
dinamica7 = np.loadtxt("./dinamica/dinamica7.txt")
dinamica8 = np.loadtxt("./dinamica/dinamica8.txt")
dinamica9 = np.loadtxt("./dinamica/dinamica9.txt")
dinamica10 = np.loadtxt("./dinamica/dinamica10.txt")


#Sample data with given frequency
w = 1  # (Modifica questo valore a tuo piacimento)
sample3 = dinamica3[::w]
sample4 = dinamica4[::w]
sample5 = dinamica5[::w]
sample6 = dinamica6[::w]
sample7 = dinamica7[::w]
sample8 = dinamica8[::w]
sample9 = dinamica9[::w]
sample10 = dinamica10[::w]

#Name columns
t = sample9[:,0]
mzt3 = sample3[:,1]
mzt4 = sample4[:,1]
mzt5 = sample5[:,1]
mzt6 = sample6[:,1]
mzt7 = sample7[:,1]
mzt8 = sample8[:,1]
mzt9 = sample9[:,1]
mzt10 = sample10[:,1]

#Plottings
#plt.errorbar(t, mzt3,fmt=' ', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
#plt.errorbar(t, mzt4,fmt=' ', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
#plt.errorbar(t, mzt5,fmt=' ', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt6,fmt=' ', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
#plt.errorbar(t, mzt7,fmt=' ', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt8,fmt=' ', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt9,fmt=' ', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt10,fmt=' ', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('t')
plt.ylabel('Mz(t)')
plt.title('Quench g = 1.0')
plt.grid(True)
plt.legend()
#plt.ylim(-0.25,0.2)
#plt.xlim(0.0, 200.0)
plt.show()


#RISCALING TIME
#(for pasting: βγνχΔ)
b = 1./8.
nu = 1.
gc = 1.0
gamma = 7./4.

t3 = t/3
t4 = t/4
t5 = t/5
t6 = t/6
t7 = t/7
t8 = t/8
t9 = t/9
t10 = t/10

mzt3 *= 3**(b/nu)
mzt4 *= 4**(b/nu)
mzt5 *= 5**(b/nu)
mzt6 *= 6**(b/nu)
mzt7 *= 7**(b/nu)
mzt8 *= 8**(b/nu)
mzt9 *= 9**(b/nu)
mzt10 *= 10**(b/nu)

#Plottings
#plt.errorbar(t3, mzt3,fmt=' ', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
#plt.errorbar(t4, mzt4,fmt=' ', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
#plt.errorbar(t5, mzt5,fmt=' ', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t6, mzt6,fmt=' ', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
#plt.errorbar(t7, mzt7,fmt=' ', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t8, mzt8,fmt=' ', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t9, mzt9,fmt=' ', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t10, mzt10,fmt=' ', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('t N^(-1)')
plt.ylabel('Mz N^(β/ν)')
plt.title('Dinamica di Mz per tempi riscalati')
plt.grid(True)
plt.legend()
#plt.ylim(-0.25,0.2)
#plt.xlim(0.0, 35.0)
plt.show()
'''

#QUENCH (hQuench = 0.01, g = 1.0)
#Load data from text file
dinamica3 = np.loadtxt("./dinamica/hdinamica3.txt")
dinamica4 = np.loadtxt("./dinamica/hdinamica4.txt")
dinamica5 = np.loadtxt("./dinamica/hdinamica5.txt")
dinamica6 = np.loadtxt("./dinamica/hdinamica6.txt")
dinamica7 = np.loadtxt("./dinamica/hdinamica7.txt")
dinamica8 = np.loadtxt("./dinamica/hdinamica8.txt")
dinamica9 = np.loadtxt("./dinamica/hdinamica9.txt")
dinamica10 = np.loadtxt("./dinamica/hdinamica10.txt")


#Sample data with given frequency
w = 1  # (Modifica questo valore a tuo piacimento)
sample3 = dinamica3[::w]
sample4 = dinamica4[::w]
sample5 = dinamica5[::w]
sample6 = dinamica6[::w]
sample7 = dinamica7[::w]
sample8 = dinamica8[::w]
sample9 = dinamica9[::w]
sample10 = dinamica10[::w]

#Name columns
t = sample9[:,0]
mzt3 = sample3[:,1]
mzt4 = sample4[:,1]
mzt5 = sample5[:,1]
mzt6 = sample6[:,1]
mzt7 = sample7[:,1]
mzt8 = sample8[:,1]
mzt9 = sample9[:,1]
mzt10 = sample10[:,1]

#Plottings
plt.errorbar(t, mzt3,fmt=' ', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt4,fmt=' ', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt5,fmt=' ', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt6,fmt=' ', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt7,fmt=' ', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt8,fmt=' ', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt9,fmt=' ', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t, mzt10,fmt=' ', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('t')
plt.ylabel('Mz(t)')
plt.title('Quench su h (g = 1)')
plt.grid(True)
plt.legend(loc = 'lower right')
#plt.ylim(-0.25,0.2)
#plt.xlim(0.0, 200.0)
plt.show()


#RISCALING TIME
#(for pasting: βγνχΔ)
b = 1./8.
nu = 1.
gc = 1.0
gamma = 7./4.

t3 = t/3
t4 = t/4
t5 = t/5
t6 = t/6
t7 = t/7
t8 = t/8
t9 = t/9
t10 = t/10

mzt3 *= 3**(b/nu)
mzt4 *= 4**(b/nu)
mzt5 *= 5**(b/nu)
mzt6 *= 6**(b/nu)
mzt7 *= 7**(b/nu)
mzt8 *= 8**(b/nu)
mzt9 *= 9**(b/nu)
mzt10 *= 10**(b/nu)

#Plottings
plt.errorbar(t3, mzt3,fmt=' ', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t4, mzt4,fmt=' ', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t5, mzt5,fmt=' ', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t6, mzt6,fmt=' ', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t7, mzt7,fmt=' ', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t8, mzt8,fmt=' ', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t9, mzt9,fmt=' ', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t10, mzt10,fmt=' ', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('t N^(-1)')
plt.ylabel('Mz N^(β/ν)')
plt.title('Dinamica di Mz per tempi riscalati')
plt.grid(True)
plt.legend()
#plt.ylim(-0.25,0.2)
#plt.xlim(0.0, 35.0)
plt.show()

'''
#SUSC AS PARTIAL DERIVATIVE
#Load data from text file
suscTrasv3 = np.loadtxt("./dinamica/suscTrasv3.txt")
suscTrasv4 = np.loadtxt("./dinamica/suscTrasv4.txt")
suscTrasv5 = np.loadtxt("./dinamica/suscTrasv5.txt")
suscTrasv6 = np.loadtxt("./dinamica/suscTrasv6.txt")
suscTrasv7 = np.loadtxt("./dinamica/suscTrasv7.txt")
suscTrasv8 = np.loadtxt("./dinamica/suscTrasv8.txt")
suscTrasv9 = np.loadtxt("./dinamica/suscTrasv9.txt")
suscTrasv10 = np.loadtxt("./dinamica/suscTrasv10.txt")

#Name columns
g = suscTrasv3[:,0]
cv3 = suscTrasv3[:,1]*g
cv4 = suscTrasv4[:,1]*g
cv5 = suscTrasv5[:,1]*g
cv6 = suscTrasv6[:,1]*g
cv7 = suscTrasv7[:,1]*g
cv8 = suscTrasv8[:,1]*g
cv9 = suscTrasv9[:,1]*g
cv10 = suscTrasv10[:,1]*g

#PLOTTING CALORE SPECIFICO IN FUNZIONE DI g (h = 0.00)
plt.errorbar(g, cv3,fmt='.', label = 'N=3', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, cv4,fmt='.', label = 'N=4', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, cv5,fmt='.', label = 'N=5', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, cv6,fmt='.', label = 'N=6', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, cv7,fmt='.', label = 'N=7', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, cv8,fmt='.', label = 'N=8', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, cv9,fmt='.', label = 'N=9', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(g, cv10,fmt='.', label = 'N=10', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('g')
plt.ylabel('Cv')
plt.title('Calore specifico (h = 0.0)')
plt.grid(True)
plt.legend()
#plt.ylim(-0.25,0.2)
#plt.xlim(0.0, 7.0)
plt.show()
'''
