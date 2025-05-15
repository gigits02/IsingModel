import numpy as np
import matplotlib.pyplot as plt
t = np.ndarray((24), dtype = np.double)
x = np.ndarray((24), dtype = np.double)
t, x = np.loadtxt("dati.txt", unpack = True)
plt.plot(t,x, color = 'r')
plt.show()