import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Dati di esempio
x_data = np.array([1, 2, 3, 4, 5])
y_data = np.array([2.5, 3.5, 4.5, 5.5, 6.5])

# Definizione di una funzione lineare (y = mx + b)
def linear_func(x, m, b):
    return m * x + b

# Fit lineare
popt_linear, pcov_linear = curve_fit(linear_func, x_data, y_data)

# Generazione dei punti per i plot delle curve
x_fit = np.linspace(min(x_data), max(x_data), 100)

# Plot dei dati e delle curve adattate
plt.scatter(x_data, y_data, label='Dati')
plt.plot(x_fit, linear_func(x_fit, *popt_linear), 'r--', label='Fit lineare: {:.2f}x + {:.2f}'.format(*popt_linear))
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()