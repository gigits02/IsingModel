import numpy as np
from pylanczos import PyLanczos

#Use of SciPy sparse matrix is recommended to take full advantage of Lanczos algorithm.

matrix = np.array([[1., 2., 3.],
                    [4., 5., 6.],
                    [7., 8., 9.]])


engine = PyLanczos(matrix, True, 2)  # Find 2 maximum eigenpairs
eigenvalues, eigenvectors = engine.run()
print('Eigenvalues:')
print(eigenvalues)
print('Eigenvectors:')
print(eigenvectors)
