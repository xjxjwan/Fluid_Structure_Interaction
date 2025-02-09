import numpy as np

a = np.ones((3, 3))
b = np.ones((3, 3))
b[0, 0] = 0
print(a)
print(b)
print(a * (b > 0))