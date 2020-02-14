#partial derivatives
import numpy as np
import numdifftools as nd

n = 3
a = np.arange(3*n+0)

print(a[:n*3].reshape((n,3)))
print(a[3*n:])