import numpy as np
import matplotlib.pyplot as plt 
from scipy import interpolate

a = np.zeros(1000)
a[40] = 20
a[42] = 20
a[41] = 25
a[43:] = 5

x = np.arange(1000)

y = x+6

f = interpolate.interp1d(y,a,kind = 'cubic',fill_value = "extrapolate")

z = f(x)

plt.plot(y,a)
plt.show()
plt.plot(x,z)
plt.show()