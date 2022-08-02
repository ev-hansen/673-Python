import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

list_x = [1, 2, 3, 5, 6, 8]
list_y = [3, 4, 6, 7, 8, 9]
list_x2 = [1.2, 4.2, 4.3, 6, 7]

xarr = np.array(list_x)
x2arr = np.array(list_x2)
yarr = np.array(list_y)

f = interpolate.interp1d(xarr, yarr)

y_new = f(x2arr)

plt.plot(xarr, yarr, '-', x2arr, y_new, 'o')

plt.show()
