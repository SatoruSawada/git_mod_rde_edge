import numpy as np
import matplotlib.pylab as plt


### (x0, y0) = (angle_fm,0)
### (x1, y1) = (angle_sl,num_ch_up-1)

num_ch_up = 36.
C = 10
h = 1.0

# a = np.arctan(C/h+np.tan())
a = np.zeros(int(num_ch_up))

for i in range(1,int(num_ch_up)):
    a[i] = np.arctan(C/h+np.tan(a[i-1]))

plt.plot(a)
plt.show()



