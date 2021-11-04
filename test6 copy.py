import numpy as np
import matplotlib.pylab as plt


### (x0, y0) = (angle_fm,0)
### (x1, y1) = (angle_sl,num_ch_up-1)

num_ch_up = 360.

### (angle_fm ~ angle_sl) ==>> (0 ~ 360)
ax = np.linspace(0,360, int(num_ch_up))
C1 = 40. / (np.pi/2.)  ### Cy に対して
C2 = 1
Cx = 3.1                          ### 可変
Cy = 20.  ### Cx に対して

a = (C1*np.arctan(C2*ax/360.*2.*np.pi-Cx) + Cy)# / 2. / np.pi * 360.
print(a)

plt.plot(ax,a)
plt.show()



