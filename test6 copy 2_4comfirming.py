import numpy as np
import matplotlib.pylab as plt

### 逆にすればいいじゃないの
### (0 ~ 360) ==>> (angle_fm ~ angle_sl)
### (x0, y0) = (angle_fm,0)
### (x1, y1) = (angle_sl,num_ch_up-1)

### ============================================
C1 = 1.  ### Cy に対して
C2 = 0.7
Cx = 1.7  ### rad以外 可変
Cy = np.pi  ### rad Cx に対して
### ============================================

y0 = 0.
angle_fm = 10.
y1 = 360.
angle_sl = 40.
num_ch_up = 100.
ax = np.linspace(0,360,int(num_ch_up))

a = (C1*np.arctan(C2*ax/360.*2.*np.pi-Cx) + Cy)# / 2. / np.pi * 360.
### (0 ~ 360) ==>> (angle_fm ~ angle_sl)
Ca = (angle_fm- angle_sl) / (a[0]-a[-1])
Cb = angle_fm - Ca * a[0]
a = Ca * a + Cb

# print(a)
plt.plot(a)
plt.show()



