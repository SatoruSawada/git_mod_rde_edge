import numpy as np
import matplotlib.pylab as plt

### 逆にすればいいじゃないの
### (0 ~ 360) ==>> (angle_fm ~ angle_sl)
### (x0, y0) = (angle_fm,0)
### (x1, y1) = (angle_sl,num_ch_up-1)

### 組み合わせる

### ============================================
C1 = 1.  ### Cy に対して
C2 = 0.7
Cx = 1.9  ### rad以外 可変
Cy = np.pi  ### rad Cx に対して
### ============================================
y0 = 0.
angle_fm = 10.
y1 = 360.
angle_sl = 40.
num_ch_up = 100.
num_ratio = 1./5.
ang_ratio = 10e-4


num_ch_up0 = int(num_ch_up*num_ratio)
num_ch_up1 = num_ch_up - num_ch_up0
a0 = np.arange(num_ch_up0)

### (0 ~ 360) ==>> (angle_fm ~ angle_sl)
Ca = (angle_fm - angle_fm*(1+ang_ratio)) / (a0[0]-a0[-1])
Cb = angle_fm - Ca * a0[0]
a0 = Ca * a0 + Cb

ax = np.linspace(angle_fm*(1+ang_ratio),360,int(num_ch_up1))

a = (C1*np.arctan(C2*ax/360.*2.*np.pi-Cx) + Cy)# / 2. / np.pi * 360.

### (0 ~ 360) ==>> (angle_fm ~ angle_sl)
Ca = (angle_fm*(1+ang_ratio) - angle_sl) / (a[0]-a[-1])
Cb = angle_fm*(1+ang_ratio) - Ca * a[0]
a = Ca * a + Cb

a0 = np.delete(a0,-1,0)

a = np.hstack((a0,a))

print(a)
plt.plot(a)
plt.show()



