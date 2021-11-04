import numpy as np
import matplotlib.pylab as plt

### 逆にすればいいじゃないの
### (0 ~ 360) ==>> (angle_fm ~ angle_sl)
### (x0, y0) = (angle_fm,0)
### (x1, y1) = (angle_sl,num_ch_up-1)

angle_fm = 10.
angle_sl = 40.

num_ch_up = 360
a = np.linspace(0, num_ch_up, int(num_ch_up))
### (0 ~ 360) ==>> (angle_fm ~ angle_sl)

Ca = (angle_fm- angle_sl) / (a[0]-a[-1])
Cb = angle_fm - Ca * a[0]
a = Ca * a + Cb

print(a)
plt.plot(a)
plt.show()



