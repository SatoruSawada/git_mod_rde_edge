import numpy as np
import matplotlib.pylab as plt

num_ch_up = 30
angle_fm = 10. / 360. * 2. * np.pi
angle_sl = 40. / 360. * 2. * np.pi
height_dw=1

C2 = np.tan(angle_sl-angle_fm)/num_ch_up
array_sample_up = np.zeros(int(num_ch_up))
array_sample_up[0] = angle_fm
for i in range(1,int(num_ch_up)):
    array_sample_up[i] = np.arctan(C2/height_dw + np.tan(array_sample_up[i-1]))
plt.plot(array_sample_up)
plt.show()
