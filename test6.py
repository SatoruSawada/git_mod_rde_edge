import numpy as np
import matplotlib.pylab as plt


### (x0, y0) = (angle_fm,0)
### (x1, y1) = (angle_sl,num_ch_up-1)

y0 = 0.
angle_fm = 10.
y1 = 360.
angle_sl = 40.
num_ch_up = 100.

### (angle_fm ~ angle_sl) ==>> (0 ~ 360)
ax = np.linspace(angle_fm,angle_sl, int(num_ch_up))
# Ca = (y0-y1) / (angle_fm-angle_sl)
# Cb = y0 - Ca * angle_fm
# a = Ca * ax + Cb

C1 = abs(angle_sl-angle_fm) / (np.pi/2.)        ### Cy に対して
C2 = 0.9 ###固定
Cx = (angle_fm+angle_sl)/2./360.*2.*np.pi       ### 可変
Cy = np.pi       ### Cx に対して
eps = 10e-10
delta = 1.0

def func_delta_arctan(Cx):
    Cy_fm = angle_fm - C1 * np.arctan(C2*angle_fm/360.*2.*np.pi-Cx)
    Cy_sl = angle_sl - C1 * np.arctan(C2*angle_sl/360.*2.*np.pi-Cx)
    delta = abs((Cy_sl-Cy_fm)/Cy_fm)
    # C1_new = (angle_fm - Cy) / np.arctan(angle_fm/360.*2.*np.pi-Cx)
    # print('angle_fm/360.*2.*np.pi-Cx', angle_fm/360.*2.*np.pi-Cx)
    # print('delta =', delta, 'C1', C1, 'C1_new', C1_new)
    # C1 = C1_new
    print(delta, ' /// ', C1)
    return delta, Cy_sl

C1_a = 0.5
C1_b = 0.4
dC1_a,Cy_sl = func_delta_arctan(C1_a)
dC1_b,Cy_sl = func_delta_arctan(C1_b)
while abs(dC1_b) > eps:
    C1_s = (C1_a * dC1_b - C1_b * dC1_a)/(dC1_b - dC1_a)
    C1_a, C1_b = C1_b, C1_s
    dC1_a = dC1_b
    dC1_b,Cy_sl = func_delta_arctan(C1_b)

C1 = C1_b
Cy = Cy_sl
a = (C1*np.arctan(ax/360.*2.*np.pi-Cx) + Cy) / 2. / np.pi * 360.
print(a)

plt.plot(ax,a)
plt.show()



