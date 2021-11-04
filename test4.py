import numpy as np
import matplotlib.pylab as plt
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)


def func_delta_M(neu_target, M, gamma=1.3):
    part1 = ((gamma+1.)/(gamma-1.))**(1./2.)
    part2 = np.arctan(((gamma-1.)/(gamma+1.)*(M**2.-1.))**(1./2.))
    part3 = np.arctan((M**2.-1.)**(1./2.))
    neu_cal = part1*part2-part3
    delta_neu = neu_target - neu_cal
    return delta_neu


def func_neu2M(neu_target,eps = 10e-6):
    M_a = 1.2
    M_b = 1.4
    dM_a = func_delta_M(neu_target,M_a)
    dM_b = func_delta_M(neu_target,M_b)
    while abs(dM_b) > eps:
        M_s = (M_a * dM_b - M_b * dM_a)/(dM_b - dM_a)
        M_a, M_b = M_b, M_s
        dM_a = dM_b
        dM_b = func_delta_M(neu_target, M_b)
    M_result = M_b
    return M_result


num_ch_up = 30
angle_fm = 0. / 360. * 2. * np.pi
angle_sl = 40. / 360. * 2. * np.pi


### C0 * x ** 4. + C1
deg1_up = angle_fm
deg2_up = angle_sl
C1_up = deg1_up
C0_up = (deg2_up-C1_up)/((num_ch_up-1)**(1))
array_sample0 = np.arange(num_ch_up)
array_sample_up = array_sample0
# print("array_sample_up ============", array_sample_up * 360. / 2./ np.pi)
array_sample_up = array_sample_up * C0_up + C1_up
array_theta_up = array_sample_up
# print("array_theta_up ============", array_theta_up * 360. / 2./ np.pi)
array_neu_up = array_sample_up - angle_fm

array_neu_up = np.tan(array_theta_up+array_neu_up)

plt.plot(array_neu_up)
plt.show()
