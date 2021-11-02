"""
目的：
gas dynamics, Vol.2 p. 200 example 17.1 を解いて
RDE with the method of characteristics の標準特性線のバリデーションを行う

問題：
１回目の corrector のlambda_o が 0.309... なのに 0.308... 
"""

## 意味ないけれども
## (__init__) -> 書いてもいいところ
## (__call__) -> 書いてはいけないところ


import os
import time
# import math
from matplotlib import colors, interactive
import numpy as np
import matplotlib.pylab as plt
from numpy.core.defchararray import array
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)


#### ================================================
#### 関数
#### ================================================


### 同じやん，いつか統一して
### こいつが一般版
def func_cross_gas_dynamics(point1, point2, lambda1, lambda2): ### ごめん，くそみたいな名前で
    x = (-(point1[1]-point2[1])+(lambda1*point1[0]-lambda2*point2[0]))/(lambda1-lambda2)
    y = point1[1] - lambda1 * point1[0] + lambda1 * x
    return x, y

### よく理解していない
## とりあえず vol.2 p.203 の例題にしたがって書いてみる
def func_MEPC_theta3_mode1(theta1, theta2, point1, point2, point4, lambda_12, lambda_o, eps=10e-10):
    delta = 1.0
    while abs(delta) >= eps:
        x, y = func_cross_gas_dynamics(point2, point4, lambda_12,lambda_o)
        theta = theta2 + (y-point2[1])/(point1[1]-point2[1]) * (theta1-theta2)
        lambda_o_new = np.tan(theta)
        delta = ((lambda_o-lambda_o_new)/lambda_o)
        lambda_o = lambda_o_new
        # print('============================================ delta =',delta)
        print('============================================ (x, y) =',(x,y))
        print("check==========================", theta*360./2./np.pi)
        print("check==========================", lambda_o)
    return x, y, theta, lambda_o

def func_MEPC_theta3_mode2(theta1, theta2,theta4, point1, point2, point4, lambda_12, lambda_o, eps=10e-10):
    delta = 1.0
    while abs(delta) >= eps:
        x, y = func_cross_gas_dynamics(point2, point4, lambda_12,lambda_o)
        theta = theta2 + (y-point2[1])/(point1[1]-point2[1]) * (theta1-theta2)
        lambda_o_new = np.tan((theta+theta4)/2.)
        delta = ((lambda_o-lambda_o_new)/lambda_o)
        lambda_o = lambda_o_new
        # print('============================================ delta =',delta)
        print('============================================ (x, y) =',(x,y))
        print("check==========================", theta*360./2./np.pi)
        print("check==========================", lambda_o)
    return x, y, theta, lambda_o

# def func_MEPC_theta3(theta1, theta2, point1, point2, point4, lambda_12, lambda_o, eps=10e-10):

#     def func_delta_lambda_o(lambda_o):
#         x, y = func_cross_gas_dynamics(point2, point4, lambda_12,lambda_o)
#         theta = theta2 + (y-point2[1])/(point1[1]-point2[1]) * (theta1-theta2)
#         lambda_o_new = np.tan(theta)
#         delta = ((lambda_o-lambda_o_new)/lambda_o)
#         return delta, x, y, theta

#     lambda_o_a = lambda_o
#     lambda_o_b = lambda_o * 1.1
#     eps_lambda_o = 10e-6
#     dlambda_o_a, x, y, theta = func_delta_lambda_o(lambda_o_a)
#     dlambda_o_b, x, y, theta = func_delta_lambda_o(lambda_o_b)
#     while abs(dlambda_o_b) > eps_lambda_o:
#         lambda_o_s = (lambda_o_a * dlambda_o_b - lambda_o_b * dlambda_o_a)/(dlambda_o_b - dlambda_o_a)
#         lambda_o_a, lambda_o_b = lambda_o_b, lambda_o_s
#         dlambda_o_a = dlambda_o_b
#         dlambda_o_b, x, y, theta = func_delta_lambda_o(lambda_o_b)
#     return x, y, theta, lambda_o_b


angle_bottom = 0. * 2. * np.pi /360.

#------------------------------------------------------------------
#### 仮定：
#------------------------------------------------------------------
### cantera_setting
import cantera as ct
from sdtoolbox.thermo import soundspeed_fr
mech = 'gri30_highT.cti'
gas = ct.Solution(mech)
ER = 1.
coef_C2H4 = 1
coef_O2 = 3
mol_ratio = coef_C2H4/coef_O2*ER
fuel = 'C2H4'
oxygen = 'O2'
q = fuel + ':' + str(mol_ratio) + ', ' + oxygen + ':1.0';


#------------------------------------------------------------------
#### 1. characteristic lines -1st
#------------------------------------------------------------------
### num_chの謎の発散の限界30くらい？
### 理由は分からぬ
num_ch_up = 5 # number of initial characteristic lines (upper side)
num_ch_down = 5 # number of initial characteristic lines (down side)

### i方向（横）にtheta-neu=const.確認
### j方向（縦）にtheta+neu=const.確認

array_zero0 = np.zeros((int(num_ch_down),int(num_ch_up-1)))
array_zero1 = np.zeros((int(num_ch_up + num_ch_down - 1),int(num_ch_up)))


### x for characteristics
array_x = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))

### y for characteristics
array_y = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))


### ============================================================================================ 20211018_sawada
### ============================================================================================ 初期値の設定がよくわからない
### ============================================================================================ とりあえずリーマン不変量で
### あくまで格子点における値
### 特性線とは違うし，あくまで区別してくれ
### theta
### neu, M, alpha
array_theta = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_neu = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_M = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_alpha = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_p = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_t = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_R = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_rho = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_a_fr = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_V = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_gamma = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
#============================================================================


### あくまで特性線における値
### 格子点とは違うし，あくまで区別してくれ
### =====

### C+ 上のパラメーター？？？
array_T_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_p_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_theta_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_V_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_rho_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_y_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_a_fr_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_M_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_alpha_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_lambda_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_Q_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_S_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))


### C- 上のパラメーター？？？
array_T_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_p_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_theta_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_V_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_rho_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_y_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_a_fr_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_M_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_alpha_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_lambda_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_Q_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_S_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))


### 真ん中の流線？？？
array_x_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_y_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_theta_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_p_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_rho_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_V_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_a_fr_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_gamma_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))

array_lambda_12 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_lambda_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_p_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_rho_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_V_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_a_fr_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_R_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_A_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_T_o1 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))
array_T_o2 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up*2)))



#============================================================================
### 値を突っ込んでいく
#============================================================================
### up (point1)
### P_b, T, R, rho, a_fr, V, gamma
array_x[num_ch_up-2][1] = 0.131460
array_y[num_ch_up-2][1] = 0.040118
array_theta[num_ch_up-2][1] = 18.191 * 2. * np.pi / 360.
array_V[num_ch_up-2][1] = 2603.5
array_p[num_ch_up-2][1] = 34042
array_rho[num_ch_up-2][1] = 0.086151
### down
### P_b, T, R, rho, a_fr, V, gamma
array_x[num_ch_up][0] = 0.135683
array_y[num_ch_up][0] = 0.037123
array_theta[num_ch_up][0] = 16.422 * 2. * np.pi / 360.
array_V[num_ch_up][0] = 2609.2
array_p[num_ch_up][0] = 32781
array_rho[num_ch_up][0] = 0.083482
#============================================================================
#============================================================================

### (a)
### state1
# gas.DP = array_rho[num_ch_up-2][1], array_p[num_ch_up-2][1]
# array_a_fr[num_ch_up-2][1] = soundspeed_fr(gas)
array_a_fr[num_ch_up-2][1] = np.sqrt((1.2*array_p[num_ch_up-2][1])/(array_rho[num_ch_up-2][1])) ### gamma = 1.2
array_M[num_ch_up-2][1] = array_V[num_ch_up-2][1] / array_a_fr[num_ch_up-2][1]
array_alpha[num_ch_up-2][1] = np.arcsin(1./array_M[num_ch_up-2][1])
### state2
# gas.DP = array_rho[num_ch_up][0], array_p[num_ch_up][0]
# array_a_fr[num_ch_up][0] = soundspeed_fr(gas)
array_a_fr[num_ch_up][0] = np.sqrt((1.2*array_p[num_ch_up][0])/(array_rho[num_ch_up][0])) ### gamma = 1.2
array_M[num_ch_up][0] = array_V[num_ch_up][0] / array_a_fr[num_ch_up][0]
array_alpha[num_ch_up][0] = np.arcsin(1./array_M[num_ch_up][0])

print("===============================================================")
print("flow properties at the initial-value points")
print("===============================================================")
print("a_fr_1      ", array_a_fr[num_ch_up-2][1])
print("M_1         ", array_M[num_ch_up-2][1])
print("alpha_1     ", array_alpha[num_ch_up-2][1]/2./np.pi*360.)
print("a_fr_2      ", array_a_fr[num_ch_up][0])
print("M_2         ", array_M[num_ch_up][0])
print("alpha_2     ", array_alpha[num_ch_up][0]/2./np.pi*360.)




### ============================

### 今まで数値の格納方法を踏襲

### C+ : [j-1][i]
### 24 : [j-1][i]

### C- : [j+1][i-1]
### 12 : [j+1][i-1]

### Co : [j][i-1]
###  3 : [j][i-1]

### ============================




### 対象は theta_3 ???
### predictor : 一つ目の比較対象の値
### corrector : 二つ目の比較対象の値
### itration  : abs(predictor - corrector) が許容誤差未満に至るまで計算（absでなくてもいいかもしれないが）




### process
### 1. predictor 計算
### 2. 偏差 1.0，許容誤差 10e-6 で仮定
### 3. corrector の while 計算を開始
### 4. while 中に abs(theta_3_predictor - theta_3_corrector) を計算
### 5. abs(theta_3_predictor - theta_3_corrector) が 0 になるまで while 繰り返し計算 




S_add = 0.


### 必要になったら随時，パラメーターの行列を追加していく方針でお願いします．

for i in range(1,2):
    for j in range(int(num_ch_up-i), int((num_ch_up-i)+1)):

        ### =====================================================================
        ### predictor
        ### =====================================================================
        # print(array_theta[j+1][i-1]/2./np.pi*360., '///', array_alpha[j+1][i-1]/2./np.pi*360.)
        # print(array_theta[j-1][i]/2./np.pi*360., '///', array_alpha[j-1][i]/2./np.pi*360.)
        # print(np.tan(array_theta[j+1][i-1]+array_alpha[j+1][i-1]))
        # print(np.tan(array_theta[j-1][i]-array_alpha[j-1][i]))
        #####################################################################################################(c)
        ### lambda_plus & lambda_minus - eq17dot47_eq17dot48 (first step predictor)
        array_lambda_plus[j+1][i-1] = np.tan(array_theta[j+1][i-1]+array_alpha[j+1][i-1])
        array_lambda_minus[j-1][i] = np.tan(array_theta[j-1][i]-array_alpha[j-1][i])
        ### Q+ & Q- - eq17dot54_eq17dot55 (first step predictor) (noting delta is removed)
        array_Q_plus[j+1][i-1] = np.sqrt((array_M[j+1][i-1])**2.-1.) / (array_rho[j+1][i-1]*array_V[j+1][i-1]**2.)
        array_Q_minus[j-1][i] = np.sqrt((array_M[j-1][i])**2.-1.) / (array_rho[j-1][i]*array_V[j-1][i]**2.)
        ### S+ & S- - eq17dot54_eq17dot55 (first step predictor) (noting delta is removed)
        array_S_plus[j+1][i-1] = np.sin(array_theta[j+1][i-1]) / \
            ((array_y[j+1][i-1]+S_add) * array_M[j+1][i-1] * np.cos(array_theta[j+1][i-1]+array_alpha[j+1][i-1]))
        array_S_minus[j-1][i] = np.sin(array_theta[j-1][i]) / \
            ((array_y[j-1][i]+S_add) * array_M[j-1][i] * np.cos(array_theta[j-1][i]-array_alpha[j-1][i]))
        
        #####################################################################################################(d)
        ### eq17dot44_eq17dot45
        array_x[j][i], array_y[j][i] = func_cross_gas_dynamics(\
            (array_x[j-1][i], array_y[j-1][i]),\
                (array_x[j+1][i-1], array_y[j+1][i-1]),\
                    array_lambda_minus[j-1][i], \
                        array_lambda_plus[j+1][i-1])
        ### T+ & T- - eq17dot52_eq17dot53 (first step predictor)
        array_T_plus[j+1][i-1] = -array_S_plus[j+1][i-1] * (array_x[j][i]-array_x[j+1][i-1]) + \
            array_Q_plus[j+1][i-1] * array_p[j+1][i-1] + array_theta[j+1][i-1]
        array_T_minus[j-1][i] = -array_S_minus[j-1][i] * (array_x[j][i]-array_x[j-1][i]) + \
            array_Q_minus[j-1][i] * array_p[j-1][i] - array_theta[j-1][i]

        #####################################################################################################(e)
        ### eq17dot49 (the modified euler predictor-corrector)
        ### eq17dot43_eq17dot46
        array_lambda_12[j-1][i] = (array_y[j+1][i-1] - array_y[j-1][i]) / (array_x[j+1][i-1] - array_x[j-1][i])
        array_lambda_o[j][i-1] = np.tan((array_theta[j-1][i] + array_theta[j+1][i-1])/2.)
        array_x_3[j][i-1],\
            array_y_3[j][i-1],\
                array_theta_3[j][i-1],\
                    array_lambda_o[j][i-1] = func_MEPC_theta3_mode1(\
            array_theta[j-1][i],\
                array_theta[j+1][i-1],\
                    (array_x[j-1][i], array_y[j-1][i]),\
                        (array_x[j+1][i-1], array_y[j+1][i-1]),\
                            (array_x[j][i], array_y[j][i]),\
                                array_lambda_12[j-1][i],\
                                    array_lambda_o[j][i-1])
        ### interpolating for the remaining flow properties gives... (p.203) 
        array_p_3[j][i-1] = array_p[j+1][i-1]+(array_y_3[j][i-1]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_p[j-1][i]-array_p[j+1][i-1])
        array_rho_3[j][i-1] = array_rho[j+1][i-1]+(array_y_3[j][i-1]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_rho[j-1][i]-array_rho[j+1][i-1])
        array_V_3[j][i-1] = array_V[j+1][i-1]+(array_y_3[j][i-1]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_V[j-1][i]-array_V[j+1][i-1])

        #####################################################################################################(f)
        # gas.DP = array_rho_3[j][i], array_p_3[j][i] ### RDE -> gas.SPX = s2, array_p_3[j][i], x2
        # array_R_o[j][i] = ct.gas_constant / gas.mean_molecular_weight
        # array_a_fr_3[j][i] = soundspeed_fr(gas)
        array_a_fr_3[j][i-1] = np.sqrt((1.2*array_p_3[j][i-1])/array_rho_3[j][i-1])
        array_R_o[j][i-1] = array_rho_3[j][i-1] * array_V_3[j][i-1]
        array_A_o[j][i-1] = array_a_fr_3[j][i-1] ** 2.
        array_T_o1[j][i-1] = array_R_o[j][i-1] * array_V_3[j][i-1] + array_p_3[j][i-1] # using state3
        array_T_o2[j][i-1] = array_p_3[j][i-1] - array_A_o[j][i-1] * array_rho_3[j][i-1] # using state3
        ### eq.(g) & (h) for calculating p4, theta4
        array_p[j][i] = (array_T_plus[j+1][i-1] + array_T_minus[j-1][i]) / (array_Q_plus[j+1][i-1] + array_Q_minus[j-1][i])
        array_theta[j][i] = array_T_plus[j+1][i-1] - array_Q_plus[j+1][i-1] * array_p[j][i]
        ### eq.(i) & (j) for calculating V4, rho4
        array_V[j][i] = (array_T_o1[j][i-1]-array_p[j][i]) / array_R_o[j][i-1]
        array_rho[j][i] = (array_p[j][i]-array_T_o2[j][i-1]) / array_A_o[j][i-1]

        print("===============================================================")
        print("predictor")
        print("===============================================================")
        print("lambda_plus ", array_lambda_plus[j+1][i-1])
        print("lambda_minus", array_lambda_minus[j-1][i])
        print("lambda_o    ", array_lambda_o[j][i-1])
        print("x_4         ", array_x[j][i])
        print("y_4         ", array_y[j][i])
        print("x_3         ", array_x_3[j][i-1])
        print("y_3         ", array_y_3[j][i-1])
        print("R_o         ", array_R_o[j][i-1])
        print("A_o         ", array_A_o[j][i-1])
        print("T_o1        ", array_T_o1[j][i-1])
        print("T_o2        ", array_T_o2[j][i-1])
        print("Q_plus      ", array_Q_plus[j+1][i-1])
        print("S_plus      ", array_S_plus[j+1][i-1])
        print("T_plus      ", array_T_plus[j+1][i-1])
        print("Q_minus     ", array_Q_minus[j-1][i])
        print("S_minus     ", array_S_minus[j-1][i])
        print("T_minus     ", array_T_minus[j-1][i])
        print("p_4         ", array_p[j][i])
        print("theta_4     ", array_theta[j][i]/2./np.pi*360.)
        print("V_4         ", array_V[j][i])
        print("rho_4       ", array_rho[j][i])


        ### set predictor
        theta_3 = array_theta_3[j][i-1]
        delta_c = 1.0
        eps_c = 10e-6
        n = 0
        ### =====================================================================
        ### corrector : 全て入れなおせているのだろうか？
        ### =====================================================================
        while abs(delta_c) >= eps_c:
            #####################################################################################################(g)
            ### along Mach line 24 (C+)
            array_p_plus[j+1][i-1] = (array_p[j+1][i-1] + array_p[j][i]) /2.
            array_theta_plus[j+1][i-1] = (array_theta[j+1][i-1] + array_theta[j][i]) /2.
            array_V_plus[j+1][i-1] = (array_V[j+1][i-1] + array_V[j][i]) /2.
            array_rho_plus[j+1][i-1] = (array_rho[j+1][i-1] + array_rho[j][i]) /2.
            array_y_plus[j+1][i-1] = (array_y[j+1][i-1] + array_y[j][i]) /2.
            array_a_fr_plus[j+1][i-1] = np.sqrt((1.2*array_p_plus[j+1][i-1])/array_rho_plus[j+1][i-1]) ### gamma = 1.2
            array_M_plus[j+1][i-1] = array_V_plus[j+1][i-1] / array_a_fr_plus[j+1][i-1]
            array_alpha_plus[j+1][i-1] = np.arcsin(1./array_M_plus[j+1][i-1])
            array_lambda_plus[j+1][i-1] = np.tan(array_theta_plus[j+1][i-1]+array_alpha_plus[j+1][i-1])
            array_Q_plus[j+1][i-1] = np.sqrt(array_M_plus[j+1][i-1]**2.-1.) / (array_rho_plus[j+1][i-1]*array_V_plus[j+1][i-1]**2.)
            array_S_plus[j+1][i-1] = np.sin(array_theta_plus[j+1][i-1]) / \
                ((array_y_plus[j+1][i-1])*array_M_plus[j+1][i-1]*np.cos(array_theta[j][i]+array_theta_plus[j+1][i-1]))
            ### along Mach line 14 (C-)
            array_p_minus[j-1][i] = (array_p[j-1][i] + array_p[j][i]) /2.
            array_theta_minus[j-1][i] = (array_theta[j-1][i] + array_theta[j][i]) /2.
            array_V_minus[j-1][i] = (array_V[j-1][i] + array_V[j][i]) /2.
            array_rho_minus[j-1][i] = (array_rho[j-1][i] + array_rho[j][i]) /2.
            array_y_minus[j-1][i] = (array_y[j-1][i] + array_y[j][i]) /2.
            array_a_fr_minus[j-1][i] = np.sqrt((1.2*array_p_minus[j-1][i])/array_rho_minus[j-1][i]) ### gamma = 1.2
            array_M_minus[j-1][i] = array_V_minus[j-1][i] / array_a_fr_minus[j-1][i]
            array_alpha_minus[j-1][i] = np.arcsin(1./array_M_minus[j-1][i])
            array_lambda_minus[j-1][i] = np.tan(array_theta_minus[j-1][i]-array_alpha_minus[j-1][i])
            array_Q_minus[j-1][i] = np.sqrt(array_M_minus[j-1][i]**2.-1.) / (array_rho_minus[j-1][i]*array_V_minus[j-1][i]**2.)
            array_S_minus[j-1][i] = np.sin(array_theta_minus[j-1][i]) / \
                ((array_y_minus[j-1][i])*array_M_minus[j-1][i]*np.cos(array_theta[j][i]-array_theta_minus[j-1][i]))

            #####################################################################################################(h)
            ### eq17dot44_eq17dot45
            array_x[j][i], array_y[j][i] = func_cross_gas_dynamics(\
                (array_x[j-1][i], array_y[j-1][i]),\
                    (array_x[j+1][i-1], array_y[j+1][i-1]),\
                        array_lambda_minus[j-1][i], \
                            array_lambda_plus[j+1][i-1])
            ### T+ & T- - eq17dot52_eq17dot53 (first step predictor
            array_T_plus[j+1][i-1] = -array_S_plus[j+1][i-1] * (array_x[j][i]-array_x[j+1][i-1]) + \
                array_Q_plus[j+1][i-1] * array_p[j+1][i-1] + array_theta[j+1][i-1]
            array_T_minus[j-1][i] = -array_S_minus[j-1][i] * (array_x[j][i]-array_x[j-1][i]) + \
                array_Q_minus[j-1][i] * array_p[j-1][i] - array_theta[j-1][i]

            #####################################################################################################(i)
            ### eq17dot49 (the modified euler predictor-corrector)
            ### eq17dot43_eq17dot46
            # array_lambda_12[j-1][i] = (array_y[j+1][i-1] - array_y[j-1][i]) / (array_x[j+1][i-1] - array_x[j-1][i])
            ### (theta3, theta4) -> corrector
            array_lambda_o[j][i-1] = np.tan((array_theta_3[j][i-1]+array_theta[j][i])/2.)
            print("check==========================", array_theta[j-1][i]*360./2./np.pi)
            print("check==========================", array_theta[j+1][i-1]*360./2./np.pi)
            print("check==========================", (array_x[j-1][i], array_y[j-1][i]))
            print("check==========================", (array_x[j+1][i-1], array_y[j+1][i-1]))
            print("check==========================", (array_x[j][i], array_y[j][i]))
            print("check==========================", array_lambda_12[j-1][i])
            print("check==========================", array_lambda_o[j][i-1])

            array_x_3[j][i-1],\
                array_y_3[j][i-1],\
                    array_theta_3[j][i-1],\
                        array_lambda_o[j][i-1] = func_MEPC_theta3_mode2(\
                array_theta[j-1][i],\
                    array_theta[j+1][i-1],\
                        array_theta[j][i],\
                        (array_x[j-1][i], array_y[j-1][i]),\
                            (array_x[j+1][i-1], array_y[j+1][i-1]),\
                                (array_x[j][i], array_y[j][i]),\
                                    array_lambda_12[j-1][i],\
                                        array_lambda_o[j][i-1])
            print("check==========================", (array_x_3[j][i-1], array_y_3[j][i-1]))
            ### interpolating for the remaining flow properties gives... (p.203) 
            array_p_3[j][i-1] = array_p[j+1][i-1]+(array_y_3[j][i-1]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_p[j-1][i]-array_p[j+1][i-1])
            array_rho_3[j][i-1] = array_rho[j+1][i-1]+(array_y_3[j][i-1]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_rho[j-1][i]-array_rho[j+1][i-1])
            array_V_3[j][i-1] = array_V[j+1][i-1]+(array_y_3[j][i-1]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_V[j-1][i]-array_V[j+1][i-1])

            #####################################################################################################(j)
            ### (start) ここで計算が「predictor」と「corrector」で異なる
            array_p_o[j][i-1] = (array_p_3[j][i-1] + array_p[j][i]) / 2.
            array_rho_o[j][i-1] = (array_rho_3[j][i-1] + array_rho[j][i]) / 2.
            array_V_o[j][i-1] = (array_V_3[j][i-1] + array_V[j][i]) / 2.
            ### (end) 
            # gas.DP = array_rho_o[j][i], array_p_o[j][i]
            # array_R_o[j][i] = ct.gas_constant / gas.mean_molecular_weight
            # array_a_fr_3[j][i] = soundspeed_fr(gas)
            array_a_fr_3[j][i-1] = np.sqrt((1.2*array_p_o[j][i-1])/array_rho_o[j][i-1]) ### gamma = 1.2
            array_R_o[j][i-1] = array_rho_o[j][i-1] * array_V_o[j][i-1]
            array_A_o[j][i-1] = array_a_fr_3[j][i-1] ** 2.
            array_T_o1[j][i-1] = array_R_o[j][i-1] * array_V_3[j][i-1] + array_p_3[j][i-1] # using state3
            array_T_o2[j][i-1] = array_p_3[j][i-1] - array_A_o[j][i-1] * array_rho_3[j][i-1] # using state3

            ### eq.(g) & (h) for calculating p4, theta4
            array_p[j][i] = (array_T_plus[j+1][i-1] + array_T_minus[j-1][i]) / (array_Q_plus[j+1][i-1] + array_Q_minus[j-1][i])
            array_theta[j][i] = array_T_plus[j+1][i-1] - array_Q_plus[j+1][i-1] * array_p[j][i]

            ### eq.(i) & (j) for calculating V4, rho4
            array_V[j][i] = (array_T_o1[j][i-1]-array_p[j][i]) / array_R_o[j][i-1]
            array_rho[j][i] = (array_p[j][i]-array_T_o2[j][i-1]) / array_A_o[j][i-1]

            ### delta_c
            theta_3_new = array_theta_3[j][i-1]
            delta_c = ((theta_3-theta_3_new)/theta_3)
            theta_3 = theta_3_new

            ### count
            n += 1
            print(n)
            

            print("===============================================================")
            print("corrector")
            print("===============================================================")
            print("lambda_plus ", array_lambda_plus[j+1][i-1])
            print("lambda_minus", array_lambda_minus[j-1][i])
            print("lambda_o    ", array_lambda_o[j][i-1])
            print("x_4         ", array_x[j][i])
            print("y_4         ", array_y[j][i])
            print("x_3         ", array_x_3[j][i-1])
            print("y_3         ", array_y_3[j][i-1])
            print("R_o         ", array_R_o[j][i-1])
            print("A_o         ", array_A_o[j][i-1])
            print("T_o1        ", array_T_o1[j][i-1])
            print("T_o2        ", array_T_o2[j][i-1])
            print("Q_plus      ", array_Q_plus[j+1][i-1])
            print("S_plus      ", array_S_plus[j+1][i-1])
            print("T_plus      ", array_T_plus[j+1][i-1])
            print("Q_minus     ", array_Q_minus[j-1][i])
            print("S_minus     ", array_S_minus[j-1][i])
            print("T_minus     ", array_T_minus[j-1][i])
            print("p_4         ", array_p[j][i])
            print("theta_4     ", array_theta[j][i]/2./np.pi*360.)
            print("V_4         ", array_V[j][i])
            print("rho_4       ", array_rho[j][i])



