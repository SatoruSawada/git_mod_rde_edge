"""
目的：
傾きの初期値を「〇〇deg.毎にする」ように上下で統一する

方法：
1. class(__init__有)にグラフの骨組みを与える
2. class(__call__)に新しいグラフを加えていく
3. とりあえずMOC
4. 等エントロピー＆全エンタルピー仮定　ー＞　(s,h,M) 2 (p,T,v,rho,tau,Reynolds_num,etc.)
5. コンター図作成
6. なんとかして軸方向-(p,T,v,rho,tau,Reynolds_num,etc.)のグラフ出せないかな？
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

class CL_graph_setting:
    def __init__(self):
        #### ==================================================================================================================
        #### setting
        #### ==================================================================================================================
        #### x軸
        self.x_label = 'radial direction [-]'
        self.x_min = -2.        #### x軸最小値
        self.x_max = 8.        #### x軸最大値，目盛りの表示の都合でここに 0.00001 % 加算し申す
        self.x_main_dis = 2.   #### x軸主目盛り間隔
        # x_sub_num = 5       #### x軸主目盛り間の小目盛りの個数
        #### y軸
        self.y_label = 'azimuthal direction [-]'
        self.y_min = -1.        #### y軸最小値
        self.y_max = 4.        #### y軸最大値，目盛りの表示の都合でここに 0.00001 % 加算し申す
        self.y_main_dis = 1.   #### y軸主目盛り間隔
        # y_sub_num = 5       #### y軸主目盛り間の小目盛りの個数
        #### 軸の大きさ・太さ
        self.major_width = 1.5    #### 軸主目盛りの線幅
        self.major_size = 10      #### 軸主目盛りの長さ
        self.minor_width = 1.0    #### 軸補助目盛りの線幅
        self.minor_size = 5      #### 軸補助目盛りの長さ
        #### 書式
        self.letter_font = "Arial" #### 文字書式
        self.letter_size = 16      #### 文字サイズ
        self.number_font = "Arial" #### 数字書式
        self.number_size = 16      #### 数字サイズ
        #### marker
        self.marker_size = 1.
        self.marker_color = "black"
        #### 枠&図全体の設定
        self.fig_x_size = 16.     #### グラフ全体のx軸方向の大きさ
        self.fig_y_size = 8.     #### グラフ全体のy軸方向の大きさ
        self.axes_width = 1       #### 枠の太さ

        #### ================================================
        #### 原点を通過する軸の表示
        #### ================================================
        #### x軸方向
        self.fig_x_zero = True
        #### y 軸方向
        self.fig_y_zero = True

class CL_graph_create(CL_graph_setting):
    def __call__(self):
        #### ==================================================================================================================
        #### graph settings
        #### ==================================================================================================================
        ### font settings
        plt.rcParams["font.family"] = "Arial"      #全体のフォントを設定
        plt.rcParams["mathtext.fontset"] = "cm" 
        ### 目盛り外内側指定（外側(defolt) -> 'out'，内側 -> 'in'）
        plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
        plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
        ### 主目盛り
        plt.rcParams["xtick.major.width"] = self.major_width              #x軸主目盛り線の線幅
        plt.rcParams["ytick.major.width"] = self.major_width              #y軸主目盛り線の線幅
        plt.rcParams["xtick.major.size"] = self.major_size                #x軸主目盛り線の長さ
        plt.rcParams["ytick.major.size"] = self.major_size                #y軸主目盛り線の長さ
        ### 補助目盛り
        plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
        plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
        plt.rcParams["xtick.minor.width"] = self.minor_width              #x軸補助目盛り線の線幅
        plt.rcParams["ytick.minor.width"] = self.minor_width              #y軸補助目盛り線の線幅
        plt.rcParams["xtick.minor.size"] = self.minor_size                 #x軸補助目盛り線の長さ
        plt.rcParams["ytick.minor.size"] = self.minor_size                 #y軸補助目盛り線の長さ
        plt.rcParams["font.size"] = self.letter_size                       #フォントの大きさ
        plt.rcParams["axes.linewidth"] = self.axes_width                 #囲みの太さ
        ### 上，右側の目盛り表示
        plt.rcParams['xtick.top'] = True
        plt.rcParams['ytick.right'] = True

        #### ================================================
        #### set graph frame
        #### ================================================
        self.fig, self.ax = plt.subplots(1,1, figsize=(self.fig_x_size, self.fig_y_size)) #### 第1,2変数はfig,axの配列の数を指定
        #### x_axis
        x_range_max_array = np.arange(0.0, self.x_max*1.00000000001, self.x_main_dis)
        x_range_min_array = np.arange(0.0, abs(self.x_min)*1.00000000001, self.x_main_dis) * (-1)
        x_range_array = np.sort(np.unique(np.hstack([x_range_max_array, x_range_min_array])))
        self.ax.set_xticks(x_range_array)
        self.ax.set_xlim(self.x_min, self.x_max*1.00000000001)
        self.ax.set_xlabel(self.x_label)
        #### y_axis
        y_range_max_array = np.arange(0.0, self.y_max*1.00000000001, self.y_main_dis)
        y_range_min_array = np.arange(0.0, abs(self.y_min)*1.00000000001, self.y_main_dis) * (-1)
        y_range_array = np.sort(np.unique(np.hstack([y_range_max_array, y_range_min_array])))
        self.ax.set_yticks(y_range_array)
        self.ax.set_ylim(self.y_min, self.y_max*1.00000000001)
        self.ax.set_ylabel(self.y_label)

    def func_rde_inlet(self):
        list_x_axis_x = [-10e1, 10e1]
        list_x_axis_y = [0., 0.]
        self.ax.plot(list_x_axis_x, list_x_axis_y, color='k')

    def func_rde_exit(self, rde_l):
        list_x_axis_x = [-10e1, 10e1]
        list_x_axis_y = [rde_l, rde_l]
        self.ax.plot(list_x_axis_x, list_x_axis_y, color='k')

    def func_graph_add(self, list_x, list_y, color=None):
        #### ==================================================================================================================
        #### graph depict
        #### ==================================================================================================================
        # colorlist = ["r", "g", "b", "c", "m", "y", "k", "w"]
        # self.ax.scatter(list_x, list_y, s=marker_size, c=marker_color)
        if color==None:
            self.ax.plot(list_x, list_y)
        else:
            self.ax.plot(list_x, list_y, c=color)

    def func_scatter_add(self,array_x,array_y):
        array_x = array_x.flatten()
        array_y = array_y.flatten()
        self.ax.scatter(array_x,array_y)

    def func_show(self):
        plt.show()



#### ================================================
#### 関数
#### ================================================

# def func_differencial_eq_10_32(x):
#     dx = 10e-6
#     # dif = ((func_eq_10_32(x+dx) - func_eq_10_32(x)) / (dx))
#     # dif = ((func_eq_10_32(x) - func_eq_10_32(x-dx)) / (dx))
#     dif = ((func_eq_10_32(x+dx/2.) - func_eq_10_32(x-dx/2.)) / (dx)) # 中間を取ります，方法の名前忘れ申した
#     return dif

def func_cross(list_a, list_b):
    a0 = list_a[0]
    a1 = list_a[1]
    b0 = list_b[0]
    b1 = list_b[1]
    cross_x = -(b0-b1)/(a0-a1)
    cross_y = a0 * cross_x + b0
    return cross_x, cross_y

def func_delta_M(neu_target, M):
    part1 = ((gamma+1.)/(gamma-1.))**(1./2.)
    part2 = np.atan(((gamma-1.)/(gamma+1.)*(M**2.-1.))**(1./2.))
    part3 = np.atan((M**2.-1.)**(1./2.))
    neu_cal = part1*part2-part3
    delta_neu = neu_target - neu_cal
    return delta_neu

def func_neu2M(neu_target):
    M_a = 1.2
    M_b = 1.4
    eps = 10e-6
    dM_a = func_delta_M(neu_target,M_a)
    dM_b = func_delta_M(neu_target,M_b)
    while abs(dM_b) > eps:
        M_s = (M_a * dM_b - M_b * dM_a)/(dM_b - dM_a)
        M_a, M_b = M_b, M_s
        dM_a = dM_b
        dM_b = func_delta_M(neu_target, M_b)
    M_result = M_b
    return M_result

def func_intercept(slope, point):
    x_A = point[0]
    y_A = point[1]
    intercept = y_A - slope * x_A
    return intercept 




### 同じやん，いつか統一して
### こいつが一般版
def func_cross_gas_dynamics(point1, point2, lambda1, lambda2): ### ごめん，くそみたいな名前で
    x = ((point1[1]-point2[1])-(lambda1*point1[0]-lambda2*point2[0]))/(lambda2-lambda1)
    y = point1[1] - lambda1 * point1[0] + lambda1 * x
    return x, y

### よく理解していない
### とりあえず vol.2 p.203 の例題にしたがって書いてみる
def func_MEPC_theta3(theta1, theta2, point1, point2, point4, lambda_12, eps=10e-6):
    lambda_o = np.tan((theta1+theta2)/2.)
    delta = 1.0
    # n=0
    while delta >= eps:
        x, y = func_cross_gas_dynamics(point2, point4, lambda_12,lambda_o)
        theta = theta2 + (y-point2[1])/(point1[1]-point2[1]) * (theta1-theta2)
        lambda_o_new = np.tan(theta)
        delta = abs((lambda_o-lambda_o_new)/lambda_o)
        lambda_o = lambda_o_new
        # n += 1
        # print(n)
    return x, y, theta, lambda_o




















#------------------------------------------------------------------
#### 0. parameters
#------------------------------------------------------------------
gamma = 1.4 # 比熱比[-]
rde_l = 3.0 # [-]: RDE's combustion chamber length normalized by injection fill height 

#------------------------------------------------------------------
#### 0. graph_prepare
#------------------------------------------------------------------
## x deto波を「1」として無次元化
## y deto波を「1」として無次元化
graph0 = CL_graph_create()
graph0() # reset
graph0.func_rde_inlet()
graph0.func_rde_exit(rde_l)

#------------------------------------------------------------------
#### 0. assumptions for detonation wave
#------------------------------------------------------------------
## deto波と燃焼室底面の接点はどの条件であろうと不変である -> 原点 （本番でもこのつもり）
## この計算では「deto_height」固定
angle_dw = 100. / 360. * 2. * np.pi # [rad]: detonation angle from horizontal axis (theta axis)
height_dw = 1.0 # [-]: injection fill height normalized by deto_height (z axis)
## deto波描画
## そういえばatanの値域って -np.pi/2. ~ +np.pi/2. だったっけか 
array_point_dw = (height_dw*np.tan(-(angle_dw-np.pi/2.)), height_dw)
graph0.func_graph_add((0., array_point_dw[0]), (0., array_point_dw[1]), color="r")

#------------------------------------------------------------------
#### 0. assumptions for fresh mixture layer
#------------------------------------------------------------------
## 「detonation wave」と「fresh mixture layer」のどちらを先に描画するべきなのか分からないからとりあえず
angle_fm = 10. / 360. * 2. * np.pi # deto_angle - np.pi/2.
slope_fm = np.tan(angle_fm)
intercept_fm = func_intercept(slope_fm, array_point_dw)
x_cross, y_cross = func_cross((0., slope_fm), (0., intercept_fm))
graph0.func_graph_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="b")

#------------------------------------------------------------------
#### 0. assumptions for slip line
#------------------------------------------------------------------
angle_sl = 20. / 360. * 2. * np.pi ### [rad]: slip line angle from horizontal axis (theta axis)
slope_sl = np.tan(angle_sl)
intercept_sl = func_intercept(slope_sl, array_point_dw)
x_cross, y_cross = func_cross((0., slope_sl), (rde_l, intercept_sl))
graph0.func_graph_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="b")

# #------------------------------------------------------------------
# #### 0. assumptions for oblique-shock
# #------------------------------------------------------------------
# angle_os = 60. / 360. * 2. * np.pi ### [rad]: slip line angle from horizontal axis (theta axis)
# slope_os = math.tan(angle_os)
# intercept_os = func_intercept(slope_os, array_point_dw)
# x_cross, y_cross = func_cross((0, slope_os), (rde_l, intercept_os))
# graph0.func_graph_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="r")

# #------------------------------------------------------------------
# #### 0. assumptions for dw (2nd)
# #------------------------------------------------------------------
# ## deto波と燃焼室底面の接点はどの条件であろうと不変である -> 原点 （本番でもこのつもり）
# ## この計算では「deto_height」固定
# ## deto波描画
# ## そういえばatanの値域って -np.pi/2. ~ +np.pi/2. だったっけか 
# array_point_dw_2nd = (2.*np.pi-1.*math.tan(angle_fm), height_dw)
# graph0.func_graph_add((2.*np.pi, array_point_dw_2nd[0]), (0., array_point_dw_2nd[1]), color="r")

# #------------------------------------------------------------------
# #### 0. assumptions for fresh-mixture (2nd)
# #------------------------------------------------------------------
# intercept_fm_2nd = func_intercept(slope_fm, array_point_dw_2nd)
# x_cross, y_cross = func_cross((0., slope_fm), (0., intercept_fm_2nd))
# graph0.func_graph_add((array_point_dw_2nd[0], x_cross), (array_point_dw_2nd[1], y_cross), color="b")

# #------------------------------------------------------------------
# #### 0. assumptions for slip line (2nd)
# #------------------------------------------------------------------
# intercept_sl_2nd = func_intercept(slope_sl, array_point_dw_2nd)
# x_cross, y_cross = func_cross((0., slope_sl), (rde_l, intercept_sl_2nd))
# graph0.func_graph_add((array_point_dw_2nd[0], x_cross), (array_point_dw_2nd[1], y_cross), color="b")

# #------------------------------------------------------------------
# #### 0. assumptions for oblique-shock (2nd)
# #------------------------------------------------------------------
# intercept_os_2nd = func_intercept(slope_os, array_point_dw_2nd)
# x_cross, y_cross = func_cross((0, slope_os), (rde_l, intercept_os_2nd))
# graph0.func_graph_add((array_point_dw_2nd[0], x_cross), (array_point_dw_2nd[1], y_cross), color="r")















angle_bottom = 0. * 2. * np.pi /360.

#------------------------------------------------------------------
#### 仮定：
#------------------------------------------------------------------
### cantera_setting
import cantera as ct
from sdtoolbox.thermo import soundspeed_fr, soundspeed_eq
mech = 'gri30_highT.cti'
gas = ct.Solution(mech)

# Layer detonation computation for gri30_highT.cti with composition C2H4:0.3333333333333333, O2:1.0
# State 1 - Initial state of reacting layer
P1 = 109278.90298573895 #(Pa)
T1 = 256.5331450748727 #(K)
rho1 = 1.5872829592437772 #(kg/m3)
a_fr1 = 305.2743367541806 #(m/s)
h1 = 379882.37074241653 #(J/kg)
s1 = 6703.726763985172 #(J/kg K)
gamma_fr1 = 1.3554605754094926

# State 2 - CJ 
CJ_speed = 2386.591152672143 #(m/s)
P2 = 4270213.871314236 #(Pa)
T2 = 3979.0578618774894 #(K)
rho2 = 2.9407834010709206 #(kg/m3)
h2 = 2398111.5517674736 #(J/kg K)
s2 = 11619.80436302621 #(J/kg K)
w2 = 1288.1586130889275 #(m/s)
u2 = 1098.4325395832154 #(m/s)
a2 = 1286.5955416302932 #(m/s)
x2 = [5.66228974e-02, 5.63248708e-02, 7.06901613e-02, 1.03138136e-01,\
 1.26319510e-01, 2.19178349e-01, 3.61837693e-04, 2.43583610e-05,\
 9.87591241e-09, 2.04220827e-09, 1.05416894e-09, 1.41773264e-10,\
 4.55493375e-10, 3.81726776e-11, 2.55874711e-01, 1.11443091e-01,\
 2.18018804e-05, 2.62082032e-07, 2.03546674e-10, 1.69701442e-11,\
 4.93273199e-12, 7.88957542e-13, 1.64881692e-12, 1.42225878e-15,\
 1.09671898e-16, 1.42327644e-19, 2.28514968e-21, 2.89655410e-11,\
 4.99422907e-12, 1.24757748e-13, 0.00000000e+00, 0.00000000e+00,\
 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,\
 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,\
 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,\
 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,\
 0.00000000e+00, 2.36527675e-29, 3.48786919e-31, 2.78780489e-15,\
 1.59215284e-16]

### total enthalpy
U2 = w2
h2_U2 = (h2 + U2**2./2.) # SDT

def func_dM(M,p):
    gas.SPX = s2, p, x2
    v2 = 2.0*(h2_U2 - gas.enthalpy_mass)
    v = np.sqrt(v2)
    dM = abs(M - v/soundspeed_fr(gas)) / M
    return dM

def func_M2P(M, eps=10e-6):
    P_a = 0.001 * 101300
    P_b = 0.002 * 101300
    dM_b = func_dM(M, P_b)
    dM_a = func_dM(M, P_a)
    while abs(dM_b) > eps:
        P_s = (P_a * dM_b - P_b * dM_a)/(dM_b - dM_a)
        P_a, P_b = P_b, P_s
        dM_a = dM_b
        dM_b = func_dM(M, P_b)
    gas.SPX = s2, P_b,x2
    T = gas.T
    R = ct.gas_constant / gas.mean_molecular_weight
    rho = gas.density
    # h = gas.enthalpy_mass
    # h_total = hwc + W ** 2. / 2.
    a_fr = soundspeed_fr(gas)
    V = M * a_fr
    gamma = a_fr**2.*rho/P_b
    return P_b, T, R, rho, a_fr, V, gamma 



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
array_x_up = np.ones((int(num_ch_up))) * array_point_dw[0]
array_x = np.flipud(np.diag(array_x_up))
array_x = np.delete(array_x,-1,0)
array_x_down = np.zeros((int(num_ch_down)))
array_x_down = np.transpose(np.array([array_x_down]))
array_x_down = np.hstack((array_x_down,array_zero0))
array_x = np.vstack((array_x,array_x_down))
array_x = np.hstack((array_x, array_zero1))
del array_x_up
del array_x_down


### y for characteristics
array_y_up = np.ones((int(num_ch_up))) * array_point_dw[1]
array_y = np.flipud(np.diag(array_y_up))
array_y = np.delete(array_y,-1,0)
array_y_down = np.zeros((int(num_ch_down)))
array_y_down = np.transpose(np.array([array_y_down]))
array_y_down = np.hstack((array_y_down,array_zero0))
array_y = np.vstack((array_y,array_y_down))
array_y = np.hstack((array_y, array_zero1))
del array_y_up
del array_y_down


### ============================================================================================ 20211018_sawada
### ============================================================================================ 初期値の設定がよくわからない
### ============================================================================================ とりあえずリーマン不変量で
### あくまで格子点における値
### 特性線とは違うし，あくまで区別してくれ
### theta
### neu, M, alpha
array_theta_up = np.linspace(angle_fm,angle_sl,num_ch_up)
array_neu_up = np.linspace(angle_fm,angle_sl,num_ch_up)
array_neu_up = array_neu_up - angle_fm
array_M_up = np.zeros((int(num_ch_up)))
array_alpha_up = np.zeros((int(num_ch_up)))
array_p_up = np.zeros((int(num_ch_up)))
array_t_up = np.zeros((int(num_ch_up)))
array_R_up = np.zeros((int(num_ch_up)))
array_rho_up = np.zeros((int(num_ch_up)))
array_a_fr_up = np.zeros((int(num_ch_up)))
array_V_up = np.zeros((int(num_ch_up)))
array_gamma_up = np.zeros((int(num_ch_up)))
### P_b, T, R, rho, a_fr, V, gamma 
for i0 in range(int(num_ch_up)):
    array_M_up[i0] = func_neu2M(array_neu_up[i0])
    array_alpha_up[i0] = np.asin(1./array_M_up[i0])
    array_p_up[i0], \
        array_t_up[i0], \
            array_R_up[i0], \
                array_rho_up[i0], \
                    array_a_fr_up[i0], \
                        array_V_up[i0], \
                            array_gamma_up[i0] = func_M2P(array_M_up[i0])

array_theta = np.flipud(np.diag(array_theta_up))
array_neu = np.flipud(np.diag(array_neu_up))
array_M = np.flipud(np.diag(array_M_up))
array_alpha = np.flipud(np.diag(array_alpha_up))
array_p = np.flipud(np.diag(array_p_up))
array_t = np.flipud(np.diag(array_t_up))
array_R = np.flipud(np.diag(array_R_up))
array_rho = np.flipud(np.diag(array_rho_up))
array_a_fr = np.flipud(np.diag(array_a_fr_up))
array_V = np.flipud(np.diag(array_V_up))
array_gamma = np.flipud(np.diag(array_gamma_up))


#============================================================================
array_theta = np.delete(array_theta,-1,0)
array_neu = np.delete(array_neu,-1,0)
array_M = np.delete(array_M,-1,0)
array_alpha = np.delete(array_alpha,-1,0)
array_p = np.delete(array_p,-1,0)
array_t = np.delete(array_t,-1,0)
array_R = np.delete(array_R,-1,0)
array_rho = np.delete(array_rho,-1,0)
array_a_fr = np.delete(array_a_fr,-1,0)
array_V = np.delete(array_V,-1,0)
array_gamma = np.delete(array_gamma,-1,0)


#============================================================================
array_theta_down = np.linspace(angle_fm,angle_bottom,num_ch_down)
array_neu_down = np.linspace(angle_fm,angle_bottom,num_ch_down)
array_neu_down = angle_fm - array_neu_down
array_M_down = np.zeros((int(num_ch_down)))
array_alpha_down = np.zeros((int(num_ch_down)))
array_p_down = np.zeros((int(num_ch_down)))
array_t_down = np.zeros((int(num_ch_down)))
array_R_down = np.zeros((int(num_ch_down)))
array_rho_down = np.zeros((int(num_ch_down)))
array_a_fr_down = np.zeros((int(num_ch_down)))
array_V_down = np.zeros((int(num_ch_down)))
array_gamma_down = np.zeros((int(num_ch_down)))
### P_b, T, R, rho, a_fr, V, gamma 
for i0 in range(int(num_ch_down)):
    array_M_down[i0] = func_neu2M(array_neu_down[i0])
    array_alpha_down[i0] = np.asin(1./array_M_down[i0])
    array_p_down[i0], \
        array_t_down[i0], \
            array_R_down[i0], \
                array_rho_down[i0], \
                    array_a_fr_down[i0], \
                        array_V_down[i0], \
                            array_gamma_down[i0] = func_M2P(array_M_down[i0])


#============================================================================
### =====
array_theta_down = np.transpose(np.array([array_theta_down]))
array_theta_down = np.hstack((array_theta_down,array_zero0))
array_theta = np.vstack((array_theta,array_theta_down))
array_theta = np.hstack((array_theta, array_zero1))
### =====
array_M_down = np.transpose(np.array([array_M_down]))
array_M_down = np.hstack((array_M_down,array_zero0))
array_M = np.vstack((array_M,array_M_down))
array_M = np.hstack((array_M, array_zero1))
### =====
array_alpha_down = np.transpose(np.array([array_alpha_down]))
array_alpha_down = np.hstack((array_alpha_down,array_zero0))
array_alpha = np.vstack((array_alpha,array_alpha_down))
array_alpha = np.hstack((array_alpha, array_zero1))
### =====
array_p_down = np.transpose(np.array([array_p_down]))
array_p_down = np.hstack((array_p_down,array_zero0))
array_p = np.vstack((array_p,array_p_down))
array_p = np.hstack((array_p, array_zero1))
### =====
array_t_down = np.transpose(np.array([array_t_down]))
array_t_down = np.hstack((array_t_down,array_zero0))
array_t = np.vstack((array_t,array_t_down))
array_t = np.hstack((array_t, array_zero1))
### =====
array_R_down = np.transpose(np.array([array_R_down]))
array_R_down = np.hstack((array_R_down,array_zero0))
array_R = np.vstack((array_R,array_R_down))
array_R = np.hstack((array_R, array_zero1))
### =====
array_rho_down = np.transpose(np.array([array_rho_down]))
array_rho_down = np.hstack((array_rho_down,array_zero0))
array_rho = np.vstack((array_rho,array_rho_down))
array_rho = np.hstack((array_rho, array_zero1))
### =====
array_a_fr_down = np.transpose(np.array([array_a_fr_down]))
array_a_fr_down = np.hstack((array_a_fr_down,array_zero0))
array_a_fr = np.vstack((array_a_fr,array_a_fr_down))
array_a_fr = np.hstack((array_a_fr, array_zero1))
### =====
array_V_down = np.transpose(np.array([array_V_down]))
array_V_down = np.hstack((array_V_down,array_zero0))
array_V = np.vstack((array_V,array_V_down))
array_V = np.hstack((array_V, array_zero1))
### =====
array_gamma_down = np.transpose(np.array([array_gamma_down]))
array_gamma_down = np.hstack((array_gamma_down,array_zero0))
array_gamma = np.vstack((array_gamma,array_gamma_down))
array_gamma = np.hstack((array_gamma, array_zero1))


#============================================================================
del array_theta_up
del array_theta_down
del array_neu_up
del array_neu_down
del array_M_up
del array_M_down
del array_alpha_up
del array_alpha_down
del array_p_up
del array_p_down
del array_t_up
del array_t_down
del array_R_up
del array_R_down
del array_rho_up
del array_rho_down
del array_a_fr_up
del array_a_fr_down
del array_V_up
del array_V_down
del array_gamma_up
del array_gamma_down
# ====











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






### 対象は theta_3 ???
### predictor : 一つ目の比較対象の値
### corrector : 二つ目の比較対象の値
### itration  : abs(predictor - corrector) が許容誤差未満に至るまで計算（absでなくてもいいかもしれないが）

### 必要になったら随時，パラメーターの行列を追加していく方針でお願いします．

for i in range(1,int(num_ch_up)):
    for j in range(int(num_ch_up-i), int((num_ch_up+num_ch_down)-2)):


        ### =====================================================================
        ### predictor
        ### =====================================================================

        ### lambda_plus & lambda_minus - eq17dot47_eq17dot48 (first step predictor)
        array_lambda_plus[j+1][i-1] = np.tan(array_theta[j+1][i-1]+array_alpha[j+1][i-1])
        array_lambda_minus[j-1][i] = np.tan(array_theta[j-1][i]-array_alpha[j-1][i])
        ### Q+ & Q- - eq17dot54_eq17dot55 (first step predictor) (noting delta is removed)
        array_Q_plus[j+1][i-1] = np.sqrt((array_M[j+1][i-1])**2.-1.) / (array_rho[j+1][i-1]*array_V[j+1][i-1])
        array_Q_minus[j-1][i] = np.sqrt((array_M[j-1][i])**2.-1.) / (array_rho[j-1][i]*array_V[j-1][i])
        ### S+ & S- - eq17dot54_eq17dot55 (first step predictor) (noting delta is removed)
        array_S_plus[j+1][i-1] = np.sin(array_theta[j+1][i-1]) / \
            (array_y[j+1][i-1] * array_M[j+1][i-1] * np.cos(array_theta[j+1][i-1]+array_alpha[j+1][i-1]))
        array_S_minus[j-1][i] = np.sin(array_theta[j-1][i]) / \
            (array_y[j-1][i] * array_M[j-1][i] * np.cos(array_theta[j-1][i]-array_alpha[j-1][i]))

        ### eq17dot44_eq17dot45
        array_x[j][i], array_y[j][i] = func_cross_gas_dynamics(\
            (array_x[j-1][i], array_y[j-1][i]),\
                (array_x[j+1][i-1], array_y[j+1][i-1]),\
                    array_lambda_minus[j-1][i], \
                        array_lambda_plus[j+1][i-1])

        ### T+ & T- - eq17dot52_eq17dot53 (first step predictor)
        array_T_plus[j+1][i-1] = -array_S_plus * (array_x[j][i]-array_x[j+1][i-1]) + \
            array_Q_plus[j+1][i-1] * array_p[j+1][i-1] + array_theta[j+1][i-1]
        array_T_minus[j-1][i] = -array_S_minus * (array_x[j][i]-array_x[j-1][i]) + \
            array_Q_minus[j-1][i] * array_p[j-1][i] - array_theta[j-1][i]
            
        ### eq17dot49 (the modified euler predictor-corrector)
        ### eq17dot43_eq17dot46
        array_lambda_12[j][i] = (array_y[j+1][i-1] - array_y[j-1][i]) / (array_x[j+1][i-1] - array_x[j-1][i])
        array_x_3[j][i],\
            array_y_3[j][i],\
                array_theta_3[j][i],\
                    array_lambda_o[j][i] = func_MEPC_theta3(\
            array_theta[j-1][i],\
                array_theta[j+1][i-1],\
                    (array_x[j-1][i], array_y[j-1][i]),\
                        (array_x[j+1][i-1], array_y[j+1][i-1]),\
                            (array_x[j][i], array_y[j][i]),\
                                array_lambda_12[j][i])

        ### interpolating for the remaining flow properties gives... (p.203) 
        array_p_3[j][i] = array_p[j+1][i-1]+(array_y_3[j][i]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_p[j-1][i]-array_p[j+1][i-1])
        array_rho_3[j][i] = array_rho[j+1][i-1]+(array_y_3[j][i]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_rho[j-1][i]-array_rho[j+1][i-1])
        array_V_3[j][i] = array_V[j+1][i-1]+(array_y_3[j][i]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_V[j-1][i]-array_V[j+1][i-1])

        ### (f)
        gas.SPX = s2, array_p_3[j][i], x2
        # array_gamma_3[j][i] = a1_fr**2.*rho1/Pcom
        array_R_o[j][i] = ct.gas_constant / gas.mean_molecular_weight
        array_a_fr_3[j][i] = soundspeed_fr(gas)
        array_A_o[j][i] = array_a_fr_3[j][i] ** 2.
        array_T_o1[j][i] = array_R_o[j][i] * array_V_3[j][i] + array_p_3[j][i] # using state3
        array_T_o2[j][i] = array_p_3[j][i] - array_A_o[j][i] * array_rho_3[j][i] # using state3

        ### eq.(g) & (h) for calculating p4, theta4
        array_p[j][i] = (array_T_plus[j][i] + array_T_minus[j][i]) / (array_Q_plus[j][i] + array_Q_minus[j][i])
        array_theta[j][i] = array_T_plus[j][i] - array_Q_plus[j][i] * array_p[j][i]

        ### eq.(i) & (j) for calculating V4, rho4
        array_V[j][i] = (array_T_o1[j][i]-array_p[j][i]) / array_R_o[j][i]
        array_rho[j][i] = (array_p[j][i]-array_T_o2[j][i]) / array_A_o[j][i]













        ### =====================================================================
        ### corrector : 全て入れなおせているのだろうか？
        ### =====================================================================
        ### along Mach line 24 (C+)
        array_p_plus[j][i] = (array_p[j+1][i-1] + array_p[j][i]) /2.
        array_theta_plus[j][i] = (array_theta[j+1][i-1] + array_theta[j][i]) /2.
        array_V_plus[j][i] = (array_V[j+1][i-1] + array_V[j][i]) /2.
        array_rho_plus[j][i] = (array_rho[j+1][i-1] + array_rho[j][i]) /2.
        array_y_plus[j][i] = (array_y[j+1][i-1] + array_y[j][i]) /2.
        gas.SPX = s2, array_p_plus[j][i], x2
        array_a_fr_plus[j][i] = soundspeed_fr(gas)
        array_M_plus[j][i] = array_V_plus[j][i] / array_a_fr_plus[j][i]
        array_alpha_plus[j][i] = np.arcsin(1./array_M_plus[j][i])
        array_lambda_plus[j][i] = np.tan(array_theta_plus[j][i]+array_alpha_plus[j][i])
        array_Q_plus[j][i] = np.sqrt(array_M_plus[j][i]**2.-1.) / (array_rho_plus[j][i]*array_V_plus[j][i]**2.)
        array_S_plus[j][i] = np.sin(array_theta_plus[j][i]) / \
            (array_y_plus[j][i]*array_M_plus[j][i]*np.cos(array_theta[j][i]+array_theta_plus[j][i]))

        ### along Mach line 14 (C-)
        array_p_minus[j][i] = (array_p[j-1][i] + array_p[j][i]) /2.
        array_theta_minus[j][i] = (array_theta[j-1][i] + array_theta[j][i]) /2.
        array_V_minus[j][i] = (array_V[j-1][i] + array_V[j][i]) /2.
        array_rho_minus[j][i] = (array_rho[j-1][i] + array_rho[j][i]) /2.
        array_y_minus[j][i] = (array_y[j-1][i] + array_y[j][i]) /2.
        gas.SPX = s2, array_p_minus[j][i], x2
        array_a_fr_minus[j][i] = soundspeed_fr(gas)
        array_M_minus[j][i] = array_V_minus[j][i] / array_a_fr_minus[j][i]
        array_alpha_minus[j][i] = np.arcsin(1./array_M_minus[j][i])
        array_lambda_minus[j][i] = np.tan(array_theta_minus[j][i]+array_alpha_minus[j][i])
        array_Q_minus[j][i] = np.sqrt(array_M_minus[j][i]**2.-1.) / (array_rho_minus[j][i]*array_V_minus[j][i]**2.)
        array_S_minus[j][i] = np.sin(array_theta_minus[j][i]) / \
            (array_y_minus[j][i]*array_M_minus[j][i]*np.cos(array_theta[j][i]-array_theta_minus[j][i]))

        ### eq17dot44_eq17dot45
        array_x[j][i], array_y[j][i] = func_cross_gas_dynamics(\
            (array_x[j-1][i], array_y[j-1][i]),\
                (array_x[j+1][i-1], array_y[j+1][i-1]),\
                    array_lambda_minus[j-1][i], \
                        array_lambda_plus[j+1][i-1])

        ### T+ & T- - eq17dot52_eq17dot53 (first step predictor)
        array_T_plus[j+1][i-1] = -array_S_plus * (array_x[j][i]-array_x[j+1][i-1]) + \
            array_Q_plus[j+1][i-1] * array_p[j+1][i-1] + array_theta[j+1][i-1]
        array_T_minus[j-1][i] = -array_S_minus * (array_x[j][i]-array_x[j-1][i]) + \
            array_Q_minus[j-1][i] * array_p[j-1][i] - array_theta[j-1][i]
            
        ### eq17dot49 (the modified euler predictor-corrector)
        ### eq17dot43_eq17dot46
        array_lambda_12[j][i] = (array_y[j+1][i-1] - array_y[j-1][i]) / (array_x[j+1][i-1] - array_x[j-1][i])
        array_x_3[j][i],\
            array_y_3[j][i],\
                array_theta_3[j][i],\
                    array_lambda_o[j][i] = func_MEPC_theta3(\
            array_theta[j-1][i],\
                array_theta[j+1][i-1],\
                    (array_x[j-1][i], array_y[j-1][i]),\
                        (array_x[j+1][i-1], array_y[j+1][i-1]),\
                            (array_x[j][i], array_y[j][i]),\
                                array_lambda_12[j][i])

        ### interpolating for the remaining flow properties gives... (p.203) 
        array_p_3[j][i] = array_p[j+1][i-1]+(array_y_3[j][i]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_p[j-1][i]-array_p[j+1][i-1])
        array_rho_3[j][i] = array_rho[j+1][i-1]+(array_y_3[j][i]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_rho[j-1][i]-array_rho[j+1][i-1])
        array_V_3[j][i] = array_V[j+1][i-1]+(array_y_3[j][i]-array_y[j+1][i-1])/(array_y[j-1][i]-array_y[j+1][i-1])*(array_V[j-1][i]-array_V[j+1][i-1])

        ### (f)
        gas.SPX = s2, array_p_3[j][i], x2
        # array_gamma_3[j][i] = a1_fr**2.*rho1/Pcom
        array_R_o[j][i] = ct.gas_constant / gas.mean_molecular_weight
        array_a_fr_3[j][i] = soundspeed_fr(gas)
        array_A_o[j][i] = array_a_fr_3[j][i] ** 2.
        array_T_o1[j][i] = array_R_o[j][i] * array_V_3[j][i] + array_p_3[j][i] # using state3
        array_T_o2[j][i] = array_p_3[j][i] - array_A_o[j][i] * array_rho_3[j][i] # using state3

        ### eq.(g) & (h) for calculating p4, theta4
        array_p[j][i] = (array_T_plus[j][i] + array_T_minus[j][i]) / (array_Q_plus[j][i] + array_Q_minus[j][i])
        array_theta[j][i] = array_T_plus[j][i] - array_Q_plus[j][i] * array_p[j][i]

        ### eq.(i) & (j) for calculating V4, rho4
        array_V[j][i] = (array_T_o1[j][i]-array_p[j][i]) / array_R_o[j][i]
        array_rho[j][i] = (array_p[j][i]-array_T_o2[j][i]) / array_A_o[j][i]





























print("/// array_x ///", array_x.shape)
print(array_x)

print("/// array_y ///", array_y.shape)
print(array_y)

print("/// array_theta ///", array_theta.shape)
print(array_theta/2./np.pi*360.)

print("/// array_M ///", array_M.shape)
print(array_M)

print("/// array_alpha ///", array_alpha.shape)
print(array_alpha/2./np.pi*360.)

print("/// array_lambda_plus ///", array_lambda_plus.shape)
print(array_lambda_plus)

print("/// array_lambda_minus ///", array_lambda_minus.shape)
print(array_lambda_minus)

print("/// array_lambda_o ///", array_lambda_o.shape)
print(array_lambda_o)

# graph0.func_scatter_add(array_x,array_y)

graph0.func_show()

















