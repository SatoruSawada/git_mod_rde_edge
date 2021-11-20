"""
20211029_sawada
未燃混合気層　ー　燃焼室底面　間の特性曲線
"""

## 意味ないけれども
## (__init__) -> 書いてもいいところ
## (__call__) -> 書いてはいけないところ


import numpy as np
import matplotlib.pylab as plt
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)

class CL_graph_setting:
    def __init__(self):
        #### ==================================================================================================================
        #### setting
        #### ==================================================================================================================
        # zoom = 0.6
        zoom = 4.
        # zoom = 1.
        #### x軸
        self.x_label = 'radial direction [-]'
        self.x_min = -0.02 * zoom        #### x軸最小値
        self.x_max = 0.14 * zoom        #### x軸最大値，目盛りの表示の都合でここに 0.00001 % 加算し申す
        self.x_main_dis = 0.02 * zoom   #### x軸主目盛り間隔
        # x_sub_num = 5       #### x軸主目盛り間の小目盛りの個数
        #### y軸
        self.y_label = 'azimuthal direction [-]'
        self.y_min = -0.01 * zoom        #### y軸最小値
        self.y_max = 0.07 * zoom        #### y軸最大値，目盛りの表示の都合でここに 0.00001 % 加算し申す
        self.y_main_dis = 0.01 * zoom   #### y軸主目盛り間隔
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

    def func_plot_add(self, list_x, list_y, color=None):
        #### ==================================================================================================================
        #### graph depict
        #### ==================================================================================================================
        # colorlist = ["r", "g", "b", "c", "m", "y", "k", "w"]
        # self.ax.scatter(list_x, list_y, s=marker_size, c=marker_color)
        if color==None:
            self.ax.plot(list_x, list_y)
        else:
            self.ax.plot(list_x, list_y, c=color)

    def func_scatter_add(self,array_x,array_y, color='k'):
        array_x = array_x.flatten()
        array_y = array_y.flatten()
        self.ax.scatter(array_x,array_y,color=color)

    def func_show(self):
        plt.show()


#### ================================================
#### 関数
#### ================================================

def func_cross(list_a, list_b):
    a0 = list_a[0]
    a1 = list_a[1]
    b0 = list_b[0]
    b1 = list_b[1]
    cross_x = -(b0-b1)/(a0-a1)
    cross_y = a0 * cross_x + b0
    return cross_x, cross_y

def func_intercept(slope, point):
    x_A = point[0]
    y_A = point[1]
    intercept = y_A - slope * x_A
    return intercept 

def func_cross_gas_dynamics(point1, point2, lambda1, lambda2): ### ごめん，くそみたいな名前で
    x = ((point1[1]-point2[1])-(lambda1*point1[0]-lambda2*point2[0]))/(lambda2-lambda1)
    y = point1[1] - lambda1 * point1[0] + lambda1 * x
    return x, y

def func_MEPC_theta3_mode1(theta1, theta2, point1, point2, point4, lambda_12, lambda_o, eps=10e-6):
    lambda_o = np.tan((theta1+theta2)/2.)
    delta = 1.0
    while delta >= eps:
        x, y = func_cross_gas_dynamics(point2, point4, lambda_12,lambda_o)
        theta = theta2 + (y-point2[1])/(point1[1]-point2[1]) * (theta1-theta2)
        lambda_o_new = np.tan(theta)
        delta = abs((lambda_o-lambda_o_new)/lambda_o)
        lambda_o = lambda_o_new
    return x, y, theta, lambda_o

def func_MEPC_theta3_mode2(theta1, theta2,theta4, point1, point2, point4, lambda_12, lambda_o, eps=10e-6):
    delta = 1.0
    while delta >= eps:
        x, y = func_cross_gas_dynamics(point2, point4, lambda_12,lambda_o)
        theta = theta2 + (y-point2[1])/(point1[1]-point2[1]) * (theta1-theta2)
        lambda_o_new = np.tan((theta+theta4)/2.)
        delta = abs((lambda_o-lambda_o_new)/lambda_o)
        lambda_o = lambda_o_new
    return x, y, theta, lambda_o


#------------------------------------------------------------------
#### 0. parameters
#------------------------------------------------------------------
rde_l = 0.06 # [-]: RDE's combustion chamber length normalized by injection fill height 
angle_fm = 10. / 360. * 2. * np.pi # deto_angle - np.pi/2.
angle_dw = angle_fm + 90. / 360. * 2. * np.pi # [rad]: detonation angle from horizontal axis (theta axis)
angle_sl = 43. / 360. * 2. * np.pi ### [rad]: slip line angle from horizontal axis (theta axis)
angle_sh = 60. / 360. * 2. * np.pi ### [rad]: oblique shock angle from horizontal axis (theta axis)
angle_bottom = 0. * 2. * np.pi /360.

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
height_dw = 0.02 # [-]: injection fill height normalized by deto_height (z axis)
array_point_dw = (height_dw*np.tan(-(angle_dw-np.pi/2.)), height_dw)
graph0.func_plot_add((0., array_point_dw[0]), (0., array_point_dw[1]), color="r")

#------------------------------------------------------------------
#### 0. assumptions for fresh mixture layer
#------------------------------------------------------------------
## 「detonation wave」と「fresh mixture layer」のどちらを先に描画するべきなのか分からないからとりあえず
slope_fm = np.tan(angle_fm)
intercept_fm = func_intercept(slope_fm, array_point_dw)
x_cross, y_cross = func_cross((0., slope_fm), (0., intercept_fm))
graph0.func_plot_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="b")

#------------------------------------------------------------------
#### 0. assumptions for slip line
#------------------------------------------------------------------
slope_sl = np.tan(angle_sl)
intercept_sl = func_intercept(slope_sl, array_point_dw)
x_cross, y_cross = func_cross((0., slope_sl), (rde_l, intercept_sl))
graph0.func_plot_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="b")

#------------------------------------------------------------------
#### 仮定：
#------------------------------------------------------------------
import cantera as ct
from sdtoolbox.thermo import soundspeed_fr, soundspeed_eq
mech = 'gri30_highT.cti'
gas = ct.Solution(mech)
ER = 1.
coef_C2H4 = 1
coef_O2 = 3
mol_ratio = coef_C2H4/coef_O2*ER
fuel = 'C2H4'
oxygen = 'O2'
q = fuel + ':' + str(mol_ratio) + ', ' + oxygen + ':1.0';

######## State 0 - plenum
P_ple = 202871.4251657753
T_ple = 300.10181107013716 #(K)
rho_ple = 2.5189165307697503 #(kg/m3)
a_fr_ple = 328.2601530165985 #(m/s)
h_ple = 425176.32528145495 #(J/kg)
s_ple = 6700.9149607432355 #(J/kg K)
gamma_fr_ple = 1.339709604023228
R_ple = 268.10468909303626

######## Layer detonation computation for gri30_highT.cti with composition C2H4:0.3333333333333333, O2:1.0
######## State 1 - Initial state of reacting layer
P_pre = 109278.90298573895 #(Pa)
T_pre = 256.5331450748727 #(K)
rho_pre = 1.5872829592437772 #(kg/m3)
a_fr_pre = 305.2743367541806 #(m/s)
h_pre = 379882.37074241653 #(J/kg)
s_pre = 6703.726763985172 #(J/kg K)
gamma_fr_pre = 1.3554605754094926
R_pre = 268.10468909303626

######## State 2 - CJ 
CJ_speed = 2386.591152672143 #(m/s)
P_post = 4270213.871314236 #(Pa)
T_post = 3979.0578618774894 #(K)
rho_post = 2.9407834010709206 #(kg/m3)
h_post = 2398111.5517674736 #(J/kg K)
s_post = 11619.80436302621 #(J/kg K)
w_post = 1288.1586130889275 #(m/s)
u_post = 1098.4325395832154 #(m/s)
a_eq_post = 1286.5955416302932 #(m/s)
gamma_eq_post = 1.1399806928419203
x_post = [5.66228974e-02, 5.63248708e-02, 7.06901613e-02, 1.03138136e-01,\
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

####### total enthalpy
U_post = w_post
h_post_U_post = (h_post + U_post**2./2.) # SDT

def func_delta_M(neu_target, M, gamma=gamma_eq_post):
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

def func_dM(M,p):
    gas.SPX = s_post, p, x_post
    v2 = 2.0*(h_post_U_post - gas.enthalpy_mass)
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
    gas.SPX = s_post, P_b,x_post
    T = gas.T
    R = ct.gas_constant / gas.mean_molecular_weight
    rho = gas.density
    a_fr = soundspeed_fr(gas)
    V = M * a_fr
    gamma = a_fr**2.*rho/P_b
    return P_b, T, R, rho, a_fr, V, gamma 



### ============================================================= 滑り面＆底面の間の特性線_start
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### ============================================================= 滑り面＆底面の間の特性線_end
from scipy import interpolate

array_x_fm = np.loadtxt('array_x_fm.csv')
array_y_fm = np.loadtxt('array_y_fm.csv')
array_theta_fm = np.loadtxt('array_theta_fm.csv')
print(array_x_fm.shape)
print(array_y_fm.shape)
print(array_theta_fm.shape)


### 滑り面：
func_fm_x_y = interpolate.PchipInterpolator(array_x_fm, array_y_fm)
func_fm_x_theta = interpolate.PchipInterpolator(array_x_fm, array_theta_fm)

### 底面：














### ============================================================= free boundaries condition reflection _example 17.5_start
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### =============================================================
### ============================================================= region between oblique shock and selip wall _example 17.5_start


### 二周目の滑り面の作成
### 仮で x = 0.45 にデトネーション波



# np.savetxt('array_x.csv', array_x, delimiter=',')
# np.savetxt('array_y.csv', array_y, delimiter=',')
# np.savetxt('array_theta.csv', array_theta/2./np.pi*360., delimiter=',')
# np.savetxt('array_p.csv', array_p, delimiter=',')
# np.savetxt('array_V.csv', array_V, delimiter=',')
# np.savetxt('array_M.csv', array_M, delimiter=',')
# np.savetxt('array_lambda_plus.csv', array_lambda_plus, delimiter=',')
# np.savetxt('array_lambda_minus.csv', array_lambda_minus, delimiter=',')










# print("/// array_x ///", array_x.shape)
# print(array_x)

# print("/// array_y ///", array_y.shape)
# print(array_y)

# print("/// array_theta ///", array_theta.shape)
# print(array_theta/2./np.pi*360.)

# print("/// array_M ///", array_M.shape)
# print(array_M)

# print("/// array_alpha ///", array_alpha.shape)
# print(array_alpha/2./np.pi*360.)

# print("/// array_V ///", array_V.shape)
# print(np.round(array_V,0))

# print("/// array_lambda_plus ///", array_lambda_plus.shape)
# print(array_lambda_plus)

# print("/// array_lambda_minus ///", array_lambda_minus.shape)
# print(array_lambda_minus)

# print("/// array_lambda_o ///", array_lambda_o.shape)
# print(array_lambda_o)

# graph0.func_scatter_add(array_x,array_y)
graph0.func_plot_add(array_x_fm,array_y_fm, color='r')

graph0.func_show()














