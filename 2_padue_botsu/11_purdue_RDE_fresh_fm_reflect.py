"""
目的：
deto波固定・斜め衝撃波無視
ー＞簡便な特性曲線を描く

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
import math
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
    part2 = math.atan(((gamma-1.)/(gamma+1.)*(M**2.-1.))**(1./2.))
    part3 = math.atan((M**2.-1.)**(1./2.))
    neu_cal = part1*part2-part3
    delta_neu = neu_target - neu_cal
    return delta_neu

def func_neu2Mach(neu_target):
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
array_point_dw = (height_dw*math.tan(-(angle_dw-np.pi/2.)), height_dw)
graph0.func_graph_add((0., array_point_dw[0]), (0., array_point_dw[1]), color="r")

#------------------------------------------------------------------
#### 0. assumptions for fresh mixture layer
#------------------------------------------------------------------
## 「detonation wave」と「fresh mixture layer」のどちらを先に描画するべきなのか分からないからとりあえず
angle_fm = 10. / 360. * 2. * np.pi # deto_angle - np.pi/2.
slope_fm = math.tan(angle_fm)
intercept_fm = func_intercept(slope_fm, array_point_dw)
x_cross, y_cross = func_cross((0., slope_fm), (0., intercept_fm))
graph0.func_graph_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="b")

#------------------------------------------------------------------
#### 0. assumptions for slip line
#------------------------------------------------------------------
angle_sl = 40. / 360. * 2. * np.pi ### [rad]: slip line angle from horizontal axis (theta axis)
slope_sl = math.tan(angle_sl)
intercept_sl = func_intercept(slope_sl, array_point_dw)
x_cross, y_cross = func_cross((0., slope_sl), (rde_l, intercept_sl))
graph0.func_graph_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="b")

#------------------------------------------------------------------
#### 0. assumptions for oblique-shock
#------------------------------------------------------------------
angle_os = 60. / 360. * 2. * np.pi ### [rad]: slip line angle from horizontal axis (theta axis)
slope_os = math.tan(angle_os)
intercept_os = func_intercept(slope_os, array_point_dw)
x_cross, y_cross = func_cross((0, slope_os), (rde_l, intercept_os))
graph0.func_graph_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="r")

#------------------------------------------------------------------
#### 0. assumptions for dw (2nd)
#------------------------------------------------------------------
## deto波と燃焼室底面の接点はどの条件であろうと不変である -> 原点 （本番でもこのつもり）
## この計算では「deto_height」固定
## deto波描画
## そういえばatanの値域って -np.pi/2. ~ +np.pi/2. だったっけか 
array_point_dw_2nd = (2.*np.pi-1.*math.tan(angle_fm), height_dw)
graph0.func_graph_add((2.*np.pi, array_point_dw_2nd[0]), (0., array_point_dw_2nd[1]), color="r")

#------------------------------------------------------------------
#### 0. assumptions for fresh-mixture (2nd)
#------------------------------------------------------------------
# intercept_fm_2nd = func_intercept(slope_fm, array_point_dw_2nd)
# x_cross, y_cross = func_cross((0., slope_fm), (0., intercept_fm_2nd))
slope_fm_2nd = slope_fm * 2
intercept_fm_2nd = func_intercept(slope_fm_2nd, array_point_dw_2nd)
x_cross, y_cross = func_cross((0., slope_fm_2nd), (0., intercept_fm_2nd))
graph0.func_graph_add((array_point_dw_2nd[0], x_cross), (array_point_dw_2nd[1], y_cross), color="b")

#------------------------------------------------------------------
#### 0. assumptions for slip line (2nd)
#------------------------------------------------------------------
intercept_sl_2nd = func_intercept(slope_sl, array_point_dw_2nd)
x_cross, y_cross = func_cross((0., slope_sl), (rde_l, intercept_sl_2nd))
graph0.func_graph_add((array_point_dw_2nd[0], x_cross), (array_point_dw_2nd[1], y_cross), color="b")

#------------------------------------------------------------------
#### 0. assumptions for oblique-shock (2nd)
#------------------------------------------------------------------
intercept_os_2nd = func_intercept(slope_os, array_point_dw_2nd)
x_cross, y_cross = func_cross((0, slope_os), (rde_l, intercept_os_2nd))
graph0.func_graph_add((array_point_dw_2nd[0], x_cross), (array_point_dw_2nd[1], y_cross), color="r")















angle_bottom = 0. * 2. * np.pi /360.








def func_fm(x):
    y = math.tan(angle_fm) * x + intercept_fm_2nd
    return y


def func_os(x):
    y = math.tan(angle_os) * x + intercept_os_2nd
    return y










#------------------------------------------------------------------
#### 1. characteristic lines -1st
#------------------------------------------------------------------
### num_chの謎の発散の限界30くらい？
### 理由は分からぬ
num_ch = 25 # number of characteristic lines

### i方向（横）にtheta-neu=const.確認
### j方向（縦）にtheta+neu=const.確認

array_zero0 = np.zeros((int(num_ch-1),int(num_ch)))
array_zero1 = np.zeros((int(2*num_ch-1),int(num_ch)))


### theta
array_theta_up = np.linspace(angle_fm,angle_sl,num_ch)
array_theta = np.flipud(np.diag(array_theta_up))
array_theta = np.delete(array_theta,-1,0)
array_theta_down = np.linspace(angle_fm,angle_bottom,num_ch)
array_theta_down = np.transpose(np.vstack((array_theta_down,array_zero0)))
array_theta = np.vstack((array_theta,array_theta_down))
array_theta = np.hstack((array_theta, array_zero1))
del array_theta_up
del array_theta_down

### neu, Mach, beta
array_neu_up = np.linspace(angle_fm,angle_sl,num_ch)
array_neu_up = array_neu_up - angle_fm
array_Mach_up = np.zeros((int(num_ch)))
array_beta_up = np.zeros((int(num_ch)))
for i0 in range(int(num_ch)):
    array_Mach_up[i0] = func_neu2Mach(array_neu_up[i0])
    array_beta_up[i0] = math.asin(1./array_Mach_up[i0])
array_neu = np.flipud(np.diag(array_neu_up))
array_Mach = np.flipud(np.diag(array_Mach_up))
array_beta = np.flipud(np.diag(array_beta_up))
#====
array_neu = np.delete(array_neu,-1,0)
array_Mach = np.delete(array_Mach,-1,0)
array_beta = np.delete(array_beta,-1,0)
#====
array_neu_down = np.linspace(angle_fm,angle_bottom,num_ch)
array_neu_down = angle_fm - array_neu_down
array_Mach_down = np.zeros((int(num_ch)))
array_beta_down = np.zeros((int(num_ch)))
for i0 in range(int(num_ch)):
    array_Mach_down[i0] = func_neu2Mach(array_neu_down[i0])
    array_beta_down[i0] = math.asin(1./array_Mach_down[i0])
array_neu_down = np.transpose(np.vstack((array_neu_down,array_zero0)))
array_Mach_down = np.transpose(np.vstack((array_Mach_down,array_zero0)))
array_beta_down = np.transpose(np.vstack((array_beta_down,array_zero0)))
array_neu = np.vstack((array_neu,array_neu_down))
array_Mach = np.vstack((array_Mach,array_Mach_down))
array_beta = np.vstack((array_beta,array_beta_down))
### array_neu&Mach&beta のすべり面要素は未知のため，そのままarray_zero1を与える
array_neu = np.hstack((array_neu, array_zero1))
array_Mach = np.hstack((array_Mach, array_zero1))
array_beta = np.hstack((array_beta, array_zero1))
del array_neu_up
del array_neu_down
del array_Mach_up
del array_Mach_down
del array_beta_up
del array_beta_down

### alpha_plus, alpha_minus
array_alpha_plus_up = np.zeros((int(num_ch)))
array_alpha_plus = np.flipud(np.diag(array_alpha_plus_up))
array_alpha_plus = np.delete(array_alpha_plus,-1,0)
array_alpha_plus_down = np.zeros((int(num_ch)))
for i0 in range(1, int(num_ch)):
    array_alpha_plus_down[i0] = (array_theta[int(num_ch-2)+i0][0]+array_beta[int(num_ch-2)+i0][0])/2. + \
        (array_theta[int(num_ch-2)+i0+1][0]+array_beta[int(num_ch-2)+i0+1][0])/2.
array_alpha_plus_down = np.transpose(np.vstack((array_alpha_plus_down,array_zero0)))
array_alpha_plus = np.vstack((array_alpha_plus,array_alpha_plus_down))
### array_alpha_plus のすべり面要素は未知のため，そのままarray_zero1を与える
array_alpha_plus = np.hstack((array_alpha_plus,array_zero1))
del array_alpha_plus_up
del array_alpha_plus_down

array_alpha_minus_up = np.zeros((int(num_ch)))
for i0 in range(1,int(num_ch)):
    array_alpha_minus_up[i0] = (array_theta[int(num_ch)-i0][i0-1]-array_beta[int(num_ch)-i0][i0-1])/2. + \
        (array_theta[int(num_ch)-i0-1][i0]-array_beta[int(num_ch)-i0-1][i0])/2.
array_alpha_minus = np.flipud(np.diag(array_alpha_minus_up))
array_alpha_minus = np.delete(array_alpha_minus,-1,0)
array_alpha_minus_down = np.zeros((int(num_ch)))
array_alpha_minus_down = np.transpose(np.vstack((array_alpha_minus_down,array_zero0)))
array_alpha_minus = np.vstack((array_alpha_minus,array_alpha_minus_down))
### array_alpha_minus のすべり面要素は未知のため，そのままarray_zero1を与える
array_alpha_minus = np.hstack((array_alpha_minus,array_zero1))
del array_alpha_minus_up
del array_alpha_minus_down


### intercept_plus, intercept_minus
array_intercept_plus_up = np.zeros((int(num_ch)))
array_intercept_plus = np.flipud(np.diag(array_intercept_plus_up))
array_intercept_plus = np.delete(array_intercept_plus,-1,0)
array_intercept_plus_down = np.zeros((int(num_ch)))
for i0 in range(1, int(num_ch)):
    array_intercept_plus_up[i0] = func_intercept(math.tan(array_alpha_plus[int(num_ch-2)+i0+1][0]), (0.,0.)) ### ただの原点
array_intercept_plus_down = np.transpose(np.vstack((array_intercept_plus_down,array_zero0)))
array_intercept_plus = np.vstack((array_intercept_plus,array_intercept_plus_down))
### array_intercept_plus のすべり面要素は未知のため，そのままarray_zero1を与える
array_intercept_plus = np.hstack((array_intercept_plus,array_zero1))
del array_intercept_plus_up
del array_intercept_plus_down

array_intercept_minus_up = np.zeros((int(num_ch)))
for i0 in range(int(num_ch)):
    array_intercept_minus_up[i0] = func_intercept(math.tan(array_alpha_minus[int(num_ch)-i0-1][i0]), array_point_dw)
array_intercept_minus = np.flipud(np.diag(array_intercept_minus_up))
array_intercept_minus = np.delete(array_intercept_minus,-1,0)
array_intercept_minus_down = np.zeros((int(num_ch)))
array_intercept_minus_down = np.transpose(np.vstack((array_intercept_minus_down,array_zero0)))
array_intercept_minus = np.vstack((array_intercept_minus,array_intercept_minus_down))
### array_intercept_minus のすべり面要素は未知のため，そのままarray_zero1を与える
array_intercept_minus = np.hstack((array_intercept_minus,array_zero1))
del array_intercept_minus_up
del array_intercept_minus_down

### x
array_x_up = np.ones((int(num_ch))) * array_point_dw[0]
array_x = np.flipud(np.diag(array_x_up))
array_x = np.delete(array_x,-1,0)
array_x_down = np.zeros((int(num_ch)))
array_x_down = np.transpose(np.vstack((array_x_down,array_zero0)))
array_x = np.vstack((array_x,array_x_down))
### array_x のすべり面要素は未知のため，そのままarray_zero1を与える
array_x = np.hstack((array_x,array_zero1))
del array_x_up
del array_x_down

### y
array_y_up = np.ones((int(num_ch))) * array_point_dw[1]
array_y = np.flipud(np.diag(array_y_up))
array_y = np.delete(array_y,-1,0)
array_y_down = np.zeros((int(num_ch)))
array_y_down = np.transpose(np.vstack((array_y_down,array_zero0)))
array_y = np.vstack((array_y,array_y_down))
### array_y のすべり面要素は未知のため，そのままarray_zero1を与える
array_y = np.hstack((array_y, array_zero1))
del array_y_up
del array_y_down

del array_zero0
del array_zero1

### array_bool
array_bool = np.zeros((int(num_ch*2-1),(int(2*num_ch))))


### 
for i1 in range(1,int(num_ch)):

    for j1 in range(int(num_ch-i1), int(2.*num_ch-2)):
        array_theta[j1][i1] = (array_neu[j1-1][i1]-array_neu[j1+1][i1-1])/2. + \
             (array_theta[j1-1][i1]+array_theta[j1+1][i1-1])/2.
        array_neu[j1][i1] = (array_neu[j1-1][i1]+array_neu[j1+1][i1-1])/2. + \
             (array_theta[j1-1][i1]-array_theta[j1+1][i1-1])/2.
        array_Mach[j1][i1] = func_neu2Mach(array_neu[j1][i1])
        array_beta[j1][i1] = math.asin(1./array_Mach[j1][i1])
        array_alpha_plus[j1][i1] = (array_theta[j1-1][i1]+array_beta[j1-1][i1])/2. + (array_theta[j1][i1]+array_beta[j1][i1])/2.
        array_alpha_minus[j1][i1] = (array_theta[j1+1][i1-1]-array_beta[j1+1][i1-1])/2. + (array_theta[j1][i1]-array_beta[j1][i1])/2.

        if j1 <= (2*num_ch-3): # array_x&y[-2(2*num_ch-3)]以降の範囲は array_intercept_plus が存在しないので計算できない
            array_x[j1][i1], array_y[j1][i1] = func_cross((math.tan(array_alpha_plus[j1+1][i1-1]),math.tan(array_alpha_minus[j1-1][i1])), \
                (array_intercept_plus[j1+1][i1-1], array_intercept_minus[j1-1][i1]))
            array_intercept_plus[j1][i1] = func_intercept(math.tan(array_alpha_plus[j1][i1]), (array_x[j1][i1], array_y[j1][i1]))
            array_intercept_minus[j1][i1] = func_intercept(math.tan(array_alpha_minus[j1][i1]), (array_x[j1][i1], array_y[j1][i1]))





        ### =====================================================================================================
        ### fm を下回ったとき
        ### =====================================================================================================
        ### 計算した (x, y) が fm の直線を下回ったとき，計算し終えた (x, y) を計算し直す
        ### 初めて特性線がfmに侵入するとき，theta, neu, Mach, beta, alpha_plus&minus は計算し直す必要なし？
        ### 初めて特性線がfmに侵入するとき，x, y, intercept_plus&minus は計算し直す必要あり
        ### 初めて特性線がfmに侵入するとき，array[-1]がかならずfmにぶつかる
        ### （つづき）array[-1]以外はfmにぶつからないようになり
        ### （つづき）array[-1]以外のif文処理は必要ないのかもしれない
        ### 要素数に変化はない？
        if array_y[j1][i1] < func_fm(array_x[j1][i1]):
            array_bool[j1][i1] = 1

            array_theta[j1][i1] = angle_fm
            array_neu[j1][i1] = array_neu[j1-1][i1] + array_theta[j1-1][i1] - array_theta[j1][i1] ### C-横断しかない
            array_Mach[j1][i1] = func_neu2Mach(array_neu[j1][i1])
            array_beta[j1][i1] = math.asin(1./array_Mach[j1][i1])
            array_alpha_plus[j1][i1] = (array_theta[j1-1][i1]+array_beta[j1-1][i1])/2. + (array_theta[j1][i1]+array_beta[j1][i1])/2.
            array_alpha_minus[j1][i1] = 0.

            if j1 <= (2*num_ch-3): # array_x&y[-2(2*num_ch-3)]以降の範囲は array_intercept_plus が存在しないので計算できない
                array_x[j1][i1], array_y[j1][i1] = func_cross((math.tan(angle_fm),math.tan(array_alpha_minus[j1-1][i1])), \
                    (intercept_fm_2nd, array_intercept_minus[j1-1][i1]))
                array_intercept_plus[j1][i1] = func_intercept(math.tan(array_alpha_plus[j1][i1]), (array_x[j1][i1], array_y[j1][i1]))
                array_intercept_minus[j1][i1] = 0.





    ### 下側反射
    array_theta[-1][i1] = angle_bottom
    array_neu[-1][i1] = array_neu[-2][i1]+array_theta[-2][i1] - array_theta[-1][i1]
    array_Mach[-1][i1] = func_neu2Mach(array_neu[-1][i1])
    array_beta[-1][i1] = math.asin(1./array_Mach[-1][i1])
    array_alpha_plus[-1][i1] = (array_theta[-2][i1]+array_beta[-2][i1])/2. + (array_theta[-1][i1]+array_beta[-1][i1])/2.
    array_alpha_minus[-1][i1] = 0.
    array_x[-1][i1], array_y[-1][i1] = func_cross((math.tan(array_alpha_minus[-2][i1]), math.tan(angle_bottom)), \
        (array_intercept_minus[-2][i1], 0.))
    array_intercept_plus[-1][i1] = func_intercept(math.tan(array_alpha_plus[-1][i1]), (array_x[-1][i1], array_y[-1][i1]))
    array_intercept_minus[-1][i1] = 0.


    ### =====================================================================================================
    if array_y[-1][i1] < func_fm(array_x[-1][i1]):
        array_bool[-1][i1] = 1
        array_theta[-1][i1] = angle_fm
        array_neu[-1][i1] = array_neu[-2][i1]+array_theta[-2][i1] - array_theta[-1][i1]
        array_Mach[-1][i1] = func_neu2Mach(array_neu[-1][i1])
        array_beta[-1][i1] = math.asin(1./array_Mach[-1][i1])
        array_alpha_plus[-1][i1] = (array_theta[-2][i1]+array_beta[-2][i1])/2. + (array_theta[-1][i1]+array_beta[-1][i1])/2.
        array_alpha_minus[-1][i1] = 0.
        array_x[-1][i1], array_y[-1][i1] = func_cross((math.tan(array_alpha_minus[-2][i1]), math.tan(angle_fm)), \
            (array_intercept_minus[-2][i1], intercept_fm_2nd))
        array_intercept_plus[-1][i1] = func_intercept(math.tan(array_alpha_plus[-1][i1]), (array_x[-1][i1], array_y[-1][i1]))
        array_intercept_minus[-1][i1] = 0.




 
for i1 in range(int(num_ch),int(num_ch*2)):
    ### 上側反射（上側反射の計算があるせいでこれまでのと合わせて一般化することは難しい）
    array_theta[0][i1] = angle_sl
    array_neu[0][i1] = array_neu[1][i1-1]-array_theta[1][i1-1] + array_theta[0][i1]
    array_Mach[0][i1] = func_neu2Mach(array_neu[0][i1])
    array_beta[0][i1] = math.asin(1./array_Mach[0][i1])
    array_alpha_plus[0][i1] = 0.
    array_alpha_minus[0][i1] = (array_theta[1][i1-1]-array_beta[1][i1-1])/2. + (array_theta[0][i1]-array_beta[0][i1])/2.
    array_x[0][i1], array_y[0][i1] = func_cross((math.tan(array_alpha_plus[1][i1-1]), math.tan(angle_sl)), \
        (array_intercept_plus[1][i1-1], intercept_sl))
    array_intercept_plus[0][i1] = 0.
    array_intercept_minus[0][i1] = func_intercept(math.tan(array_alpha_minus[0][i1]), (array_x[0][i1], array_y[0][i1]))


    # for j1 in range(1, int(2.*num_ch-2)):
    for j1 in range(1, int(2.*num_ch-2+(-i1+1))):
        array_theta[j1][i1] = (array_neu[j1-1][i1]-array_neu[j1+1][i1-1])/2. + \
             (array_theta[j1-1][i1]+array_theta[j1+1][i1-1])/2.
        array_neu[j1][i1] = (array_neu[j1-1][i1]+array_neu[j1+1][i1-1])/2. + \
             (array_theta[j1-1][i1]-array_theta[j1+1][i1-1])/2.
        array_Mach[j1][i1] = func_neu2Mach(array_neu[j1][i1])
        array_beta[j1][i1] = math.asin(1./array_Mach[j1][i1])
        array_alpha_plus[j1][i1] = (array_theta[j1-1][i1]+array_beta[j1-1][i1])/2. + (array_theta[j1][i1]+array_beta[j1][i1])/2.
        array_alpha_minus[j1][i1] = (array_theta[j1+1][i1-1]-array_beta[j1+1][i1-1])/2. + (array_theta[j1][i1]-array_beta[j1][i1])/2.

        if j1 <= (2*num_ch-3): # array_x&y[-2(2*num_ch-3)]以降の範囲は array_intercept_plus が存在しないので計算できない
            array_x[j1][i1], array_y[j1][i1] = func_cross((math.tan(array_alpha_plus[j1+1][i1-1]),math.tan(array_alpha_minus[j1-1][i1])), \
                (array_intercept_plus[j1+1][i1-1], array_intercept_minus[j1-1][i1]))
            array_intercept_plus[j1][i1] = func_intercept(math.tan(array_alpha_plus[j1][i1]), (array_x[j1][i1], array_y[j1][i1]))
            array_intercept_minus[j1][i1] = func_intercept(math.tan(array_alpha_minus[j1][i1]), (array_x[j1][i1], array_y[j1][i1]))




print("///array_theta///", array_theta.shape)
print(array_theta*360./2./np.pi)
print("///array_neu///", array_neu.shape)
print(array_neu*360./2./np.pi)

# array_plus = np.round((array_theta + array_neu),2)
# array_minus = np.round((array_theta - array_neu),2)
# print("///theta + neu///", array_theta.shape)
# print(array_plus*360./2./np.pi)
# print("///theta - neu///", array_neu.shape)
# print(array_minus*360./2./np.pi)

print("///array_Mach///", array_Mach.shape)
print(array_Mach)
print("///array_beta///", array_beta.shape)
print(array_beta*360./2./np.pi)

print("///array_alpha_plus///", array_alpha_plus.shape)
print(array_alpha_plus)
print("///array_alpha_minus///", array_alpha_minus.shape)
print(array_alpha_minus)

print("///array_intercept_plus///", array_intercept_plus.shape)
# print(np.round(array_intercept_plus,2))
print(array_intercept_plus)
print("///array_intercept_minus///", array_intercept_minus.shape)
# print(np.round(array_intercept_minus,2))
print(array_intercept_minus)


print("///array_x///", array_x.shape)
# print(np.round(array_x,2))
print(array_x)
print("///array_y///", array_y.shape)
print(np.round(array_y,2))

print("///array_bool///", array_bool.shape)
print(np.round(array_bool,2))



# print("so far ok_20210922")



graph0.func_scatter_add(array_x,array_y)

graph0.func_show()

















