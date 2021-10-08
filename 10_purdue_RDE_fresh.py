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
np.set_printoptions(precision=3)

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
dw_angle = 100. / 360. * 2. * np.pi # [rad]: detonation angle from horizontal axis (theta axis)
dw_height = 1.0 # [-]: injection fill height normalized by deto_height (z axis)
## deto波描画
## そういえばatanの値域って -np.pi/2. ~ +np.pi/2. だったっけか 
array_point_dw = (dw_height*math.atan(-(dw_angle-np.pi/2.)), dw_height)
graph0.func_graph_add((0., array_point_dw[0]), (0., array_point_dw[1]), color="r")

#------------------------------------------------------------------
#### 0. assumptions for fresh mixture layer
#------------------------------------------------------------------
## 「detonation wave」と「fresh mixture layer」のどちらを先に描画するべきなのか分からないからとりあえず
fm_angle = 20. / 360. * 2. * np.pi # deto_angle - np.pi/2.
slope_fm = math.atan(fm_angle)
intercept_fm = func_intercept(slope_fm, array_point_dw)
x_cross, y_cross = func_cross((0., slope_fm), (0., intercept_fm))
graph0.func_graph_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="b")

#------------------------------------------------------------------
#### 0. assumptions for slip line
#------------------------------------------------------------------
sl_angle = 30. / 360. * 2. * np.pi ### [rad]: slip line angle from horizontal axis (theta axis)
slope_sl = math.atan(sl_angle)
intercept_sl = func_intercept(slope_sl, array_point_dw)
x_cross, y_cross = func_cross((0., slope_sl), (rde_l, intercept_sl))
graph0.func_graph_add((array_point_dw[0], x_cross), (array_point_dw[1], y_cross), color="b")

#------------------------------------------------------------------
#### 0. assumptions for oblique-shock
#------------------------------------------------------------------
os_angle = 60. / 360. * 2. * np.pi ### [rad]: slip line angle from horizontal axis (theta axis)
slope_os = math.atan(os_angle)
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
array_point_dw_2nd = (2.*np.pi-1.*math.tan(fm_angle), dw_height)
graph0.func_graph_add((2.*np.pi, array_point_dw_2nd[0]), (0., array_point_dw_2nd[1]), color="r")

#------------------------------------------------------------------
#### 0. assumptions for fresh-mixture (2nd)
#------------------------------------------------------------------
intercept_fm_2nd = func_intercept(slope_fm, array_point_dw_2nd)
x_cross, y_cross = func_cross((0., slope_fm), (0., intercept_fm_2nd))
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



#------------------------------------------------------------------
#### 1. characteristic lines -1st
#------------------------------------------------------------------
num_ch = 5 # number of characteristic lines




## デトネーション波の角度とすべり面の角度が入ってしまうので追加して削除
# もっと賢い方法教えてください




array_theta = np.zeros((num_ch, num_ch))
array_neu = np.zeros((num_ch, num_ch))
array_beta = np.zeros((num_ch, num_ch))
array_Mach = np.zeros((num_ch, num_ch))


## for C-
array_theta[0] = np.linspace((dw_angle-np.pi/2.),sl_angle,num_ch)
array_neu[0] = array_theta[0]-(dw_angle-np.pi/2.)


## for C+
array_theta = np.transpose(array_theta)
array_theta[0] = np.linspace((dw_angle-np.pi/2.),0.,num_ch)
array_theta = np.transpose(array_theta)
# print(array_theta*360./2./np.pi)
array_neu = np.transpose(array_neu)
array_neu[0] = np.linspace(0.,(dw_angle-np.pi/2.),num_ch)
array_neu = np.transpose(array_neu)
# print(array_neu*360./2./np.pi)
### i方向（横）にtheta-neu=const.確認
### j方向（縦）にtheta+neu=const.確認

for i0 in range(int(num_ch)):
    array_Mach[0][i0] = func_neu2Mach(array_neu[0][i0])
    array_Mach[i0][0] = func_neu2Mach(array_neu[i0][0])
    array_beta[0][i0] = math.asin(1./array_Mach[0][i0])
    array_beta[i0][0] = math.asin(1./array_Mach[i0][0])

for j1 in range(1, int(num_ch)):
    for i1 in range(1, int(num_ch)):
        array_theta[j1][i1] = (array_neu[j1-1][i1]-array_neu[j1][i1-1])/2. \
            +(array_theta[j1-1][i1]+array_theta[j1][i1-1])/2.
        array_neu[j1][i1] = (array_neu[j1-1][i1]+array_neu[j1][i1-1])/2. \
            +(array_theta[j1-1][i1]-array_theta[j1][i1-1])/2.
        array_Mach[j1][i1] = func_neu2Mach(array_neu[j1][i1])
        array_beta[j1][i1] = math.asin(1./array_Mach[j1][i1])

# print("///array_theta///")
# print(array_theta*360./2./np.pi)
# print("///array_neu///")
# print(array_neu*360./2./np.pi)
# print("///array_Mach///")
# print(array_Mach)
# print("///array_beta///")
# print(array_beta*360./2./np.pi)

array_alpha_plus = np.zeros((num_ch, num_ch))
array_alpha_minus = np.zeros((num_ch, num_ch))
array_intercept_plus = np.zeros((num_ch, num_ch))
array_intercept_minus = np.zeros((num_ch, num_ch))
array_point = np.zeros((num_ch+1, num_ch+1, 2))

### array_point 境界条件
array_point[0][:] = [array_point_dw]*(num_ch+1)
array_point[0][0] = None
array_point[0][-1] = None

for j1 in range(1,int(num_ch)):
    for i1 in range(int(num_ch-1)):
        array_alpha_plus[j1][i1] = ((array_theta[j1-1][i1]+array_beta[j1-1][i1]) \
            + (array_theta[j1][i1]+array_beta[j1][i1]))/2.
for j1 in range(int(num_ch-1)):
    for i1 in range(1,int(num_ch)):
        array_alpha_minus[j1][i1] = ((array_theta[j1][i1-1]-array_beta[j1][i1-1]) \
            + (array_theta[j1][i1]-array_beta[j1][i1]))/2.

for j1 in range(int(num_ch-1)):
    array_point[j1+1][0] = np.zeros((1,2))
    for i1 in range(int(num_ch-1)):
        array_intercept_minus[j1][i1+1] = func_intercept(math.tan(array_alpha_minus[j1][i1+1]),array_point[j1][i1+1])
        array_intercept_plus[j1+1][i1] = func_intercept(math.tan(array_alpha_plus[j1+1][i1]),array_point[j1+1][i1])
        x_cross, y_cross = func_cross([math.tan(array_alpha_minus[j1][i1+1]),math.tan(array_alpha_plus[j1+1][i1])], \
            [array_intercept_minus[j1][i1+1],array_intercept_plus[j1+1][i1]])
        array_point[j1+1][i1+1] = np.array([x_cross, y_cross])

# print("///array_alpha_plus///")
# print(array_alpha_plus*360./2./np.pi)
# print("///array_alpha_minus///")
# print(array_alpha_minus*360./2./np.pi)
# print("///array_intercept_plus///")
# print(array_intercept_plus)
# print("///array_intercept_minus///")
# print(array_intercept_minus)
# print("///array_point_A///")
# print(array_point)

for j1 in range(int(num_ch-1)):
    for i1 in range(int(num_ch-1)):
        graph0.func_graph_add((array_point[j1][i1+1][0], array_point[j1+1][i1+1][0]), \
            (array_point[j1][i1+1][1], array_point[j1+1][i1+1][1]))
        graph0.func_graph_add((array_point[j1+1][i1][0], array_point[j1+1][i1+1][0]), \
            (array_point[j1+1][i1][1], array_point[j1+1][i1+1][1]))
        

### 本当であれば一般化したかったが
### 無理矢理燃焼器壁面＆すべり面の境界を代入し申す
### 燃焼器底部との交点
for i1 in range(1, int(num_ch)):
    x_cross, y_cross = func_cross([0.,math.tan(array_alpha_minus[-2][i1])], \
        [0.,array_intercept_minus[-2][i1]]) # array_alpha(or intercept)_minusの格納方向がくそなので
    array_point[-1][i1] = np.array([x_cross, y_cross])
### すべり面との交点
for j1 in range(1, int(num_ch)):
    x_cross, y_cross = func_cross([slope_sl,math.tan(array_alpha_plus[j1][-2])], \
        [intercept_sl,array_intercept_plus[j1][-2]])
    array_point[j1][-1] = np.array([x_cross, y_cross])
### 描画
for i1 in range(int(num_ch)):
    graph0.func_graph_add((array_point[-2][i1][0], array_point[-1][i1][0]), \
        (array_point[-2][i1][1], array_point[-1][i1][1]))
for j1 in range(int(num_ch)):
    graph0.func_graph_add((array_point[j1][-2][0], array_point[j1][-1][0]), \
        (array_point[j1][-2][1], array_point[j1][-1][1]))











graph0.func_show()















