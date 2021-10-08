"""
目的：
絵の描き方を調べる

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
from matplotlib import colors
import numpy as np
import matplotlib.pylab as plt

class CL_graph_setting:
    def __init__(self):
        #### ==================================================================================================================
        #### setting
        #### ==================================================================================================================
        #### x軸
        self.x_label = 'time from ignition, $\it{t}$, s'
        self.x_min = -0.4        #### x軸最小値
        self.x_max = 1.2        #### x軸最大値，目盛りの表示の都合でここに 0.00001 % 加算し申す
        self.x_main_dis = 0.5   #### x軸主目盛り間隔
        # x_sub_num = 5       #### x軸主目盛り間の小目盛りの個数
        #### y軸
        self.y_label = 'self-luminescence angle''\n'r'ratio, $\it{ε}$, deg./360deg.'
        self.y_min = -0.1        #### y軸最小値
        self.y_max = 1.5        #### y軸最大値，目盛りの表示の都合でここに 0.00001 % 加算し申す
        self.y_main_dis = 0.2   #### y軸主目盛り間隔
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
        self.fig_x_size = 8.     #### グラフ全体のx軸方向の大きさ
        self.fig_y_size = 4.     #### グラフ全体のy軸方向の大きさ
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
        self.i = 0
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

    def func_x_axis(self):
        list_x_axis_x = [-10e5, 10e5]
        list_x_axis_y = [0., 0.]
        self.ax.plot(list_x_axis_x, list_x_axis_y, color='k')

    def func_graph_add(self, list_x, list_y):
        #### ==================================================================================================================
        #### graph depict
        #### ==================================================================================================================
        # self.ax.scatter(list_x, list_y, s=marker_size, c=marker_color)
        self.ax.plot(list_x, list_y)
        self.i = self.i + 1
        print('for comfirming iteration number {}'.format(self.i))

    def func_show(self):
        plt.show()







#### ================================================
#### 計算
#### ================================================



#------------------------------------------------------------------
#### 0. parameters
#------------------------------------------------------------------
gamma = 1.32 # 比熱比[-]



#------------------------------------------------------------------
#### 1. Me
#------------------------------------------------------------------
Me = 3. # 出口Mach数[-]
print("(1) Me = {0} [-]".format(Me))



#------------------------------------------------------------------
#### 2. y0
#------------------------------------------------------------------
ye_over_y0 = 1/Me * (((gamma-1.)*Me**2.+2.)/(gamma+1.)) ** ((gamma+1.)/(2*(gamma-1.)))
y0 = 1.
ye = ye_over_y0 * y0 # 無次元ノズル出口半径[-]
print("(2) y0 = {0} [-] /// ye = {1} [-]".format(y0, ye))



#------------------------------------------------------------------
#### 3. inlet section (subsonic flow, before throat)
#------------------------------------------------------------------
## 必ずスロート部で音速になるためここの部分を描画するしても意味がない
## 亜音速部と任意の軸の描画の練習
yi = 1.2 * y0 #（任意）入口マッハ数[-]
x_4_yi = - 0.2 # yi の存在する x 座標（膨張部の長さを1とする）[-]
coef_4_func_subsonic = (yi-y0)/(x_4_yi**2.) # 亜音速部を二次関数(ax**2+b)としたときの a
print("(3) yi = {0} [-] /// x_4_yi = {1} [-] /// coef_4_func_subsonic = {2}".format(yi, x_4_yi, coef_4_func_subsonic))
## 描画する亜音速部の関数
func_subsonic = lambda x: coef_4_func_subsonic * x ** 2. + y0 # 亜音速部の二次関数(ax**2+b)
number_func_subsonic = 10 # 亜音速部の関数を描画する点数
list_x_func_subsonic = [x_4_yi + (yi-y0) / number_func_subsonic * i0 \
    for i0 in range(int(number_func_subsonic+1))] # 描画するためのx座標
list_y_func_subsonic = [func_subsonic(i1) \
    for i1 in list_x_func_subsonic] # 描画するためのy座標

graph0 = CL_graph_create()
graph0() # reset
graph0.func_graph_add(list_x_func_subsonic, list_y_func_subsonic)



#------------------------------------------------------------------
#### 4. expansion section
#------------------------------------------------------------------
number_func_eq_10_32 = 10. # x0 - xa までを分割する数
xa = 0.2 # 膨張部最短部のx座標[-]
theta_max_deg = 20. # 膨張最短部の水平軸に対する角度[deg]
theta_max_rad = theta_max_deg * 2 * np.pi / 360.

def func_eq_10_32(x):
    eq_10_32_right = math.tan(theta_max_rad*((x/xa)**2.-(1./3.)*(x/xa)**3.))
    y = xa * eq_10_32_right + y0
    return y

list_x_func_eq_10_32 = [(xa) / number_func_eq_10_32 * i2 \
    for i2 in range(int(number_func_eq_10_32+1))] # 描画するためのx座標
list_y_func_eq_10_32 = [func_eq_10_32(i3) \
    for i3 in list_x_func_eq_10_32] # 描画するためのy座標

print(list_x_func_eq_10_32)
print(list_y_func_eq_10_32)

graph0.func_graph_add(list_x_func_eq_10_32, list_y_func_eq_10_32) # 膨張部描画
graph0.func_x_axis()

"""
(パデューのRDE解析をする上でコメント)
今回のノズル　ー＞　膨張部が複数個所で偏向してイクスパンションファンが発生する
パデュー解析　ー＞　一点で複数のイクスパンションを描く必要がある？？？
"""



#------------------------------------------------------------------
#### 5. contraction section
#------------------------------------------------------------------
## 無理矢理微分を計算し申す
number_expansion_fan = 10 # number_func_eq_10_32 とそんなに変わらないけれども

def func_differencial_eq_10_32(x):
    dx = 10e-6
    # dif = ((func_eq_10_32(x+dx) - func_eq_10_32(x)) / (dx))
    # dif = ((func_eq_10_32(x) - func_eq_10_32(x-dx)) / (dx))
    dif = ((func_eq_10_32(x+dx/2.) - func_eq_10_32(x-dx/2.)) / (dx)) # 中間を取ります，方法の名前忘れ申した
    return dif




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


array_x_expansion = np.empty((0,1))
array_y_expansion = np.empty((0,1))

array_theta = np.empty((0,1))

array_neu = np.empty((0,1))

array_alpha_plus = np.empty((0,1))
array_alpha_minus = np.empty((0,1))

array_beta = np.empty((0,1))

array_slope_plus = np.empty((0,1))
array_slope_minus = np.empty((0,1))

array_intercept_plus = np.empty((0,1))
array_intercept_minus = np.empty((0,1))


start = 0
end = 1


## 一巡目
## minus -> 下向き
for i4 in range(number_expansion_fan): ### ここの「1」を「+1」ずつfor構文で変化させていくつもり
    x_i4_a = xa/number_expansion_fan*(i4)
    x_i4_b = xa/number_expansion_fan*(i4+1)
    array_x_expansion = np.append(array_x_expansion, [x_i4_a])
    array_y_expansion = np.append(array_y_expansion, [func_eq_10_32(x_i4_a)])
    array_theta = np.append(array_theta, [math.atan(func_differencial_eq_10_32(x_i4_b))])
    array_neu = np.append(array_neu, [math.atan(func_differencial_eq_10_32(x_i4_b))])
    array_alpha_minus = None
    array_beta = np.append(array_beta, [math.asin(1./func_neu2Mach(array_neu[-1]))])
    array_slope_minus = np.append(array_slope_minus, [-math.tan(np.pi/2.-array_theta[int(-1)])])
    array_intercept_minus = np.append(array_intercept_minus, [func_eq_10_32(x_i4_b)-array_slope_minus[-1]*x_i4_b])

print("intercept_minus")
print(array_intercept_minus)

## from wall
list_x_0 = [0,0] # reset
list_y_0 = [0,0] # reset
list_x_0[int(start)] = array_x_expansion[0]
list_y_0[int(start)] = array_y_expansion[0]
list_x_0[int(end)] = - array_intercept_minus[0] / array_slope_minus[0]
list_y_0[int(end)] = 0.

graph0.func_graph_add(list_x_0, list_y_0) # 膨張部描画








## 二巡目
# 対称面での反射

array_slope_plus = np.append(array_slope_plus, [-array_slope_minus[int(0)]]) # plus for alpha_plus
array_intercept_plus = np.append(array_intercept_plus, \
    [-array_slope_plus[int(0)]*(- array_intercept_minus[0] / array_slope_minus[0])])

cross_x0, cross_y0 = func_cross([0., array_slope_minus[int(0)]],\
    [0., array_intercept_minus[int(0)]])
cross_x1, cross_y1 = func_cross([array_slope_minus[int(1)], array_slope_plus[int(0)]],\
    [array_intercept_minus[int(1)], array_intercept_plus[int(0)]])

list_x_0[int(start)] = cross_x0
list_y_0[int(start)] = cross_y0
list_x_0[int(end)] = cross_x1
list_y_0[int(end)] = cross_y1
graph0.func_graph_add(list_x_0, list_y_0)


## from wall
cross_x0, cross_y0 = func_cross([array_slope_minus[int(1)], array_slope_plus[int(0)]],\
    [array_intercept_minus[int(1)], array_intercept_plus[int(0)]])

list_x_0[int(start)] = array_x_expansion[1]
list_y_0[int(start)] = array_y_expansion[1]
list_x_0[int(end)] = cross_x0
list_y_0[int(end)] = cross_y0
graph0.func_graph_add(list_x_0, list_y_0)










## 交点2対称面(plus)
## alpha_minus
array_theta_new = np.empty((0,1))
array_neu_new = np.empty((0,1))

# wall
array_theta_new = np.append(array_theta_new, [0.])
array_neu_new = np.append(array_neu_new, [2.*array_theta[0]])


# not wall
neu_1 = array_neu[0]                           ## inclination 候補
theta_1 = array_theta[0]                       ## inclination 候補
beta_1 = array_beta[0]
print('neu_1 = {0} [rad] >>> {1}'.format(neu_1, func_neu2Mach(neu_1)))
print('theta_1 = {0} /// beta_1 = {1}'.format(theta_1, beta_1))

array_beta = np.empty((0,1))
neu_2 = array_neu_new[0]                       ## inclination 候補
theta_2 = array_theta_new[0]                   ## inclination 候補
beta_2 = math.asin(1./func_neu2Mach(neu_2))
array_beta = np.append(array_beta, [beta_2])
print('neu_2 = {0} [rad] >>> {1}'.format(neu_2, func_neu2Mach(neu_2)))
print('theta_2 = {0} /// beta_2 = {1}'.format(theta_2, beta_2))

neu_3 = (neu_1+neu_2)/2. + (theta_1-theta_2)/2.
theta_3 = (neu_1-neu_2)/2. + (theta_1+theta_2)/2.
beta_3 = math.asin(1./func_neu2Mach(neu_3))
array_beta = np.append(array_beta, [beta_3])
print('neu_3 = {0} [rad] >>> {1}'.format(neu_3, func_neu2Mach(neu_3)))
print('theta_3 = {0} /// beta_3 = {1}'.format(theta_3, beta_3))

array_theta_new = np.append(array_theta_new, [theta_3])
array_neu_new = np.append(array_neu_new, [neu_3])

alpha_plus = ((theta_1+beta_1)+(theta_3+beta_3))/2.
alpha_minus = ((theta_2-beta_2)+(theta_3-beta_3))/2.

array_slope_minus_new = np.empty((0,1))
array_intercept_minus_new = np.empty((0,1))
array_slope_minus_new = np.append(array_slope_minus_new, [math.tan(alpha_minus)]) # 下向き
array_intercept_minus_new = np.append(array_intercept_minus_new, [cross_y0-array_slope_minus_new[-1]*cross_x0]) # 下向き
print("alpha_minus = {} [rad]".format(alpha_minus))
print("tan(alpha_minus) = {} ".format(math.tan(alpha_minus)))

cross_x1, cross_y1 = func_cross([0., array_slope_minus_new[int(-1)]],\
    [0., array_intercept_minus_new[int(-1)]])

list_x_0[int(start)] = cross_x0
list_y_0[int(start)] = cross_y0
list_x_0[int(end)] = cross_x1
list_y_0[int(end)] = cross_y1
# graph0.func_graph_add(list_x_0, list_y_0)



# 交点2交点(plus)
# alpha_plus
array_slope_plus = np.append(array_slope_plus, [math.tan(alpha_plus)]) # 下向き
array_intercept_plus = np.append(array_intercept_plus, [cross_y0-array_slope_plus[0]*cross_x0]) # 下向き
print(array_slope_plus)
print(array_intercept_plus)
print("=================")
print(array_slope_minus[int(0)])
print(array_intercept_minus[int(0)])

cross_x2, cross_y2 = func_cross([array_slope_minus[int(2)], array_slope_plus[int(-1)]],\
    [array_intercept_minus[int(2)], array_intercept_plus[int(-1)]])

print(cross_x2, cross_y2)

list_x_0[int(start)] = cross_x0
list_y_0[int(start)] = cross_y0
list_x_0[int(end)] = cross_x2
list_y_0[int(end)] = cross_y2
# graph0.func_graph_add(list_x_0, list_y_0)













graph0.func_show() # plt.show って一度しか使えないの？？？



