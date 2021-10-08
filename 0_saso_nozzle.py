"""
目的：
絵の描き方を調べる
Mach_number でコンター図が描けるのかどうかを調べる

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

import math
# import cantera as ct
import numpy as np
import matplotlib.pylab as plt
np.set_printoptions(precision=0)

### =================
### === cantera
### =================


class CL_graph_setting:
    def __init__(self):
        #### ==================================================================================================================
        #### setting
        #### ==================================================================================================================
        #### x軸
        self.x_label = 'radial direction [-]'
        self.x_min = -0.4        #### x軸最小値
        self.x_max = 3.        #### x軸最大値，目盛りの表示の都合でここに 0.00001 % 加算し申す
        self.x_main_dis = 0.5   #### x軸主目盛り間隔
        # x_sub_num = 5       #### x軸主目盛り間の小目盛りの個数
        #### y軸
        self.y_label = 'axial direction [-]'
        self.y_min = -2.        #### y軸最小値
        self.y_max = 2.        #### y軸最大値，目盛りの表示の都合でここに 0.00001 % 加算し申す
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
    
    def func_fill_add(self, list_x, list_y, color_value):
        #### ==================================================================================================================
        #### filling space surrounded by graphs
        #### ==================================================================================================================
        ###   list_x  : 四つ角のx座標
        ###   list_y  : 四つ角のy座標
        ### color(str): コンターにするための色（0.0 ~ 1.0 のグレースケール）
        self.ax.fill(list_x, list_y, color=str(color_value))

    def func_show(self):
        plt.show()

#### ================================================
#### 計算
#### ================================================


#------------------------------------------------------------------
#### 0. parameters
#------------------------------------------------------------------
gamma = 1.4 # 比熱比[-]



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
xa = 0.8 # 膨張部最短部のx座標[-]
theta_max_deg = 10. # 膨張最短部の水平軸に対する角度[deg]
theta_max_rad = theta_max_deg * 2 * np.pi / 360.

def func_eq_10_32(x):
    eq_10_32_right = math.tan(theta_max_rad)*((x/xa)**2.-(1./3.)*(x/xa)**3.)
    y = xa * eq_10_32_right + y0
    return y

list_x_func_eq_10_32 = [(xa) / number_func_eq_10_32 * i2 \
    for i2 in range(int(number_func_eq_10_32+1))] # 描画するためのx座標
list_y_func_eq_10_32 = [func_eq_10_32(i3) \
    for i3 in list_x_func_eq_10_32] # 描画するためのy座標

graph0.func_graph_add(list_x_func_eq_10_32, list_y_func_eq_10_32) # 膨張部描画
graph0.func_x_axis()

"""
(パデューのRDE解析をする上でコメント)
今回のノズル　ー＞　膨張部が複数個所で偏向してイクスパンションファンが発生する
パデュー解析　ー＞　一点で複数のイクスパンションを描く必要がある？？？
"""



#------------------------------------------------------------------
#### 5_preparation. contraction section
#------------------------------------------------------------------
## 無理矢理微分を計算し申す
number_expansion_fan = 8 # number_func_eq_10_32 とそんなに変わらないけれども

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
    # print(M_result)
    return M_result

def func_intercept(slope, list_point_A):
    x_A = list_point_A[0]
    y_A = list_point_A[1]
    intercept = y_A - slope * x_A
    return intercept 



array_theta = np.zeros((number_expansion_fan+1, number_expansion_fan+1+1))
array_neu = np.zeros((number_expansion_fan+1, number_expansion_fan+1+1))
array_Mach = np.zeros((number_expansion_fan+1, number_expansion_fan+1+1))
array_beta = np.zeros((number_expansion_fan+1, number_expansion_fan+1+1))



#------------------------------------------------------------------
#### 5. contraction section
#------------------------------------------------------------------

array_theta_base = np.empty((0, 1))
for i4 in range(number_expansion_fan+1):
    x_i4 = xa/number_expansion_fan*(i4)
    array_theta_base = np.append(array_theta_base, [math.atan(func_differencial_eq_10_32(x_i4))])
# print(array_theta_base*360/2./np.pi)
# print(math.atan(func_differencial_eq_10_32(0.8))/2./np.pi*360.)

#------------------------------------------------------------------
#### 6. create theta & neu through eq.(10.24) & (10.25)
#------------------------------------------------------------------
## 横にi並んでいく
## 縦にj並んでいく
for j6 in range(number_expansion_fan+1):
    for i6 in range(j6, number_expansion_fan+1): ### 重複OK
        theta_new = array_theta_base[i6] - array_theta_base[j6]
        neu_new = array_theta_base[i6] + array_theta_base[j6]
        Mach_new = func_neu2Mach(neu_new)
        beta_new = math.asin(1./Mach_new)
        array_theta[j6][i6] = theta_new
        array_neu[j6][i6] = neu_new
        array_Mach[j6][i6] = Mach_new
        array_beta[j6][i6] = beta_new
# print(array_theta)

#------------------------------------------------------------------
#### 7. create section to kill reflection wave through eq.(10.9)
#------------------------------------------------------------------
## eq.(10.9) からしてこれって theta_edge -> 違う
## theta_max_2_theta(n+1) の範囲の直線はなんでもいい -> 違うかも
## theta_max_2_theta(n+1) の範囲の直線は theta_max の直線？ 
## -> あくまで最短部でイクスパンションファンがでているのはそこがまだ滑らかに変化しているため

## 直接array_theta(or neu)に追加する
theta_kill = theta_max_rad
array_theta[0][-1] = theta_max_rad
neu_kill = array_neu[0][-2] - array_theta[0][-1] + theta_kill
array_neu[0][-1] = neu_kill
Mach_kill = func_neu2Mach(neu_kill)
array_Mach[0][-1] = Mach_kill
beta_kill = math.asin(1./Mach_kill)
array_beta[0][-1] = beta_kill

# for i7 in range(1, number_expansion_fan-1):
# 縦に初期値をすでに計算して与えているため
# 横にtheta(or neu)_kill を含むため1つ多い
# 縦横の数に限らず，全て縦方向の数に依存する
# 各行の最後の数字 -> kill
# 各行の最後の数字より左の数字 -> edge
# つまり，[ , , , ..., edge, kill, 0, ..., 0]
for i7 in range(1, number_expansion_fan):
    theta_kill = (array_neu[i7-1][-1]-array_neu[i7-1][-2])/2. + (array_theta[i7-1][-1]+array_theta[i7-1][-2])/2.
    neu_kill = (array_neu[i7-1][-1]+array_neu[i7-1][-2])/2. + (array_theta[i7-1][-1]-array_theta[i7-1][-2])/2.
    Mach_kill = func_neu2Mach(neu_kill)
    beta_kill = math.asin(1./Mach_kill)
    array_theta[i7][-1] = theta_kill
    array_neu[i7][-1] = neu_kill
    array_Mach[i7][-1] = Mach_kill
    array_beta[i7][-1] = beta_kill

print("==========================array_theta")
print(array_theta.shape)
print(array_theta*360/2./np.pi)

print("==========================array_neu")
print(array_neu.shape)
print(array_neu*360/2./np.pi)

print("==========================array_Mach")
print(array_Mach.shape)
print(array_Mach)

print("==========================array_beta")
print(array_beta.shape)
print(array_beta*360/2./np.pi)

#------------------------------------------------------------------
#### 8. depict characteristic lines through eq.(10.5) & (10.6)
#------------------------------------------------------------------

array_alpha_plus = np.zeros((number_expansion_fan+1,number_expansion_fan+1))
array_alpha_minus = np.zeros((number_expansion_fan+1,number_expansion_fan+1))
array_slope_plus = np.zeros((number_expansion_fan+1,number_expansion_fan+1))
array_slope_minus = np.zeros((number_expansion_fan+1,number_expansion_fan+1))
array_intercept_plus = np.zeros((number_expansion_fan+1,number_expansion_fan+1))
array_intercept_minus = np.zeros((number_expansion_fan+1,number_expansion_fan+1))

## 下上　下上　下上　下上
## 通る点と傾き
## 通る点を上書き（その代わり，徐々に上書きされない交点が残る，でも使わないからいい）
array_x = np.zeros((number_expansion_fan))
array_y = np.zeros((number_expansion_fan))

for j8 in range(number_expansion_fan):
    x_j8 = xa / number_expansion_fan * j8
    array_x[j8] = x_j8
    array_y[j8] = func_eq_10_32(x_j8)


## 横にi並んでいく
## 縦にj並んでいく
for j8 in range(number_expansion_fan):
    slope_plus = 0. # それぞれのはじめはx軸との交点
    intercept_plus = 0. # それぞれのはじめはx軸との交点
    for i8 in range(j8, number_expansion_fan): ### 重複OK
        alpha_minus = ((array_theta[j8][i8]-array_beta[j8][i8])+(array_theta[j8][i8+1]-array_beta[j8][i8+1]))/2. #+0.15*j8 # 下向き
        array_alpha_minus[j8][i8] = alpha_minus
        slope_minus = math.tan(alpha_minus)
        array_slope_minus[j8][i8] = slope_minus
        intercept_minus = func_intercept(slope_minus, (array_x[i8], array_y[i8]))
        array_intercept_minus[j8][i8] = intercept_minus
        cross_x, cross_y = func_cross((slope_plus, slope_minus), (intercept_plus, intercept_minus))
        graph0.func_graph_add((cross_x,array_x[i8]), (cross_y,array_y[i8])) ### 描いた
        array_x[i8], array_y[i8] = cross_x, cross_y

        alpha_plus = ((array_theta[j8][i8+1]+array_beta[j8][i8+1])+(array_theta[j8+1][i8+1]+array_beta[j8+1][i8+1]))/2. # 上向き
        array_alpha_plus[j8][i8] = alpha_plus
        slope_plus = math.tan(alpha_plus)
        array_slope_plus[j8][i8] = slope_plus
        intercept_plus = func_intercept(slope_plus, (array_x[i8], array_y[i8]))
        array_intercept_plus[j8][i8] = intercept_plus
    
    for k8 in range(j8, number_expansion_fan-1):
        graph0.func_graph_add((array_x[k8],array_x[k8+1]), (array_y[k8],array_y[k8+1])) ### 描いた

# print("==============array_alpha_plus")
# print(array_alpha_plus.shape)
# print(array_alpha_plus*360/2./np.pi)

# print("==============array_alpha_minus")
# print(array_alpha_minus.shape)
# print(array_alpha_minus*360/2./np.pi)

# print("============================array_slope_plus")
# print(array_slope_plus.shape)
# print(array_slope_plus)

# print("============================array_slope_minus")
# print(array_slope_minus.shape)
# print(array_slope_minus)

# print("==========================================array_intercept_plus")
# print(array_intercept_plus.shape)
# print(array_intercept_plus)

# print("==========================================array_intercept_minus")
# print(array_intercept_minus.shape)
# print(array_intercept_minus)


#### testing filling color_function














































































slope_kill = math.tan(theta_max_rad)
list_point_A = (xa, func_eq_10_32(xa))
### 相殺部を描画
# for i9 in range(number_expansion_fan):
#     intercept_kill = func_intercept(slope_kill,list_point_A)
#     x_cross, y_cross = func_cross((slope_kill, array_slope_plus[i9][-1]), (intercept_kill, array_intercept_plus[i9][-1]))
#     print(x_cross,y_cross)
#     graph0.func_graph_add((list_point_A[0], x_cross), (list_point_A[1], y_cross))
#     slope_kill = math.tan(array_theta[i9][-1])
#     x_cross, y_cross = func_cross((slope_kill, array_slope_plus[i9][-1]), (intercept_kill, array_intercept_plus[i9][-1]))
#     list_point_A = (x_cross, y_cross)

### 20210728_comment
### あとはここの相殺部を計算してい描画するだけ
### 「x_cross」，「y_cross」と「x_cross_new」，「y_cross_new」を組み合わせて更新していくだけ














graph0.func_show() # plt.show って一度しか使えないの？？？












