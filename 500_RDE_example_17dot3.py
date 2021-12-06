"""
特性線を追加する

1. 上からの特性線が底面に到達するまで計算
2. 特性線追加
3. 下からの特性線が上側境界条件に到達するまでの計算

特性線の状態量をしっかりと計算
追加の特性線の状態に関してarray_minus は同じ値にすべきなのか？
"""

## 意味ないけれども
## (__init__) -> 書いてもいいところ
## (__call__) -> 書いてはいけないところ


import numpy as np
import matplotlib.pylab as plt
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)

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


### point 1
x_1 = 5.495
y_1 = 26.020
V_1 = 1728.9
theta_1 = 24.090 /360.*2.*np.pi
p_1 = 11.545 * 10e4
rho_1 = 1.6240
### ===========
a_fr_1 = np.sqrt(1.2*p_1/rho_1)

### point 2
### none

### point 3
x_3 = 5.085
y_3 = 26.080
V_3 = 1726.8
theta_3 = 24.000 /360.*2.*np.pi
p_3 = 11.604 * 10e4
rho_3 = 1.6309
### ===========
a_fr_3 = np.sqrt(1.2*p_3/rho_3)

### point 4
x_4 = 5.283
y_4 = 26.170
theta_4 = 25.000 /360.*2.*np.pi


#------------------------------------------------------------------
#### 1. characteristic lines -1st
#------------------------------------------------------------------
### num_ch_up & num_ch_down が小さすぎても問題（num_ch_up & num_ch_down >= 7）
### i方向（横）にtheta-neu=const.確認
### j方向（縦）にtheta+neu=const.確認
num_ch_up = 20 # number of initial characteristic lines (upper side)
num_ch_down = 10 # number of initial characteristic lines (down side)
num_ch_add = 5

array_zero0 = np.zeros((int(num_ch_down),int(num_ch_up-1)))
array_zero1 = np.zeros((int(num_ch_up + num_ch_down - 1),int(num_ch_down + num_ch_add)))

### x for characteristics
array_x_up = np.ones((int(num_ch_up))) * x_3 ### point 3
array_x = np.flipud(np.diag(array_x_up))
array_x_down = np.ones((int(num_ch_down))) * x_1 ### point 1
array_x_down = np.transpose(np.array([array_x_down]))
array_x_down = np.hstack((array_x_down,array_zero0))
array_x_down = np.delete(array_x_down,0,0)
array_x = np.vstack((array_x,array_x_down))
array_x = np.hstack((array_x, array_zero1))
del array_x_up
del array_x_down

### y for characteristics
array_y_up = np.ones((int(num_ch_up))) * y_3 ### point 3
array_y = np.flipud(np.diag(array_y_up))
array_y_down = np.ones((int(num_ch_down))) * y_1 ### point 1
array_y_down = np.transpose(np.array([array_y_down]))
array_y_down = np.hstack((array_y_down,array_zero0))
array_y_down = np.delete(array_y_down,0,0)
array_y = np.vstack((array_y,array_y_down))
array_y = np.hstack((array_y, array_zero1))
del array_y_up
del array_y_down

### ============================================================================================ 20211018_sawada
### ============================================================================================ 初期値の設定がよくわからない
### ============================================================================================ とりあえずリーマン不変量で
array_M_up = np.zeros((int(num_ch_up)))
array_neu_up = np.zeros((int(num_ch_up)))
array_alpha_up = np.zeros((int(num_ch_up)))
array_p_up = np.ones((int(num_ch_up))) * p_3 ### point 3
array_t_up = np.zeros((int(num_ch_up)))
array_R_up = np.zeros((int(num_ch_up)))
array_rho_up = np.ones((int(num_ch_up))) * rho_3 ### point3
array_a_fr_up = np.zeros((int(num_ch_up))) * a_fr_3 ### point3
array_V_up = np.ones((int(num_ch_up))) * V_3 ### point 3
array_gamma_up = np.zeros((int(num_ch_up)))
array_theta_up = np.ones((int(num_ch_up))) * theta_3 ### point 3

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
array_neu = np.delete(array_neu,-1,0)

#============================================================================
array_M_down = np.zeros((int(num_ch_down)))
array_neu_down = np.zeros((int(num_ch_down)))
array_alpha_down = np.zeros((int(num_ch_down)))
array_p_down = np.ones((int(num_ch_down))) * p_1 ### point 1
array_t_down = np.zeros((int(num_ch_down)))
array_R_down = np.zeros((int(num_ch_down)))
array_rho_down = np.ones((int(num_ch_down))) * rho_1 ### point 1
array_a_fr_down = np.zeros((int(num_ch_down))) * a_fr_1 ### point 1
array_V_down = np.ones((int(num_ch_down))) * V_1 ### point 1
array_gamma_down = np.zeros((int(num_ch_down)))
array_theta_down = np.ones((int(num_ch_down))) * theta_1 ### point 1

#============================================================================
### =====
array_theta_down = np.transpose(np.array([array_theta_down]))
array_theta_down = np.hstack((array_theta_down,array_zero0))
array_theta_down = np.delete(array_theta_down,0,0)
array_theta = np.vstack((array_theta,array_theta_down))
array_theta = np.hstack((array_theta, array_zero1))
### =====
array_M_down = np.transpose(np.array([array_M_down]))
array_M_down = np.hstack((array_M_down,array_zero0))
array_M_down = np.delete(array_M_down,0,0)
array_M = np.vstack((array_M,array_M_down))
array_M = np.hstack((array_M, array_zero1))
### =====
array_alpha_down = np.transpose(np.array([array_alpha_down]))
array_alpha_down = np.hstack((array_alpha_down,array_zero0))
array_alpha_down = np.delete(array_alpha_down,0,0)
array_alpha = np.vstack((array_alpha,array_alpha_down))
array_alpha = np.hstack((array_alpha, array_zero1))
### =====
array_p_down = np.transpose(np.array([array_p_down]))
array_p_down = np.hstack((array_p_down,array_zero0))
array_p_down = np.delete(array_p_down,0,0)
array_p = np.vstack((array_p,array_p_down))
array_p = np.hstack((array_p, array_zero1))
### =====
array_t_down = np.transpose(np.array([array_t_down]))
array_t_down = np.hstack((array_t_down,array_zero0))
array_t_down = np.delete(array_t_down,0,0)
array_t = np.vstack((array_t,array_t_down))
array_t = np.hstack((array_t, array_zero1))
### =====
array_R_down = np.transpose(np.array([array_R_down]))
array_R_down = np.hstack((array_R_down,array_zero0))
array_R_down = np.delete(array_R_down,0,0)
array_R = np.vstack((array_R,array_R_down))
array_R = np.hstack((array_R, array_zero1))
### =====
array_rho_down = np.transpose(np.array([array_rho_down]))
array_rho_down = np.hstack((array_rho_down,array_zero0))
array_rho_down = np.delete(array_rho_down,0,0)
array_rho = np.vstack((array_rho,array_rho_down))
array_rho = np.hstack((array_rho, array_zero1))
### =====
array_a_fr_down = np.transpose(np.array([array_a_fr_down]))
array_a_fr_down = np.hstack((array_a_fr_down,array_zero0))
array_a_fr_down = np.delete(array_a_fr_down,0,0)
array_a_fr = np.vstack((array_a_fr,array_a_fr_down))
array_a_fr = np.hstack((array_a_fr, array_zero1))
### =====
array_V_down = np.transpose(np.array([array_V_down]))
array_V_down = np.hstack((array_V_down,array_zero0))
array_V_down = np.delete(array_V_down,0,0)
array_V = np.vstack((array_V,array_V_down))
array_V = np.hstack((array_V, array_zero1))
### =====
array_gamma_down = np.transpose(np.array([array_gamma_down]))
array_gamma_down = np.hstack((array_gamma_down,array_zero0))
array_gamma_down = np.delete(array_gamma_down,0,0)
array_gamma = np.vstack((array_gamma,array_gamma_down))
array_gamma = np.hstack((array_gamma, array_zero1))

### ============================================================================
del array_theta_up; del array_theta_down; del array_neu_up; del array_neu_down
del array_M_up; del array_M_down; del array_alpha_up; del array_alpha_down
del array_p_up; del array_p_down; del array_t_up; del array_t_down
del array_R_up; del array_R_down; del array_rho_up; del array_rho_down
del array_a_fr_up; del array_a_fr_down; del array_V_up; del array_V_down
del array_gamma_up; del array_gamma_down
### ====

### C+ 上のパラメーター？？？
array_T_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_p_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_theta_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_V_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_rho_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_y_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_a_fr_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_M_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_alpha_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_lambda_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_Q_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_S_plus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))

### C- 上のパラメーター？？？
array_T_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_p_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_theta_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_V_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_rho_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_y_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_a_fr_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_M_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_alpha_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_lambda_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_Q_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_S_minus = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))

### point3
array_x_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_y_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_theta_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_p_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_rho_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_V_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_a_fr_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_gamma_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))

### point_o
array_lambda_12 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_lambda_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_p_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_rho_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_V_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_a_fr_3 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_R_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_A_o = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_T_o1 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))
array_T_o2 = np.zeros((int(num_ch_up+num_ch_down-1), int(num_ch_up + num_ch_down + num_ch_add)))



### point 4 の設定
# array_x[num_ch_up-1][1] = x_4
# array_y[num_ch_up-1][1] = y_4
# array_theta[num_ch_up-1][1] = theta_4

### lambda_minus の計算
array_lambda_minus[num_ch_up-2][1] = (array_y[num_ch_up-2][1]-array_y[num_ch_up][0]) / \
    (array_x[num_ch_up-2][1]-array_x[num_ch_up][0])

print('array_x[num_ch_up-2][1] =',array_x[num_ch_up-2][1])
print('array_y[num_ch_up-2][1] =',array_y[num_ch_up-2][1])
print('array_x[num_ch_up][0] =', array_x[num_ch_up][0])
print('array_y[num_ch_up][0] =', array_y[num_ch_up][0])
print('lambda_minus =',array_lambda_minus[num_ch_up-2][1])



mode_top = 0
mode_bottom = 0
j_bottom = num_ch_up+num_ch_down-2
eps_c = 10e-10

### point 2 って保存しておくべきなのだろうか？とりあえず保存しない
### point 4 って保存しておくべきなのだろうか？とりあえず保存しない

for i in range(1,int(2)):### 20211022_sawada : 次の列の計算をしていないためにエラーが起きている
    for j in range(int(num_ch_up-i), int(num_ch_up-i+1)):
        if j < j_bottom:
            ### a_interior
            ### =====================================================================
            ### predictor
            ### =====================================================================
            
            ### Mach_line 24 -> _plus_2
            
            ### (c) 
            p_2 = (array_p[j][i-1]+array_p[j+1][i-1])/2.
            theta_2 = (array_theta[j][i-1]+array_theta[j+1][i-1])/2.
            V_2 = (array_V[j][i-1]+array_V[j+1][i-1])/2.
            rho_2 = (array_rho[j][i-1]+array_rho[j+1][i-1])/2.
            a_fr_2 = np.sqrt(1.2*p_2/rho_2)
            M_2 = V_2/a_fr_2
            alpha_2 = np.arcsin(1./M_2)
            theta_plus_24 = (theta_2+theta_4)/2.
            lambda_plus_24 = np.tan(theta_plus_24+alpha_2)
            ### （point 4 の傾きと座標）& 特性線 C- 1-3
            ### point 1 & lambda_minus
            ### point 4 & theta_4
            x_2, y_2 = func_cross_gas_dynamics(\
                (x_4, y_4),\
                    (array_x[j+1][i-1], array_y[j+1][i-1]),\
                        lambda_plus_24, \
                            array_lambda_minus[j-1][i])
            ### point 2 interpolated between 1 and 3
            
            n = 0
            while n <= 3:
                n+=1
                
                p_2 = array_p[j+1][i-1] + (x_2-array_x[j+1][i-1])/(array_x[j][i-1]-array_x[j+1][i-1])*(array_p[j][i-1]-array_p[j+1][i-1])
                theta_2 = array_theta[j+1][i-1] + (x_2-array_x[j+1][i-1])/(array_x[j][i-1]-array_x[j+1][i-1])*(array_theta[j][i-1]-array_theta[j+1][i-1])
                V_2 = array_V[j+1][i-1] + (x_2-array_x[j+1][i-1])/(array_x[j][i-1]-array_x[j+1][i-1])*(array_V[j][i-1]-array_V[j+1][i-1])
                rho_2 = array_rho[j+1][i-1] + (x_2-array_x[j+1][i-1])/(array_x[j][i-1]-array_x[j+1][i-1])*(array_rho[j][i-1]-array_rho[j+1][i-1])

                a_fr_2 = np.sqrt(1.2*p_2/rho_2)
                M_2 = V_2 /a_fr_2
                alpha_2 = np.arcsin(1./M_2)
                theta_plus_24 = (theta_2+theta_4)/2.
                lambda_plus_24 = np.tan(theta_plus_24+alpha_2)
                ### （point 4 の傾きと座標）& 特性線 C- 1-3
                ### point 1 & lambda_minus
                ### point 4 & theta_4
                x_2, y_2 = func_cross_gas_dynamics(\
                    (x_4, y_4),\
                        (array_x[j+1][i-1], array_y[j+1][i-1]),\
                            lambda_plus_24, \
                                array_lambda_minus[j-1][i])
                print('x_2 =', x_2, '/// y_2 =',y_2, '/// lambda_plus_24 =',lambda_plus_24)
            
            
            ### (c)
            theta_plus_24 = (theta_2+theta_4)/2.
            y_plus_2 = (y_2+y_4)/2.
            a_fr_2 = np.sqrt(1.2*p_2/rho_2)
            M_2 = V_2 / a_fr_2
            alpha_2 = np.arcsin(1./M_2)
            lambda_plus_24 = np.tan(theta_plus_24 + alpha_2)
            
            Q_plus_2 = np.sqrt(M_2**2.-1.)/(rho_2*V_2**2.)
            S_plus_2 = np.sin(theta_plus_24)/\
                (y_plus_2*M_2*np.cos(theta_plus_24+alpha_2))
            T_plus_2 = -S_plus_2*(x_4-x_2)+Q_plus_2*p_2+theta_2
            
            ### (d) - stream line 34
            array_a_fr[j][i-1] = np.sqrt(1.2*array_p[j][i-1]/array_rho[j][i-1])
            array_R_o[j][i-1] = array_rho[j][i-1]*array_V[j][i-1]
            array_A_o[j][i-1] = array_a_fr[j][i-1]**2.
            array_T_o1[j][i-1] = array_R_o[j][i-1]*array_V[j][i-1]+array_p[j][i-1]
            array_T_o2[j][i-1] = array_p[j][i-1]-array_A_o[j][i-1]*array_rho[j][i-1]
            
            ### (e)
            p_4 = (T_plus_2 - theta_4) / Q_plus_2
            V_4 = (array_T_o1[j][i-1]-p_4)/array_R_o[j][i-1]
            rho_4 = -(array_T_o2[j][i-1]-p_4)/array_A_o[j][i-1]
            
            print("===============================================================")
            print("predictor")
            print("===============================================================")
            print("lambda_plus_24 ", lambda_plus_24)
            print("x_2           ", x_2)
            print("y_2           ", y_2)
            print("R_o           ", array_R_o[j][i-1])
            print("A_o           ", array_A_o[j][i-1])
            print("T_o1          ", array_T_o1[j][i-1])
            print("T_o2          ", array_T_o2[j][i-1])
            print("Q_plus_2      ", Q_plus_2)
            print("S_plus_2      ", S_plus_2)
            print("T_plus_2      ", T_plus_2)
            print("p_4           ", p_4)
            print("V_4           ", V_4)
            print("rho_4         ", rho_4)




            ### set predictor
            theta = p_4
            delta_c = 1.0
            ### =====================================================================
            ### corrector : 全て入れなおせているのだろうか？
            ### ここで stream line 34 の間の値を計算して利用
            ### ここで 特性線 C+ 24 の間の値を計算して利用（average ）
            ### =====================================================================
            while delta_c >= eps_c:
                    
                ### Mach_line 24 -> _plus_2

                n = 0
                while n <= 3:
                    n+=1

                    p_2 = array_p[j+1][i-1] + (x_2-array_x[j+1][i-1])/(array_x[j][i-1]-array_x[j+1][i-1])*(array_p[j][i-1]-array_p[j+1][i-1])
                    theta_2 = array_theta[j+1][i-1] + (x_2-array_x[j+1][i-1])/(array_x[j][i-1]-array_x[j+1][i-1])*(array_theta[j][i-1]-array_theta[j+1][i-1])
                    V_2 = array_V[j+1][i-1] + (x_2-array_x[j+1][i-1])/(array_x[j][i-1]-array_x[j+1][i-1])*(array_V[j][i-1]-array_V[j+1][i-1])
                    rho_2 = array_rho[j+1][i-1] + (x_2-array_x[j+1][i-1])/(array_x[j][i-1]-array_x[j+1][i-1])*(array_rho[j][i-1]-array_rho[j+1][i-1])
                    a_fr_2 = np.sqrt(1.2*p_2/rho_2)
                    M_2 = V_2 /a_fr_2
                    
                    ### （point 4 の傾きと座標）& 特性線 C- 1-3
                    ### point 1 & lambda_minus
                    ### point 4 & theta_4
                    x_2, y_2 = func_cross_gas_dynamics(\
                        (x_4, y_4),\
                            (array_x[j+1][i-1], array_y[j+1][i-1]),\
                                lambda_plus_24, \
                                    array_lambda_minus[j-1][i])
                    print('x_2 =', x_2, '/// y_2 =',y_2, '/// lambda_plus_24 =',lambda_plus_24)
                
                    ### (c)
                    theta_plus_24 = (theta_2+theta_4)/2.
                    y_24 = (y_2+y_4)/2.
                    p_24 = (p_2+p_4)/2.
                    rho_24 = (rho_2+rho_4)/2.
                    V_24 = (V_2+V_4)/2.
                    
                    a_fr_24 = np.sqrt(1.2*p_24/rho_24)
                    M_24 = V_24 / a_fr_24
                    alpha_24 = np.arcsin(1./M_24)
                    lambda_plus_24 = np.tan(theta_plus_24 + alpha_24)
                
                Q_plus_2 = np.sqrt(M_24**2.-1.)/(rho_24*V_24**2.)
                S_plus_2 = np.sin(theta_plus_24)/\
                    (y_24*M_24*np.cos(theta_plus_24+alpha_24))
                T_plus_2 = -S_plus_2*(x_4-x_2)+Q_plus_2*p_2+theta_2
                
                
                ### (d) - stream line 34
                p_34 = (array_p[j][i-1]+p_4)/2.
                rho_34 = (array_rho[j][i-1]+rho_4)/2.
                V_34 = (array_V[j][i-1]+V_4)/2.
                a_fr_34 = np.sqrt(1.2*p_34/rho_34)
                array_R_o[j][i-1] = rho_34*V_34
                array_A_o[j][i-1] = a_fr_34**2.
                array_T_o1[j][i-1] = array_R_o[j][i-1]*V_34+p_34
                array_T_o2[j][i-1] = p_34-array_A_o[j][i-1]*rho_34
                
                ### (e)
                p_4 = (T_plus_2 - theta_4) / Q_plus_2
                V_4 = (array_T_o1[j][i-1]-p_4)/array_R_o[j][i-1]
                rho_4 = -(array_T_o2[j][i-1]-p_4)/array_A_o[j][i-1]

                ### delta_c
                theta_new = p_4
                delta_c = abs((theta-theta_new)/theta)
                theta = theta_new


                print("===============================================================")
                print("corrector")
                print("===============================================================")
                print("lambda_plus_24 ", lambda_plus_24)
                print("x_2           ", x_2)
                print("y_2           ", y_2)
                print("R_o           ", array_R_o[j][i-1])
                print("A_o           ", array_A_o[j][i-1])
                print("T_o1          ", array_T_o1[j][i-1])
                print("T_o2          ", array_T_o2[j][i-1])
                print("Q_plus_2      ", Q_plus_2)
                print("S_plus_2      ", S_plus_2)
                print("T_plus_2      ", T_plus_2)
                print("p_4           ", p_4)
                print("V_4           ", V_4)
                print("rho_4         ", rho_4)










