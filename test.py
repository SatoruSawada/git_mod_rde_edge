import numpy as np
import math

### 同じやん，いつか統一して
### こいつが一般版
def func_cross_gas_dynamics(point1, point2, lambda1, lambda2): ### ごめん，くそみたいな名前で
    x = ((point1[1]-point2[1])-(lambda2*point1[0]-lambda1*point2[0]))/(lambda1-lambda2)
    y = point1[1] - lambda2 * point1[0] + lambda2 * x
    return x, y

### よく理解していない
### とりあえず vol.2 p.203 の例題にしたがって書いてみる
def func_func_MEPC_theta3(theta1, theta2, point1, point2, lambda_12, eps=10e-6):
    lambda_o = math.tan((theta1+theta2)/2.)
    delta = 1.0
    while delta >= eps:
        x, y = func_cross_gas_dynamics(point1, point2, lambda_o,lambda_12)
        theta = theta2 + (y-point2[1])/(point1[1]-point2[1]) * (theta1-theta2)
        lambda_o_new = math.tan(theta)
        delta = abs((lambda_o-lambda_o_new)/lambda_o)
        lambda_o = lambda_o_new
    return x, y, theta


theta1 = 18.191 / 360. *2.*np.pi
theta2 = 16.422 / 360. *2.*np.pi
x1 = 0.13146
y1 = 0.040118
x2 = 0.135683
y2 = 0.037123
lambda2 = -0.70921
x4 = 0.141335
y4 = 0.04061






