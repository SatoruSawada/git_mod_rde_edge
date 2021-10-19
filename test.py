import numpy as np
import math

### 同じやん，いつか統一して
### こいつが一般版
def func_cross_gas_dynamics(point1, point2, lambda1, lambda2): ### ごめん，くそみたいな名前で
    x = ((point1[1]-point2[1])-(lambda1*point1[0]-lambda2*point2[0]))/(lambda2-lambda1)
    y = point1[1] - lambda1 * point1[0] + lambda1 * x
    return x, y

### よく理解していない
### とりあえず vol.2 p.203 の例題にしたがって書いてみる
def func_MEPC_theta3(theta1, theta2, point1, point2, point4, lambda_12, eps=10e-6):
    lambda_o = math.tan((theta1+theta2)/2.)
    delta = 1.0
    n=0
    while delta >= eps:
        x, y = func_cross_gas_dynamics(point2, point4, lambda_12,lambda_o)
        theta = theta2 + (y-point2[1])/(point1[1]-point2[1]) * (theta1-theta2)
        lambda_o_new = math.tan(theta)
        delta = abs((lambda_o-lambda_o_new)/lambda_o)
        lambda_o = lambda_o_new
        n += 1
        print(n)
    return x, y, theta


# theta1 = 18.191 / 360. *2.*np.pi
# theta2 = 16.422 / 360. *2.*np.pi
# x1 = 0.13146
# y1 = 0.040118
# x2 = 0.135683
# y2 = 0.037123
# lambda2 = -0.70921
# x4 = 0.141335
# y4 = 0.04061

theta1 = 0.22
theta2 = 0.17
x1 = -0.18
y1 = 1.
x2 = 0.
y2 = 0.
lambda12 = -5.67
x4 = 0.23
y4 = 0.56

x, y, theta = func_MEPC_theta3(theta1,theta2,(x1,y1),(x2,y2),(x4,y4),lambda12)
print(x, "///", y, "///", theta*360/2/np.pi)










