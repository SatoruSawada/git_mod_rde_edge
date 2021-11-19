import numpy as np


#------------------------------------------------------------------
#### 1. characteristic lines -1st
#------------------------------------------------------------------
### num_ch_up & num_ch_down が小さすぎても問題（num_ch_up & num_ch_down >= 7）
### i方向（横）にtheta-neu=const.確認
### j方向（縦）にtheta+neu=const.確認
num_ch_up = 20 # number of initial characteristic lines (upper side)
num_ch_down = 10 # number of initial characteristic lines (down side)
num_ch_add = 5
# init_theta_delta = 10e-11
inflow_distance = 0.
array_x_fm = np.empty(0)
array_y_fm = np.empty(0)
array_theta_fm = np.empty(0)
array_x_sl = np.empty(0)
array_y_sl = np.empty(0)
array_zero0 = np.zeros((int(num_ch_down),int(num_ch_up-1)))
array_zero1 = np.zeros((int(num_ch_up + num_ch_down - 1),int(num_ch_down + num_ch_add)))

### x for characteristics
array_x_up = np.ones((int(num_ch_up))) * 1
array_x = np.flipud(np.diag(array_x_up))
array_x = np.delete(array_x,-1,0)
array_x_down = np.zeros((int(num_ch_down)))
array_x_down = np.transpose(np.array([array_x_down]))
array_x_down = np.hstack((array_x_down,array_zero0))
array_x = np.vstack((array_x,array_x_down))
array_x = np.hstack((array_x, array_zero1))
del array_x_up
del array_x_down
array_x[0][int(num_ch_up-1)] = 10.
array_x[1][int(num_ch_up-1)] = 1.

num_ch_add = 5

def func_add_ch(array_target):
    target_sl = array_target[0][int(num_ch_up-1)]
    print(target_sl)
    array_target = np.delete(array_target,0,0)
    array_target_add0 = np.linspace(target_sl,array_target[0][int(num_ch_up-1)],int(num_ch_add+1))
    print(array_target_add0)
    array_target_add0 = np.array([np.delete(array_target_add0, -1, 0)])
    print(array_target_add0)
    array_target_add1 = np.zeros((int(num_ch_add), int(num_ch_up-1)))
    array_target_add1 = np.hstack((array_target_add1, np.transpose(array_target_add0)))
    array_target_add1 = np.hstack((array_target_add1, np.zeros((int(num_ch_add),int(num_ch_down+num_ch_add)))))
    array_target = np.vstack((array_target_add1,array_target))
    return array_target


array_x = func_add_ch(array_x)

# print(array_x)
np.savetxt('array_x.csv', array_x, delimiter=',')
