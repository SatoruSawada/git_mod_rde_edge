import numpy as np

num_ch = 5
angle_fm = 10.
angle_sl = 30.
angle_bottom = 0.

# array0 = np.linspace(angle_fm,angle_sl,num_ch)
# print("======== array0")
# print(array0)


# # array1_sum = np.sum(array1)
# array1 = np.linspace(1,num_ch,num_ch)
# array1 = array1 * array1
# array1_delta = angle_sl / array1[-1]
# array1 = array1 * array1_delta
# print(array1)
# # print(array1_sum)
# print(array1_delta)
# # print(array1)



# # array1_sum = np.sum(array1)
# array2 = np.linspace(1,num_ch,num_ch)
# array2 = array2 * array2
# array2_delta = angle_fm / array2[-1]
# array2 = angle_fm - array2 * array2_delta

# print(array2)
# # print(array1_sum)
# print(array2_delta)
# print(array2)



# array1_sum = np.sum(array1)
array1 = np.linspace(angle_fm,angle_sl,num_ch)
array1 = array1 * array1
array1_delta = (angle_sl-angle_fm) / np.sum(array1)
array1 = array1 * array1_delta + angle_fm
print(array1)
# print(array1_sum)
print(array1_delta)
# print(array1)


