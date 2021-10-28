import numpy as np


judge = 0
judge_new = 1

for i in range(10):
    a = 0
    while judge != judge_new or a == 0:
        if i <= 3:
            print(i, 'state1')
            judge = judge_new
            judge = 1
            a += 1
        elif 3 < i <= 6:
            print(i, 'state2')
            judge = judge_new
            judge_new = 2
            a += 1
        else:
            print(i, 'state3')
            judge = judge_new
            judge_new = 3
            a += 1




