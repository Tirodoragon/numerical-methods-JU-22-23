import numpy as np

f = 1.0  # 0.5

x = 1  # 0.5
x4 = x * x * x * x  # 3.5

exp1 = np.exp(-1)  # 30.5
exp4 = exp1 * exp1 * exp1 * exp1  # 3.5

current_exp1 = exp1  # 0.5
current_exp4 = exp4  # 0.5

cos1 = np.cos(x4)  # 30.5
sin1 = np.sin(x4)  # 30.5

current_sin = sin1  # 0.5
current_cos = cos1  # 0.5

# 101.5

for i in range(1, 24, 1):  # (4.5 + 1.5 + 1.5 + 0.5 + 3.5 + 3.5) x23 = 345
    f += current_sin * current_sin * current_exp1 + current_cos * current_exp4  # 4.5

    current_exp1 *= exp1  # 1.5
    current_exp4 *= exp4  # 1.5

    # sin(a+b) = sin(a)cos(b) + sin(b)cos(a)
    sin_cpy = current_sin  # 0.5
    current_sin = current_sin * cos1 + sin1 * current_cos  # 3.5

    # cos(a+b) = cos(a)cos(b) = sin(a)sin(b)
    current_cos = current_cos * cos1 - sin_cpy * sin1  # 3.5

print(f'{f:.10f}')  # 0

# 446.5
