# 47982
import numpy as np
import matplotlib.pyplot as plt

N = 1000  # 0.5


def thomas(d):  #  15979.5 * 2 = 31959
    n = N - 1  # 1.5
    beta = np.empty(n)  # 999.5
    gamma = np.empty(n)  # 999.5
    u = np.empty(n)  # 999.5
    beta[0] = 1 / 3  # 1.5
    gamma[0] = d[0] / 3  # 1.5
    for i in range(1, n - 1):  # 997 * 7.5 = 7477.5
        w = 4 - beta[i - 1]  # 2.5
        beta[i] = 1 / w  # 1.5
        gamma[i] = (d[i] - gamma[i - 1]) / w  # 3.5
    gamma[n - 1] = (d[n - 1] - gamma[n - 2]) / (3 - beta[n - 2])  # 7.5
    u[n - 1] = gamma[n - 1]  # 2.5
    for i in range(n - 1, 0, -1):  # 998 * 5.5 = 5489
        u[i - 1] = gamma[i - 1] - beta[i - 1] * u[i]  # 5.5
    return u


h = 2 / N  # 1.5
f = [(1 / (1 + 25 * x * x)) for x in np.arange(-1, h, h)]  # 2004

g = np.empty(N - 1)  # 999.5

if (N / 2 == int(N / 2)):  # 5
    r = int(N / 2) - 1  # 3.5
    g[r] = 2 * f[r] - 2 * f[r + 1]  # 4.5
else:
    r = int(N / 2)  # 2.5

for i in range(r):  # 500 * 12 = 6000
    g[i] = f[i] - 2 * f[i + 1] + f[i + 2]  # 5.5
    g[-1 - i] = f[i] - 2 * f[i + 1] + f[i + 2]  # 6.5

g *= (6 / (h * h))  # 999 * 3.5 = 3496.5

g_prim = np.zeros(N - 1)  # 999.5
g_prim[0] = 1.0  # 0.5
g_prim[-1] = 1.0  # 0.5
z = thomas(g)  # 1.5
q = thomas(g_prim)  # 1.5
factor = (z[0] + z[-1]) / (1 + q[0] + q[-1])  # 4.5
for i in range(N - 1):  # 999 * 2 = 1998
    z[i] -= factor * q[i] # 2

plt.plot([x for x in np.arange(-1 + h, 1, h)], z)  # 0.5 * 999 = 499.5
plt.title('Wykres zaleznosci (xk, uk)')  # 0
plt.xlabel('xk')  # 0
plt.ylabel('uk')  # 0
plt.show()  # 0
