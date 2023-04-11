import numpy as np
import matplotlib.pyplot as plt

N = 1001
h_sq = 1 / 1e6


def TDMA(a, b, c, d):
    n = len(d)
    w = np.zeros(n - 1, float)
    g = np.zeros(n, float)
    p = np.zeros(n, float)

    w[0] = c[0] / b[0]
    g[0] = d[0] / b[0]

    for i in range(1, n - 1):
        tmp = (b[i] - a[i - 1] * w[i - 1])
        w[i] = c[i] / tmp
        g[i] = (d[i] - a[i - 1] * g[i - 1]) / tmp
    g[-1] = (d[-1] - a[-1] * g[-1]) / (b[-1] - a[-1] * w[-1])
    p[-1] = g[-1]
    for i in range(n - 1, 0, -1):
        p[i - 1] = g[i - 1] - w[i - 1] * p[i]
    return p


L = np.full((N - 1), -1)
L[-1] = 0
D = np.full(N, h_sq + 2)
D[0] = D[-1] = h_sq
U = np.full(N - 1, -1)
U[0] = 0
b = np.zeros(N)
b[0] = b[-1] = h_sq
x = np.zeros(N)
x[0] = x[-1] = 1

plt.plot(np.arange(0, 1 + 1 / 1e3, 1 / 1e3), TDMA(L, D, U, b))
plt.show()
