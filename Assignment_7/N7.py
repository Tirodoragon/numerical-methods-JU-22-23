import numpy as np


def LinearCG():
    iter = 0
    xk = x0 = np.ones(1001)
    rk = np.empty(len(x0))
    rk[-1] = rk[0] = 0
    h_sq = 1 / 1e6
    for i in range(1, len(x0) - 1, 1):
        rk[i] = -xk[i - 1] + (h_sq + 2) * xk[i] - xk[i + 1]

    pk = -rk
    rk_norm = np.linalg.norm(rk)

    while rk_norm > 1e-10:
        apk = np.empty(len(pk))
        apk[0] = h_sq * pk[0]
        for i in range(1, len(pk) - 1, 1):
            apk[i] = -pk[i - 1] + (h_sq + 2) * pk[i] - pk[i + 1]
        apk[-1] = h_sq * pk[-1]
        rkrk = np.dot(rk, rk)

        alpha = rkrk / np.dot(pk, apk)
        xk = xk + alpha * pk
        rk = rk + alpha * apk
        beta = np.dot(rk, rk) / rkrk
        pk = -rk + beta * pk

        rk_norm = np.linalg.norm(rk)
        iter += 1

    print(f'Rozwiązanie: \t x = {xk} \t w {iter} iteracjach')


LinearCG()
