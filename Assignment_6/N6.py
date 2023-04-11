import numpy as np
import copy

N = 1001
h_sq = 1 / 1e6


def rel(x, gamma, b, L, D, U, xp):
    iterations = 0
    N = len(x) - 1
    while xp is None:
        xp = copy.deepcopy(x)
        x[0] = x[0] + gamma * (b[0] - D[0] * x[0] - U[0] * x[1])
        for i in range(1, N, 1):
            x[i] = x[i] + gamma * (b[i] - L[i - 1] * x[i - 1] - D[i] * x[i] - U[i] * x[i + 1])
        x[N] = x[N] + gamma * (b[N] - L[N - 1] * x[N - 1] - D[N] * x[N])
        iterations += 1
        for i in range(len(x)):
            if abs(x[i] - xp[i]) > 1e-10:
                xp = None
                break
    return iterations, x


def jac(x, b, L, D, U, xp):
    iterations = 0
    N = len(x) - 1
    while xp is None:
        xp = copy.deepcopy(x)
        x[0] = (b[0] - U[0] * xp[1]) / D[0]
        for i in range(1, N, 1):
            x[i] = (b[i] - L[i - 1] * xp[i - 1] - U[i] * xp[i + 1]) / D[i]
        x[N] = (b[N] - L[N - 1] * xp[N - 1]) / D[N]
        iterations += 1
        for i in range(len(x)):
            if abs(x[i] - xp[i]) > 1e-10:
                xp = None
                break
    return iterations, x


def gau(x, b, L, D, U, xp):
    iterations = 0
    N = len(x) - 1
    while xp is None:
        xp = copy.deepcopy(x)
        x[0] = (b[0] - U[0] * x[1]) / D[0]
        for i in range(1, N, 1):
            x[i] = (b[i] - L[i - 1] * x[i - 1] - U[i] * x[i + 1]) / D[i]
        x[N] = (b[N] - L[N - 1] * x[N - 1]) / D[N]
        iterations += 1
        for i in range(len(x)):
            if abs(x[i] - xp[i]) > 1e-10:
                xp = None
                break
    return iterations, x


def sor(x, omega, b, L, D, U, xp):
    iterations = 0
    N = len(x) - 1
    while xp is None:
        xp = copy.deepcopy(x)
        x[0] = (1 - omega) * x[0] + omega * (b[0] - U[0] * x[1]) / D[0]
        for i in range(1, N, 1):
            x[i] = (1 - omega) * x[i] + omega * (b[i] - L[i - 1] * x[i - 1] - U[i] * x[i + 1]) / D[i]
        x[N] = (1 - omega) * x[N] + omega * (b[N] - L[N - 1] * x[N - 1]) / D[N]
        iterations += 1
        for i in range(len(x)):
            if abs(x[i] - xp[i]) > 1e-10:
                xp = None
                break
    return iterations, x


# prostsza macierz - 80 punktów

iter1 = rel(np.array([0., 0., 0.]), 0.25878, np.array([1., 0., 1.]), np.array([-1., -1.]), np.array([4., 4., 4.]),
            np.array([-1., -1.]), None)
print(iter1)
iter2 = jac(np.array([0., 0., 0.]), np.array([1., 0., 1.]), np.array([-1., -1.]), np.array([4., 4., 4.]),
            np.array([-1., -1.]), None)
print(iter2)
iter3 = gau(np.array([0., 0., 0.]), np.array([1., 0., 1.]), np.array([-1., -1.]), np.array([4., 4., 4.]),
            np.array([-1., -1.]), None)
print(iter3)
iter4 = sor(np.array([0., 0., 0.]), 1.03511, np.array([1., 0., 1.]), np.array([-1., -1.]), np.array([4., 4., 4.]),
            np.array([-1., -1.]), None)
print(iter4)

# macierz N5 - 100 punktów

L = np.full((N - 1), -1)
L[-1] = 0
D = np.full(N, h_sq + 2)
D[0] = D[-1] = h_sq
U = np.full(N - 1, -1)
U[0] = 0
b = np.zeros(N)
b[0] = b[-1] = h_sq
x = np.ones(N)

iter5 = rel(copy.deepcopy(x), 0.978, b, L, D, U, None)
print(iter5)
iter6 = jac(copy.deepcopy(x), b, L, D, U, None)
print(iter6)
iter7 = gau(copy.deepcopy(x), b, L, D, U, None)
print(iter7)
iter8 = sor(copy.deepcopy(x), 1.941, b, L, D, U, None)
print(iter8)
