import numpy as np
import matplotlib.pyplot as plt


def solve(A, B, C, D):
    for i in range(N - 3):
        w = A[i] / B[i]
        B[i + 1] -= w * C[i]
        D[i + 1] -= w * D[i]
    u = np.zeros(N - 2)
    u[-1] = D[-1] / B[-3]
    for i in range(N - 4, -1, -1):
        u[i] = (D[i] - C[i] * u[i + 1]) / B[i]
    return u


def spline(x_x, x, y, ksi, step):
    f_x = np.zeros_like(x_x)
    for i in range(x_x.size - 1):
        j = int((x_x[i] - x[0]) / step)
        Ax = (x[j + 1] - x_x[i]) / step
        Bx = (x_x[i] - x[j]) / step
        Cx = (1 / 6) * (Ax * Ax * Ax - Ax) * step * step
        Dx = (1 / 6) * (Bx * Bx * Bx - Bx) * step * step
        f_x[i] = Ax * y[j] + Bx * y[j + 1] + Cx * ksi[j] + Dx * ksi[j + 1]
    Ax = ((x[-1] - x_x[-1]) / step)
    Bx = ((x_x[-1] - x[-2]) / step)
    f_x[-1] = Ax * y[-2] + Bx * y[-1] + ((1 / 6) * (Ax * Ax * Ax - Ax) * step * step) * ksi[-2] + ((1 / 6) * (Bx * Bx * Bx - Bx) * step * step) * ksi[-1]
    return f_x


cmp = 1
N = 3
while 0.24684427850543922 <= cmp:
    # 1 / (1 + x * x)
    step = 10 / (N - 1)
    x = np.array(np.arange(-5, 5 + step, step))
    y = 1 / (1 + x * x)
    A = np.full(N - 1, 1.)
    B = np.full(N, 4.)
    C = np.full(N - 1, 1.)
    D = np.array([y[i] - 2 * y[i + 1] + y[i + 2] for i in range(0, N - 2)])
    ksi = solve(A, B, C, (6 / (step * step)) * D)
    ksi = np.insert(ksi, (0, N - 2), 0.)

    plot_x = np.linspace(x[0], x[-1], 200)
    plotf_x = spline(plot_x, x, y, ksi, step)

    cmp = np.amax(abs((1 / (1 + plot_x * plot_x)) - plotf_x))
    if 0.24684427850543922 >= cmp:
        print(f"Uzyskany wynik dla N={N}, uzyskany błąd: {cmp}")
        plt.scatter(x, y, c="red")
        plt.plot(plot_x, plotf_x)
        plt.show()
    N += 1

cmp = 1
N = 3
while 0.28372299145685076 <= cmp:
    # exp(-x * x)
    step = 10 / (N - 1)
    x = np.array(np.arange(-5, 5 + step, step))
    y = np.exp(-x * x)
    A = np.full(N - 1, 1.)
    B = np.full(N, 4.)
    C = np.full(N - 1, 1.)
    D = np.array([y[i] - 2 * y[i + 1] + y[i + 2] for i in range(0, N - 2)])
    ksi = solve(A, B, C, (6 / (step * step)) * D)
    ksi = np.insert(ksi, (0, N - 2), 0.)

    plot_x = np.linspace(x[0], x[-1], 200)
    plotf_x = spline(plot_x, x, y, ksi, step)

    cmp = np.amax(abs(np.exp(-plot_x * plot_x) - plotf_x))
    if 0.28372299145685076 >= cmp:
        print(f"Uzyskany wynik dla N={N}, uzyskany błąd: {cmp}")
        plt.scatter(x, y, c="red")
        plt.plot(plot_x, plotf_x)
        plt.show()
    N += 1

cmp = 1
N = 3
while 0.24684427850543922 <= cmp:
    # 1 / (1 + x * x)
    step = 10 / (N - 1)
    x = np.array(np.arange(-5, 5 + step, step))
    y = 1 / (1 + x * x)

    plot_x = np.linspace(x[0], x[-1], 200)
    sinc_m = np.tile(plot_x, (len(x), 1)) - np.tile(x[:, np.newaxis], (1, len(plot_x)))
    plotf_x = np.dot(y, np.sinc(sinc_m / step))

    cmp = np.amax(abs((1 / (1 + plot_x * plot_x)) - plotf_x))
    if 0.24684427850543922 >= cmp:
        print(f"Uzyskany wynik dla N={N}, uzyskany błąd: {cmp}")
        plt.scatter(x, y, c="red")
        plt.plot(plot_x, plotf_x)
        plt.show()
    N += 1

cmp = 1
N = 3
while 0.28372299145685076 <= cmp:
    # exp(-x * x)
    step = 10 / (N - 1)
    x = np.array(np.arange(-5, 5 + step, step))
    y = np.exp(-x * x)

    plot_x = np.linspace(x[0], x[-1], 200)
    sinc_m = np.tile(plot_x, (len(x), 1)) - np.tile(x[:, np.newaxis], (1, len(plot_x)))
    plotf_x = np.dot(y, np.sinc(sinc_m / step))

    cmp = np.amax(abs(np.exp(-plot_x * plot_x) - plotf_x))
    if 0.28372299145685076 >= cmp:
        print(f"Uzyskany wynik dla N={N}, uzyskany błąd: {cmp}")
        plt.scatter(x, y, c="red")
        plt.plot(plot_x, plotf_x)
        plt.show()
    N += 1
