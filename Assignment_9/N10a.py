import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")


def lagrange(j, x, N):
    xj = np.full(N - 1, x[j])
    xn = np.concatenate((x[:j], x[j + 1:]))
    poly = np.prod(-xn)
    prod = np.prod(xj - xn)
    return poly / prod


def count_ak(fjk, L, j, a):
    a[j] = np.sum(L * fjk)


def newfunction(y, x, k, a):
    return (y - a[k]) / x


def function_test(opt):
    y = 0
    N = 2

    a = np.zeros(N)

    step = 10 / (N - 1)

    if opt == 1:
        x = np.array(np.arange(-5, 5 + (step / 2), step))
        y = 1 / (1 + x * x)
    elif opt == 2:
        x = np.array(np.arange(-5, 5 + (step / 2), step))
        y = np.exp(-x * x)
    elif opt == 3:
        x = 5 * np.cos((np.arange(0, N, 1) * np.pi) / (N - 1))
        y = 1 / (1 + x * x)
    else:
        x = 5 * np.cos((np.arange(0, N, 1) * np.pi) / (N - 1))
        y = np.exp(-x * x)

    L = np.zeros(N)
    for i in range(N):
        L[i] = lagrange(i, x, N)

    for i in range(N):
        count_ak(y, L, i, a)
        y = newfunction(y, x, i, a)

    plot_x = np.linspace(np.amin(x), np.amax(x), 200)
    plot_y = np.zeros(len(plot_x))
    for i in range(len(plot_x)):
        y = 0
        for j in range(N):
            y += plot_x[i] ** j * a[j]

        plot_y[i] = y

    if opt == 1 or opt == 3:
        cur = np.amax(abs((1 / (1 + plot_x * plot_x)) - plot_y))
    else:
        cur = np.amax(abs((np.exp(-plot_x * plot_x) - plot_y)))
    prev = cur + 1

    BestN = 0
    Bestcur = 10000

    while N < 40:
        N += 1
        prev = cur
        a = np.zeros(N)

        step = 10 / (N - 1)
        x = np.array(np.arange(-5, 5 + (step / 2), step))
        if opt == 1:
            x = np.array(np.arange(-5, 5 + (step / 2), step))
            y = 1 / (1 + x * x)
        elif opt == 2:
            x = np.array(np.arange(-5, 5 + (step / 2), step))
            y = np.exp(-x * x)
        elif opt == 3:
            x = 5 * np.cos((np.arange(0, N, 1) * np.pi) / (N - 1))
            y = 1 / (1 + x * x)
        else:
            x = 5 * np.cos((np.arange(0, N, 1) * np.pi) / (N - 1))
            y = np.exp(-x * x)

        L = np.zeros(N)
        for i in range(N):
            L[i] = lagrange(i, x, N)

        for i in range(N):
            count_ak(y, L, i, a)
            y = newfunction(y, x, i, a)

        plot_x = np.linspace(np.amin(x), np.amax(x), 200)
        plot_y = np.zeros(len(plot_x))

        for i in range(len(plot_x)):
            y = 0
            for j in range(N):
                y += plot_x[i] ** j * a[j]

            plot_y[i] = y

        if opt == 1 or opt == 3:
            cur = np.amax(abs((1 / (1 + plot_x * plot_x)) - plot_y))
        else:
            cur = np.amax(abs((np.exp(-plot_x * plot_x) - plot_y)))

        if cur < Bestcur:
            BestN = N
            Bestcur = cur

    print(f'Najlepszą dokładność uzyskano dla N={BestN}, wyniosła ona max(|f(x) - faprox(x)|)={Bestcur}')

    step = 10 / (BestN - 1)
    x = np.array(np.arange(-5, 5 + (step / 2), step))
    if opt == 1:
        x = np.array(np.arange(-5, 5 + (step / 2), step))
        y = 1 / (1 + x * x)
    elif opt == 2:
        x = np.array(np.arange(-5, 5 + (step / 2), step))
        y = np.exp(-x * x)
    elif opt == 3:
        x = 5 * np.cos((np.arange(0, BestN, 1) * np.pi) / (BestN - 1))
        y = 1 / (1 + x * x)
    else:
        x = 5 * np.cos((np.arange(0, BestN, 1) * np.pi) / (BestN - 1))
        y = np.exp(-x * x)
    plt.scatter(x, y, c='red')

    a = np.zeros(BestN)

    L = np.zeros(BestN)
    for i in range(BestN):
        L[i] = lagrange(i, x, BestN)

    for i in range(BestN):
        count_ak(y, L, i, a)
        y = newfunction(y, x, i, a)

    plot_x = np.linspace(np.amin(x), np.amax(x), 200)
    plot_y = np.zeros(len(plot_x))
    for i in range(len(plot_x)):
        y = 0
        for j in range(BestN):
            y += plot_x[i] ** j * a[j]

        plot_y[i] = y

    plt.plot(plot_x, plot_y)

    plt.show()


function_test(1)
function_test(2)
function_test(3)
function_test(4)
