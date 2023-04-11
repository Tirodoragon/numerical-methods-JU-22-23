import numpy as np
import matplotlib.pyplot as plt


def plot(x, y1, y2, y3):
    plt.xscale('log')
    plt.yscale('log')
    plt.scatter(x, y1, s=16, c='purple', label='Dx1')
    plt.scatter(x, y2, s=8, c='orange', label='Dx2')
    plt.scatter(x, y3, s=1, c='green', label='Dx3')
    plt.legend(loc="lower left")
    plt.show()


def printer(m, number, dx, h0):
    print('Plot ', number, ': ', 'min = ', m, 'h = ', h0[dx.index(m)])


def dx1(x, h1):
    return (np.sin(x + h1) - np.sin(x)) / h1


def dx2(x, h2):
    return (np.sin(x + h2) - np.sin(x - h2)) / (2.0 * h2)


def dx3(x, h3):
    return (-np.sin(x + (2.0 * h3)) + 8.0 * np.sin(x + h3) - 8.0 * np.sin(x - h3) + np.sin(x - (2.0 * h3))) / (12.0 * h3)


a = 1e-16
b = 1e-16
h = []
i = 0
while a < 1:
    h.append(a)
    a += 0.001 * b
    i += 1
    if i == 9000:
        b *= 10
        i = 0

h = np.array(h)

DX11 = [abs(dx1(1, i) - np.cos(1)) for i in h]
DX21 = [abs(dx2(1, i) - np.cos(1)) for i in h]
DX31 = [abs(dx3(1, i) - np.cos(1)) for i in h]
DX12 = [abs(dx1(np.pi / 2, i) - np.cos(np.pi / 2)) for i in h]
DX22 = [abs(dx2(np.pi / 2, i) - np.cos(np.pi / 2)) for i in h]
DX32 = [abs(dx3(np.pi / 2, i) - np.cos(np.pi / 2)) for i in h]

min1 = np.amin(DX11)
min2 = np.amin(DX21)
min3 = np.amin(DX31)
min4 = np.amin(DX12)
min5 = np.amin(DX22)
min6 = np.amin(DX32)

DX1 = (np.array(DX11) + np.array(DX12)) / 2
DX2 = (np.array(DX21) + np.array(DX22)) / 2
DX3 = (np.array(DX31) + np.array(DX32)) / 2

plot(h, DX11, DX21, DX31)

plot(h, DX12, DX22, DX32)

printer(min1, 1, DX11, h)
printer(min2, 2, DX21, h)
printer(min3, 3, DX31, h)
printer(min4, 4, DX12, h)
printer(min5, 5, DX22, h)
printer(min6, 6, DX32, h)

DX1 = list(DX1)
DX2 = list(DX2)
DX3 = list(DX3)

print('Optimal h for Dx1 = ', h[DX1.index(np.amin(DX1))])
print('Optimal h for Dx2 = ', h[DX2.index(np.amin(DX2))])
print('Optimal h for Dx3 = ', h[DX3.index(np.amin(DX3))])
