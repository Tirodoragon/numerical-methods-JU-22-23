import numpy as np
from multipledispatch import dispatch


def perpendicular(x):
    if x[0] != 0:
        return np.array([-(x[1] + x[2]) / x[0], 1, 1])
    if x[1] != 0:
        return np.array([1, -(x[0] + x[2]) / x[1], 1])
    if x[2] != 0:
        return np.array([1, 1, -(x[0] + x[1]) / x[2]])


@dispatch(np.ndarray, np.ndarray)
def powmet(A, yk):
    yp = np.ones(len(yk))
    iter = 0
    while abs((np.dot(A, yk)[0] / yk[0]) - np.dot(A, yp)[0] / yp[0]) > 1e-8:
        iter += 1
        yp = yk
        zk = np.dot(A, yk)
        yk = zk / np.linalg.norm(zk)
    return np.dot(A, yk)[0] / yk[0], yk, iter


@dispatch(np.ndarray, np.ndarray, np.ndarray)
def powmet(A, yk, e1):
    yp = np.ones(len(yk))
    iter = 0
    while abs((np.dot(A, yk)[0] / yk[0]) - np.dot(A, yp)[0] / yp[0]) > 1e-8:
        iter += 1
        yp = yk
        zk = np.dot(A, yk)
        zk -= e1 * (e1.dot(zk))
        yk = zk / np.linalg.norm(zk)
    return np.dot(A, yk)[0] / yk[0], yk, iter


@dispatch(np.ndarray, np.ndarray, np.ndarray, np.ndarray)
def powmet(A, yk, e1, e2):
    yp = np.ones(len(yk))
    iter = 0
    while abs((np.dot(A, yk)[0] / yk[0]) - np.dot(A, yp)[0] / yp[0]) > 1e-8:
        iter += 1
        yp = yk
        zk = np.dot(A, yk)
        zk -= (e1 * (e1.dot(zk)) + e2 * (e2.dot(zk)))
        yk = zk / np.linalg.norm(zk)
    return np.dot(A, yk)[0] / yk[0], yk, iter


def powvec(A, y):
    eig = np.zeros(len(y))
    eig[0], e1, iter1 = powmet(A, y)
    y2 = perpendicular(e1)
    eig[1], e2, iter2 = powmet(A, y2, e1)
    y3 = perpendicular(e2)
    eig[2], e3, iter3 = powmet(A, y3, e1, e2)
    return eig, iter1 + iter2 + iter3


@dispatch(np.ndarray, np.ndarray)
def ray(A, yk):
    rp = 0
    rk = 1
    iter = 0
    while abs(rk - rp) > 1e-8:
        iter += 1
        rp = rk
        zk = np.dot(A, yk)
        yk = zk / np.linalg.norm(zk)
        rk = yk.T.dot(A).dot(yk) / yk.T.dot(yk)
    return rk, yk, iter


@dispatch(np.ndarray, np.ndarray, np.ndarray)
def ray(A, yk, e1):
    rp = 0
    rk = 1
    iter = 0
    while abs(rk - rp) > 1e-8:
        iter += 1
        rp = rk
        zk = np.dot(A, yk)
        zk -= e1 * (e1.dot(zk))
        yk = zk / np.linalg.norm(zk)
        rk = yk.T.dot(A).dot(yk) / yk.T.dot(yk)
    return rk, yk, iter


@dispatch(np.ndarray, np.ndarray, np.ndarray, np.ndarray)
def ray(A, yk, e1, e2):
    rp = 0
    rk = 1
    iter = 0
    while abs(rk - rp) > 1e-8:
        iter += 1
        rp = rk
        zk = np.dot(A, yk)
        zk -= (e1 * (e1.dot(zk)) + e2 * (e2.dot(zk)))
        yk = zk / np.linalg.norm(zk)
        rk = yk.T.dot(A).dot(yk) / yk.T.dot(yk)
    return rk, yk, iter


def rayvec(A, y):
    eig = np.zeros(len(y))
    eig[0], e1, iter1 = ray(A, y)
    y2 = perpendicular(e1)
    eig[1], e2, iter2 = ray(A, y2, e1)
    y3 = perpendicular(e2)
    eig[2], e3, iter3 = ray(A, y3, e1, e2)
    return eig, iter1 + iter2 + iter3


def maketridiag(matrix):
    # Householder
    v = np.array(matrix[0][1:3])
    I = np.identity(v.size)
    e = np.zeros(v.size)
    e[0] = 1.
    u = v - (e * np.linalg.norm(v))
    Au = np.zeros((v.size, v.size))
    Au[0] += u
    Au = np.dot(Au.T, Au)
    norm = np.linalg.norm(u)
    H = I - 2. / (norm * norm) * Au

    # Householder on matrix
    TT = np.identity(3)
    for i in range(1, 3, 1):
        for j in range(1, 3, 1):
            TT[i, j] = H[i - 1, j - 1]
    matrix = TT.dot(matrix).dot(TT.T)
    return matrix


def givens(i, A):
    I = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., -1.]])
    v = np.array(A[i])
    sq = np.sqrt(v[i] * v[i] + v[i + 1] * v[i + 1])
    c = v[i] / sq
    s = v[i + 1] / sq
    I[i, i] = I[i + 1, i + 1] = c
    I[i + 1, i] = -s
    I[i, i + 1] = s
    return I


def iterQR(A):
    Ap = A - 1
    iter = 0
    while not np.allclose(A, Ap):
        iter += 1
        Ap = A
        for i in range(A[0].size - 1):
            G = givens(i, A)
            A = G.dot(A).dot(G.T)
    return np.array([A[0, 0], A[1, 1], A[2, 2]]), iter


matrix = np.array([[1., 2., 3.], [2., 4., 5.], [3., 5., -1.]])
y = np.array([1., 0., 0.])
pv, pviter = powvec(matrix, y)
print('Metoda potęgowa: ', pv, ' w ', pviter, ' iteracjach')
rv, rviter = rayvec(matrix,y)
print('Metoda potęgowa z Rayleighem: ', rv, ' w ', rviter, ' iteracjach')
qr, qriter = iterQR(maketridiag(matrix))
print('Metoda iteracyjna QR: ', qr, ' w ', qriter, ' iteracjach')
print('Domyślna metoda numpy: ', np.linalg.eigvals(matrix))
