import numpy as np


def function_x(x):
    return np.exp(np.cos(x) * np.cos(x))


def double_points(points, indices, to_add):
    for i in range(1, len(points) + to_add, 2):
        points = np.insert(points, i, function_x((indices[i - 1] + indices[i]) / 2))
        indices = np.insert(indices, i, (indices[i - 1] + indices[i]) / 2)
    return points, indices


def integral_trap(points, indices, to_add=1):
    old = 0
    new = ((indices[1] - indices[0]) / 2) * (points[0] + points[1])
    points[0] /= 2
    points[1] /= 2
    iter = 0
    while abs(old - new) > 1e-6:
        old = new
        points, indices = double_points(points, indices, to_add)
        to_add *= 2
        new = (indices[1] - indices[0]) * np.sum(points)
        iter += 1
    return new, iter


def integral_simp(points, indices, to_add=2):
    old = 0
    new = ((indices[1] - indices[0]) / 6) * (points[0] + 4 * points[1] + points[2])
    s1 = slice(1, -1, 2)
    s2 = slice(2, -1, 2)
    iter = 0
    while abs(old - new) > 1e-6:
        old = new
        points, indices = double_points(points, indices, to_add)
        to_add *= 2
        new = ((indices[1] - indices[0]) / 3) * (points[0] + 4 * np.sum(points[s1]) + 2 * np.sum(points[s2]) + points[-1])
        iter += 1
    return new, iter


a = 0
b = np.pi
points = np.array([function_x(a), function_x(b)])
indices = np.array([a, b])
integ_trap, iter = integral_trap(points, indices)
print(f"Wynik uzyskany za pomocą metody trapezów: {integ_trap} w {iter} iteracjach")
points = np.array([function_x(a), function_x((a + b) / 2), function_x(b)])
indices = np.array([a, (a + b) / 2, b])
integ_simp, iter = integral_simp(points, indices)
print(f"Wynik uzyskany za pomocą metody Simpsona: {integ_simp} w {iter} iteracjach")
