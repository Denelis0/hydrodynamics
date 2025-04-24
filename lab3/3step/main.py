import math
import numpy as np
import matplotlib.pyplot as plt
import tabulate


def boundary_condition(i, j, scheme):
    x, y = i * scheme.h, j * scheme.k
    n2, m4 = scheme.n // 2, scheme.m // 4
    if i == 0 and j >= m4:
        return y**2
    if i == scheme.n:
        return y**2 + 4
    if j == m4 and i <= n2:
        return x**2 + 0.0625
    if j == scheme.m:
        return x**2 + 1
    if i == n2 and j < m4:
        return y**2 + 1
    if j == 0 and i >= n2:
        return x**2
    if j < scheme.m // 4 and i < scheme.n // 2:
        return np.nan
    return 0

def rhs_function(i, j, scheme):
    return -4

def exact_solution(x, y):
    return x**2 + y**2

class Scheme:
    def __init__(self, w, h, n, m, rhs_function, boundary_condition, exact_solution=None):

        self.w, self.h_, self.n, self.m = w, h, n, m
        self.h, self.k = w / n, h / m
        self.inv_h2 = 1 / self.h**2
        self.inv_k2 = 1 / self.k**2
        self.a_star = -2 * (self.inv_h2 + self.inv_k2)
        self.rhs_function = rhs_function
        self.boundary_condition = boundary_condition
        self.exact_solution = exact_solution
        self.grid = np.zeros((n + 1, m + 1))
        self.x = np.linspace(0, w, n + 1)
        self.y = np.linspace(0, h, m + 1)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')

        if exact_solution:
            self.uv = exact_solution(self.X, self.Y)

        for i in range(n + 1):
            for j in range(m + 1):
                self.grid[i, j] = boundary_condition(i, j, self)

    def solve(self, eps, max_iter):
        def update_value(i, j):
            prev = self.grid[i, j]
            self.grid[i, j] = (
                                      - self.rhs_function(i, j, self)
                                      - self.inv_h2 * (self.grid[i - 1, j] + self.grid[i + 1, j])
                                      - self.inv_k2 * (self.grid[i, j - 1] + self.grid[i, j + 1])
                              ) / self.a_star
            return abs(prev - self.grid[i, j]), i, j

        accuracy = []
        delta, it = eps + 1, 0

        while it < max_iter and delta > eps:
            delta = 0
            max_i = max_j = -1

            # Нижняя правая часть
            for i in range(self.n // 2 + 1, self.n):
                for j in range(1, self.m // 4 + 1):
                    d, ii, jj = update_value(i, j)
                    if d > delta:
                        delta, max_i, max_j = d, ii, jj

            # Верхняя часть
            for i in range(1, self.n):
                for j in range(self.m // 4 + 1, self.m):
                    d, ii, jj = update_value(i, j)
                    if d > delta:
                        delta, max_i, max_j = d, ii, jj

            accuracy.append(delta)
            print(f"Итерация {it + 1}: E_max = {delta:.3e} в точке (i={max_i}, j={max_j})")
            it += 1

        return it, delta

    def plot(self):
        fig, ax = plt.subplots()
        Z = self.grid
        cs = ax.contourf(self.X, self.Y, Z, cmap='viridis')
        fig.colorbar(cs)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.xaxis.set_major_locator(plt.MultipleLocator(self.h))
        ax.yaxis.set_major_locator(plt.MultipleLocator(self.k))
        plt.grid(True)
        rect = plt.Rectangle((0, 0), 1, 0.25, linewidth=1, edgecolor='none', facecolor='white')
        ax.add_patch(rect)
        plt.title("2D решение")
        plt.show()

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(self.X, self.Y, Z, cmap='viridis', edgecolor='none')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('U(x, y)')
        ax.set_title("3D решение")
        plt.show()

if __name__ == '__main__':
    s = Scheme(2, 1, 20, 20, rhs_function, boundary_condition, exact_solution)
    iters, acc = s.solve(1e-14, 10000)
    print("-----------------------------------------------------------------------------")
    print("Этап 3 лабораторной работы №3: уравнение Пуассона Δu(x,y) = 4 с граничными условиями Дирихле в области с прямоугольными границами.")
    print("Граничные условия:\n"
          "u(0,y)=y^2 y принадлежит [1/4,1]\n"
          "u(2,y)=y^2 + 4 y принадлежит [0,1]\n"
          "u(1,y)=y^2 + 1 y принадлежит [0,1/4]\n"
          "u(x,1/4)=x^2 + 1/16 x принадлежит [0,1]\n"
          "u(x,0)=x^2 x принадлежит [1,2]\n"
          "u(x,1)=x^2 + 1 x принадлежит [0,2]\n")
    print("Критерии остановки:\n"
          "Число шагов max_iter=10000\n"
          "Точность на шаге eps=10^-14\n")
    print("Точное решение представлено в виде u(x,y)=x^2 + y^2\n"
          "Решение задачи находилось на основе метода Зенделя.")
    print("-----------------------------------------------------------------------------")
    print("Сравнение точного и численного решения в узлах сетки:")
    data = []
    max_diff = 0
    for i in range(s.n + 1):
        for j in range(s.m + 1):
            x, y = s.x[i], s.y[j]
            numerical = s.grid[i, j]
            exact = s.uv[i, j]
            diff = abs(numerical - exact)
            max_diff = max(max_diff, diff)
            data.append([f"{x:.2f}", f"{y:.2f}", f"{numerical:.6f}", f"{exact:.6f}", f"{diff:.2e}"])

    print(tabulate.tabulate(data, headers=["x", "y", "Численное", "Точное", "Погрешность"], tablefmt='simple_grid'))
    print("-----------------------------------------------------------------------------")
    print(f"Максимальное отличие точного и численного решений: {max_diff:.6e}\n"
          f"Количество итераций: {iters}/10000\n"
          f"Точность на последнем шаге: {acc}")
    s.plot()
