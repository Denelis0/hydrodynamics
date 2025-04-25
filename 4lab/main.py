from math import sin, pi
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

def create_area(length, time_line, n, m, U0, boundary_condition):
    grid = []
    h = length / n # Шаг по пространству
    t = time_line / m # Шаг по времени
    for j in range(m):
        grid.append([0] * n) # Все зануляем
    for i in range(n):
        grid[0][i] = U0(i * h) # При времени t=0 загружаем начальные условия
    for j in range(m):
        grid[j][0] = boundary_condition(j * t) # При x=0 ставим граничные условия
    return grid, h, t

def rectangle(x): # Прямоугольный импульс
    return float(x >= 1 and x <= 2)

def triangle(x): # Треугольный импульс
    if x < 1 or x > 2:
        return 0
    elif x >= 1 and x < 1.5:
        return 2 * (x - 1)
    else:
        return 1 - 2 * (x - 1.5)

def impulse(x): # Гладкий импульс
    if x < 1 or x > 2:
        return 0
    else:
        return 0.5 * (1 + sin(2 * pi * (x - 1) - pi / 2))

def boundary_condition(t): # Граничное условие слева u(0,t)
    return 0

def asymmetric(grid, a, h, r):
    for t in range(1, len(grid)):
        for x in range(1, len(grid[t])):
            grid[t][x] = grid[t - 1][x] - a * r * ((grid[t - 1][x] - grid[t - 1][x - 1]) / h)
    return grid

def symmetrical(grid, a, h, r):
    for t in range(1, len(grid)):
        for x in range(1, len(grid[t]) - 1):
            grid[t][x] = grid[t - 1][x] - a * r * ((grid[t - 1][x + 1] - grid[t - 1][x - 1]) / (2 * h))
    return grid

def run_through(grid, a, h, r):
    for t in range(1, len(grid)):
        prev = grid[t - 1]
        curr = grid[t] # Массив в определенный t по всем х
        T = 1 / r
        mu1 = curr[0] # Точка при х=0
        mu2 = 0
        kappa_first = 0
        kappa_end = 0
        ai = [kappa_first]
        bi = [mu1]
        C = -T
        A = -a / (4 * h)
        B = a / (4 * h)

        for x in range(1, len(curr) - 1):
            phi = prev[x] * T - a / (4 * h) * (prev[x + 1] - prev[x - 1])
            alpha = B / (C - A * ai[x - 1])
            beta = (-phi + A * bi[x - 1]) / (C - A * ai[x - 1])
            ai.append(alpha)
            bi.append(beta)

        curr[-1] = (kappa_end * bi[-1] + mu2) / (1 - kappa_end * ai[-1])
        for x in range(len(curr) - 2, -1, -1):
            curr[x] = ai[x] * curr[x + 1] + bi[x]
    return grid

def generate_precise(U0, a, i, h, n, r):
    grid = []
    for x in range(n):
        grid.append(U0(x * h - i * r * a))
    return grid

a = 2 # Скорость
w = 15
h = 15 # Время графика
n = w * 10
m = h * 20

if __name__ == '__main__':
    type = impulse # Тип импульса
    area1, hx, r = create_area(w, h, n, m, type, boundary_condition)
    area2, hx, r = create_area(w, h, n, m, type, boundary_condition)
    area3, hx, r = create_area(w, h, n, m, type, boundary_condition)
    c = abs(a) * r / hx
    if c > 1:
        print('Не сходится: c = ', c)
        print(r, hx)
        exit(1)
    area_new_1 = asymmetric(area1, a, hx, r)
    area_new_2 = symmetrical(area2, a, hx, r)
    area_new_3 = run_through(area3, a, hx, r)

    fig, ax = plt.subplots()
    x = np.linspace(0, w, n)
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    print("-----------------------------------------------------------------------------")
    print(
        "Лабораторная работа №4: уравнение переноса ∂u/∂t+a∂=*∂u/∂x=0 с скоростью а=2.\n"
        "Начальные условия: Прямоугольный (rectangle), треугольный (triangle) и гладкий (impulse) импульс")
    print("Граничные условия: u(0,t)=0\n")
    print("Точное решение представлено в виде u(x,t)=(x-at)\n"
          "В коде использованы 3 метода: ассиметричная разностная схема, симметричная разностная схема и 3 метод (с использованием прогонки).")
    print("-----------------------------------------------------------------------------")
    def update(frame):
        plt.cla()
        #ax.plot(x, area_new_1[frame], '-', label='метод 1')
        #ax.plot(x,area_new_2[frame], '-', label = 'метод 2')
        #ax.plot(x, area_new_3[frame], '-', label='метод 3')
        ax.plot(x, generate_precise(type, a, frame, hx, n, r), '-', label='Точное решение')
        plt.legend()
    ani = animation.FuncAnimation(fig, update, frames=m, interval=30)
    plt.show()
