from math import sin, pi
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import time
from tabulate import tabulate

has_printed_steps = False

def create_area(length, time_line, n, m, U0, boundary_condition):
    global has_printed_steps
    grid = []
    h = length / n # Шаг по пространству
    t = time_line / m # Шаг по времени
    if not has_printed_steps:
        print(f"Шаг по пространству (h): {h}")
        print(f"Шаг по времени (t): {t}")
        has_printed_steps = True
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
h = 10 # Время графика
n = w * 10
m = h * 20

def plot():
    # 2D график
    fig, ax = plt.subplots()
    x = np.linspace(0, w, n)
    ax.plot(x, area_new_1[selected_frame], label='Метод 1')
    ax.plot(x, area_new_2[selected_frame], label='Метод 2')
    ax.plot(x, area_new_3[selected_frame], label='Метод 3')
    ax.plot(x, generate_precise(type, a, selected_frame, hx, n, r), label='Точное решение')
    ax.set_xlabel('x')
    ax.set_ylabel('u(x, t)')
    ax.set_title(f'Решения на временном слое t = {selected_frame}')
    ax.legend()
    plt.grid(True)
    plt.show()




if __name__ == '__main__':
    selected_frame = 10  # Номер временного слоя
    type = impulse  # Тип импульса
    area1, hx, r = create_area(w, h, n, m, type, boundary_condition)
    area2, hx, r = create_area(w, h, n, m, type, boundary_condition)
    area3, hx, r = create_area(w, h, n, m, type, boundary_condition)
    c = abs(a) * r / hx
    if c > 1:
        print('Не сходится: c = ', c)
        print(r, hx)
        exit(1)

    # Засекаем время для каждого метода
    start = time.time()
    area_new_1 = asymmetric(area1, a, hx, r)
    time_asymmetric = time.time() - start

    start = time.time()
    area_new_2 = symmetrical(area2, a, hx, r)
    time_symmetrical = time.time() - start

    start = time.time()
    area_new_3 = run_through(area3, a, hx, r)
    time_run_through = time.time() - start

    # Печать общей информации
    print("-----------------------------------------------------------------------------")
    print(
        "Лабораторная работа №4: уравнение переноса ∂u/∂t + a*∂u/∂x = 0 с a=2.\n"
        "Начальные условия: Прямоугольный (rectangle), треугольный (triangle) или гладкий (impulse) импульс."
    )
    print("Граничные условия: u(0,t) = 0\n")
    print("Точное решение: u(x, t) = u0(x - a*t)")
    print("Методы: асимметричная схема, симметричная схема, метод прогонки.")
    print("-----------------------------------------------------------------------------")
    # Инициализация переменных для максимальных погрешностей
    diff_max1 = diff_max2 = diff_max3 = 0

    # Вывод таблицы в терминал
    print("Таблица результатов:")
    data = []
    for i in range(n):
        x_val = i * hx
        num1 = area_new_1[selected_frame][i]
        num2 = area_new_2[selected_frame][i]
        num3 = area_new_3[selected_frame][i]
        exact = generate_precise(type, a, selected_frame, hx, n, r)[i]
        diff1 = abs(num1 - exact)
        diff2 = abs(num2 - exact)
        diff3 = abs(num3 - exact)

        # Обновление максимальных погрешностей
        diff_max1 = max(diff_max1, diff1)
        diff_max2 = max(diff_max2, diff2)
        diff_max3 = max(diff_max3, diff3)

        data.append([
            f"{x_val:.2f}",
            f"{num1:.6f}", f"{diff1:.2e}",
            f"{num2:.6f}", f"{diff2:.2e}",
            f"{num3:.6f}", f"{diff3:.2e}",
            f"{exact:.6f}"
        ])

    headers = ["x", "Метод 1", "Погрешность 1", "Метод 2", "Погрешность 2", "Метод 3", "Погрешность 3", "Точное"]
    print(tabulate(data, headers=headers, tablefmt='grid'))

    # Вывод максимальных погрешностей
    print(f"\nМаксимальная погрешность на первом методе: {diff_max1:.2e}")
    print(f"Максимальная погрешность на втором методе: {diff_max2:.2e}")
    print(f"Максимальная погрешность на третьем методе: {diff_max3:.2e}")
    plot()
