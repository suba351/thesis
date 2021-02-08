import numpy as np
from sympy import symbols, Matrix
from sympy import Symbol, re, sqrt
from sympy.solvers import solve
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import os
os.chdir(r"C:\Users\kirik\PycharmProjects\Diser\3d sem freq")
t = Symbol("T")
a, b, E, mu, D, k, k_rez, rho, h, p0_2, b_a = symbols('a b E mu D k k_rez rho h p0_2 b_a', real=True)
M_rol, A_len, l_len, E_len, h_r = symbols('M_rol A_len l_len E_len h_r', real=True)
kappa1, kappa2, kappa3, kappa4 = symbols('kappa1 kappa2 kappa3 kappa4', real=True)
x, x1, x2, xi1, xi2, xi1__, xi2__ = symbols('x x1 x2 xi1 xi2 xi1__ xi2__', real=True)
"""
    a - длина пластины
    b - ширина пластины
    E, mu, D - модуль Юнга, коэффициент Пуассона, цилиндрическая жесткость пластины
    k - жесткость пружин Ролика и Пластины (конусная бабка моделируется пружиной)
    k_rez - коэффициент резания = предел текучести материала пластины * ширину ленты
    rho - плотность материала пластины
    h - толщина пластины
    
    M_rol - масса ролика
    A_len - площадь сечения ленты
    l_len - длина ленты (расстояние от шарнира до ролика)
    E_len - модуль Юнга ленты
    h_r - толшина резания (толщина снимаемого слоя материала)
    
"""
# считываем характеристики материалов и прочее из файла data.txt
with open('data.txt') as f:
    values = {}
    for line in f:
        line = line.rstrip('\r\n')
        number, variable = line.split('#')
        values[symbols(variable, real=True)] = float(eval(number))

values[D] = values[E] * values[h]**3 / (12 * (1 - values[mu]**2))
kappa = {}
values[p0_2] = values[D] / (values[rho] * values[h] * values[a]**4)
kappa[kappa1] = (values[k] + values[E_len] * values[A_len] / values[l_len]) / (values[M_rol] * values[p0_2])
kappa[kappa2] = (values[k_rez]) / (values[M_rol] * values[p0_2])
kappa[kappa3] = (values[k] * values[a]**2) / (values[D])
kappa[kappa4] = (values[k_rez] * values[a]**2) / (values[D])

# загружаем матрицы масс и жесткости
M = Matrix(np.load("M_Matrix.npy", allow_pickle=True))
C = Matrix(np.load("C_Matrix.npy", allow_pickle=True))

# подставляем значения параметров в матрицы (кроме координат контакта)
for x in kappa:
    M = M.subs(x, kappa[x])
    C = C.subs(x, kappa[x])

M = M.subs(mu, values[mu])
C = C.subs(mu, values[mu])

# разбиваем пластину сеткой (точки контакта)
mesh1 = np.linspace(0.0, 1.0, 41)
mesh2 = np.linspace(0., float(values[b] / values[a]), 81)
os.chdir(r"C:\Users\kirik\PycharmProjects\Diser\3d sem freq\figures")
i = 0

for xi1 in mesh1:
    i += 1
    p1 = []
    p2 = []
    p3 = []
    for xi2 in mesh2:
        sol = solve((C.subs([(xi1__, xi1), (xi2__, xi2)]) - t*M).det(), t)
        p1.append(sqrt(re(sol[0])))
        p2.append(sqrt(re(sol[1])))
        p3.append(sqrt(re(sol[2])))
    fig, axs = plt.subplots(3, 1, constrained_layout=True, figsize=(6, 9))
    axs[0].plot(mesh2, p1, 'g', linewidth=2)
    axs[0].set_title('p1')
    axs[0].set_xlabel('xi_2')
    axs[0].set_ylabel('frequent')
    axs[0].grid(True)
    axs[0].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    fig.suptitle('xi1 = ' + str('{:.3f}'.format(xi1)), fontsize=16)

    axs[1].plot(mesh2, p2, 'r', linewidth=2)
    axs[1].set_title('p2')
    axs[1].set_xlabel('xi_2')
    axs[1].set_ylabel('frequent')
    axs[1].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    axs[1].grid(True)

    axs[2].plot(mesh2, p3, 'b', linewidth=2)
    axs[2].set_title('p2')
    axs[2].set_xlabel('xi_2')
    axs[2].set_ylabel('frequent')
    axs[2].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    axs[2].grid(True)
    plt.savefig(str(int(i)) + '. xi1 = ' + str('{:.3f}'.format(xi1)) + '.jpg')
    plt.close()
    

