import matplotlib.pyplot as plt
import numpy as np
import os
from sympy import symbols, Function, diff, simplify, nsolve, sin, sinh
from calc_coeff_3d import a_coeff, integrate_coeff
from plot_func import plot_graph, plot_3d
os.chdir(r"C:\Users\kirik\PycharmProjects\Diser\3d sem freq")
# Символьные переменные для характеристик пластины
a, b, E, mu, D, k, k_rez, P, rho, h = symbols('a b E mu D k k_rez P rho h', real=True)
# Символьные переменные для характеристик ролика (он же лента)
M_rol, A_len, l_len, E_len, k_len, h_r = symbols('M_rol A_len l_len E_len k_len h_r', real=True)
# Вспомогательные переменные
x, t, xi1, xi2, xi1__, xi2__ = symbols('x t xi1 xi2 xi1__ xi2__', real=True)
# Безразмерные жесткости
kappa0, kappa1, kappa2, kappa3, kappa4 = symbols('kappa0 kappa1 kappa2 kappa3 kappa4', real=True)

# Запись физических параметров материалов из файла data.txt
with open('data.txt') as f:
    values = {}
    for line in f:
        line = line.rstrip('\r\n')
        number, variable = line.split('#')
        values[symbols(variable, real=True)] = float(eval(number))

b_a = float(values[b] / values[a])
alpha_, lambda_ = symbols('alpha_ lambda_', real=True)

# функции Крылова
K4 = 1 / 2 * (sinh(x) - sin(x))
K3 = diff(K4, x)
K2 = diff(K3, x)
K1 = diff(K2, x)

# Определение собственных значений для случая заделка / свободный край
# С1 = С2 = 0
Eq1 = simplify(K1 ** 2 - K4 * K2)

# x1_lim = [0, 6]  # Границы по x1
# plot_graph(Eq1, x1_lim, graph_name='Determine 1', x_axis_name=r'$\alpha$',
#            title_name='Определитель матрицы для "Балки 1"')  # строим график детерминанта
# plt.xlim(x1_lim[0], x1_lim[1])
# plt.ylim(-9, 6)

alfa1 = nsolve(Eq1, x, 1.2)
alfa2 = nsolve(Eq1, x, 4.5)

# Определение собственных значений для случая свободный край / свободный край
# C3 = C4 = 0
Eq2 = simplify(K3 ** 2 - K2 * K4)

# x2_lim = [0, 10]  # Границы по x2
# plot_graph(Eq2, x2_lim, graph_name='Determine 2', x_axis_name=r'$\lambda$',
#            title_name='Определитель матрицы для "Балки 2"')  # строим график детерминанта
# plt.xlim(x2_lim[0], x2_lim[1])
# plt.ylim(-210, 110)

lambda1 = nsolve(Eq2, x, 5) / b_a
lambda2 = nsolve(Eq2, x, 8) / b_a

# Уравнение функций форм в общем виде
fi = K3.subs(x, xi1 * alpha_) - K1.subs(x, alpha_) / K2.subs(x, alpha_) * K4.subs(x, xi1 * alpha_)
psi = K1.subs(x, xi2 * lambda_) - K3.subs(x, lambda_ * b_a) / K4.subs(x, lambda_ * b_a) * K2.subs(x, xi2 * lambda_)

# Выражения функций форм для полученных собственных значений
fi1 = fi.subs(alpha_, alfa1)
fi2 = fi.subs(alpha_, alfa2)
psi1 = psi.subs(lambda_, lambda1)
psi2 = psi.subs(lambda_, lambda2)

# # Построение графики функций форм для "балок"
# plot_graph(fi1, [0, 1], graph_name='fi1', variable=xi1, x_axis_name=r'$\xi_1$',
#            y_axis_name=r'$\varphi (\alpha_1 \xi_1)$', title_name='Первая форма для "Балки 1"')
# plot_graph(fi2, [0, 1], graph_name='fi2', variable=xi1, x_axis_name=r'$\xi_1$',
#            y_axis_name=r'$\varphi (\alpha_2 \xi_1)$', title_name='Вторая форма для "Балки 1"')
# plot_graph(psi1, [0, b_a], graph_name='psi1', variable=xi2, x_axis_name=r'$\xi_2$',
#            y_axis_name=r'$\psi (\lambda_1 \xi_2)$', title_name='Первая форма для "Балки 2"')
# plot_graph(psi2, [0, b_a], graph_name='psi2', variable=xi2, x_axis_name=r'$\xi_2$',
#            y_axis_name=r'$\psi (\lambda_2 \xi_2)$', title_name='Вторая форма для "Балки 2"')

# # Строим графики функций форм для пластины
# plot_3d(fi1*psi1, 1, b_a, graph_name='fi1_psi1')
# plot_3d(fi1*psi2, 1, b_a, graph_name='fi1_psi2')
# plot_3d(fi2*psi1, 1, b_a, graph_name='fi2_psi1')
# plot_3d(fi2*psi2, 1, b_a, graph_name='fi2_psi2')
# plt.show()

# Коэффициенты при обобщённых координатах в выражении потенциальной энергии при сжатии пружины:
u11_spring = (fi1.subs(xi1, 1) * psi1.subs(xi2, 0.5 * b_a)) ** 2 * kappa3
u22_spring = (fi2.subs(xi1, 1) * psi2.subs(xi2, 0.5 * b_a)) ** 2 * kappa3
u12_spring = fi1.subs(xi1, 1) * psi1.subs(xi2, 0.5 * b_a) * fi2.subs(xi1, 1) * psi2.subs(xi2, 0.5 * b_a) * kappa3

# Коэффициент при обобщённых координатах в выражении потенциальной энергии при деформации пластины:
A = a_coeff(fi1, psi1, fi2, psi2, b_a)
a_11 = A[0] + u11_spring
a_22 = A[1] + u22_spring
a_12 = A[2] + u12_spring
a_21 = a_12

# Создание матриц коэффициентов
# M * ddz + C * z = 0
m11 = integrate_coeff(fi1**2 * psi1**2, b_a)
m12 = integrate_coeff(fi1 * fi2 * psi1 * psi2, b_a)
m21 = m12
m22 = integrate_coeff(fi2**2 * psi2**2, b_a)

M = np.array([
    [1*xi1/xi1, 0*xi1, 0*xi1],
    [0*xi1, m11, m12],
    [0*xi1, m21, m22]
])

# Значение перемещения пластины в точке касания ролика
u1 = fi1.subs(xi1, xi1__)*psi1.subs(xi2, xi2__)
u2 = fi2.subs(xi1, xi1__)*psi2.subs(xi2, xi2__)

c11 = a_11 + 2 * kappa4 * u1 ** 2
c12 = a_12 + 2 * kappa4 * u1 * u2
c21 = a_21 + 2 * kappa4 * u1 * u2
c22 = a_22 + 2 * kappa4 * u2 ** 2

C = np.array([
    [kappa1 + kappa2, -kappa2 * u1, -kappa2 * u2],
    [-kappa4 * u1, c11, c12],
    [-kappa4 * u2, c21, c22]
])

np.save('M_Matrix', M)
np.save('C_Matrix', C)
