import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, lambdify


def plot_graph(func, x_limit=None, variable=symbols('x', real=True),
               graph_name=None, x_axis_name=None, y_axis_name=None, title_name=None):
    """

    функция строит график функции

    """
    if x_limit is None:
        x_limit = [0, 1]
    x_left, x_right = x_limit[0], x_limit[1]
    x_space = np.linspace(x_left, x_right, 101)
    y_space = lambdify(variable, func, 'numpy')
    plt.figure(graph_name)
    plt.grid(True)
    plt.xlabel(x_axis_name, fontsize=12)
    plt.ylabel(y_axis_name, fontsize=12)
    plt.title(title_name, fontsize=12)
    plt.plot(x_space, y_space(x_space), 'g', linewidth=2)
    plt.plot([x_left, x_right], [0, 0], 'r', linewidth=1)
    return ()


def plot_3d(func, x1_lim, x2_lim, graph_name=None, x_axis_name=None, y_axis_name=None):
    xi1, xi2 = symbols('xi1 xi2', real=True)
    x1_space, x2_space = np.mgrid[0:x1_lim:101j, 0:x2_lim:101j]
    x3_space = lambdify((xi1, xi2), func, 'numpy')
    x3_space = x3_space(x1_space, x2_space)
    fig = plt.figure(graph_name)
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(x_axis_name, fontsize=12)
    plt.ylabel(y_axis_name, fontsize=12)
    ax.plot_wireframe(x1_space, x2_space, x3_space, rstride=10, cstride=10)
    return ()
