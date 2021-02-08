import sympy
a, b, D, mu, k, P, xi1, xi2 = sympy.symbols('a b D mu k P xi1 xi2', real=True)


def diff_coeff(func1, func2, deriv_ord1=2):
    """
    Функция считает вторую производную произведения двух функций func1 и func2 каждая из которых зависит от одной из
    переменных, либо вторая производная по xi1, либо вторая производная по xi2, либо первая смешанная производная
    по xi1 и xi2 одновременно

    :param func1: первая функция
    :param func2: вторая функция
    :param deriv_ord1: Порядок производной.
                    Если =2 -> берется вторая производная по переменной первой функции
                    Если =1 -> берется смешанная производная по переменной первой функции и переменной второй функции
                    Если =0 -> берется вторая производная по переменной второй функции
    :return:
    """
    symbol1 = list(func1.free_symbols)[0]
    symbol2 = list(func2.free_symbols)[0]
    if deriv_ord1 == 2:
        value = sympy.diff(func1, symbol1, deriv_ord1) * func2
    elif deriv_ord1 == 1:
        value = sympy.diff(func1, symbol1, deriv_ord1) * sympy.diff(func2, symbol2, 2 - deriv_ord1)
    return value


def integrate_coeff(u, b_a):
    """
    Считается двойной интеграл в пределах [0, 1] от функции u
    :param u: интегрируемая функция
    :return:func1
    """
    value = sympy.integrate(u, (xi1, 0, 1), (xi2, 0, b_a))
    return value


def a_coeff(fi1, psi1, fi2, psi2, b_a):
    """
       Вычисляет коэффициент, стоящий перед множителем f1 * f2
    """
    u1_11 = diff_coeff(fi1, psi1)
    u2_11 = diff_coeff(fi2, psi2)
    u1_22 = diff_coeff(psi1, fi1)
    u2_22 = diff_coeff(psi2, fi2)
    a_11 = integrate_coeff(u1_11**2, b_a) + 2*integrate_coeff(u1_11*u1_22, b_a) + integrate_coeff(u1_22**2, b_a) - \
        2 * (1 - mu) * (integrate_coeff(u1_11 * u1_22, b_a) - integrate_coeff(diff_coeff(fi1, psi1, 1) ** 2, b_a))

    a_22 = integrate_coeff(u2_11 ** 2, b_a) + 2 * integrate_coeff(u2_11 * u2_22, b_a) + integrate_coeff(u2_22 ** 2, b_a) - \
           2 * (1 - mu) * (integrate_coeff(u2_11 * u2_22, b_a) - integrate_coeff(diff_coeff(psi2, fi2, 1) ** 2, b_a))

    a_12 = integrate_coeff(u1_11 * u2_11, b_a) + integrate_coeff(u1_11 * u2_22 + u2_11 * u1_22, b_a) + integrate_coeff(u1_22 * u2_22, b_a) - \
           (1 - mu) * (integrate_coeff(u1_11 * u2_22 + u2_11 * u1_22, b_a) - 2 * integrate_coeff(diff_coeff(fi1, psi1, 1) * diff_coeff(fi2, psi2, 1), b_a))
    value = [a_11, a_22, a_12]
    return value
