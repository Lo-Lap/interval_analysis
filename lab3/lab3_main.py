import copy

import struct
import numpy as np
import intvalpy as ip
import matplotlib.pyplot as plt
from functools import cmp_to_key

from data_preparation import GetData

ip.precision.extendedPrecisionQ = False


def union_intervals(x, y):
    return ip.Interval(min(x.a, y.a), max(x.b, y.b))


def mode(X):
    print("Calculate mode")
    if X is None:
        return None

    # InterSec = X[0]
    # for el in X[1:]:
    #     InterSec = ip.intersection(InterSec, el)
    #
    # if not (np.isnan(InterSec.a) and np.isnan(InterSec.b)):
    #     return InterSec

    Y_g = []
    for el in X:
        Y_g.append(el.a)
        Y_g.append(el.b)

    Y_g.sort()

    Z = [ip.Interval(Y_g[i], Y_g[i + 1]) for i in range(len(Y_g) - 1)]

    mu = [sum(1 for x_i in X if z_i in x_i) for z_i in Z]

    max_mu = max(mu)
    K = [index for index, element in enumerate(mu) if element == max_mu]

    m = [Z[k] for k in K]
    mode_ = []

    current_interval = m[0]

    for next_interval in m[1:]:
        print(current_interval, next_interval)
        res_inter = ip.intersection(current_interval, next_interval)
        if not (np.isnan(res_inter.a) and np.isnan(res_inter.b)):
            current_interval = union_intervals(current_interval, next_interval)
        else:
            mode_.append(current_interval)
            current_interval = next_interval

    mode_.append(current_interval)
    print(mode_)
    return np.array(mode_)


def med_K(X_data):
    c_inf = [ip.inf(el) for el in X_data]
    c_sup = [ip.sup(el) for el in X_data]

    return ip.Interval(np.median(c_inf), np.median(c_sup))


def med_P(X_data):
    x = sorted(X_data, key=cmp_to_key(lambda x, y: (x.a + x.b) / 2 - (y.a + y.b) / 2))

    index_med = len(x) // 2

    if len(x) % 2 == 0:
        return (x[index_med - 1] + x[index_med]) / 2

    return x[index_med]


def coefficient_Jakkard(X_data, Y_data=None):
    if Y_data is None:
        x_inf = [ip.inf(x) for x in X_data]
        x_sup = [ip.sup(x) for x in X_data]
        return (min(x_sup) - max(x_inf)) / (max(x_sup) - min(x_inf))

    # if isinstance(X_data, ip.ClassicalArithmetic) and isinstance(Y_data, ip.ClassicalArithmetic):
    #     return (min(ip.sup(X_data), ip.sup(Y_data)) - max(ip.inf(X_data), ip.inf(Y_data))) / \
    #         (max(ip.sup(X_data), ip.sup(Y_data)) - min(ip.inf(X_data), ip.inf(Y_data)))
    #
    # # jakkard_v = []
    # # for x, y in zip(X_data, Y_data):
    # #     coeff = (min(ip.sup(x), ip.sup(y)) - max(ip.inf(x), ip.inf(y))) / (
    # #                 max(ip.sup(x), ip.sup(y)) - min(ip.inf(x), ip.inf(y)))
    # #     jakkard_v.append(coeff)
    #
    # x_inf = [ip.inf(x) for x in X_data]
    # x_sup = [ip.sup(x) for x in X_data]
    # y_inf = [ip.inf(x) for x in Y_data]
    # y_sup = [ip.sup(x) for x in Y_data]
    # coeff = (min(max(x_sup), max(y_sup)) - max(min(x_inf), min(y_inf))) / \
    #         (max(max(x_sup), max(y_sup)) - min(min(x_inf), min(y_inf)))
    #
    # return coeff


def argmaxF(f, a, b, eps):
    lmbd = a + (3 - 5 ** 0.5) * (b - a) / 2
    mu = b - (3 - 5 ** 0.5) * (b - a) / 2
    f_lambda = f(lmbd)
    f_mu = f(mu)

    while 1:
        if f_lambda <= f_mu:
            a = lmbd
            if eps > b - a:
                break
            lmbd = mu
            f_lambda = f_mu
            mu = b - (3 - 5 ** 0.5) * (b - a) / 2
            f_mu = f(mu)
        else:
            b = mu
            if eps > b - a:
                break
            mu = lmbd
            f_mu = f_lambda
            lmbd = a + (3 - 5 ** 0.5) * (b - a) / 2
            f_lambda = f(lmbd)

        # print(a)
        # print(b)

    return (a + b) / 2


def func_a(a):
    new_X = X + a
    XY = np.concatenate((new_X, Y))
    return coefficient_Jakkard(XY)


def func_t(t):
    new_X = X*t
    XY = np.concatenate((new_X, Y))
    return coefficient_Jakkard(XY)


def func_mode_a(a):
    new_X = X + a
    XY = np.concatenate((mode(new_X), Y_mode))
    # XY = np.concatenate((X_mode+a, Y_mode))
    return coefficient_Jakkard(XY)


def func_mode_t(t):
    new_X = X*t
    XY = np.concatenate((mode(new_X), Y_mode))
    # XY = np.concatenate((X_mode*t, Y_mode))
    return coefficient_Jakkard(XY)


def func_med_p_a(a):
    new_X = X + a
    XY = np.array([med_P(new_X), med_P(Y)])
    return coefficient_Jakkard(XY)


def func_med_p_t(t):
    new_X = X * t
    XY = np.array([med_P(new_X), med_P(Y)])
    return coefficient_Jakkard(XY)


def func_med_k_a(a):
    new_X = X + a
    XY = np.array([med_K(new_X), med_K(Y)])
    return coefficient_Jakkard(XY)


def func_med_k_t(t):
    new_X = X * t
    XY = np.array([med_K(new_X), med_K(Y)])
    return coefficient_Jakkard(XY)


def draw_func(f, a, b, parametr: str, func=""):
    if parametr == "a":
        # X_linsp = np.linspace(a, b, 175)
        X_linsp = np.linspace(a, b, 50)
    else:
        # X_linsp = np.linspace(a, b, 300)
        X_linsp = np.linspace(a, b, 175)
    y = [f(x) for x in X_linsp]
    y_max = max(y)
    y_min = min(y)
    ind_max = y.index(y_max)
    ind_min = y.index(y_min)
    x_max = X_linsp[ind_max]
    print(x_max, y_max)

    from scipy.interpolate import interp1d
    f1 = interp1d(np.array(y[0:ind_max + 1]), np.array(X_linsp[0:ind_max + 1]), kind='linear')
    f2 = interp1d(np.array(y[ind_max:]), np.array(X_linsp[ind_max:]), kind='linear')

    if y_max < 0:
        int_line = (3 * y[ind_max] + y[ind_min]) / 4
        ext_line = (y[ind_max] + 3 * y[ind_min]) / 4
    else:
        ext_line = 0.5
        int_line = 0.8

    # print(f"internal: {int_line}, x=[{f1(int_line)},{f2(int_line)}]")
    # print(f"external: {ext_line}, x=[{f1(ext_line)},{f2(ext_line)}]")

    plt.figure(figsize=(12, 9))
    plt.plot(X_linsp, y, color='b')
    plt.xlabel(f"{parametr}, {parametr}_max={round(x_max, 5)}")
    plt.ylabel(f"Ji({parametr}, {func}(X), {func}(Y))")
    plt.axvline(x=x_max, linestyle='--', color='black')

    # plt.hlines(int_line, f1(int_line), f2(int_line), linestyles='dashed', label='internal', colors='green')
    # plt.hlines(ext_line, f1(ext_line), f2(ext_line), linestyles='dashed', label='external', colors='red')
    # plt.legend()
    # plt.text(x_max + 0.15, (y_max + y_min) / 2, f"{parametr} = {round(x_max, 4)}", color='red',
    #          ha='center')  # Positioned near the line
    plt.title("Jaccard Index")
    plt.savefig(f"Jaccadrd-{parametr}-{func}")
    plt.show()


def draw_func_all(i, f, a, b, parametr: str, func=""):
    colors = ["#EF476F", "#F78C6B", "#FFD166", "#83D483", "#06D6A0", "#0CB0A9", "#118AB2", "#073B4C"]
    X_linsp = np.linspace(a, b, 100)
    y = np.array([f(x) for x in X_linsp])
    plt.plot(X_linsp, y, color=colors[i], label=f"Ji({parametr}, {func}(X), {func}(Y))", alpha=0.7)

    plt.xlabel(f"{parametr}")
    plt.ylabel(f"Ji({parametr}, {func}(X), {func}(Y))")
    plt.title("Jaccard Index")
    plt.show()
    plt.savefig(f"Jaccadrd-{parametr}-{func}")


if __name__ == "__main__":
    X, Y = GetData()
    # Функционал = Ji(const, X, Y)
    draw_func(func_a, 0.32, 0.37, "a")
    # a_f = argmaxF(func_a, 0, 1, 1e-3)
    # print(a_f, func_a(a_f))

    draw_func(func_t, -1.06, -0.98, "t")
    # t_f = argmaxF(func_t, -4, 0, 1e-3)
    # print(t_f, func_t(t_f))

    # # Функционал = Ji(const,mode(X), mode(Y))

    # X_mode = []
    # Y_mode = []

    # X_mode = mode(X)
    Y_mode = mode(Y)

    draw_func(func_mode_a, 0.34692, 0.34693, "a", "mode")
    # a_f_mode = argmaxF(func_mode_a, 0.3425, 0.3445, 1e-5)
    # print(a_f_mode, func_mode_a(a_f_mode))
    draw_func(func_mode_t, -1.0398, -1.0394, "t", "mode")
    # t_f_mode = argmaxF(func_mode_t, -1.017, -1.011, 1e-5)
    # print(t_f_mode, func_mode_t(t_f_mode))

    # Функционал = Ji(const,med_K(X), med_K(Y))
    draw_func(func_med_k_a, 0.3425, 0.3445, "a", "med_K")
    # # a_f_med_k = argmaxF(func_med_k_a, 0, 1, 1e-3)
    # # print(a_f_med_k, func_med_k_a(a_f_med_k))
    draw_func(func_med_k_t, -1.017, -1.011, "t", "med_K")
    # # t_f_med_k = argmaxF(func_med_k_t, -4, 0, 1e-3)
    # # print(t_f_med_k, func_med_k_t(t_f_med_k))

    # Функционал = Ji(const,med_р(X), med_р(Y))
    draw_func(func_med_p_a, 0.3425, 0.3445, "a", "med_p")
    # # a_f_med_p = argmaxF(func_med_p_a, 0, 1, 1e-3)
    # # print(a_f_med_p, func_med_p_a(a_f_med_p))
    draw_func(func_med_p_t, -1.017, -1.011, "t", "med_p")
    # # t_f_med_p = argmaxF(func_med_p_t, -4, 0, 1e-3)
    # # print(t_f_med_p, func_med_p_t(t_f_med_p))

    # funcs = [func_a, func_t, func_a_med_k, func_t_med_k, func_a_med_p, func_t_med_p]
    # funcs_str = ["", "", "med_k", "med_k", "med_p", "med_p"]
    # bounds = [[0, 1], [-4, 0], [0, 1], [-4, 0], [0, 1], [-4, 0]]
    # params = ["a", "t", "a", "t", "a", "t"]
    #
    # for i in range(1, len(funcs)+1, 2):
    #     draw_func_all(i, funcs[i], bounds[i][0], bounds[i][1], params[i], funcs_str[i])
    # plt.xlabel(f"const")
    # plt.ylabel(f"Ji(const, func(X), func(Y))")
    # plt.title("Jaccard Index")
    # plt.legend()
    # plt.savefig(f"Jaccadrd-all-in-one-T")
    # plt.show()
    #
    # for i in range(0, len(funcs), 2):
    #     draw_func_all(i, funcs[i], bounds[i][0], bounds[i][1], params[i], funcs_str[i])
    # plt.xlabel(f"const")
    # plt.ylabel(f"Ji(const, func(X), func(Y))")
    # plt.title("Jaccard Index")
    # plt.legend()
    # plt.savefig(f"Jaccadrd-all-in-one-A")
    # plt.show()
