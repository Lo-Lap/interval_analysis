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

    Y = []
    for el in X:
        Y.append(el.a)
        Y.append(el.b)

    Y.sort()

    Z = [ip.Interval(Y[i], Y[i + 1]) for i in range(len(Y) - 1)]

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

    return mode_


def med_K(X):
    c_inf = [ip.inf(el) for el in X]
    c_sup = [ip.sup(el) for el in X]

    return ip.Interval(np.median(c_inf), np.median(c_sup))


def med_P(X):
    x = sorted(X, key=cmp_to_key(lambda x, y: (x.a + x.b) / 2 - (y.a + y.b) / 2))

    index_med = len(x) // 2

    if len(x) % 2 == 0:
        return (x[index_med - 1] + x[index_med]) / 2

    return x[index_med]


def coefficient_Jakkard(X_data, Y_data=None):
    if Y_data is None:
        x_inf = [ip.inf(x) for x in X_data]
        x_sup = [ip.sup(x) for x in X_data]
        return (min(x_sup) - max(x_inf)) / (max(x_sup) - min(x_inf))

    if isinstance(X_data, ip.ClassicalArithmetic) and isinstance(Y_data, ip.ClassicalArithmetic):
        return (min(ip.sup(X_data), ip.sup(Y_data)) - max(ip.inf(X_data), ip.inf(Y_data))) / \
            (max(ip.sup(X_data), ip.sup(Y_data)) - min(ip.inf(X_data), ip.inf(Y_data)))

    jakkard_v = []
    for x, y in zip(X_data, Y_data):
        coeff = (min(ip.sup(x), ip.sup(y)) - max(ip.inf(x), ip.inf(y))) / (max(ip.sup(x), ip.sup(y)) - min(ip.inf(x), ip.inf(y)))
        jakkard_v.append(coeff)

    return jakkard_v


def argmaxF(f, a, b, eps):
    lmbd = a + (3 - 5 ** 0.5) * (b - a)/2
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
            lmbd = a + (3 - 5 ** 0.5) * (b - a)/2
            f_lambda = f(lmbd)

        print(a)
        print(b)

    return (a+b) / 2


def func_a(a):
    return np.mean(coefficient_Jakkard(X + a, Y))


def func_t(t):
    return np.mean(coefficient_Jakkard(X * t, Y))


def func_mode_a(a):
    return np.mean(coefficient_Jakkard(mode(X + a), mode(Y)))


def func_mode_t(t):
    return np.mean(coefficient_Jakkard(mode(X * t), mode(Y)))


def func_med_p_a(a):
    return np.mean(coefficient_Jakkard(med_P(X + a), med_P(Y)))


def func_med_p_t(t):
    return np.mean(coefficient_Jakkard(med_P(X * t), med_P(Y)))


def func_med_k_a(a):
    return np.mean(coefficient_Jakkard(med_K(X + a), med_K(Y)))


def func_med_k_t(t):
    return np.mean(coefficient_Jakkard(med_K(X * t), med_K(Y)))


def draw_func(f, a, b, parametr: str, func=""):
    X_linsp = np.linspace(a, b, 100)
    y = [f(x) for x in X_linsp]
    y_max = max(y)
    y_min = min(y)
    ind_max = y.index(y_max)
    x_max = X_linsp[ind_max]

    plt.plot(X_linsp, y, color='b')
    plt.xlabel(f"{parametr}")
    plt.ylabel(f"Ji({parametr}, {func}(X), {func}(Y))")
    plt.axvline(x=x_max, linestyle='--', color='r')
    plt.text(x_max+0.15, (y_max + y_min)/2, f"x = {round(x_max, 4)}", color='red',
             ha='center')  # Positioned near the line

    plt.title("Jaccard Index")
    plt.savefig(f"Jaccadrd-{parametr}-{func}")
    plt.show()


def draw_func_all(i, f, a, b, parametr: str, func=""):
    colors = ["#EF476F", "#F78C6B", "#FFD166", "#83D483", "#06D6A0", "#0CB0A9", "#118AB2", "#073B4C"]
    X_linsp = np.linspace(a, b, 100)
    y = np.array([f(x) for x in X_linsp])
    plt.plot(X_linsp, y, color=colors[i], label=f"Ji({parametr}, {func}(X), {func}(Y))", alpha=0.7)

    # plt.xlabel(f"{parametr}")
    # plt.ylabel(f"Ji({parametr}, {func}(X), {func}(Y))")
    # plt.title("Jaccard Index")
    # plt.show()
    # plt.savefig(f"Jaccadrd-{parametr}-{func}")


if __name__ == "__main__":
    X, Y = GetData()

    # # Функционал = Ji(const, X, Y)
    draw_func(func_a, 0, 1, "a")
    # a_f = argmaxF(func_a, 0, 1, 1e-3)
    # print(a_f, func_a(a_f))
    draw_func(func_t, -4, 0, "t")
    # t_f = argmaxF(func_t, -4, 0, 1e-3)
    # print(t_f, func_t(t_f))

    # # Функционал = Ji(const,mode(X), mode(Y))
    # draw_func(func_mode_a, 0, 1)
    # # a_f_mode = argmaxF(func_mode_a, 0, 1, 1e-3)
    # # print(a_f_mode, func_mode_a(a_f_mode))
    # draw_func(func_mode_t, -4, 0)
    # t_f_mode = argmaxF(func_mode_t, -4, 0, 1e-3)
    # print(t_f_mode, func_mode_t(t_f_mode))

    # # Функционал = Ji(const,med_K(X), med_K(Y))
    draw_func(func_med_k_a, 0, 1, "a", "med_K")
    # a_f_med_k = argmaxF(func_med_k_a, 0, 1, 1e-3)
    # print(a_f_med_k, func_med_k_a(a_f_med_k))
    draw_func(func_med_k_t, -4, 0, "t", "med_K")
    # t_f_med_k = argmaxF(func_med_k_t, -4, 0, 1e-3)
    # print(t_f_med_k, func_med_k_t(t_f_med_k))

    # # Функционал = Ji(const,med_р(X), med_р(Y))
    draw_func(func_med_p_a, 0, 1, "a", "med_p")
    # a_f_med_p = argmaxF(func_med_p_a, 0, 1, 1e-3)
    # print(a_f_med_p, func_med_p_a(a_f_med_p))
    draw_func(func_med_p_t, -4, 0, "t", "med_p")
    # t_f_med_p = argmaxF(func_med_p_t, -4, 0, 1e-3)
    # print(t_f_med_p, func_med_p_t(t_f_med_p))

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
