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


if __name__ == "__main__":
    X, Y = GetData()

    # Функционал = Ji(const,X,Y)
    a_f = argmaxF(func_a, 0, 1, 1e-3)
    print(a_f, func_a(a_f))

    t_f = argmaxF(func_t, -4, 0, 1e-3)
    print(t_f, func_t(t_f))
    #
    # Функционал = Ji(const,mode(X), mode(Y))
    a_f_mode = argmaxF(func_mode_a, 0, 1, 1e-3)
    print(a_f_mode, func_mode_a(a_f_mode))

    t_f_mode = argmaxF(func_mode_t, -4, 0, 1e-3)
    print(t_f_mode, func_mode_t(t_f_mode))

    # Функционал = Ji(const,med_K(X), med_K(Y))
    a_f_med_k = argmaxF(func_med_k_a, 0, 1, 1e-3)
    print(a_f_med_k, func_med_k_a(a_f_med_k))

    t_f_med_k = argmaxF(func_med_k_t, -4, 0, 1e-3)
    print(t_f_med_k, func_med_k_t(t_f_med_k))

    # Функционал = Ji(const,mode(X), mode(Y))
    a_f_med_p = argmaxF(func_med_p_a, 0, 1, 1e-3)
    print(a_f_med_p, func_med_p_a(a_f_med_p))

    t_f_med_p = argmaxF(func_med_p_t, -4, 0, 1e-3)
    print(t_f_med_p, func_med_p_t(t_f_med_p))
