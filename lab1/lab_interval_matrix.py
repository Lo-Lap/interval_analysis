import intvalpy as ip
from sympy import *


def determinant_2(A):
    return ip.Interval(A_eps[0][0])*ip.Interval(A_eps[1][1]) - ip.Interval(A_eps[0][1])*ip.Interval(A_eps[1][0])


def intersection_int(A, B):
    return [max(min(A), min(B)), min(max(A), max(B))]


eps_list = [e/10000 for e in range(0, 10000)]
eps_suitable = []
det_in_0 = []
for eps in eps_list:
    A_eps = [
        [[1.05-eps, 1.05+eps], [0.95-eps, 0.95+eps]],
        [[1-eps, 1+eps], [1-eps, 1+eps]]
    ]

    det = determinant_2(A_eps)
    if 0 in det:
        eps_suitable.append(eps)
        det_in_0.append(det)

print("подходящие eps_min=", min(eps_suitable), " eps_max", max(eps_suitable))
print("соответсвующие определители det_min", det_in_0[0], " det_max", det_in_0[-1])

print()
eps = 0.0251
print("eps=", eps)
A_eps = [
        [[1.05-eps, 1.05+eps], [0.95-eps, 0.95+eps]],
        [[1-eps, 1+eps], [1-eps, 1+eps]]
    ]
print(f"det A({eps}) =", ip.Interval(A_eps[0][0])*ip.Interval(A_eps[1][1]) - ip.Interval(A_eps[0][1])*ip.Interval(A_eps[1][0]))
print("A_eps = ", A_eps)
print()
print("A_11 cap A_21", intersection_int(A_eps[0][0], A_eps[1][0]))
print("A_12 cap A_22", intersection_int(A_eps[0][1], A_eps[1][1]))
print("det A'=", 1.0251*0.9751-0.9751*1.0251)
