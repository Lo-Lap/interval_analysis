import intvalpy as ip
from sympy import *


def determinant_2(A):
    return ip.Interval(A_eps[0][0])*ip.Interval(A_eps[1][1]) - ip.Interval(A_eps[0][1])*ip.Interval(A_eps[1][0])


eps_list = [e/10000 for e in range(1, 10000)]
min_eps = 2
min_det = ip.Interval([0, 0])
for eps in eps_list:
    A_eps = [
        [[1-eps, 1+eps], [0.9-eps, 0.9+eps]],
        [[1-eps, 1+eps], [1.1-eps, 1.1+eps]]
    ]

    det = determinant_2(A_eps)
    if 0 in det and eps < min_eps:
        min_eps = eps
        min_det = det

print("min_eps=", min_eps)
print(min_det)
# inter = ip.Interval([1, 5])
print()
eps = 0.0
print("eps=", eps)
A_eps = [
        [[1-eps, 1+eps], [0.9-eps, 0.9+eps]],
        [[1-eps, 1+eps], [1.1-eps, 1.1+eps]]
    ]
# print(A_eps[1][1])
print(ip.Interval(A_eps[0][0])*ip.Interval(A_eps[1][1]) - ip.Interval(A_eps[0][1])*ip.Interval(A_eps[1][0]))
