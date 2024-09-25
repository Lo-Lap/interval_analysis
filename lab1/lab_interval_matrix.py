import intvalpy as ip
from sympy import *


def determinant_2(A):
    return ip.Interval(A_eps[0][0])*ip.Interval(A_eps[1][1]) - ip.Interval(A_eps[0][1])*ip.Interval(A_eps[1][0])


def intersection_int(A, B):
    return [max(min(A), min(B)), min(max(A), max(B))]


def A_epsilon(eps):
    return [
        [[1.05 - eps, 1.05 + eps], [0.95 - eps, 0.95 + eps]],
        [[1 - eps, 1 + eps], [1 - eps, 1 + eps]]
    ]


h = 0.01
eps = 0
A_eps = A_epsilon(eps)
det = determinant_2(A_eps)
iter_not_zero = 1
while not 0 in det:
    print("eps=", eps)
    print(f"det A({eps}) =", determinant_2(A_eps))
    print("A_eps = ", A_eps)
    print()
    eps += h
    A_eps = A_epsilon(eps)
    det = determinant_2(A_eps)
    iter_not_zero += 1
print("eps with h=0.01: ", eps)
print(f"det A({eps}) =", det)
print("A_eps = ", A_eps)

h = 0.0001
min_eps = eps
iter_in = 0
while 0 in det:
    if eps < min_eps:
        min_eps = eps
        iter_in += 1
    eps -= h
    A_eps = A_epsilon(eps)
    det = determinant_2(A_eps)

print()
print("min eps=", min_eps)
print("A = ", A_epsilon(min_eps))
print("det = ", determinant_2(min_eps))
print("Iteration with h=0.01 = ", iter_not_zero)
print("Iteration with h=0.0001 = ", iter_in)

# Поведение при увеличении epsilon
print()
for eps in [1, 5, 10, 15]:
    print("eps=", eps)
    A_eps = A_epsilon(eps)
    print(f"det A({eps}) =", determinant_2(A_eps))
    print("A_eps = ", A_eps)
    print()
print()

# Для поиска точечной матрицы
eps = 0.025
print("eps=", eps)
A_eps = A_epsilon(eps)
print(f"det A({eps}) =", determinant_2(A_eps))
print("A_eps = ", A_eps)
print()
print("A_11 cap A_21", intersection_int(A_eps[0][0], A_eps[1][0]))
print("A_12 cap A_22", intersection_int(A_eps[0][1], A_eps[1][1]))
print("det A'=", 1.0251*0.9751-0.9751*1.0251)

# eps_list = [e * 0.001 for e in range(0, 1000)]
# eps_suitable = []
# det_in_0 = []
# for eps in eps_list:
#     A_eps = [
#         [[1.05-eps, 1.05+eps], [0.95-eps, 0.95+eps]],
#         [[1-eps, 1+eps], [1-eps, 1+eps]]
#     ]
#
#     det = determinant_2(A_eps)
#     iter_ += 1
#     if 0 in det:
#         if eps < min_eps:
#             min_eps = eps
#             iteration_for_min = iter_
#         eps_suitable.append(eps)
#         det_in_0.append(det)

# print("Итерации для min = ", iteration_for_min)
# print("подходящие eps_min=", min(eps_suitable), " eps_max", max(eps_suitable))
# print("соответсвующие определители det_min", det_in_0[0], " det_max", det_in_0[-1])

