import sys
import os

sys.path.append(os.path.abspath('C:/Users/andre/Desktop/repo_2/src/main'))

from assignment_2 import neville_int, newton_forward_difference, hermite_divided_difference, print_hermite_matrix, cubic_spline

import numpy as np

#1
def test_neville():
    x = [3.6, 3.8, 3.9]
    val = [1.675, 1.436, 1.318]
    w = 3.7
    interpolated_value = neville_int(x, val, w)
    print(f"2nd degree interpolation at x={w}: {interpolated_value:.7f}")

#2
def test_newton():
    print("\n")
    xi = [7.2, 7.4, 7.5, 7.6]
    fxi = [23.5492, 25.3913, 26.8224, 27.4589]
    newton_forward_difference(xi, fxi)

#3
def test_hermite():
    print("\n")
    x = [3.6, 3.8, 3.9]
    f = [1.675, 1.436, 1.318]
    f_prime = [-1.195, -1.188, -1.182]
    div_diff = hermite_divided_difference(x, f, f_prime)
    print_hermite_matrix(div_diff)

#4
def test_cubic_spline():
    x = [2, 5, 8, 10]
    f_x = [3, 5, 7, 9]
    A, b, c, a, b_spline, d = cubic_spline(x, f_x)

if __name__ == "__main__":
    test_neville()
    test_newton()
    test_hermite()
    test_cubic_spline()