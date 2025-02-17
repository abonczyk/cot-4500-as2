import numpy as np
#1

x = [3.6, 3.8, 3.9]
val = [1.675, 1.436, 1.318]
w = 3.7
n = len(x)
i, j = 0, 0
term1, term2 = 0, 0
neville = [[0] * n for i in range(n)]

for i in range(n):
  neville[i][0] = val[i]
for i in range(1, n):
  for j in range(1, i + 1):
    term1 = (w - x[i - j]) * neville[i][j - 1]
    term2 = (w - x[i]) * neville[i - 1][j - 1]
    neville[i][j] = (term1 - term2) / (x[i] - x[i - j])

print(f"2nd degree interpolation at x={w}: {neville[2][2]:.7f}")

#2

print("\n")
i = 0
xi = [7.2, 7.4, 7.5, 7.6]
fxi = [23.5492, 25.3913, 26.8224, 27.4589]
def newton_forward_difference(xi, fxi):
    n = len(xi)
    diffs = [[0] * n for i in range(n)]


    for i in range(n):
        diffs[i][0] = fxi[i]
    for i in range(1, n):
        for j in range(i, n):
            diffs[j][i] = (diffs[j][i-1] - diffs[j-1][i-1]) / (xi[j] - xi[j-i])

    print("Divided difference table:")
    for i in range(n):
        for j in range(i + 1):
            print(f"{diffs[i][j]:.7f}", end=" ")
        print()


    print("\nGeneral forms of polynomial approximations:")
    for degree in range(1, 4):
        poly_str = f"P{degree}(x) = {fxi[0]:.7f}"
        for i in range(degree):
            term = f" + {diffs[i+1][i+1]:.7f}"
            for j in range(i+1):
                term += f"(x - {xi[j]})"
            poly_str += term
        print(f"Degree {degree}: {poly_str}")

    print("\nEvaluations at x = 7.3:")
    x = 7.3
    for degree in range(1, 4):
        result = fxi[0]
        product = 1.0
        for i in range(degree):
            product *= (x - xi[i])
            result += diffs[i+1][i+1] * product
        print(f"P{degree}({x}) = {result:.7f}")

newton_forward_difference(xi, fxi)

print("\n")

#3

def hermite_divided_difference(x_vals, f_vals, f_prime_vals):

    n = len(x_vals)

    div_diff = np.zeros((2 * n, 2 * n))

    for i in range(n):
        div_diff[2 * i][0] = f_vals[i]     
        div_diff[2 * i + 1][0] = f_prime_vals[i]  

    for j in range(1, 2 * n):
        for i in range(2 * n - j):
            if x_vals[i // 2] != x_vals[(i + j) // 2]:
                div_diff[i][j] = (div_diff[i + 1][j - 1] - div_diff[i][j - 1]) / (x_vals[(i + j) // 2] - x_vals[i // 2])
            else:
                div_diff[i][j] = div_diff[i + 1][j - 1]

    return div_diff

def print_hermite_matrix(div_diff):
    
    print("Hermite Polynomial Approximation Matrix:")
    print(div_diff)

x = [3.6, 3.8, 3.9] 
f = [1.675, 1.436, -1.188] 
f_prime = [3.9, 1.318, -1.182]

div_diff = hermite_divided_difference(x, f, f_prime)

print_hermite_matrix(div_diff)

#4


x = [2, 5, 8, 10]
f_x = [3, 5, 7, 9]

def cubic_spline_interpolation(x, f_x):
    n = len(x)

    h = np.diff(x)

    A = np.zeros((n, n))
    b = np.zeros(n)

    A[0, 0] = 1
    A[-1, -1] = 1


    for i in range(1, n-1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
        b[i] = 3 * ((f_x[i+1] - f_x[i]) / h[i] - (f_x[i] - f_x[i-1]) / h[i-1])

    c = np.linalg.solve(A, b)


    a = np.zeros(n-1)
    b_spline = np.zeros(n-1)
    d = np.zeros(n-1)

    for i in range(n-1):
        a[i] = (f_x[i+1] - f_x[i]) / h[i] - h[i] * (2*c[i] + c[i+1]) / 3
        b_spline[i] = c[i]
        d[i] = (c[i+1] - c[i]) / (3 * h[i])


    print("\nMatrix A:")
    print(A)
    print("\nVector b:")
    print(b)
    print("\nVector c:")
    print(c)
    
    return A, b, c, a, b_spline, d

A, b, c, a, b_spline, d = cubic_spline_interpolation(x, f_x)