"""
given
A_i ∈ ℝ^(m × n)  
b_i ∈ ℝ^m  
`x₀` ∈ ℝ^n  

min_(x ∈ ℝ^n) ∑_i ‖A_i x + b_i‖ + (1/2)‖x-`x₀`‖²
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_276ResultType:
    def __init__( self, ret):
        self.ret = ret


def convex_optimization_276(A, b, x0):
    A = np.asarray(A, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    x0 = np.asarray(x0, dtype=np.float64)

    dim_0 = A.shape[0]
    m = A.shape[1]
    n = A.shape[2]
    assert A.shape == (dim_0, m, n)
    assert b.shape == (dim_0, m, )
    assert x0.shape == (n,)

    def target_0(x):
        sum_0 = 0
        for i in range(1, len(A)+1):
            sum_0 += np.linalg.norm(A[i-1] @ x + b[i-1], 2)
        return sum_0 + (1 / 2) * np.power(np.linalg.norm(x - x0, 2), 2)
    ret = minimize(target_0, np.zeros(n)).fun
    return convex_optimization_276ResultType(ret)


def generateRandomData():
    dim_0 = np.random.randint(10)
    m = np.random.randint(10)
    n = np.random.randint(10)
    A = np.random.randn(dim_0, m, n)
    b = np.random.randn(dim_0, m, )
    x0 = np.random.randn(n)
    return A, b, x0


if __name__ == '__main__':
    A, b, x0 = generateRandomData()
    print("A:", A)
    print("b:", b)
    print("x0:", x0)
    func_value = convex_optimization_276(A, b, x0)
    print("return value: ", func_value.ret)