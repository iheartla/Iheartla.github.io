"""
given
A ∈ ℝ^(k×k) 
B ∈ ℝ^(k×m) 
C ∈ ℝ^(m×m) 

S = C - BᵀA⁻¹B
[A⁻¹+A⁻¹BS⁻¹BᵀA⁻¹   -A⁻¹BS⁻¹
   -S⁻¹BᵀA⁻¹           S⁻¹]
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_650ResultType:
    def __init__( self, S, ret):
        self.S = S
        self.ret = ret


def convex_optimization_650(A, B, C):
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    C = np.asarray(C, dtype=np.float64)

    k = A.shape[1]
    m = B.shape[1]
    assert A.shape == (k, k)
    assert B.shape == (k, m)
    assert C.shape == (m, m)

    S = C - B.T @ np.linalg.solve(A, B)
    ret_0 = np.block([[np.linalg.inv(A) + np.linalg.solve(A, B) @ np.linalg.solve(S, B.T) @ np.linalg.inv(A), -np.linalg.solve(A, B) @ np.linalg.inv(S)], [-np.linalg.solve(S, B.T) @ np.linalg.inv(A), np.linalg.inv(S)]])
    ret = ret_0
    return convex_optimization_650ResultType(S, ret)


def generateRandomData():
    k = np.random.randint(10)
    m = np.random.randint(10)
    A = np.random.randn(k, k)
    B = np.random.randn(k, m)
    C = np.random.randn(m, m)
    return A, B, C


if __name__ == '__main__':
    A, B, C = generateRandomData()
    print("A:", A)
    print("B:", B)
    print("C:", C)
    func_value = convex_optimization_650(A, B, C)
    print("return value: ", func_value.ret)