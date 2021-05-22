"""
[`P₁`  0
  0   `P₃`][      L          0
            `P₃`ᵀC`P₂`ᵀU⁻¹  -L̃][U   L⁻¹`P₁`ᵀB
                                0       Ũ    ][`P₂`  0
                                                0   I_n]

where

`P₁` ∈ ℝ^(m×m) 
`P₂` ∈ ℝ^(m×m) 
`P₃` ∈ ℝ^(n×n) 
  B ∈ ℝ^(m×n) 
  C ∈ ℝ^(n×m) 
  L ∈ ℝ^(m×m) 
  L̃ ∈ ℝ^(n×n) 
  U ∈ ℝ^(m×m) 
  Ũ ∈ ℝ^(n×n)
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_680ResultType:
    def __init__( self, ret):
        self.ret = ret


def convex_optimization_680(P1, P2, P3, B, C, L, L̃, U, Ũ):
    P1 = np.asarray(P1, dtype=np.float64)
    P2 = np.asarray(P2, dtype=np.float64)
    P3 = np.asarray(P3, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    C = np.asarray(C, dtype=np.float64)
    L = np.asarray(L, dtype=np.float64)
    L̃ = np.asarray(L̃, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)
    Ũ = np.asarray(Ũ, dtype=np.float64)

    m = P1.shape[1]
    n = P3.shape[1]
    assert P1.shape == (m, m)
    assert P2.shape == (m, m)
    assert P3.shape == (n, n)
    assert B.shape == (m, n)
    assert C.shape == (n, m)
    assert L.shape == (m, m)
    assert L̃.shape == (n, n)
    assert U.shape == (m, m)
    assert Ũ.shape == (n, n)

    ret_0 = np.block([[P1, np.zeros((m, n))], [np.zeros((n, m)), P3]])
    ret_1 = np.block([[L, np.zeros((m, n))], [P3.T @ C @ P2.T @ np.linalg.inv(U), -L̃]])
    ret_2 = np.block([[U, np.linalg.solve(L, P1.T) @ B], [np.zeros((n, m)), Ũ]])
    ret_3 = np.block([[P2, np.zeros((m, n))], [np.zeros((n, m)), np.identity(n)]])
    ret = ret_0 @ ret_1 @ ret_2 @ ret_3
    return convex_optimization_680ResultType(ret)


def generateRandomData():
    m = np.random.randint(10)
    n = np.random.randint(10)
    P1 = np.random.randn(m, m)
    P2 = np.random.randn(m, m)
    P3 = np.random.randn(n, n)
    B = np.random.randn(m, n)
    C = np.random.randn(n, m)
    L = np.random.randn(m, m)
    L̃ = np.random.randn(n, n)
    U = np.random.randn(m, m)
    Ũ = np.random.randn(n, n)
    return P1, P2, P3, B, C, L, L̃, U, Ũ


if __name__ == '__main__':
    P1, P2, P3, B, C, L, L̃, U, Ũ = generateRandomData()
    print("P1:", P1)
    print("P2:", P2)
    print("P3:", P3)
    print("B:", B)
    print("C:", C)
    print("L:", L)
    print("L̃:", L̃)
    print("U:", U)
    print("Ũ:", Ũ)
    func_value = convex_optimization_680(P1, P2, P3, B, C, L, L̃, U, Ũ)
    print("return value: ", func_value.ret)