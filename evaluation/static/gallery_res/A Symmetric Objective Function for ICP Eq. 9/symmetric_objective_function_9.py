"""
from trigonometry: tan, cos

t̃ = t/cos(θ)
ã = a tan(θ)
∑_i cos²(θ)((p_i - q_i)⋅n_i +((p_i+q_i)×n_i)⋅ã+n_i⋅t̃)² 

where
a ∈ ℝ³ : axis of rotation
θ ∈ ℝ  : angle of rotation
p_i ∈ ℝ³
q_i ∈ ℝ³
n_i ∈ ℝ³
t ∈ ℝ³
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class symmetric_objective_function_9ResultType:
    def __init__( self, t̃, ã, ret):
        self.t̃ = t̃
        self.ã = ã
        self.ret = ret


def symmetric_objective_function_9(a, θ, p, q, n, t):
    """
    :param :a : axis of rotation
    :param :θ : angle of rotation
    """
    a = np.asarray(a, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)

    dim_0 = p.shape[0]
    assert a.shape == (3,)
    assert np.ndim(θ) == 0
    assert p.shape == (dim_0, 3, )
    assert q.shape == (dim_0, 3, )
    assert n.shape == (dim_0, 3, )
    assert t.shape == (3,)

    t̃ = t / np.cos(θ)
    ã = a * np.tan(θ)
    sum_0 = 0
    for i in range(1, len(p)+1):
        sum_0 += np.power(np.cos(θ), 2) * np.power((np.dot(((p[i-1] - q[i-1])).ravel(), (n[i-1]).ravel()) + np.dot(((np.cross((p[i-1] + q[i-1]), n[i-1]))).ravel(), (ã).ravel()) + np.dot((n[i-1]).ravel(), (t̃).ravel())), 2)
    ret = sum_0
    return symmetric_objective_function_9ResultType(t̃, ã, ret)


def generateRandomData():
    θ = np.random.randn()
    dim_0 = np.random.randint(10)
    a = np.random.randn(3)
    p = np.random.randn(dim_0, 3, )
    q = np.random.randn(dim_0, 3, )
    n = np.random.randn(dim_0, 3, )
    t = np.random.randn(3)
    return a, θ, p, q, n, t


if __name__ == '__main__':
    a, θ, p, q, n, t = generateRandomData()
    print("a:", a)
    print("θ:", θ)
    print("p:", p)
    print("q:", q)
    print("n:", n)
    print("t:", t)
    func_value = symmetric_objective_function_9(a, θ, p, q, n, t)
    print("return value: ", func_value.ret)