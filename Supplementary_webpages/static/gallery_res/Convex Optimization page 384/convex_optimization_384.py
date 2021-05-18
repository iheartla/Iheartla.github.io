"""
given
a_i ∈ ℝ^n : the measurement vectors  
x ∈ ℝ^n   : original vector 
w_i ∈ ℝ   : measurement noise 

y_i = a_iᵀ x + w_i
x̂ = (∑_i a_i a_iᵀ)⁻¹ ∑_i y_i a_i
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_384ResultType:
    def __init__( self, y, x̂):
        self.y = y
        self.x̂ = x̂


def convex_optimization_384(a, x, w):
    """
    :param :a : the measurement vectors  
    :param :x : original vector 
    :param :w : measurement noise 
    """
    a = np.asarray(a, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)

    dim_0 = w.shape[0]
    n = a.shape[1]
    assert a.shape == (dim_0, n, )
    assert x.shape == (n,)
    assert w.shape == (dim_0,)

    y = np.zeros(dim_0)
    for i in range(1, dim_0+1):
        y[i-1] = (a[i-1].T.reshape(1, n) @ x).item() + w[i-1]
    sum_0 = np.zeros((n, n))
    for i in range(1, len(a)+1):
        sum_0 += (a[i-1]).reshape(n, 1) @ a[i-1].T.reshape(1, n)
    sum_1 = np.zeros((n, ))
    for i in range(1, len(a)+1):
        sum_1 += y[i-1] * a[i-1]
    x̂ = np.linalg.solve((sum_0), sum_1)
    return convex_optimization_384ResultType(y, x̂)


def generateRandomData():
    dim_0 = np.random.randint(10)
    n = np.random.randint(10)
    a = np.random.randn(dim_0, n, )
    x = np.random.randn(n)
    w = np.random.randn(dim_0)
    return a, x, w


if __name__ == '__main__':
    a, x, w = generateRandomData()
    print("a:", a)
    print("x:", x)
    print("w:", w)
    func_value = convex_optimization_384(a, x, w)
    print("return value: ", func_value.x̂)