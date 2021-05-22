"""
given
α_T : ℝ
n_T : ℝ³

`n(v)` = ( ∑_(T for T ∈ `N₁`_v) α_T n_T )/‖ ∑_(T for T ∈ `N₁`_v) α_T n_T ‖

where
v ∈ ℤ
`N₁`_i ∈ {ℤ}
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class pmp_42ResultType:
    def __init__( self, n_left_parenthesis_v_right_parenthesis):
        self.n_left_parenthesis_v_right_parenthesis = n_left_parenthesis_v_right_parenthesis


def pmp_42(α, n, v, N1):
    """
    :param :α : ℝ
    :param :n : ℝ³
    """
    α = np.asarray(α, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)

    dim_0 = α.shape[0]
    dim_1 = N1.shape[0]
    assert α.shape == (dim_0,)
    assert n.shape == (dim_0, 3, )
    assert np.ndim(v) == 0

    sum_0 = np.zeros((3, ))
    for T in range(1, len(n)+1):
        if((T) in N1[v-1]):
            sum_0 += α[T-1] * n[T-1]
    sum_1 = np.zeros((3, ))
    for T in range(1, len(n)+1):
        if((T) in N1[v-1]):
            sum_1 += α[T-1] * n[T-1]
    n_left_parenthesis_v_right_parenthesis = (sum_0) / np.linalg.norm(sum_1, 2)
    return pmp_42ResultType(n_left_parenthesis_v_right_parenthesis)


def generateRandomData():
    v = np.random.randint(10)
    dim_0 = np.random.randint(10)
    dim_1 = np.random.randint(10)
    α = np.random.randn(dim_0)
    n = np.random.randn(dim_0, 3, )
    N1 = []
    for i in range(dim_1):
        N1_tmp = []
        dim_2 = np.random.randint(1, 10)
        for j in range(dim_2):
            N1_tmp.append((np.random.randint(10)))
        N1.append(N1_tmp)
    N1 = np.asarray(N1)
    return α, n, v, N1


if __name__ == '__main__':
    α, n, v, N1 = generateRandomData()
    print("α:", α)
    print("n:", n)
    print("v:", v)
    print("N1:", N1)
    func_value = pmp_42(α, n, v, N1)
    print("return value: ", func_value.n_left_parenthesis_v_right_parenthesis)