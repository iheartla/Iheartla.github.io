/*
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
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_680ResultType {
    Eigen::MatrixXd ret;
    convex_optimization_680ResultType(const Eigen::MatrixXd & ret)
    : ret(ret)
    {}
};

convex_optimization_680ResultType convex_optimization_680(
    const Eigen::MatrixXd & P₁,
    const Eigen::MatrixXd & P₂,
    const Eigen::MatrixXd & P₃,
    const Eigen::MatrixXd & B,
    const Eigen::MatrixXd & C,
    const Eigen::MatrixXd & L,
    const Eigen::MatrixXd & L̃,
    const Eigen::MatrixXd & U,
    const Eigen::MatrixXd & Ũ)
{
    const long m = P₁.cols();
    const long n = P₃.cols();
    assert( P₁.rows() == m );
    assert( P₂.rows() == m );
    assert( P₂.cols() == m );
    assert( P₃.rows() == n );
    assert( B.rows() == m );
    assert( B.cols() == n );
    assert( C.rows() == n );
    assert( C.cols() == m );
    assert( L.rows() == m );
    assert( L.cols() == m );
    assert( L̃.rows() == n );
    assert( L̃.cols() == n );
    assert( U.rows() == m );
    assert( U.cols() == m );
    assert( Ũ.rows() == n );
    assert( Ũ.cols() == n );

    Eigen::MatrixXd ret_0(m+n, m+n);
    ret_0 << P₁, Eigen::MatrixXd::Zero(m, n),
    Eigen::MatrixXd::Zero(n, m), P₃;
    Eigen::MatrixXd ret_1(m+n, m+n);
    ret_1 << L, Eigen::MatrixXd::Zero(m, n),
    P₃.transpose() * C * P₂.transpose() * U.inverse(), -L̃;
    Eigen::MatrixXd ret_2(m+n, m+n);
    ret_2 << U, L.colPivHouseholderQr().solve(P₁.transpose()) * B,
    Eigen::MatrixXd::Zero(n, m), Ũ;
    Eigen::MatrixXd ret_3(m+n, m+n);
    ret_3 << P₂, Eigen::MatrixXd::Zero(m, n),
    Eigen::MatrixXd::Zero(n, m), Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd ret = ret_0 * ret_1 * ret_2 * ret_3;
    return convex_optimization_680ResultType(ret);
}


void generateRandomData(Eigen::MatrixXd & P₁,
    Eigen::MatrixXd & P₂,
    Eigen::MatrixXd & P₃,
    Eigen::MatrixXd & B,
    Eigen::MatrixXd & C,
    Eigen::MatrixXd & L,
    Eigen::MatrixXd & L̃,
    Eigen::MatrixXd & U,
    Eigen::MatrixXd & Ũ)
{
    const int m = rand()%10;
    const int n = rand()%10;
    P₁ = Eigen::MatrixXd::Random(m, m);
    P₂ = Eigen::MatrixXd::Random(m, m);
    P₃ = Eigen::MatrixXd::Random(n, n);
    B = Eigen::MatrixXd::Random(m, n);
    C = Eigen::MatrixXd::Random(n, m);
    L = Eigen::MatrixXd::Random(m, m);
    L̃ = Eigen::MatrixXd::Random(n, n);
    U = Eigen::MatrixXd::Random(m, m);
    Ũ = Eigen::MatrixXd::Random(n, n);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::MatrixXd P₁;
    Eigen::MatrixXd P₂;
    Eigen::MatrixXd P₃;
    Eigen::MatrixXd B;
    Eigen::MatrixXd C;
    Eigen::MatrixXd L;
    Eigen::MatrixXd L̃;
    Eigen::MatrixXd U;
    Eigen::MatrixXd Ũ;
    generateRandomData(P₁, P₂, P₃, B, C, L, L̃, U, Ũ);
    convex_optimization_680ResultType func_value = convex_optimization_680(P₁, P₂, P₃, B, C, L, L̃, U, Ũ);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}