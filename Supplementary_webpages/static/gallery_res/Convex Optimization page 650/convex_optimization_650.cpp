/*
given
A ∈ ℝ^(k×k) 
B ∈ ℝ^(k×m) 
C ∈ ℝ^(m×m) 

S = C - BᵀA⁻¹B
[A⁻¹+A⁻¹BS⁻¹BᵀA⁻¹   -A⁻¹BS⁻¹
   -S⁻¹BᵀA⁻¹           S⁻¹]
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_650ResultType {
    Eigen::MatrixXd S;
    Eigen::MatrixXd ret;
    convex_optimization_650ResultType(const Eigen::MatrixXd & S,
               const Eigen::MatrixXd & ret)
    : S(S),
    ret(ret)
    {}
};

convex_optimization_650ResultType convex_optimization_650(
    const Eigen::MatrixXd & A,
    const Eigen::MatrixXd & B,
    const Eigen::MatrixXd & C)
{
    const long k = A.cols();
    const long m = B.cols();
    assert( A.rows() == k );
    assert( B.rows() == k );
    assert( C.rows() == m );
    assert( C.cols() == m );

    Eigen::MatrixXd S = C - B.transpose() * A.colPivHouseholderQr().solve(B);

    Eigen::MatrixXd ret_0(k+m, k+m);
    ret_0 << A.inverse() + A.colPivHouseholderQr().solve(B) * S.colPivHouseholderQr().solve(B.transpose()) * A.inverse(), -A.colPivHouseholderQr().solve(B) * S.inverse(),
    -S.colPivHouseholderQr().solve(B.transpose()) * A.inverse(), S.inverse();
    Eigen::MatrixXd ret = ret_0;
    return convex_optimization_650ResultType(S, ret);
}


void generateRandomData(Eigen::MatrixXd & A,
    Eigen::MatrixXd & B,
    Eigen::MatrixXd & C)
{
    const int k = rand()%10;
    const int m = rand()%10;
    A = Eigen::MatrixXd::Random(k, k);
    B = Eigen::MatrixXd::Random(k, m);
    C = Eigen::MatrixXd::Random(m, m);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd C;
    generateRandomData(A, B, C);
    convex_optimization_650ResultType func_value = convex_optimization_650(A, B, C);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}