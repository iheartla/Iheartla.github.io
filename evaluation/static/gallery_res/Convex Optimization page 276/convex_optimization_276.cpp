/*
given
A_i ∈ ℝ^(m_i × n)  
b_i ∈ ℝ^m_i 
`x₀` ∈ ℝ^n  

min_(x ∈ ℝ^n) ∑_i ‖A_i x + b_i‖ + (1/2)‖x-`x₀`‖²
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_276ResultType {
    double ret;
    convex_optimization_276ResultType(const double & ret)
    : ret(ret)
    {}
};

convex_optimization_276ResultType convex_optimization_276(
    const std::vector<Eigen::MatrixXd> & A,
    const std::vector<Eigen::VectorXd> & b,
    const Eigen::VectorXd & x₀)
{
    const long dim_0 = A.size();
    const long n = A[0].cols();
    assert( A.size() == dim_0 );
    assert( b.size() == dim_0 );
    assert( x₀.size() == n );

    double ret = ;
    return convex_optimization_276ResultType(ret);
}


void generateRandomData(std::vector<Eigen::MatrixXd> & A,
    std::vector<Eigen::VectorXd> & b,
    Eigen::VectorXd & x₀)
{
    const int dim_0 = rand()%10;
    const int n = rand()%10;
    A.resize(dim_0);
    b.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        int m_1 = rand()%10;
        A[i] = Eigen::MatrixXd::Random(m_1, n);
        b[i] = Eigen::VectorXd::Random(m_1);
    }
    x₀ = Eigen::VectorXd::Random(n);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::MatrixXd> A;
    std::vector<Eigen::VectorXd> b;
    Eigen::VectorXd x₀;
    generateRandomData(A, b, x₀);
    convex_optimization_276ResultType func_value = convex_optimization_276(A, b, x₀);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}