/*
given
f ∈ ℝ^(n)
p ∈ ℝ^(n)

∑_i f_i²p_i - (∑_i f_i p_i)²
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_154ResultType {
    double ret;
    convex_optimization_154ResultType(const double & ret)
    : ret(ret)
    {}
};

convex_optimization_154ResultType convex_optimization_154(
    const Eigen::VectorXd & f,
    const Eigen::VectorXd & p)
{
    const long n = f.size();
    assert( p.size() == n );

    double sum_0 = 0;
    for(int i=1; i<=f.size(); i++){
        sum_0 += pow(f[i-1], 2) * p[i-1];
    }
    double sum_1 = 0;
    for(int i=1; i<=f.size(); i++){
        sum_1 += f[i-1] * p[i-1];
    }
    double ret = sum_0 - pow((sum_1), 2);
    return convex_optimization_154ResultType(ret);
}


void generateRandomData(Eigen::VectorXd & f,
    Eigen::VectorXd & p)
{
    const int n = rand()%10;
    f = Eigen::VectorXd::Random(n);
    p = Eigen::VectorXd::Random(n);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::VectorXd f;
    Eigen::VectorXd p;
    generateRandomData(f, p);
    convex_optimization_154ResultType func_value = convex_optimization_154(f, p);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}