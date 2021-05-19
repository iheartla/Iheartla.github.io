/*
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
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct symmetric_objective_function_9ResultType {
    Eigen::Matrix<double, 3, 1> t̃;
    Eigen::Matrix<double, 3, 1> ã;
    double ret;
    symmetric_objective_function_9ResultType(const Eigen::Matrix<double, 3, 1> & t̃,
               const Eigen::Matrix<double, 3, 1> & ã,
               const double & ret)
    : t̃(t̃),
    ã(ã),
    ret(ret)
    {}
};

/**
 * symmetric_objective_function_9
 *
 * @param a  axis of rotation
 * @param θ  angle of rotation
 * @return ret
 */
symmetric_objective_function_9ResultType symmetric_objective_function_9(
    const Eigen::Matrix<double, 3, 1> & a,
    const double & θ,
    const std::vector<Eigen::Matrix<double, 3, 1>> & p,
    const std::vector<Eigen::Matrix<double, 3, 1>> & q,
    const std::vector<Eigen::Matrix<double, 3, 1>> & n,
    const Eigen::Matrix<double, 3, 1> & t)
{
    const long dim_0 = p.size();
    assert( q.size() == dim_0 );
    assert( n.size() == dim_0 );

    Eigen::Matrix<double, 3, 1> t̃ = t / double(cos(θ));

    Eigen::Matrix<double, 3, 1> ã = a * tan(θ);

    double sum_0 = 0;
    for(int i=1; i<=p.size(); i++){
        sum_0 += pow(cos(θ), 2) * pow((((p.at(i-1) - q.at(i-1))).dot(n.at(i-1)) + ((((p.at(i-1) + q.at(i-1))).cross(n.at(i-1)))).dot(ã) + (n.at(i-1)).dot(t̃)), 2);
    }
    double ret = sum_0;
    return symmetric_objective_function_9ResultType(t̃, ã, ret);
}


void generateRandomData(Eigen::Matrix<double, 3, 1> & a,
    double & θ,
    std::vector<Eigen::Matrix<double, 3, 1>> & p,
    std::vector<Eigen::Matrix<double, 3, 1>> & q,
    std::vector<Eigen::Matrix<double, 3, 1>> & n,
    Eigen::Matrix<double, 3, 1> & t)
{
    θ = rand() % 10;
    const int dim_0 = rand()%10;
    a = Eigen::VectorXd::Random(3);
    p.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        p[i] = Eigen::VectorXd::Random(3);
    }
    q.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        q[i] = Eigen::VectorXd::Random(3);
    }
    n.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        n[i] = Eigen::VectorXd::Random(3);
    }
    t = Eigen::VectorXd::Random(3);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::Matrix<double, 3, 1> a;
    double θ;
    std::vector<Eigen::Matrix<double, 3, 1>> p;
    std::vector<Eigen::Matrix<double, 3, 1>> q;
    std::vector<Eigen::Matrix<double, 3, 1>> n;
    Eigen::Matrix<double, 3, 1> t;
    generateRandomData(a, θ, p, q, n, t);
    symmetric_objective_function_9ResultType func_value = symmetric_objective_function_9(a, θ, p, q, n, t);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}