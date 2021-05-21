/*
given
α_T : ℝ
n_T : ℝ³

`n(v)` = ( ∑_(T for T ∈ `N₁`_v) α_T n_T )/‖ ∑_(T for T ∈ `N₁`_v) α_T n_T ‖

where
v ∈ ℤ
`N₁`_i ∈ {ℤ}
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct pmp_42ResultType {
    Eigen::Matrix<double, 3, 1> n_left_parenthesis_v_right_parenthesis;
    pmp_42ResultType(const Eigen::Matrix<double, 3, 1> & n_left_parenthesis_v_right_parenthesis)
    : n_left_parenthesis_v_right_parenthesis(n_left_parenthesis_v_right_parenthesis)
    {}
};

/**
 * pmp_42
 *
 * @param α  ℝ
 * @param n  ℝ³
 * @return n_left_parenthesis_v_right_parenthesis
 */
pmp_42ResultType pmp_42(
    const std::vector<double> & α,
    const std::vector<Eigen::Matrix<double, 3, 1>> & n,
    const int & v,
    const std::vector<std::set<int >> & N₁)
{
    const long dim_0 = α.size();
    const long dim_1 = N₁.size();
    assert( n.size() == dim_0 );

    Eigen::MatrixXd sum_0 = Eigen::MatrixXd::Zero(3, 1);
    for(int T=1; T<=α.size(); T++){
        std::set<int > set_0 = N₁.at(v-1);
        if(set_0.find(int(T)) != set_0.end()){
            sum_0 += α.at(T-1) * n.at(T-1);
        }
    }
    Eigen::MatrixXd sum_1 = Eigen::MatrixXd::Zero(3, 1);
    for(int T=1; T<=α.size(); T++){
        std::set<int > set_1 = N₁.at(v-1);
        if(set_1.find(int(T)) != set_1.end()){
            sum_1 += α.at(T-1) * n.at(T-1);
        }
    }
    Eigen::Matrix<double, 3, 1> n_left_parenthesis_v_right_parenthesis = (sum_0) / double((sum_1).lpNorm<2>());

    return pmp_42ResultType(n_left_parenthesis_v_right_parenthesis);
}


void generateRandomData(std::vector<double> & α,
    std::vector<Eigen::Matrix<double, 3, 1>> & n,
    int & v,
    std::vector<std::set<int >> & N₁)
{
    v = rand() % 10;
    const int dim_0 = rand()%10;
    const int dim_1 = rand()%10;
    α.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        α[i] = rand() % 10;
    }
    n.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        n[i] = Eigen::VectorXd::Random(3);
    }
    N₁.resize(dim_1);
    for(int i=0; i<dim_1; i++){
        const int dim_3 = rand()%10;
        for(int j=0; j<dim_3; j++){
            N₁[i].insert(rand()%10);
        }
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<double> α;
    std::vector<Eigen::Matrix<double, 3, 1>> n;
    int v;
    std::vector<std::set<int >> N₁;
    generateRandomData(α, n, v, N₁);
    pmp_42ResultType func_value = pmp_42(α, n, v, N₁);
    std::cout<<"return value:\n"<<func_value.n_left_parenthesis_v_right_parenthesis<<std::endl;
    return 0;
}