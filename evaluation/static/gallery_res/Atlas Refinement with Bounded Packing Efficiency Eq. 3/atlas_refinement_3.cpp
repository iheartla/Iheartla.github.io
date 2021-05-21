/*
`G_σ(s^k_i)` = ∑_j l_j exp(-dist(`bᵢ`, b_j)²/(2σ²)) `s^k`_j

where
l_j ∈ ℝ : the length of bj
dist: ℝ², ℝ² → ℝ : measures the geodesic distance 
`bᵢ` ∈ ℝ²
b_j ∈ ℝ²
σ ∈ ℝ
`s^k`_j ∈ ℝ² : direction vector
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct atlas_refinement_3ResultType {
    Eigen::Matrix<double, 2, 1> G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis;
    atlas_refinement_3ResultType(const Eigen::Matrix<double, 2, 1> & G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis)
    : G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis(G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis)
    {}
};

/**
 * atlas_refinement_3
 *
 * @param l  the length of bj
 * @param dist  ℝ², ℝ² → ℝ : measures the geodesic distance 
 * @param s_circumflex_accent_k  direction vector
 * @return G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis
 */
atlas_refinement_3ResultType atlas_refinement_3(
    const std::vector<double> & l,
    const std::function<double(Eigen::Matrix<double, 2, 1>, Eigen::Matrix<double, 2, 1>)> & dist,
    const Eigen::Matrix<double, 2, 1> & bᵢ,
    const std::vector<Eigen::Matrix<double, 2, 1>> & b,
    const double & σ,
    const std::vector<Eigen::Matrix<double, 2, 1>> & s_circumflex_accent_k)
{
    const long dim_0 = l.size();
    assert( b.size() == dim_0 );
    assert( s_circumflex_accent_k.size() == dim_0 );

    Eigen::MatrixXd sum_0 = Eigen::MatrixXd::Zero(2, 1);
    for(int j=1; j<=b.size(); j++){
        sum_0 += l.at(j-1) * exp(-pow(dist(bᵢ, b.at(j-1)), 2) / double((2 * pow(σ, 2)))) * s_circumflex_accent_k.at(j-1);
    }
    Eigen::Matrix<double, 2, 1> G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis = sum_0;

    return atlas_refinement_3ResultType(G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis);
}


void generateRandomData(std::vector<double> & l,
    std::function<double(Eigen::Matrix<double, 2, 1>, Eigen::Matrix<double, 2, 1>)> & dist,
    Eigen::Matrix<double, 2, 1> & bᵢ,
    std::vector<Eigen::Matrix<double, 2, 1>> & b,
    double & σ,
    std::vector<Eigen::Matrix<double, 2, 1>> & s_circumflex_accent_k)
{
    σ = rand() % 10;
    const int dim_0 = rand()%10;
    l.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        l[i] = rand() % 10;
    }
    dist = [](Eigen::Matrix<double, 2, 1>, Eigen::Matrix<double, 2, 1>)->double{
        return rand() % 10;
    };
    bᵢ = Eigen::VectorXd::Random(2);
    b.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        b[i] = Eigen::VectorXd::Random(2);
    }
    s_circumflex_accent_k.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        s_circumflex_accent_k[i] = Eigen::VectorXd::Random(2);
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<double> l;
    std::function<double(Eigen::Matrix<double, 2, 1>, Eigen::Matrix<double, 2, 1>)> dist;
    Eigen::Matrix<double, 2, 1> bᵢ;
    std::vector<Eigen::Matrix<double, 2, 1>> b;
    double σ;
    std::vector<Eigen::Matrix<double, 2, 1>> s_circumflex_accent_k;
    generateRandomData(l, dist, bᵢ, b, σ, s_circumflex_accent_k);
    atlas_refinement_3ResultType func_value = atlas_refinement_3(l, dist, bᵢ, b, σ, s_circumflex_accent_k);
    std::cout<<"return value:\n"<<func_value.G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis<<std::endl;
    return 0;
}