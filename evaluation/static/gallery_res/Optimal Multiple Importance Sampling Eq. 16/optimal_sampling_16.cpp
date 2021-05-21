/*
 ∑_i α_i + 1/M ∑_i ∑_j (f(X_i,j)/`p_c`(X_i,j) - (∑_k α_k p_k(X_i,j))/`p_c`(X_i,j))

where

α ∈ ℝ^N
p_j ∈ ℝ → ℝ 
X_i ∈ ℝ^(n_i) 
M ∈ ℝ
f: ℝ → ℝ 
`p_c`: ℝ → ℝ
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct optimal_sampling_16ResultType {
    double ret;
    optimal_sampling_16ResultType(const double & ret)
    : ret(ret)
    {}
};

/**
 * optimal_sampling_16
 *
 * @param f  ℝ → ℝ
 * @param p_c  ℝ → ℝ
 * @return ret
 */
optimal_sampling_16ResultType optimal_sampling_16(
    const Eigen::VectorXd & α,
    const std::vector<std::function<double(double)>> & p,
    const std::vector<Eigen::VectorXd> & X,
    const double & M,
    const std::function<double(double)> & f,
    const std::function<double(double)> & p_c)
{
    const long N = α.size();
    const long dim_0 = X.size();
    const long dim_1 = p.size();
    assert( dim_1 == N );

    double sum_0 = 0;
    for(int i=1; i<=α.size(); i++){
        sum_0 += α[i-1];
    }
    double sum_1 = 0;
    for(int i=1; i<=X.size(); i++){
        double sum_2 = 0;
        for(int j=1; j<=X.at(i-1).rows(); j++){
            double sum_3 = 0;
            for(int k=1; k<=α.size(); k++){
                sum_3 += α[k-1] * p.at(k-1)(X.at(i-1)[j-1]);
            }
            sum_2 += (f(X.at(i-1)[j-1]) / double(p_c(X.at(i-1)[j-1])) - (sum_3) / double(p_c(X.at(i-1)[j-1])));
        }
        sum_1 += sum_2;
    }
    double ret = sum_0 + 1 / double(M) * sum_1;
    return optimal_sampling_16ResultType(ret);
}


void generateRandomData(Eigen::VectorXd & α,
    std::vector<std::function<double(double)>> & p,
    std::vector<Eigen::VectorXd> & X,
    double & M,
    std::function<double(double)> & f,
    std::function<double(double)> & p_c)
{
    M = rand() % 10;
    const int N = rand()%10;
    const int dim_1 = N;
    const int dim_0 = rand()%10;
    α = Eigen::VectorXd::Random(N);
    p.resize(dim_1);
    for(int i=0; i<dim_1; i++){
        p[i] = [](double)->double{
            return rand() % 10;
        };
    }
    X.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        X[i] = Eigen::VectorXd::Random(rand()%10);
    }
    f = [](double)->double{
        return rand() % 10;
    };
    p_c = [](double)->double{
        return rand() % 10;
    };
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::VectorXd α;
    std::vector<std::function<double(double)>> p;
    std::vector<Eigen::VectorXd> X;
    double M;
    std::function<double(double)> f;
    std::function<double(double)> p_c;
    generateRandomData(α, p, X, M, f, p_c);
    optimal_sampling_16ResultType func_value = optimal_sampling_16(α, p, X, M, f, p_c);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}