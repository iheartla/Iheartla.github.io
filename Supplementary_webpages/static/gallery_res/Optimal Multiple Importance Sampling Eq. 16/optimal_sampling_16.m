function output = optimal_sampling_16(alpha, p, X, M, f, p_c)
% output = optimal_sampling_16(α, p, X, M, f, `p_c`)
%
%     ∑_i α_i + 1/M ∑_i ∑_j (f(X_i,j)/`p_c`(X_i,j) - (∑_k α_k p_k(X_i,j))/`p_c`(X_i,j))
%    
%    where
%    
%    α ∈ ℝ^N
%    p_j ∈ ℝ → ℝ 
%    X_i ∈ ℝ^(n_i) 
%    M ∈ ℝ
%    f: ℝ → ℝ 
%    `p_c`: ℝ → ℝ
    if nargin==0
        warning('generating random input data');
        [alpha, p, X, M, f, p_c] = generateRandomData();
    end
    function [alpha, p, X, M, f, p_c] = generateRandomData()
        M = randn();
        N = randi(10);
        dim_1 = N;
        dim_0 = randi(10);
        alpha = randn(N,1);
        p = {};
        for i = 1:dim_1
            p_f = @(p0) randn();
            p{end+1,1} = p_f;
        end
        X = randn(dim_0,randi(10));
        f = @fFunc;
        rseed = randi(2^32);
        function tmp =  fFunc(p0)
            rng(rseed);
            tmp = randn();
        end

        p_c = @p_cFunc;
        rseed = randi(2^32);
        function tmp =  p_cFunc(p0)
            rng(rseed);
            tmp = randn();
        end

    end

    alpha = reshape(alpha,[],1);

    N = size(alpha, 1);
    dim_0 = size(X, 1);
    dim_1 = size(p, 1);
    assert( numel(alpha) == N );
    assert(numel(M) == 1);
    assert( N == dim_1 );

    sum_0 = 0;
    for i = 1:size(alpha,1)
        sum_0 = sum_0 + alpha(i);
    end
    sum_1 = 0;
    for i = 1:size(size(X, 1),1)
        sum_2 = 0;
        for j = 1:size(size(X(i), 1),1)
            sum_3 = 0;
            for k = 1:size(alpha,1)
                sum_3 = sum_3 + alpha(k) * p{k}(X(i,j));
            end
            sum_2 = sum_2 + (f(X(i,j)) / p_c(X(i,j)) - (sum_3) / p_c(X(i,j)));
        end
        sum_1 = sum_1 + sum_2;
    end
    ret = sum_0 + 1 / M * sum_1;
    output.ret = ret;
end
