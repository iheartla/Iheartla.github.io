function output = convex_optimization_208(x, p)
% output = convex_optimization_208(x, p)
%
%    `I(X;Y)` = ∑_i ∑_j x_j p_i,j log₂(p_i,j/∑_k x_k p_i,k)
%    
%    where
%    
%    x ∈ ℝ^n
%    p ∈ ℝ^(m×n)
    if nargin==0
        warning('generating random input data');
        [x, p] = generateRandomData();
    end
    function [x, p] = generateRandomData()
        n = randi(10);
        m = randi(10);
        x = randn(n,1);
        p = randn(m, n);
    end

    x = reshape(x,[],1);

    n = size(x, 1);
    m = size(p, 1);
    assert( numel(x) == n );
    assert( isequal(size(p), [m, n]) );

    sum_0 = 0;
    for i = 1:size(p,1)
        sum_1 = 0;
        for j = 1:size(x,1)
            sum_2 = 0;
            for k = 1:size(p,2)
                sum_2 = sum_2 + x(k) * p(i, k);
            end
            sum_1 = sum_1 + x(j) * p(i, j) * log2(p(i, j) / sum_2);
        end
        sum_0 = sum_0 + sum_1;
    end
    I_X_semicolon_Y = sum_0;
    output.I_X_semicolon_Y = I_X_semicolon_Y;
end
