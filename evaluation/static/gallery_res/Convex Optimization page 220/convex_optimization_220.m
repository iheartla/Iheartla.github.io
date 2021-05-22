function output = convex_optimization_220(x, W, nu)
% output = convex_optimization_220(x, W, ν)
%
%    `L(x,ν)` = xᵀWx + ∑_i ν_i(x_i²-1)
%    
%    where
%    
%    x ∈ ℝ^n
%    W ∈ ℝ^(n×n)
%    ν ∈ ℝ^n
    if nargin==0
        warning('generating random input data');
        [x, W, nu] = generateRandomData();
    end
    function [x, W, nu] = generateRandomData()
        n = randi(10);
        x = randn(n,1);
        W = randn(n, n);
        nu = randn(n,1);
    end

    x = reshape(x,[],1);
    nu = reshape(nu,[],1);

    n = size(x, 1);
    assert( numel(x) == n );
    assert( isequal(size(W), [n, n]) );
    assert( numel(nu) == n );

    sum_0 = 0;
    for i = 1:size(x,1)
        sum_0 = sum_0 + nu(i) * (x(i).^2 - 1);
    end
    L_x_nu = x' * W * x + sum_0;
    output.L_x_nu = L_x_nu;
end
