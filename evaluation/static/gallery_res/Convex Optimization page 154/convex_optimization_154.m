function output = convex_optimization_154(f, p)
% output = convex_optimization_154(f, p)
%
%    given
%    f ∈ ℝ^(n)
%    p ∈ ℝ^(n)
%    
%    ∑_i f_i²p_i - (∑_i f_i p_i)²
    if nargin==0
        warning('generating random input data');
        [f, p] = generateRandomData();
    end
    function [f, p] = generateRandomData()
        n = randi(10);
        f = randn(n,1);
        p = randn(n,1);
    end

    f = reshape(f,[],1);
    p = reshape(p,[],1);

    n = size(f, 1);
    assert( numel(f) == n );
    assert( numel(p) == n );

    sum_0 = 0;
    for i = 1:size(p,1)
        sum_0 = sum_0 + f(i).^2 * p(i);
    end
    sum_1 = 0;
    for i = 1:size(p,1)
        sum_1 = sum_1 + f(i) * p(i);
    end
    ret = sum_0 - (sum_1).^2;
    output.ret = ret;
end
