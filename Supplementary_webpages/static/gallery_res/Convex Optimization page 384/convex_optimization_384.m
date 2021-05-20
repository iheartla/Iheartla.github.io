function output = convex_optimization_384(a, x, w)
% output = convex_optimization_384(a, x, w)
%
%    given
%    a_i ∈ ℝ^n : the measurement vectors  
%    x ∈ ℝ^n   : original vector 
%    w_i ∈ ℝ   : measurement noise 
%    
%    y_i = a_iᵀ x + w_i
%    x̂ = (∑_i a_i a_iᵀ)⁻¹ ∑_i y_i a_i
    if nargin==0
        warning('generating random input data');
        [a, x, w] = generateRandomData();
    end
    function [a, x, w] = generateRandomData()
        dim_0 = randi(10);
        n = randi(10);
        a = randn(dim_0,n);
        x = randn(n,1);
        w = randn(dim_0,1);
    end

    x = reshape(x,[],1);
    w = reshape(w,[],1);

    dim_0 = size(w, 1);
    n = size(a, 2);
    assert( isequal(size(a), [dim_0, n]) );
    assert( numel(x) == n );
    assert( size(w,1) == dim_0 );

    y = zeros(dim_0,1);
    for i = 1:dim_0
        y(i) = a(i,:)'' * x + w(i);
    end
    sum_0 = zeros(n, n);
    for i = 1:size(a, 1)
        sum_0 = sum_0 + reshape(a(i,:)', [n, 1]) * a(i,:)'';
    end
    sum_1 = zeros(n,1);
    for i = 1:size(y,1)
        sum_1 = sum_1 + y(i) * a(i,:)';
    end
    x_hat = ((sum_0)\sum_1);
    output.y = y;

    output.x_hat = x_hat;
end
