function output = convex_optimization_276(A, b, x0)
% output = convex_optimization_276(A, b, `x₀`)
%
%    given
%    A_i ∈ ℝ^(m_i × n)  
%    b_i ∈ ℝ^m_i 
%    `x₀` ∈ ℝ^n  
%    
%    min_(x ∈ ℝ^n) ∑_i ‖A_i x + b_i‖ + (1/2)‖x-`x₀`‖²
    if nargin==0
        warning('generating random input data');
        [A, b, x0] = generateRandomData();
    end
    function [A, b, x0] = generateRandomData()
        dim_0 = randi(10);
        n = randi(10);
        A = {};
        b = {};
        for i = 1:dim_0
            m_2 = randi(10);
            A = [A; randn(m_2, n)];
            b = [b; randn(m_2)];
        end
        x0 = randn(n,1);
    end

    x0 = reshape(x0,[],1);

    dim_0 = size(A, 1);
    n = size(A{1}, 2);
    assert( numel(x0) == n );

    function ret = target_1(x)
        sum_0 = 0;
        for i = 1:size(A, 1)
            sum_0 = sum_0 + norm(A{i} * x + b{i}, 2);
        end
        ret = sum_0 + (1 / 2) * norm(x - x0, 2).^2;
    end
    [~,optimize_0] = fminunc(@target_1,zeros(n,1));
    ret = optimize_0;
    output.ret = ret;
end
