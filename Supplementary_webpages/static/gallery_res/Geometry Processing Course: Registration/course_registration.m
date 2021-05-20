function output = course_registration(x, n_hat, p)
% output = course_registration(x, n̂, p)
%
%    min_(u ∈ ℝ⁶) uᵀ(∑_i [x_i×n̂_i
%                           n̂_i  ][(x_i×n̂_i)ᵀ n̂_iᵀ])u - 2uᵀ(∑_i [x_i×n̂_i
%                                                                   n̂_i  ]n̂_iᵀ(p_i-x_i)) + ∑_i(p_i-x_i)ᵀn̂_i n̂_iᵀ(p_i-x_i)
%    
%    where
%    
%    x_i ∈ ℝ³
%    n̂_i ∈ ℝ³ 
%    p_i ∈ ℝ³
    if nargin==0
        warning('generating random input data');
        [x, n_hat, p] = generateRandomData();
    end
    function [x, n_hat, p] = generateRandomData()
        dim_0 = randi(10);
        x = randn(dim_0,3);
        n_hat = randn(dim_0,3);
        p = randn(dim_0,3);
    end

    dim_0 = size(x, 1);
    assert( isequal(size(x), [dim_0, 3]) );
    assert( isequal(size(n_hat), [dim_0, 3]) );
    assert( isequal(size(p), [dim_0, 3]) );

    function ret = target_1(u)
        sum_0 = zeros(6, 6);
        for i = 1:size(x, 1)
            ret_0 = [[reshape(cross(x(i,:)', n_hat(i,:)'), [3, 1])]; [reshape(n_hat(i,:)', [3, 1])]];
            ret_1 = [[(cross(x(i,:)', n_hat(i,:)'))', n_hat(i,:)'']];
            sum_0 = sum_0 + ret_0 * ret_1;
        end
        sum_1 = zeros(6,1);
        for i = 1:size(p, 1)
            ret_2 = [[reshape(cross(x(i,:)', n_hat(i,:)'), [3, 1])]; [reshape(n_hat(i,:)', [3, 1])]];
            sum_1 = sum_1 + ret_2 * n_hat(i,:)'' * (p(i,:)' - x(i,:)');
        end
        sum_2 = 0;
        for i = 1:size(n_hat, 1)
            sum_2 = sum_2 + (p(i,:)' - x(i,:)')' * n_hat(i,:)' * n_hat(i,:)'' * (p(i,:)' - x(i,:)');
        end
        ret = u' * (sum_0) * u - 2 * u' * (sum_1) + sum_2;
    end
    [~,optimize_0] = fminunc(@target_1,zeros(6,1));
    ret = optimize_0;
    output.ret = ret;
end
