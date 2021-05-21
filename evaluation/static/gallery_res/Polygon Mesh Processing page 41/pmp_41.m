function output = pmp_41(T)
% output = pmp_41(T)
%
%    `xᵢ` = T_*,1
%    `xⱼ` = T_*,2
%    `xₖ` = T_*,3
%    `n(T)` = (`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)/‖(`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)‖
%    
%    where
%     
%    T ∈ ℝ^(3×3)
    if nargin==0
        warning('generating random input data');
        [T] = generateRandomData();
    end
    function [T] = generateRandomData()
        T = randn(3, 3);
    end

    assert( isequal(size(T), [3, 3]) );

    x_i = T(:, 1);
    x_j = T(:, 2);
    x_k = T(:, 3);
    n_T = cross((x_j - x_i), (x_k - x_i)) / norm(cross((x_j - x_i), (x_k - x_i)), 2);
    output.x_i = x_i;

    output.x_j = x_j;

    output.x_k = x_k;

    output.n_T = n_T;
end
