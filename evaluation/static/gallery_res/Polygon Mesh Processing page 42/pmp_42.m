function output = pmp_42(alpha, n, v, N1)
% output = pmp_42(α, n, v, `N₁`)
%
%    given
%    α_T : ℝ
%    n_T : ℝ³
%    
%    `n(v)` = ( ∑_(T for T ∈ `N₁`_v) α_T n_T )/‖ ∑_(T for T ∈ `N₁`_v) α_T n_T ‖
%    
%    where
%    v ∈ ℤ
%    `N₁`_i ∈ {ℤ}
    if nargin==0
        warning('generating random input data');
        [alpha, n, v, N1] = generateRandomData();
    end
    function [alpha, n, v, N1] = generateRandomData()
        v = randi(10);
        dim_0 = randi(10);
        dim_1 = randi(10);
        alpha = randn(dim_0,1);
        n = randn(dim_0,3);
        N1 = {};
        for i = 1:dim_1
            N1_tmp = [];
            dim_4 = randi(10);
            for j = 1:dim_4 
                N1_tmp = [N1_tmp;randi(10)];
            end
            N1 = [N1, N1_tmp];
        end
    end

    alpha = reshape(alpha,[],1);

    dim_0 = size(alpha, 1);
    dim_1 = size(N1, 1);
    assert( size(alpha,1) == dim_0 );
    assert( isequal(size(n), [dim_0, 3]) );
    assert(numel(v) == 1);

    sum_0 = zeros(3,1);
    for T = 1:size(alpha, 1)
        if ismember([T],N1{v},'rows')
          sum_0 = sum_0 + alpha(T) * n(T,:)';
        end
    end
    sum_1 = zeros(3,1);
    for T = 1:size(alpha, 1)
        if ismember([T],N1{v},'rows')
          sum_1 = sum_1 + alpha(T) * n(T,:)';
        end
    end
    n_v = (sum_0) / norm(sum_1, 2);
    output.n_v = n_v;
end
