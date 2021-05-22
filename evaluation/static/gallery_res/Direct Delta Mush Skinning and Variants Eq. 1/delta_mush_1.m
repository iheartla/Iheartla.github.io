function output = delta_mush_1(w, M, u)
% output = delta_mush_1(w, M, u)
%
%    v_i = ∑_j w_i,j M_j u_i
%    
%    where
%    
%    w ∈ ℝ^(n×m)
%    M_j ∈ ℝ^(4×4)
%    u_i ∈ ℝ^4
    if nargin==0
        warning('generating random input data');
        [w, M, u] = generateRandomData();
    end
    function [w, M, u] = generateRandomData()
        n = randi(10);
        dim_1 = n;
        m = randi(10);
        dim_0 = m;
        w = randn(n, m);
        M = randn(dim_0,4,4);
        u = randn(dim_1,4);
    end

    n = size(w, 1);
    m = size(w, 2);
    dim_0 = size(M, 1);
    dim_1 = size(u, 1);
    assert( isequal(size(w), [n, m]) );
    assert( isequal(size(M), [dim_0, 4, 4]) );
    assert( isequal(size(u), [dim_1, 4]) );
    assert( dim_0 == m );
    assert( dim_1 == n );

    v = zeros(dim_1, 4);
    for i = 1:dim_1
        sum_0 = zeros(4,1);
        for j = 1:size(M, 1)
            sum_0 = sum_0 + w(i, j) * squeeze(M(j,:,:)) * u(i,:)';
        end
        v(i,:) = (sum_0)';
    end
    output.v = v;
end
