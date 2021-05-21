function output = course_parameterization(w, E)
% output = course_parameterization(w, E)
%
%    L_i,j = { w_i,j if (i,j) ∈ E
%              0 otherwise
%    L_i,i = -∑_(ℓ for ℓ ≠ i) L_i,ℓ
%    
%    where
%    L ∈ ℝ^(n×n)
%    w ∈ ℝ^(n×n): edge weight matrix
%    E ∈ {ℤ²} index: edges
    if nargin==0
        warning('generating random input data');
        [w, E] = generateRandomData();
    end
    function [w, E] = generateRandomData()
        n = randi(10);
        w = randn(n, n);
        E = [];
        dim_2 = randi(10);
        for i = 1:dim_2 
            E = [E;randi(10), randi(10)];
        end
    end

    n = size(w, 2);
    assert( isequal(size(w), [n, n]) );
    assert(size(E,2) == 2)

    Lij_0 = zeros(2,0);
    Lvals_0 = zeros(1,0);
    for i = 1:n
        for j = 1:n
            if ismember([i, j],E,'rows')
                Lij_0(1:2,end+1) = [i;j];
                Lvals_0(end+1) = w(i, j);
            end
        end
    end
    sparse_0 = sparse(Lij_0(1,:),Lij_0(2,:),Lvals_0,n,n);
    L = sparse_0;
    for i = 1:n
        sum_0 = 0;
        for ell = 1:size(L,2)
            if ell ~= i
              sum_0 = sum_0 + L(i, ell);
            end
        end
        Lij_0(1:2,end+1) = [i;i];
        Lvals_0(end+1) = -sum_0;
    end
    L = sparse(Lij_0(1,:),Lij_0(2,:),Lvals_0,n,n);
    output.L = L;
end
