function output = analytic_eigensystems_13(U, V)
% output = analytic_eigensystems_13(U, V)
%
%    `T₁` = 1/√2 U[0 0 0
%                  0 0 -1
%                  0 1 0]Vᵀ
%    
%    where 
%    
%    U ∈ ℝ^(3×3) 
%    V ∈ ℝ^(3×3)
    if nargin==0
        warning('generating random input data');
        [U, V] = generateRandomData();
    end
    function [U, V] = generateRandomData()
        U = randn(3, 3);
        V = randn(3, 3);
    end

    assert( isequal(size(U), [3, 3]) );
    assert( isequal(size(V), [3, 3]) );

    T1_0 = zeros(3, 3);
    T1_0(1,:) = [0, 0, 0];
    T1_0(2,:) = [0, 0, -1];
    T1_0(3,:) = [0, 1, 0];
    T1 = 1 / sqrt(2) * U * T1_0 * V';
    output.T1 = T1;
end
