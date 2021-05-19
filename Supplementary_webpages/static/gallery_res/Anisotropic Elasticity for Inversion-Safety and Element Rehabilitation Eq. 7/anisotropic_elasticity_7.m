function output = anisotropic_elasticity_7(A)
% output = anisotropic_elasticity_7(A)
%
%    `∂²I₅/∂f²` = 2[A₁,₁I₃  A₁,₂I₃  A₁,₃I₃
%                   A₂,₁I₃  A₂,₂I₃  A₂,₃I₃
%                   A₃,₁I₃  A₃,₂I₃  A₃,₃I₃] 
%    
%    where
%    
%    A ∈ ℝ^(3×3)
    if nargin==0
        warning('generating random input data');
        [A] = generateRandomData();
    end
    function [A] = generateRandomData()
        A = randn(3, 3);
    end

    assert( isequal(size(A), [3, 3]) );

    partial_differential_2I5_solidus_partial_differential_f2_0 = [[A(1, 1) * speye(3), A(1, 2) * speye(3), A(1, 3) * speye(3)]; [A(2, 1) * speye(3), A(2, 2) * speye(3), A(2, 3) * speye(3)]; [A(3, 1) * speye(3), A(3, 2) * speye(3), A(3, 3) * speye(3)]];
    partial_differential_2I5_solidus_partial_differential_f2 = 2 * partial_differential_2I5_solidus_partial_differential_f2_0;
    output.partial_differential_2I5_solidus_partial_differential_f2 = partial_differential_2I5_solidus_partial_differential_f2;
end
