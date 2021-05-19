function output = anisotropic_elasticity_47(D_m, v)
% output = anisotropic_elasticity_47(`Dₘ`, v)
%
%    from linearalgebra: tr
%    
%    `J₃` = 1₃,₃
%    `κ_angle(Dₘ)` = 3(√2 v)^(2/3)(7/4‖`Dₘ`‖_F^2-1/4tr(`J₃``Dₘ`ᵀ`Dₘ`))⁻¹
%    
%    where
%    
%    `Dₘ` ∈ ℝ^(3×3)  
%    v ∈ ℝ
    if nargin==0
        warning('generating random input data');
        [D_m, v] = generateRandomData();
    end
    function [D_m, v] = generateRandomData()
        v = randn();
        D_m = randn(3, 3);
    end

    assert( isequal(size(D_m), [3, 3]) );
    assert(numel(v) == 1);

    J3 = ones(3, 3);
    kappa_angle_D_m = 3 * (sqrt(2) * v).^(2 / 3) * 1 / ((7 / 4 * norm(D_m, 'fro').^2 - 1 / 4 * trace(J3 * D_m' * D_m)));
    output.J3 = J3;

    output.kappa_angle_D_m = kappa_angle_D_m;
end
