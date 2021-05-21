function output = course_curvature(p, k_n)
% output = course_curvature(p, `kₙ`)
%
%    `H(p)` = 1/(2π)∫_[0, 2π] `kₙ`(φ, p) ∂φ
%    
%    where 
%    
%    p ∈ ℝ^3 : point on the surface
%    `kₙ`: ℝ,ℝ^3 → ℝ : normal curvature
    if nargin==0
        warning('generating random input data');
        [p, k_n] = generateRandomData();
    end
    function [p, k_n] = generateRandomData()
        p = randn(3,1);
        k_n = @k_nFunc;
        rseed = randi(2^32);
        function tmp =  k_nFunc(p0, p1)
            rng(rseed);
            tmp = randn();
        end

    end

    p = reshape(p,[],1);

    assert( numel(p) == 3 );

    H_p = 1 / (2 * pi) * integral(@(phi) k_n(phi, p), 0, 2 * pi,'ArrayValued',true);
    output.H_p = H_p;
end
