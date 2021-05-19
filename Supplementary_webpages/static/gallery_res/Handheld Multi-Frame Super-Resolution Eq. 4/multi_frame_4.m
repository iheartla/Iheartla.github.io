function output = multi_frame_4(k1, k2, e1, e2)
% output = multi_frame_4(`k₁`, `k₂`, `e₁`, `e₂`)
%
%    Ω = [`e₁` `e₂`][`k₁`   0
%                     0    `k₂`] [`e₁`ᵀ
%                                 `e₂`ᵀ]
%    
%    where
%    `k₁` ∈ ℝ  
%    `k₂` ∈ ℝ 
%    `e₁` ∈ ℝ²
%    `e₂` ∈ ℝ²
    if nargin==0
        warning('generating random input data');
        [k1, k2, e1, e2] = generateRandomData();
    end
    function [k1, k2, e1, e2] = generateRandomData()
        k1 = randn();
        k2 = randn();
        e1 = randn(2,1);
        e2 = randn(2,1);
    end

    e1 = reshape(e1,[],1);
    e2 = reshape(e2,[],1);

    assert(numel(k1) == 1);
    assert(numel(k2) == 1);
    assert( numel(e1) == 2 );
    assert( numel(e2) == 2 );

    Omega_0 = [[reshape(e1, [2, 1]), reshape(e2, [2, 1])]];
    Omega_1 = zeros(2, 2);
    Omega_1(1,:) = [k1, 0];
    Omega_1(2,:) = [0, k2];
    Omega_2 = [[e1']; [e2']];
    Omega = Omega_0 * Omega_1 * Omega_2;
    output.Omega = Omega;
end
