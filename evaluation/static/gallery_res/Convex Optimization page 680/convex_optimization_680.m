function output = convex_optimization_680(P1, P2, P3, B, C, L, L_tilde, U, U_tilde)
% output = convex_optimization_680(`P₁`, `P₂`, `P₃`, B, C, L, L̃, U, Ũ)
%
%    [`P₁`  0
%      0   `P₃`][      L          0
%                `P₃`ᵀC`P₂`ᵀU⁻¹  -L̃][U   L⁻¹`P₁`ᵀB
%                                    0       Ũ    ][`P₂`  0
%                                                    0   I_n]
%    
%    where
%    
%    `P₁` ∈ ℝ^(m×m) 
%    `P₂` ∈ ℝ^(m×m) 
%    `P₃` ∈ ℝ^(n×n) 
%      B ∈ ℝ^(m×n) 
%      C ∈ ℝ^(n×m) 
%      L ∈ ℝ^(m×m) 
%      L̃ ∈ ℝ^(n×n) 
%      U ∈ ℝ^(m×m) 
%      Ũ ∈ ℝ^(n×n)
    if nargin==0
        warning('generating random input data');
        [P1, P2, P3, B, C, L, L_tilde, U, U_tilde] = generateRandomData();
    end
    function [P1, P2, P3, B, C, L, L_tilde, U, U_tilde] = generateRandomData()
        m = randi(10);
        n = randi(10);
        P1 = randn(m, m);
        P2 = randn(m, m);
        P3 = randn(n, n);
        B = randn(m, n);
        C = randn(n, m);
        L = randn(m, m);
        L_tilde = randn(n, n);
        U = randn(m, m);
        U_tilde = randn(n, n);
    end

    m = size(P1, 2);
    n = size(P3, 2);
    assert( isequal(size(P1), [m, m]) );
    assert( isequal(size(P2), [m, m]) );
    assert( isequal(size(P3), [n, n]) );
    assert( isequal(size(B), [m, n]) );
    assert( isequal(size(C), [n, m]) );
    assert( isequal(size(L), [m, m]) );
    assert( isequal(size(L_tilde), [n, n]) );
    assert( isequal(size(U), [m, m]) );
    assert( isequal(size(U_tilde), [n, n]) );

    ret_0 = [[P1, zeros(m, n)]; [zeros(n, m), P3]];
    ret_1 = [[L, zeros(m, n)]; [P3' * C * P2' * inv(U), -L_tilde]];
    ret_2 = [[U, (L\P1') * B]; [zeros(n, m), U_tilde]];
    ret_3 = [[P2, zeros(m, n)]; [zeros(n, m), speye(n)]];
    ret = ret_0 * ret_1 * ret_2 * ret_3;
    output.ret = ret;
end
