function output = convex_optimization_650(A, B, C)
% output = convex_optimization_650(A, B, C)
%
%    given
%    A ∈ ℝ^(k×k) 
%    B ∈ ℝ^(k×m) 
%    C ∈ ℝ^(m×m) 
%    
%    S = C - BᵀA⁻¹B
%    [A⁻¹+A⁻¹BS⁻¹BᵀA⁻¹   -A⁻¹BS⁻¹
%       -S⁻¹BᵀA⁻¹           S⁻¹]
    if nargin==0
        warning('generating random input data');
        [A, B, C] = generateRandomData();
    end
    function [A, B, C] = generateRandomData()
        k = randi(10);
        m = randi(10);
        A = randn(k, k);
        B = randn(k, m);
        C = randn(m, m);
    end

    k = size(A, 2);
    m = size(B, 2);
    assert( isequal(size(A), [k, k]) );
    assert( isequal(size(B), [k, m]) );
    assert( isequal(size(C), [m, m]) );

    S = C - B' * (A\B);
    ret_0 = [[inv(A) + (A\B) * (S\B') * inv(A), -(A\B) * inv(S)]; [-(S\B') * inv(A), inv(S)]];
    ret = ret_0;
    output.S = S;

    output.ret = ret;
end
