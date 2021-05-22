function output = pmp_74(v, u, M, A)
% output = pmp_74(v, u, M, A)
%
%    `E_LSCM` = ∑_T A_T‖M_T v_T - [0 -1
%                                  1  0] M_T u_T‖²
%    where
%     
%    v_T ∈ ℝ^3
%    u_T ∈ ℝ^3
%    M_T ∈ ℝ^(2×3)
%    A_T ∈ ℝ
    if nargin==0
        warning('generating random input data');
        [v, u, M, A] = generateRandomData();
    end
    function [v, u, M, A] = generateRandomData()
        dim_0 = randi(10);
        v = randn(dim_0,3);
        u = randn(dim_0,3);
        M = randn(dim_0,2,3);
        A = randn(dim_0,1);
    end

    A = reshape(A,[],1);

    dim_0 = size(A, 1);
    assert( isequal(size(v), [dim_0, 3]) );
    assert( isequal(size(u), [dim_0, 3]) );
    assert( isequal(size(M), [dim_0, 2, 3]) );
    assert( size(A,1) == dim_0 );

    sum_0 = 0;
    for T = 1:size(M, 1)
        E_LSCM_0 = zeros(2, 2);
        E_LSCM_0(1,:) = [0, -1];
        E_LSCM_0(2,:) = [1, 0];
        sum_0 = sum_0 + A(T) * norm(squeeze(M(T,:,:)) * v(T,:)' - E_LSCM_0 * squeeze(M(T,:,:)) * u(T,:)', 2).^2;
    end
    E_LSCM = sum_0;
    output.E_LSCM = E_LSCM;
end
