function output = morphable_model_5(sigma_N, E_I, alpha, beta, sigma_S, sigma_T, rho, rho_bar, sigma_rho, latin_small_letter_a_with_macron)
% output = morphable_model_5(`σ_N`, `E_I`, α, β, `σ_S`, `σ_T`, ρ, ρ̄, `σ_ρ`, ā)
%
%    E = 1/`σ_N`²`E_I` + ∑_j α_j²/`σ_S`_j² + ∑_j β_j²/`σ_T`_j²  + ∑_j (ρ_j-ρ̄_j)²/`σ_ρ`_j²
%    
%    where
%    
%    `σ_N` ∈ ℝ 
%    `E_I` ∈ ℝ
%    α_i ∈ ℝ
%    β_i ∈ ℝ
%    `σ_S`_i ∈ ℝ 
%    `σ_T`_i ∈ ℝ 
%    ρ_j ∈ ℝ 
%    ρ̄_j ∈ ℝ 
%    `σ_ρ`_j ∈ ℝ 
%    ā_i ∈ ℝ 
    if nargin==0
        warning('generating random input data');
        [sigma_N, E_I, alpha, beta, sigma_S, sigma_T, rho, rho_bar, sigma_rho, latin_small_letter_a_with_macron] = generateRandomData();
    end
    function [sigma_N, E_I, alpha, beta, sigma_S, sigma_T, rho, rho_bar, sigma_rho, latin_small_letter_a_with_macron] = generateRandomData()
        sigma_N = randn();
        E_I = randn();
        dim_0 = randi(10);
        dim_1 = randi(10);
        alpha = randn(dim_0,1);
        beta = randn(dim_0,1);
        sigma_S = randn(dim_0,1);
        sigma_T = randn(dim_0,1);
        rho = randn(dim_1,1);
        rho_bar = randn(dim_1,1);
        sigma_rho = randn(dim_1,1);
        latin_small_letter_a_with_macron = randn(dim_0,1);
    end

    alpha = reshape(alpha,[],1);
    beta = reshape(beta,[],1);
    sigma_S = reshape(sigma_S,[],1);
    sigma_T = reshape(sigma_T,[],1);
    rho = reshape(rho,[],1);
    rho_bar = reshape(rho_bar,[],1);
    sigma_rho = reshape(sigma_rho,[],1);
    latin_small_letter_a_with_macron = reshape(latin_small_letter_a_with_macron,[],1);

    dim_0 = size(alpha, 1);
    dim_1 = size(rho, 1);
    assert(numel(sigma_N) == 1);
    assert(numel(E_I) == 1);
    assert( size(alpha,1) == dim_0 );
    assert( size(beta,1) == dim_0 );
    assert( size(sigma_S,1) == dim_0 );
    assert( size(sigma_T,1) == dim_0 );
    assert( size(rho,1) == dim_1 );
    assert( size(rho_bar,1) == dim_1 );
    assert( size(sigma_rho,1) == dim_1 );
    assert( size(latin_small_letter_a_with_macron,1) == dim_0 );

    sum_0 = 0;
    for j = 1:size(size(alpha, 1),1)
        sum_0 = sum_0 + alpha(j).^2 / sigma_S(j).^2;
    end
    sum_1 = 0;
    for j = 1:size(size(beta, 1),1)
        sum_1 = sum_1 + beta(j).^2 / sigma_T(j).^2;
    end
    sum_2 = 0;
    for j = 1:size(size(rho_bar, 1),1)
        sum_2 = sum_2 + (rho(j) - rho_bar(j)).^2 / sigma_rho(j).^2;
    end
    E = 1 / sigma_N.^2 * E_I + sum_0 + sum_1 + sum_2;
    output.E = E;
end
