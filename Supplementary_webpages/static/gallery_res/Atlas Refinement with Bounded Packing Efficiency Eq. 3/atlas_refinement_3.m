function output = atlas_refinement_3(l, dist, b_i, b, sigma, s_k)
% output = atlas_refinement_3(l, dist, `bᵢ`, b, σ, `s^k`)
%
%    `G_σ(s^k_i)` = ∑_j l_j exp(-dist(`bᵢ`, b_j)²/(2σ²)) `s^k`_j
%    
%    where
%    l_j ∈ ℝ : the length of bj
%    dist: ℝ², ℝ² → ℝ : measures the geodesic distance 
%    `bᵢ` ∈ ℝ²
%    b_j ∈ ℝ²
%    σ ∈ ℝ
%    `s^k`_j ∈ ℝ² : direction vector
    if nargin==0
        warning('generating random input data');
        [l, dist, b_i, b, sigma, s_k] = generateRandomData();
    end
    function [l, dist, b_i, b, sigma, s_k] = generateRandomData()
        sigma = randn();
        dim_0 = randi(10);
        l = randn(dim_0,1);
        dist = @distFunc;
        rseed = randi(2^32);
        function tmp =  distFunc(p0, p1)
            rng(rseed);
            tmp = randn();
        end

        b_i = randn(2,1);
        b = randn(dim_0,2);
        s_k = randn(dim_0,2);
    end

    l = reshape(l,[],1);
    b_i = reshape(b_i,[],1);

    dim_0 = size(l, 1);
    assert( size(l,1) == dim_0 );
    assert( numel(b_i) == 2 );
    assert( isequal(size(b), [dim_0, 2]) );
    assert(numel(sigma) == 1);
    assert( isequal(size(s_k), [dim_0, 2]) );

    sum_0 = zeros(2,1);
    for j = 1:size(size(b, 1),1)
        sum_0 = sum_0 + l(j) * exp(-dist(b_i, b(j,:)').^2 / (2 * sigma.^2)) * s_k(j,:)';
    end
    G_sigma_s_k_i = sum_0;
    output.G_sigma_s_k_i = G_sigma_s_k_i;
end
