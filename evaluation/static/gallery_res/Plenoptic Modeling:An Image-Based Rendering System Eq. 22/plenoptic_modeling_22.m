function output = plenoptic_modeling_22(v_bar, o_bar, u_bar, V_bar, C_bar_a, theta, v, D_A, delta)
% output = plenoptic_modeling_22(v̄, ō, ū, V̄, `C̄ₐ`, θ, v, `D_A`, δ)
%
%    r̄ = v̄×ō
%    s̄ = ō×ū
%    n̄ = ū×v̄
%    
%    `kᵣ` = r̄⋅(`C̄ₐ`-V̄)
%    `kₛ` = s̄⋅(`C̄ₐ`-V̄)
%    `kₙ` = n̄⋅(`C̄ₐ`-V̄)
%    
%    `x(θ,v)` =  (r̄⋅`D_A`(θ, v)+`kᵣ`δ(θ, v))/(n̄⋅`D_A`(θ, v)+`kₙ`δ(θ, v))
%    `y(θ,v)` =  (s̄⋅`D_A`(θ, v)+`kₛ`δ(θ, v))/(n̄⋅`D_A`(θ, v)+`kₙ`δ(θ, v))
%    
%    where
%    
%    v̄ ∈ ℝ^3
%    ō ∈ ℝ^3
%    ū ∈ ℝ^3
%    V̄ ∈ ℝ^3
%    `C̄ₐ` ∈ ℝ^3
%    θ ∈ ℝ 
%    v ∈ ℝ 
%    `D_A`: ℝ,ℝ → ℝ^3
%    δ: ℝ,ℝ → ℝ 
    if nargin==0
        warning('generating random input data');
        [v_bar, o_bar, u_bar, V_bar, C_bar_a, theta, v, D_A, delta] = generateRandomData();
    end
    function [v_bar, o_bar, u_bar, V_bar, C_bar_a, theta, v, D_A, delta] = generateRandomData()
        theta = randn();
        v = randn();
        v_bar = randn(3,1);
        o_bar = randn(3,1);
        u_bar = randn(3,1);
        V_bar = randn(3,1);
        C_bar_a = randn(3,1);
        D_A = @D_AFunc;
        rseed = randi(2^32);
        function tmp =  D_AFunc(p0, p1)
            rng(rseed);
            tmp = randn(3,1);
        end

        delta = @deltaFunc;
        rseed = randi(2^32);
        function tmp =  deltaFunc(p0, p1)
            rng(rseed);
            tmp = randn();
        end

    end

    v_bar = reshape(v_bar,[],1);
    o_bar = reshape(o_bar,[],1);
    u_bar = reshape(u_bar,[],1);
    V_bar = reshape(V_bar,[],1);
    C_bar_a = reshape(C_bar_a,[],1);

    assert( numel(v_bar) == 3 );
    assert( numel(o_bar) == 3 );
    assert( numel(u_bar) == 3 );
    assert( numel(V_bar) == 3 );
    assert( numel(C_bar_a) == 3 );
    assert(numel(theta) == 1);
    assert(numel(v) == 1);

    r_bar = cross(v_bar, o_bar);
    s_bar = cross(o_bar, u_bar);
    n_bar = cross(u_bar, v_bar);
    k_r = dot(r_bar,(C_bar_a - V_bar));
    k_s = dot(s_bar,(C_bar_a - V_bar));
    k_n = dot(n_bar,(C_bar_a - V_bar));
    x_theta_v = (dot(r_bar,D_A(theta, v)) + k_r * delta(theta, v)) / (dot(n_bar,D_A(theta, v)) + k_n * delta(theta, v));
    y_theta_v = (dot(s_bar,D_A(theta, v)) + k_s * delta(theta, v)) / (dot(n_bar,D_A(theta, v)) + k_n * delta(theta, v));
    output.r_bar = r_bar;

    output.s_bar = s_bar;

    output.n_bar = n_bar;

    output.k_r = k_r;

    output.k_s = k_s;

    output.k_n = k_n;

    output.x_theta_v = x_theta_v;

    output.y_theta_v = y_theta_v;
end
