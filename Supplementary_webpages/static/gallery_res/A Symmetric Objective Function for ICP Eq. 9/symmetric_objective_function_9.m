function output = symmetric_objective_function_9(a, theta, p, q, n, t)
% output = symmetric_objective_function_9(a, θ, p, q, n, t)
%
%    from trigonometry: tan, cos
%    
%    t̃ = t/cos(θ)
%    ã = a tan(θ)
%    ∑_i cos²(θ)((p_i - q_i)⋅n_i +((p_i+q_i)×n_i)⋅ã+n_i⋅t̃)² 
%    
%    where
%    a ∈ ℝ³ : axis of rotation
%    θ ∈ ℝ  : angle of rotation
%    p_i ∈ ℝ³
%    q_i ∈ ℝ³
%    n_i ∈ ℝ³
%    t ∈ ℝ³
    if nargin==0
        warning('generating random input data');
        [a, theta, p, q, n, t] = generateRandomData();
    end
    function [a, theta, p, q, n, t] = generateRandomData()
        theta = randn();
        dim_0 = randi(10);
        a = randn(3,1);
        p = randn(dim_0,3);
        q = randn(dim_0,3);
        n = randn(dim_0,3);
        t = randn(3,1);
    end

    a = reshape(a,[],1);
    t = reshape(t,[],1);

    dim_0 = size(p, 1);
    assert( numel(a) == 3 );
    assert(numel(theta) == 1);
    assert( isequal(size(p), [dim_0, 3]) );
    assert( isequal(size(q), [dim_0, 3]) );
    assert( isequal(size(n), [dim_0, 3]) );
    assert( numel(t) == 3 );

    t_tilde = t / cos(theta);
    a_tilde = a * tan(theta);
    sum_0 = 0;
    for i = 1:size(size(p, 1),1)
        sum_0 = sum_0 + cos(theta).^2 * (dot((p(i,:)' - q(i,:)'),n(i,:)') + dot((cross((p(i,:)' + q(i,:)'), n(i,:)')),a_tilde) + dot(n(i,:)',t_tilde)).^2;
    end
    ret = sum_0;
    output.t_tilde = t_tilde;

    output.a_tilde = a_tilde;

    output.ret = ret;
end
