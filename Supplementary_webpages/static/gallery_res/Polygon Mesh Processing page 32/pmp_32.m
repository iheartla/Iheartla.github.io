function output = pmp_32(phi, theta, R)
% output = pmp_32(ϕ, θ, R)
%
%    from trigonometry: sin, cos
%    `x(θ, ϕ)` = [Rcos(θ)cos(ϕ)
%                 Rsin(θ)cos(ϕ)
%                 Rsin(ϕ)]
%    
%    where
%    
%    ϕ ∈ ℝ : angle between 0 and 2π
%    θ ∈ ℝ : angle between -π/2 and π/2
%    R ∈ ℝ : the radius of the sphere
    if nargin==0
        warning('generating random input data');
        [phi, theta, R] = generateRandomData();
    end
    function [phi, theta, R] = generateRandomData()
        phi = randn();
        theta = randn();
        R = randn();
    end

    assert(numel(phi) == 1);
    assert(numel(theta) == 1);
    assert(numel(R) == 1);

    x_theta_phi_0 = zeros(3, 1);
    x_theta_phi_0(1,:) = [R * cos(theta) * cos(phi)];
    x_theta_phi_0(2,:) = [R * sin(theta) * cos(phi)];
    x_theta_phi_0(3,:) = [R * sin(phi)];
    x_theta_phi = x_theta_phi_0;
    output.x_theta_phi = x_theta_phi;
end
