function output = hand_modeling_3(x, R)
% output = hand_modeling_3(x, R)
%
%    min_(C ∈ ℝ^3) ∑_i ‖x_i + (R_i - I₃)C‖²
%    
%    where
%    
%    x_i ∈ ℝ^3
%    R_i ∈ ℝ^(3×3)
    if nargin==0
        warning('generating random input data');
        [x, R] = generateRandomData();
    end
    function [x, R] = generateRandomData()
        dim_0 = randi(10);
        x = randn(dim_0,3);
        R = randn(dim_0,3,3);
    end

    dim_0 = size(x, 1);
    assert( isequal(size(x), [dim_0, 3]) );
    assert( isequal(size(R), [dim_0, 3, 3]) );

    function ret = target_1(C)
        sum_0 = 0;
        for i = 1:size(R, 1)
            sum_0 = sum_0 + norm(x(i,:)' + (squeeze(R(i,:,:)) - speye(3)) * C, 2).^2;
        end
        ret = sum_0;
    end
    [~,optimize_0] = fminunc(@target_1,zeros(3,1));
    ret = optimize_0;
    output.ret = ret;
end
