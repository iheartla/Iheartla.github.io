function output = multi_frame_1(c, w, R_hat)
% output = multi_frame_1(c, w, R̂)
%
%    `C(x,y)` = (∑_n ∑_i c_n,i ⋅ w_n,i ⋅ R̂_n) / (∑_n ∑_i w_n,i ⋅ R̂_n)
%    
%    where
%    
%    c ∈ ℝ^(f×s) : the value of the Bayer pixel
%    w ∈ ℝ^(f×s) : the local sample weight
%    R̂ ∈ ℝ^f     : the local robustness
    if nargin==0
        warning('generating random input data');
        [c, w, R_hat] = generateRandomData();
    end
    function [c, w, R_hat] = generateRandomData()
        f = randi(10);
        s = randi(10);
        c = randn(f, s);
        w = randn(f, s);
        R_hat = randn(f,1);
    end

    R_hat = reshape(R_hat,[],1);

    f = size(c, 1);
    s = size(c, 2);
    assert( isequal(size(c), [f, s]) );
    assert( isequal(size(w), [f, s]) );
    assert( numel(R_hat) == f );

    sum_0 = 0;
    for n = 1:size(R_hat,1)
        sum_1 = 0;
        for i = 1:size(w,2)
            sum_1 = sum_1 + c(n, i) * w(n, i) * R_hat(n);
        end
        sum_0 = sum_0 + sum_1;
    end
    sum_2 = 0;
    for n = 1:size(R_hat,1)
        sum_3 = 0;
        for i = 1:size(w,2)
            sum_3 = sum_3 + w(n, i) * R_hat(n);
        end
        sum_2 = sum_2 + sum_3;
    end
    C_x_y = (sum_0) / (sum_2);
    output.C_x_y = C_x_y;
end
