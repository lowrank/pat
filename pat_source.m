function [h] = pat_source()
    K = 12;
    h = cell(K, 1);
    dec = 2.0;
    theta = 2 * pi / K;
    for d = 1 : K
        c = cos(theta * d);
        s = sin(theta * d);
        h{d} = @(x)(exp(dec * (c * x(:, 1) + s * x(:, 2))));
    end
end

