function [h] = pat_source()
    % ot source are geometry related. Now on a circle.
    K = 36;
    h = cell(K, 1);
    dec = 10.0;
    theta = 2 * pi / K;
    for d = 1 : K
        c = cos(theta * d);
        s = sin(theta * d);
        h{d} = @(x)(1e2 * exp(-dec * ((x(1, :) - c).^2 + (x(2, :) - s).^2)));
    end
end

