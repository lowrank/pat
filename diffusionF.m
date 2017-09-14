function [ret] = diffusionF(x)
    % 0.04 - 0.044
    ret = 0.04;
    K = 4;
    dec = 18.0;
    theta = 2 * pi / K;
    phi = rand();
    for id = 1:K
        s = 0.6 * cos(id * theta + phi);
        c = 0.6 * sin(id * theta + phi);
        rsq = ((x(1,:) - s).^2 + (x(2,:) - c).^2);
        ret = ret + 0.004 * exp(-dec * rsq);
    end

end

