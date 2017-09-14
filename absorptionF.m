function [ret] = absorptionF(x)
    % 0.1 - 0.2
    ret = 0.1;
    K = 5;
    dec = 18.0;
    theta = 2 * pi / K;
    for id = 1:K
        s = 0.6 * cos(id * theta);
        c = 0.6 * sin(id * theta);
        rsq = ((x(1,:) - s).^2 + (x(2,:) - c).^2);
        ret = ret + 0.1 * exp(-dec * rsq);
    end
end

