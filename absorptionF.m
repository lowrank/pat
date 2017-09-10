function [ret] = absorptionF(x)
    % 0.1 - 0.2
    ret = 0.1;
    dec = 4.0;
    rsq = (x(1,:) - 0.2).^2 + (x(2,:) - 0.2).^2;
    ret = ret + 0.1 * exp(-dec * rsq);
end

