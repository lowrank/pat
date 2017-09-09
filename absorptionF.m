function [ret] = absorptionF(x)
    % 0.1 - 0.2
    r = sqrt((x(1,:) - 0.2).^2 + (x(2,:) - 0.2).^2);
    ret = 0.2 + 0.05 * (1 + cos(r * 2 * pi/0.6)) .* (r < 0.3);
end

