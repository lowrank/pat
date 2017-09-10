function [ ret ] = gruneisenF(x)
    ret = 1.0;
    r = sqrt((x(1,:) - 0.2).^2 + (x(2,:) - 0.2).^2);
    ret = ret + 0.1 * (1 + cos(r * 2 * pi/0.6)) .* (r < 0.3);
    r = sqrt((x(1,:) - 0.3).^2 + (x(2,:) + 0.5).^2);
    ret = ret + 0.1 * (1 + cos(r * 2 * pi/0.6)) .* (r < 0.3);
    r = sqrt((x(1,:) + 0.5).^2 + (x(2,:) - 0.3).^2);
    ret = ret + 0.1 * (1 + cos(r * 2 * pi/0.6)) .* (r < 0.3);
%     ret = ones(size(x(1,:)));
end

