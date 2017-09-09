function [ret] = diffusionF(x)
    % 0.02 - 0.022
    r = sqrt((x(1,:) - 0.2).^2 + (x(2,:) - 0.2).^2);
    ret = 0.02* (1 + 0.1 * cos(r * pi/0.8).^3.* (r < 0.4)).^2 ;

end

