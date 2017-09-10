function [ret] = diffusionF(x)
    % 0.04 - 0.044
    ret = 0.04;
    rsq = ((x(1,:) - 0).^2 + (x(2,:) + 0.5).^2);
    dec = 12.0;
    ret = ret + 0.004 * exp(-dec * rsq);

    rsq = (x(1,:) + 0.5).^2 + (x(2,:) - 0).^2;
    ret = ret + 0.004 * exp(-dec * rsq);
    
    rsq = (x(1,:) + 0.2).^2 + (x(2,:) - 0.5).^2;
    ret = ret + 0.004 * exp(-dec * rsq);
    
    rsq =(x(1,:) - 0.5).^2 + (x(2,:) + 0.2).^2;
    ret = ret + 0.004 * exp(-dec * rsq);
    
end

