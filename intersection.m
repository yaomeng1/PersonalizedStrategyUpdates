function [bcr]=intersection(x_array, y_array, y0)
%% Calculate the numerical critical ratio C* by computing the intersection between 
% the curve of numerical fixation probability (rho_c) and horizontal line
% (rho_c=1/n)
% Input: x_array (benefit-to-cost ratio, b/c)
% Input: y_array (numerical fixation probability, rho_c)
% Input: y0 (y0=1/n)
% Output: intersection point on x axis (which is the numerical critical ratio C*)
    x_array = double(x_array);
    for i =1:length(x_array)
        if y_array(i) <= y0 && y_array(i + 1) >= y0
            deltax = x_array(i + 1) - x_array(i);
            deltay = y_array(i + 1) - y_array(i);
            bcr = x_array(i) + deltax * (y0 - y_array(i)) / deltay;
        end   
    end

end
