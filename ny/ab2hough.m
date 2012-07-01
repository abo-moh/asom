function [ theta, rho ] = ab2hough( a, b )
%AB2HOUGH Convert y = ax + b to hough line coordinates
    theta = atand(-1 / a);
    x0 = -b / (a + (1/a));
    y0 = a*x0+b;
    rho = sqrt(x0^2+y0^2);
end

