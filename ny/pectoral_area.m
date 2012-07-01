function [ rho, theta, A ] = pectoral_area( rho, theta, resize_param )
%PECTORAL_AREA Summary of this function goes here
%   Detailed explanation goes here

    width = resize_param{4}(2);
    height = resize_param{4}(1);
    
    % Calculate end points in the original image
    % points = [left_x, left_y, right_x, right_y]
    points = hough2points(theta,...
                            rho,...
                            resize_param{4},...
                            'resize', resize_param);
    left_x = points(:,1);
    left_y = points(:,2);
    right_x = points(:,3);
    right_y = points(:,4);
    
    %figure, patch([1 1 width width],[1 height height 1],4);
    %hold on, axis ij;
    %plot([left_x, right_x], [left_y, right_y]);
    %plot(1 + cosd(theta).*rho,1 + sind(theta).*rho, 'Xr');
    
    
    % Calculate area of triangle:
    A_triangle = 1/2 * abs(right_x - left_x) .* abs(right_y - left_y);
    
    % If left_x > 1 there exists a rectangle to the
    % left of the triangle that we have to calculate.
    % The area above the right endpoint is not included here:
    mask = left_x > 1;
    A_left = mask .* ((left_x - 1) .* (left_y - right_y));
    
    % If right_y > 1 there exists a rectangle above the
    % triangle as well that we have to calculate.
    % This includes the potential area above A_left:
    mask = right_y > 1;
    A_above = mask .* (right_x - 1) .* (right_y - 1);
    
    % Sum up all intermediate area calculations
    A = A_triangle + A_left + A_above;
end

