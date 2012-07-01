function [ points ] = hough2points( theta, rho, img_size, varargin )
%HOUGH2LINE Calculate end points of line defined by (theta, rho)
%   Converts from the definition theta, rho to ax + b and calculates
%   endpoints in a image of the given image size img_size
    resize = false;
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'resize'
                resize_param = varargin{i+1};
                resize = true;
            otherwise
        end
    end

    % Calculate line parameters a and b
    % If resize is set the slope a and intercept b should
    % bet set according to the resize parameters
    if resize
        [a b] = hough2xy(rho, theta, 'resize', resize_param);
    else
        [a b] = hough2xy(rho, theta);
    end
    
    % Extract width and height
    width = img_size(2);
    height = img_size(1);
    
    %% Calculate points in a coordinate system with orego (0,0)
    % Calculate width and height difference from matrix orego (1,1)
    h_diff = height - 1;
    w_diff = width - 1;
    
    % Compute endpoints as (x y):
    left_x = 0;
    left_y = b;
    right_x = -b ./ a;
    right_y = 0;
    
    % Cutt off edges to fit image
    % Left endpoint
    mask = left_y > h_diff;
    left_x = (mask .* (h_diff - b) ./ a) + (~mask .* left_x);
    left_y(mask) = h_diff;
    
    % Right endpoint
    mask = right_x > w_diff;
    right_y = (mask .* (a* w_diff + b)) + (~mask .* right_y);
    right_x(mask) = w_diff;
    
    %% Switch orego to (1,1)
    left_x = left_x + 1;
    left_y = left_y + 1;
    right_x = right_x + 1;
    right_y = right_y + 1;
    
    % Put points in output as (l_x, l_y, r_x, r_y)
    points(:,1) = left_x;
    points(:,2) = left_y;
    points(:,3) = right_x;
    points(:,4) = right_y;
end

