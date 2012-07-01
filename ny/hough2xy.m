function [ a b ] = hough2xy( rho, theta, varargin)
%HOUGH2XY Calculate (x,y)-slope a and -intercept b of hough point
    % Check for optional arguments
    resize = false;
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'resize'
                resize_param = varargin{i+1};
                resize = true;
            otherwise
        end
    end

    % Calculate (x0,y0) defined by (theta,rho)
    x0 = cosd(theta) .* rho;
    y0 = sind(theta) .* rho;
    
    if resize
        % Resize back from roi to original image
        [x_org y_org] = roi2org(x0,y0, {[0 0; 0 0], resize_param{2},resize_param{3},resize_param{4}});
    else
        x_org = x0;
        y_org = y0;
    end
    
    % Calculate slope a and start value b.
    a = - x_org ./ y_org;
    b = (y_org - a .* x_org);
    if resize
        b = b + a * resize_param{1}(2,1) + resize_param{1}(1,1);
    end
end