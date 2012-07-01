function [ NH ] = hough_normalize( H, theta, rho, gmag_size, normalize )
%HOUGH_NORMALIZE Normalizes the Hough plane
%   The normalization is done by computing a matrix of
%   lengths of the lines in the 

    % Define bounds that the pectoral line has to lie within
    theta_bound = [10 45];
    rho_min = 0;
    
    % Filter out points in the Houghplane
    % that aren't within the bounds
    rho_values = rho > rho_min;
    theta_values = theta > theta_bound(1) & theta < theta_bound(2);
    theta_mask = ones(numel(rho),1) * theta_values;
    rho_mask = rho_values' * ones(1,numel(theta));
    filter_mask = theta_mask & rho_mask;
    
    NH = filter_mask .* H;
    
    if normalize
        % Compute length of each line not filtered:
        ind = find(filter_mask);
        [row,col] = ind2sub(size(filter_mask), ind);

        % Compute end points of each line in the image
        points = hough2points(theta(col), rho(row), gmag_size);

        % Compute distance between endpoints
        length = sqrt((points(:,1) - points(:,3)).^2 + (points(:,2) - points(:,4)) .^2);

        % Multiply filtered values by the reciproc of their lengths
        % if length(rho,theta) > height/10 and otherwise multiply
        % by 10 / height of gradient image
        mask = length > (gmag_size(1)/10);
        NF = mask .* (1 ./ sqrt(length)) + ~mask * 10/gmag_size(1);

        NH(ind) = NH(ind) .* NF;
    end
end

