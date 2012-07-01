function [ mx my ] = center_of_mass ( img )
% CENTER_OF_MASS Compute center of mass of image

    % Get size of the image
    [rows, cols] = size(img);
    
    % Create matrices to weight every element with their position
    % in the x and y directions
    x = ones(rows, 1) * [1 : cols]; % Columns
    y = [1 : rows]' * ones(1, cols); % Rows

    % Calculate total sum of intensities
    area = sum(sum(img));
    
    % Calculate center of mass in the x and y direction
    mx = sum(sum(double(img).* x)) / area;
    my = sum(sum(double(img).* y)) / area;
end