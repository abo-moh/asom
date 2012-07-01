function [ H, theta, rho ] = hough_transform( gmag )
%HOUGH_TRANSFORM Summary of this function goes here
%   Detailed explanation goes here

% 702
% Calculate area of entire histogram of gradient magnitudes, A_hist
A_hist = histc(gmag(:),[0 Inf]);
A_hist = A_hist(1);

% 704: Create function Wf to use gmag and A_hist as default parameters
% The function Wf(x,y) will return the weight of g(x,y)
    function weight = Wf(x,y)
        weight = gradientWeight(gmag, A_hist, x, y);
    end

% 706
[height, width] = size(gmag)    ;
% Calculate diagonal of gradient-image
% The diagonal is the maximal value of rho:
diagonal = norm([width height]);

% Set discrete values of rho:
rho = 0:diagonal; % Discrete values of theta
H_height = numel(rho);

% Set discrete values of theta:
%{
thetaStep = pi/height; % Theta step value
theta = -pi/2:thetaStep:pi/2; % Discrete values of theta
theta_count = numel(theta); % Number of theta values
%}
%{
thetaStep = 180 / height;
theta = -90:thetaStep:90;
theta_count = numel(theta);
%}
thetaStep = 90 / height;
theta = 0:thetaStep:90;
theta_count = numel(theta);

% Allocate hough space:
H = zeros(H_height, theta_count);

% Find coordinates in magnitude plane where g(x,y) > 0.
% These are the specific pixels, that we have to map to the
% Hough plane
[mag_map_y, mag_map_x] = find(gmag);

% Number of pixels that we have to map:
mag_map_count = numel(mag_map_x);

% Allocate accumulator array. One row of values for each pixel we
% have to map. Each row has theta_count columns; one for each
% discrete value of angle theta.
acc = zeros(mag_map_count, theta_count);


% Preallocate cosine and sine arrays for calculating r later on.
% Transpose to row vectors to match mag_map_x and mag_map_y.
cos_theta = (0:width-1)' * cosd(theta);
sin_theta = (0:height-1)' * sind(theta);

% Calculate r values in accumulator cells for each (x,y):
acc(1:mag_map_count,:) = cos_theta(mag_map_x,:) + sin_theta(mag_map_y,:);

% Create vector of weights for each pixel
accW = Wf(mag_map_x, mag_map_y);

for i = 1:theta_count
    for j = 1:H_height
        % Find pixels to put into the rho bin
        [y,x] = find(acc(:,i) > (rho(j) - 0.5) & acc(:,i) < (rho(j) + 0.5));
        % Weigh each pixel and put it into the Hough plane
        if numel(y) > 0
            H(j,i) = sum(accW(y));
        end
    end
end

end