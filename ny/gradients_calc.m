function [ gradients, g_directions ] = gradients_calc( img )
%GRADIENTS_CALC Calculate gradient magnitudes and directions
% Initial things to do
img = im2double(img);

% High and low bound
bound_high = 10;
bound_low = 5;

%% Preprocessing

% Perform gaussian filtering (by convoluting)
gauss = fspecial('gaussian', [3 3], 0.5);
img_gauss = imfilter(img, gauss, 'conv');
img_gauss = img;

%% Gradient magnitudes

% Get sobel 3x3 operator for finding gradients
sobel_operator = fspecial('sobel');
%disp(sobel_operator)
%disp(rot90(sobel_operator,1));

% Apply the sobel operator to get horizontal and vertical gradient
% magnitudes
gx = imfilter(img_gauss, rot90(sobel_operator), 'conv', 'replicate');
gy = imfilter(img_gauss, sobel_operator, 'conv', 'replicate');

% Compute combined magnitudes
gradients = sqrt(gx.^2 + gy.^2);

gradients(1:2,:) = 0; gradients(end-1:end,:) = 0;
gradients(:,1:2) = 0; gradients(:,end-1:end) = 0;

%% Gradient directions
g_directions = zeros(size(gradients));
g_d_temp = g_directions;
% Handle 4 cases:
for i = 1:4
    g_d_temp(:,:) = 0;
    if i == 1
        % gx, gy >= 0:
        mask = gx >= 0 & gy >= 0;
        gt_gx = gx .* mask;
        gt_gy = gy .* mask;
        g_d_temp = atan(abs(gt_gy) ./ abs(gt_gx));
        
    elseif i == 2
        % gx < 0, gy > 0:
        mask = gx < 0 & gy > 0;
        gt_gx = gx .* mask;
        gt_gy = gy .* mask;
        g_d_temp = atan(abs(gt_gx) ./ abs(gt_gy)) + pi/2;
        
    elseif i == 3
        % gx, gy <= 0;
        mask = gx <= 0 & gy <= 0;
        gt_gx = gx .* mask;
        gt_gy = gy .* mask;
        g_d_temp = atan(abs(gt_gy) ./ abs(gt_gx)) + pi;
        
    elseif i == 4
        % gx > 0, gy < 0:
        mask = gx > 0 & gy < 0;
        gt_gx = gx .* mask;
        gt_gy = gy .* mask;
        g_d_temp = atan(abs(gt_gx) ./ abs(gt_gy)) + 3/2*pi;
        
    end
    nans = isnan(g_d_temp);
    g_d_temp(nans) = 0;
    g_directions = g_directions + g_d_temp;
end


%% display the image gradient flow for debugging and examples
test = false;
if test
    figure;clf;imagesc(img);colormap(gray);axis image;
    hold on;
    quiver(gx, gy);
end

end

