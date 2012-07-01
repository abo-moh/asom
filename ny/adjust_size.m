function [ img_out, flip, crop, scale ] = adjust_size( img, resize)
%ADJUST_SIZE Summary of this function goes here
%   Detailed explanation goes here
% Default values for optional arguments

img_out = img;

% Compute middle of mass to determine the orientation of the breast
[mx my] = center_of_mass(img_out);

% Flip image vertically if middle of mass is to the right of the image
flip = mx > size(img_out,2)/2;
if flip
    img_out = flipdim(img_out,2);
end

% Crop [y_left y_right;x_left x_right]
crop = [26 26;4 99];

% Crop away rows and columns
width = size(img_out,2) - sum(crop(2,:));
height = size(img_out,1) - sum(crop(1,:));
img_out = imcrop(img_out, [crop(2,1) crop(1,1) width height]);

% Create mask
mask = true(size(img_out));

% Set the first 25 rows to 0
mask(1:25,:) = false;

% Set the 100 right-most columns to 0
mask(:,end-99:end) = false;

% Apply mask to output image
%img_out(~mask) = 0;

% Resize image
if resize
    img_out = imresize(img_out, [NaN 128]);
end
scale_cols = size(img_out,1) / height;
scale_rows = size(img_out,2) / width;
% Scale
scale = [scale_rows scale_cols];
end