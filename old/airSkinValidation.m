function [value] = airSkinValidation(mask,im_crop)
%AIRSKINVALIDATION gives a measure of the air-skin segmentation.
% if value < 0.8, the air-skin segmentation is good.

if nargin < 2
    im_crop = [0 0 0 0];
end

% Crop image
mask = mask(im_crop(1)+1:end-im_crop(2),im_crop(3)+1:end-im_crop(4));

% Convert to logical
mask_BW = logical(mask);

% sum of circumference
O = sum(sum(bwperim(mask_BW, 8)));

% area
%A = sum(sum(mask_BW));
A = bwarea(mask_BW);

value = log(O^2/(4*pi*A));