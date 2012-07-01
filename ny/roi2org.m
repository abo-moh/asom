function [ x y ] = roi2org( x,y, resize_param )
%ROI2ORG Summary of this function goes here
%   Detailed explanation goes here
    adj_crop = resize_param{1};
    adj_scale = resize_param{2};
    roi_crop = resize_param{3};

    x = (x + roi_crop(1)) ./ adj_scale(1) + adj_crop(2,1);
    y = (y + roi_crop(2)) ./ adj_scale(2) + adj_crop(1,1);
end

