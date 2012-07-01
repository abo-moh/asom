function [ img_out ] = im_capframe( img, flip, crop )
%IM_CAPFRAME Summary of this function goes here
%   Detailed explanation goes here
    crop_cap = crop;
    img_out = img;
    if flip
        crop_cap(2,1) = crop(2,2);
        crop_cap(2,2) = crop(2,1);
    else
        crop_cap = crop{1};
    end
    
    if crop_cap(1,1) ~= 0
        img_out(1:crop_cap(1,1),:) = 0;
    end
    if crop_cap(1,2) ~= 0
        img_out(end-crop_cap(1,2)+1:end,:) = 0;
    end
    if crop_cap(2,1) ~= 0
        img_out(:,1:crop_cap(2,1)) = 0;
    end
    if crop_cap(2,2) ~= 0
        img_out(:,end-crop_cap(2,2)+1:end) = 0;
    end

end

