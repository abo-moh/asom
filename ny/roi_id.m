function [ roi, mx, my, crop ] = roi_id( img, canny )
%ROI_ID is used to find the region of interest
% Currently just static cropping:
    roi = img;
    
    % Set upper-left cropping-values [columns rows]
    crop = [0 0];
    if canny
        crop(1) = 0;
    end

    % Define lower threshold using Otsu's method
    % and use on img to create a mask
    threshold = graythresh(img);
    mask = im2bw(img,threshold);

    % Apply lower threshold mask
    img_thresh = mask & img;
   
    % Calculate center of mass (lower right-most point in roi)
    [mx my] = center_of_mass(img_thresh);
    
    % Crop image to upper-left part compared to center of mass
    roi = imcrop(roi, [(crop(1)+1) (crop(2)+1) (mx-crop(1)) (my-crop(1))]);

end

