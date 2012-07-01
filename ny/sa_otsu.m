function breast = sa_otsu( img, flip, resize_param )
%SA_OTSU Perform Otsu's method for optimal image
%   Detailed explanation goes here
    % Premodify image
    % Gauss filter
    FG = fspecial('gaussian', [9 9],1.5);
    img_gauss = imfilter(img, FG);

    % Initialize variables by performing otsu's method once:
    threshold = graythresh(img_gauss);
    threshold_t = threshold;
    img_current = im2bw(img_gauss,threshold);
    current_bg = ~img_current;
    img_t = img_current;

    % Perform Otsu's method until the difference between current background
    % and next background gets too big or next background equals 0
    while sum(~img_t(:)) > 0.4 * sum(current_bg(:)) && sum(~img_t(:)) > 0
        % Update current masks and threshold
        threshold = threshold_t;
        img_current = img_current | img_t;
        current_bg = ~img_current;
    
        % Find next threshold
        threshold_t = graythresh(img_gauss .* current_bg);
        img_t = im2bw(img_gauss, threshold_t);
    end 
    % Use found foreground as image, cap "cropped" rows and cols to 0
    breast = img_current;
    breast = im_capframe(breast, flip, resize_param{1});
    
    % Structure pixels in regions
    comps = bwconncomp(breast,8);
    
    % Get area of each region
    areas = regionprops(comps, 'Area');
    
    % Get size of biggest region = size of breast
    breast_area = max(cell2mat(struct2cell(areas)));
    
    % Remove all areas except the biggest
    breast = bwareaopen(breast,breast_area);
    
    % Remove 1s surrounded by 0s
    breast = bwmorph(breast,'clean');
    
    % Set pixels to 1 if at least 5 neighbours are 1
    breast = bwmorph(breast,'majority');
    
    % Perform erosion followed by dilation using a disk of size 10
    SE = strel('disk',10);
    breast = imerode(breast,SE);
    breast = imdilate(breast,SE);
    
    % Fill all holes = points that can't be reach from the background
    breast = imfill(breast,'holes');
    
    % Perform removal of small regions again of previous modifications gave
    % arise to new regions
    comps = bwconncomp(breast,8);
    areas = regionprops(comps, 'Area');
    breast_area = max(cell2mat(struct2cell(areas)));
    breast = bwareaopen(breast,breast_area);
    
    % Cap cropped-away cols and rows to 0 again after having dilated etc.
    breast = im_capframe(breast, flip, resize_param{1});
    
    %figure, imshow(img_org_crop);
    %figure, imshow(img_org_crop .* img_org, [min(img_org(:)) max(img_org(:))]);
    
    %background = ~sa_grow_breast(~background, [2,2]);
    
    %{
    figure;
    subplot(2,3,1), imshow(img, [min(img(:)) max(img(:))]);
    subplot(2,3,2), imshow(~background);
    subplot(2,3,3), imshow(img .* ~background, [min(img(:)) max(img(:))]);   
    subplot(2,3,4), imshow(img_org, [min(img_org(:)) max(img_org(:))]);
    subplot(2,3,5), imshow(img_org_crop);
    subplot(2,3,6), imshow(img_org .* img_org_crop, [min(img_org(:)) max(img_org(:))]);
    %}
end

