function [mask,pectoral,nipple] = breast_segmentation(image,airskin_iter,hough_det_sens,fig,im_crop,im_direction)
% BREAST_SEGMENTATION finds the pectoral muscle in the image and returns a 
% mask (0 background, 1 breast, 2 pectoral), two points on the pectoral
% line, and the nipple point.
%
% Input:
% airskin_iter   - Number of iterations in Naga's air-skin segmentation
% hough_det_sens - Sensitivity of line detection in the Hough transform.
%                  The range of the value is from 0 to 1, open ends. The 
%                  smaller the value, the more features in the image will 
%                  be considered as lines.
% fig            - Figure output. Number of seconds before next image is 
%                  processed. (fig = 0 => no output)
% im_crop        - Number of pixels to crop from borders
%                  [top, bottom, left, right]
% im_direction   - 'left' or 'right' direction of the breast.
%                  'left'  - pectoral on the right, nipple on the left.
%                  'right' - pectoral on the left, nipple on the right.
%                  Will be automatically determined if not specified.
%
% Output:
% mask           - 0 background, 1 breast, 2 pectoral
% pectoral       - Two points [I1,J1;I2,J2] on the pectoral line
% nipple         - One point [I1,J1]
%
% Andreas Eilschou, 2010

im = image;

if nargin < 2
    airskin_iter = 4;
end
if nargin < 3
    hough_det_sens = 0.08;
end
if nargin < 4
    fig = 0;
end
if nargin < 5
    im_crop = [0 0 0 0];
end
if nargin < 6
    % Automatically determine the direction of the breast (left or right)
    numberPixels = numel(im);
    if (sum(im(1:round(numberPixels/2))) < ...
            sum(im(round(numberPixels/2):numberPixels)))
        im_direction = 'left';
    else
        im_direction = 'right';
    end
end


% Initial smoothing
im = conv2(double(im),fspecial('gaussian'),'same');

% Crop image
% im_crop = [top,bottom,left,rigt]
im = im(im_crop(1)+1:end-im_crop(2),im_crop(3)+1:end-im_crop(4));

% Flip image
if strcmp('left',im_direction)
    im = fliplr(im);
end

% Air-skin segmentation
BS_mask = SaveBSRegion_simplified(im,airskin_iter);

% Pectoral segmentation
[mask,lineprm,nipple] = segmentPectoral(im,BS_mask,hough_det_sens);

% Convert line parameters (rho,theta) to (x1,y1,x2,y2)
[N,M] = size(im);
if numel(lineprm) > 0
    [x1,y1,x2,y2] = polar2borderpoints([N M],lineprm(1),lineprm(2));
end

% Re-flip mask and horizontal coordinates of points
if strcmp('left',im_direction)
    mask = fliplr(mask);
    if numel(lineprm) > 0
        y1 = M-y1+1;
        y2 = M-y2+1;
        nipple(2) = M-nipple(2)+1;
    end
end

mask(find(mask==0))=3; %Change the label of the background; SB

% Zero pad cropped mask
mask = padarray(mask,[im_crop(1) im_crop(3)],0,'pre');
mask = padarray(mask,[im_crop(2) im_crop(4)],0,'post');

% Adjust points
if numel(lineprm) > 0
    x1 = x1 + im_crop(1);
    x2 = x2 + im_crop(1);
    nipple(1) = nipple(1) + im_crop(1);
    
    y1 = y1 + im_crop(3);
    y2 = y2 + im_crop(3);
    nipple(2) = nipple(2) + im_crop(3);

    pectoral = [x1,y1;x2,y2];
else
    pectoral = [];
end

% Figure
if fig > 0
    figure(1);
    subplot(1,3,1); imagesc(image); colormap('gray'); axis image; title('image');
    if numel(lineprm) > 0
        hold on; plot(nipple(2),nipple(1),'rs'); plot([y1,y2],[x1,x2]);
        plot([y1,y2],[x1,x2],'rs'); hold off;
    end
    subplot(1,3,2); imagesc(mask); colormap('gray'); axis image; title('mask');
    if numel(lineprm) > 0
        hold on; plot(nipple(2),nipple(1),'rs'); hold off;
    end
    pause(fig)
end