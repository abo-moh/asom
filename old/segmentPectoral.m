function [mask,lineprm,nipple] = segmentPectoral(im,BS_mask,hough_detsens)
% SEGMENTPECTORAL finds the pectoral muscle in the image and returns a mask
% (0 background, 1 breast, 2 pectoral), polar coordinates for the pectoral
% line (rho,theta), and the nipple location.
%
% im            - image to be processed
% BS_mask       - mask specifying air and skin
% hough_detsens - Line detection sensitivity. Between 0 and 1
%                 The lower, the more sensitive.
%
% Andreas Eilschou, 2010

%%%% PARAMETERS %%%%
resize_height = 200;

% figure output
houghfigure = false;
pectfigure = false;
%%%%%%%%%%%%%%%%%%%%

% resize images for faster computation
resize_im = imresize(im,[resize_height NaN]);
resize_BS_mask = imresize(BS_mask,[resize_height NaN]);
resize_factor = size(im,1)/resize_height;

% the center of gravity (CG)
%STATS = regionprops(resize_BS_mask, 'Centroid');
STATS = regionprops(double(resize_BS_mask), 'Centroid');
CG = round(STATS(1).Centroid);

% reduce the ROI to above a line through CG with slope 1.5.
% CG(1) is horizontal, CG(2) is vertical
[N,M] = size(resize_BS_mask);
for j = 1:M
    line = CG(2)+floor((CG(1)-j)*1.5);
    if (line < 1)
        line = 1;
    elseif (line > N)
        line = N;
    end
    for i = line:N
        resize_BS_mask(i,j) = 0;
    end
end

% Find edges using the 3x3 Sobel operator
SOBEL = [1 2 1; 0 0 0; -1 -2 -1];
Gy = conv2(double(resize_im),SOBEL,'same');
Gx = conv2(double(resize_im),SOBEL','same');
% Gradient magnitude


% Gradient magnitude
G = sqrt(Gx.^2+Gy.^2) .* resize_BS_mask;

% Hough transform
[accum, axis_rho, axis_theta, lineprm] = Hough_Grd(G, 8, hough_detsens);

if houghfigure
    figure; imagesc(axis_theta*(180/pi), axis_rho, accum); axis xy;
    xlabel('Theta (degree)'); ylabel('Rho (pixels)');
    title('Accumulation Array from Hough Transform');
end


mask = +BS_mask;
if numel(lineprm) > 0
    % Upscale the pectoral line to the original size.
    % lineprm contains the line parameters rho,theta.    
    % Adjust rho. Theta is unchanged.
    lineprm(1) = (lineprm(1)-0.5)*resize_factor+0.5;
    
    % Identify the pectoral in the mask
    [dummy,pectoral_mask,dummy] = pectoral_area(BS_mask,lineprm(1),lineprm(2));
    mask(pectoral_mask) = 2;
    
    % Locate the nipple
    % Move a line parallel to the pectoral line to the right and locate the
    % last breast point. This is faster than computing the point-line distances
    [N,M] = size(mask);
    for line_rho = lineprm(1):sqrt(N^2+M^2)
        for I = N:-1:1
            J = round((line_rho - (I-0.5) * sin(lineprm(2))) / ...
                cos(lineprm(2)) + 0.5);
            if J < 1
                continue;
            elseif J > M
                break;
            elseif mask(I,J) == 1
                nipple = [I,J];
                break;
            end
        end
    end
else
    nipple = [];
end

% Figure output
if pectfigure
    h = figure;
    subplot(1,2,1); imshow(im,[]); title('image');
    if numel(lineprm) > 0
        DrawLines_Polar(size(im), lineprm);
        hold on; plot(nipple(2),nipple(1),'rs'); hold off;
    end
    subplot(1,2,2); imshow(mask,[]); title('mask');
    if numel(lineprm) > 0
        DrawLines_Polar(size(im), lineprm);
        hold on; plot(nipple(2),nipple(1),'rs'); hold off;
    end
    pause
    close(h)
end

end