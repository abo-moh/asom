function [ Xmin Ymin Xmax Ymax ] = roi(img,Xmax,Ymax)
%ROI Finds the region of interest
%   Detailed explanation goes here
% Der skal måske laves en ny crop funktion. Den her virker ikke optimalt!
    % Round max coordinates to integers
    Xmax = round(Xmax);
    Ymax = round(Ymax);
    
    % Create cropped image
    imgC = imcrop(img, [1 1 Xmax Ymax]);
    
    % Find max pixel value to use as cell threshold
    cellThresh = double(max(max(imgC)));
    
    % Define threshold of columns and test columns
    threshold = cellThresh * double(size(imgC,1)) * 0.8;
    x = 1;
    while (sum(imgC(:,x)) > threshold) && x < size(imgC,2)
        x = x + 1;
    end
    
    % Define threshold of rows and test rows
    threshold = cellThresh * double(size(imgC,2));
    y = 1;
    while (sum(imgC(y,:)) > threshold) && y < size(imgC,1)
        y = y + 1;
    end
    
    % Set default X and Y minimum values
    Xmin = 1;
    Ymin = 1;
    
    % If thresholding of colums gave a plausible value
    % then set Xmin to this value
    if x < 15
        Xmin = x;
    end
    
    % If thresholding of rows gave a plausible value
    % then set Ymin to this value
    if y < 15
        Ymin = y;
    end
    
    % Define cropping width and height
    Xmax = Xmax - Xmin + 1;
    Ymax = Ymax - Ymin + 1;
end