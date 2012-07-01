% clean up
clear all
close all

% which directories we're looking a´t
mammographies = '../../nijmegen/mammographies/';
mammographies_dir = dir(mammographies);
masks = '../../nijmegen/masks/';
masks_dir = dir(masks);

% pick random images
pick = 1;
images = randi(numel(mammographies), 1, pick);

for i = images
    [pathstr, name, ext] = fileparts([mammographies mammographies_dir(i).name]);
    if ~mammographies_dir(i).isdir && ~strcmp(ext, '.ini') && name(1) == 'l'
        mask_name = [masks 'mask' mammographies_dir(i).name(2 : end)];
        if exist(mask_name, 'file') == 2
            % setup figure
            figure('name', 'Canny')
            
            % read in the original image
            mammography = imread([pathstr '/' name ext], 'tif');
            subplot(2, 2, 1);
            imshow(mammography, [min(mammography(:)) max(mammography(:))]);
            title('Mammography');
            
            % adjust size and resolution using Karssemeijer
            [small, flip, adj_crop, adj_scale] = adjust_size(mammography, true);
            % find the region of interest using Karssemeijer
            [roi, col, row, crop] = roi_id(small);
            % and plot roi
            subplot(2, 2, 2);
            imshow(roi, [min(roi(:)) max(roi(:))]);
            title('ROI');
            
            % run Canny edge detection
            canny = edge(roi, 'canny');
            % and plot
            subplot(2, 2, 3);
            imshow(canny);
            title('Canny');
            
            % perform Hough transformation
            [H, theta, rho] = hough(canny);
            peaks = houghpeaks(H, 3);
            lines = houghlines(canny, theta, rho, peaks);
            if ~isempty(lines)
                max = 0;
                for i = 1:length(lines)
                    current = norm(lines(i).point1 - lines(i).point2);
                    if current > max
                        coords = [lines(i).point1; lines(i).point2];
                        max = current;
                    end
                end
            end
        end
    end
end
