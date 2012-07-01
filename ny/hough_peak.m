function [ peak ] = hough_peak( NH, theta, rho, resize_param )
%HOUGH_PREAK Find highest ranking peak

%% 902 Select local peaks in mowing KxK window across HN(phi, theta)
% K = 5 according to Karssemeijer
K = 5;
% Value to add/subtract from window center point
delta = (K-1)/2;
% Preallocate candidate array
peak_image = zeros(size(NH));
for x = 0:numel(theta)-1
    for y = 0:numel(rho)-1
        % Minimum x and y indexes
        min_coord = [x-delta, y-delta];
        
        % Maximum x and y indexes
        max_coord = [x+delta, y+delta];

        % Cap min and max values to image size
        min_coord = min_coord .* (min_coord >= 1)+1;
        overflow = (max_coord > [numel(theta) numel(rho)]);
        max_coord = max_coord .* ~overflow + [numel(theta) numel(rho)] .* overflow;
        window = NH(min_coord(2):max_coord(2),min_coord(1):max_coord(1));
        % Find highest value
        value = max(window(:));
        % Increment peak values by 1
        if value > 0
            peak_image(min_coord(2):max_coord(2), min_coord(1):max_coord(1)) = peak_image(min_coord(2):max_coord(2), min_coord(1):max_coord(1)) + (window == value);
        end
    end
end

% Convert peak image to binary image:
peak_image = peak_image > 0;

% Make image of peak values
peak_values = NH .* peak_image;

%% 904 Check for peaks > T_L if none return null
% T_L = 350 according to Karssemeijer
T_L = 3.5;
T_H = 6.2;

low_peaks = peak_values .* (peak_values > T_L);

if sum(low_peaks(:) > 0) == 0;
    peak = [NaN NaN];
else
    %% 908 Candidate peaks > T_H? Yes = 910, no = 914
    % T_H = 550 according to Karssemeijer
    index_high = find(low_peaks > T_H);
    
    if numel(index_high) > 0
        [row, col] = ind2sub(size(low_peaks), index_high);
        % 910 Compute pectoral area A foreach peak above T_H
        [cand_A_rho, cand_A_theta, cand_A_v] = pectoral_area(rho(row), theta(col), resize_param);
        
        % 912 Largest A = highest ranking peak
        [~, index] = max(cand_A_v);
        peak = [cand_A_rho(index) cand_A_theta(index)];
    else
        % 914 Select highest peaks as highest ranking peak
        [~, index] = max(low_peaks(:));
        [index_r index_c] = ind2sub(size(peak_values), index);
        peak = [rho(index_r) theta(index_c)];
    end
end
end