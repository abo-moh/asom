function [ pec_a pec_b ] = karssemeijer(directory, name, varargin)
 %KARSSEMEIJER Perform all steps in the Karssemeijer patent
% Set default values of optional arguments
arg_showall = false;
arg_resize = true;
arg_canny = false;

% Set values of given optional arguments
c_vargs = length(varargin);
for i = 1:2:c_vargs
    switch varargin{i}
        case 'showall' %Show all intermediate calculations and images
            arg_showall = varargin{i+1};
        case 'resize' %Resize image to [NaN 128]
            arg_resize = varargin{i+1};
        case 'canny'
            arg_canny = varargin{i+1};
        otherwise
            %Invalid input!
    end
end

% Initial reading of image
img = imread([directory name]);
img = im2double(img);

%% Adjust size / resolution (402)
[img_small, flip, adj_crop, adj_scale] = adjust_size(img, arg_resize);

%% Identify ROI (404)
[img_roi, mass_col, mass_row, roi_crop] = roi_id(img_small, arg_canny);

% Pack together all resize parameters
resize_param = {adj_crop, adj_scale, roi_crop, size(img)};

%% If canny is enabled:
if arg_canny
    % Perform canny edge detection followed by Houghtransformation
    edges = edge(img_roi, 'canny');
    [H, theta, rho] = hough(edges);
    NH = hough_normalize(H, theta, rho, size(edges), false);
    % Find peaks and convert to points
    peak = houghpeaks(NH,1);
    % Define peak as theta, rho instead of theta index, rho index
    peak = [rho(peak(1)) theta(peak(2))];
else
    %% Gradient magnitudes and gradient directions (406)
    [g_mag, g_dir] = gradients_calc(img_roi);
    %[g_mag, g_dir] = gradients_calc(g_mag);

    %% Filter gradient magnitudes (408)
    [g_magF, g_dirF] = gradients_filter(g_mag, g_dir);

    %% Hough transform (410)
    [H, theta, rho] = hough_transform(g_magF);
    
    %% Normalize hough plane (412)
    NH = hough_normalize(H, theta, rho, size(img_roi), true);
    %% Highest ranking peak (414)
    peak = hough_peak(NH, theta, rho, resize_param);
end

% If peak is found, set pec_a and pec_b to their respective sizes
if ~isnan(peak)
    [pec_a, pec_b] = hough2xy(peak(1),peak(2), 'resize', resize_param);
    
    % Calculate endpoints on original (flipped) image:
    peak_coords = hough2points(peak(2), peak(1), size(img_roi));
else
    pec_a = NaN;
    pec_b = NaN;
end

%% Calculate pectoral line
if ~isnan(peak)
    peak_coords = hough2points(peak(2), peak(1), size(img_roi));
    peak_mini = hough2points(peak(2), peak(1), size(img_small), 'resize', {[0 0;0 0],[1 1], roi_crop, size(img_small)});
    peak_orgc = hough2points(peak(2), peak(1), size(img), 'resize', resize_param);
    
    if flip
        width = size(img,2);
        peak_orgc(:,1) = width - (peak_orgc(:,1)-1);
        peak_orgc(:,3) = width - (peak_orgc(:,3)-1);
    end
end

%% Create mask of area contained by pectoral line
pec_area = ab2area(size(img,1),size(img,2), pec_a, pec_b);
if flip
    pec_area = flipdim(pec_area,2);
end
pec_area = im_capframe(pec_area, flip, adj_crop);

%% Show intermediate calculations
if arg_showall
    figure('name','Original image', 'numbertitle','off');
    imshow(img, [min(img(:)) max(img(:))]);
    hold on;
    if ~isnan(peak)
        plot([peak_orgc(1,1) peak_orgc(1,3)],[peak_orgc(1,2) peak_orgc(1,4)],'LineWidth',1,'Color','red');
    end
    
    hold on;
    if ~isnan(peak)
        plot([peak_orgc(1,1) peak_orgc(1,3)],[peak_orgc(1,2) peak_orgc(1,4)],'LineWidth',1,'Color','red');
    end
    
    figure('name','Intermediate calculations');
    subplot(2,4,1)
    imshow(img_small, [min(img_small(:)) max(img_small(:))]);
    title('402');
    
    subplot(2,4,2)
    imshow(img_roi, [min(img_roi(:)) max(img_roi(:))]);
    hold on;
    plot([peak_coords(1,1) peak_coords(1,3)], [peak_coords(1,2) peak_coords(1,4)], 'LineWidth', 1, 'Color', 'red');
    title('404');
    
    if ~arg_canny
        subplot(2,4,3)
        imshow(g_mag, [min(g_mag(:)) max(g_mag(:))]);
        subplot(2,4,4);
        imshow(g_dir, [min(g_dir(:)) max(g_dir(:))]);
        title('406');

        subplot(2,4,5)
        imshow(g_magF, [min(g_magF(:)) max(g_magF(:))]);
        title('408');
    end
    
    subplot(2,4,6);
    imshow(imadjust(mat2gray(H)),[],'XData',theta,'YData',rho,...
        'InitialMagnification','fit');
    xlabel('\theta (grader)'), ylabel('\rho');
    axis on, axis normal, hold on;
    plot(peak(2),peak(1),'Xr');
    title('410');
    
    subplot(2,4,7);
    imshow(imadjust(mat2gray(NH)),[],'XData',theta,'YData',rho,...
        'InitialMagnification','fit');
    xlabel('\theta (grader)'), ylabel('\rho');
    axis on, axis normal, hold on;
    plot(peak(2),peak(1),'Xr');
    
    title('412');
    subplot(2,4,8);
    imshow(img_small, [min(img_small(:)) max(img_small(:))]), hold on;
    plot([peak_mini(:,1) peak_mini(:,3)],[peak_mini(:,2) peak_mini(:,4)], 'LineWidth',1,'color','red');
    title('414');
    
end

%% Perform otsu to find skin-air boundary
% Calculate angle theta_org and rho_org of pectoral line in original image
% coordinates:

if ~isnan(pec_a) && ~isnan(pec_b)
    % Calculate original theta and rho
    [theta_org, rho_org] = ab2hough(pec_a,pec_b);
end
[A, B, C, constants, sa_bound, sa_pol] = skin_air(img, [theta_org, rho_org], flip, resize_param);
sa_bound = (pec_area & sa_bound) * 2 + (~pec_area & sa_bound);
figure, imshow(flipdim(sa_bound,2), [0 2]);

%figure, imshow(sa_bound, [min(sa_bound(:)) max(sa_bound(:))]);
%sa_pol = sa_pol + pec_area;

% Things to save in shapes/<image_name>.mat:
% A, B, C
% constant a for upper skin-air polynomial (constants(1))
% constant a for lower skin-air polynomial (constants(2))
% constants a and b for pectoral line (pec_a, pec_b)
% boolean for flipping of image or not
% 2 images with 0 = background, 1 = breast, 2 = pectoral:
% 1 where breast is modelled by skin-air (sa_bound)
% and 1 where breast is modelled by upper and lower skin-air polynomial
% (sa_pol)
if arg_canny
    save_name = ['shapes_canny/' name '.mat'];
else
    save_name = ['shapes_karssemeijer/' name '.mat'];
end
save(save_name, 'A','B','C','constants', 'flip','pec_a','pec_b','sa_bound','sa_pol');
end