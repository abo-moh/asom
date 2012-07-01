function slope = karssemeijer(image_string, varargin)
% Set default values of optional arguments
arg_mask = '';
arg_showall = false;
arg_resize = true;
arg_otsu = false;

% Set values of given optional arguments
c_vargs = length(varargin);
for i = 1:2:c_vargs
    switch varargin{i}
        case 'mask'
            arg_mask = varargin{i+1};
        case 'showall'
            arg_showall = varargin{i+1};
        case 'resize'
            arg_resize = varargin{i+1};
        case 'otsu'
            arg_otsu = varargin{i+1};
        otherwise
            %Invalid input!
    end
end

if arg_showall
    splot_cols = 3;
    splot_rows = 2;
else
    splot_cols = 2;
    splot_rows = 1;
end

% Read image and resize to Y x 128
img = imread(image_string);

% If resizing is selected
if arg_resize
    Ienh = imresize(img, [NaN 128]);
else
    Ienh = img;
end

% Define x and y scale values
x_scale = round(size(img,2)/size(Ienh,2));
y_scale = round(size(img,1)/size(Ienh,1));

% Read image mask if defined
if ~strcmp(arg_mask, '')
    Im = imread(arg_mask);
    if arg_resize
        Imenh = imresize(Im, [NaN 128]);
    else
        Imenh = Im;
    end
end

% Calculate middle of mass
[mx my] = massemidtpunkt(Ienh);
IenhSize = size(Ienh);

% Perform vertical flip if middle of mass is on the
% right-hand side of the image
if (mx > round(IenhSize(2)/2))
    mx = IenhSize(2) - mx;
    Ienh = flipdim(Ienh, 2);
    img = flipdim(img, 2);
    % Flip mask image as well
    if ~strcmp(arg_mask, '')
        Imenh = flipdim(Imenh, 2);
        Im = flipdim(Im, 2);
    end
end

% Perform cropping of white lines and find roi
[cropXmin cropYmin cropXmax cropYmax] = roi(Ienh, mx, my);
Ienh2 = imcrop(Ienh, [cropXmin cropYmin cropXmax cropYmax]);

% Perform upscaling of cropping points
crop = [[cropXmin;cropXmax].*x_scale, [cropYmin; cropYmax].*y_scale];

% Perform canny edge detection of roi
% @edges: binary image of edges
% @thresh: high and low threshholds used by canny
[edges, thresh] = edge(Ienh2,'canny');

% Perform Hough Transformation
[H,theta,rho] = hough(edges);

%Find 5 peaks og kig kun på celler med 0.3 * den højeste værdi.
P=houghpeaks(H,20,'threshold',ceil(0.3*max(H(:))));
%Tegn punkterne ind på H:
x = theta(P(:,2));
y = rho(P(:,1));

%Find linjerne der svarer til hough peaks'ne:
lines = houghlines(edges,theta,rho,P,'FillGap',5,'MinLength',7);

%Find den længste linje og kompenser for tilpasning af billedet
max_len = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    xy(:,1) = xy(:,1)+cropXmin; %Kompensering for massemidtpunktet
    xy(:,2) = xy(:,2)+cropYmin;
    %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

    % Plot beginnings and ends of lines
    %plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    %plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

    % Determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end

% highlight the longest line segment
%plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');

% Udregn hældning a = (y2 - y1) / (x2 - x1)
pl_a = (xy_long(2,2) - xy_long(1,2)) /(xy_long(2,1) - xy_long(1,1));
if pl_a == Inf
    %Hvis hældningen er uendelig, er der en lodret linje:
    pl_b = IenhSize(1);
    slope_x = [xy_long(1,1);xy_long(1,1)];
else
    %Ikke-lodret linje, udregn x-koordinater:
    pl_b = round(xy_long(1,2) - pl_a * xy_long(1,1));
    slope_x = [0;round(pl_b/-pl_a)];
end

slope_y = [pl_b;0];
slope = [pl_a pl_b];

o_slope_x = slope_x * x_scale;
o_slope_y = slope_y * y_scale;

if arg_otsu
    % Perform Otsu's method to find threshold of image
    I_otsu = Ienh;
    otsu_c = 2; % Use the method otsu_c times
    final = false(size(I_otsu));
    for i = 1:otsu_c
        bound = graythresh(I_otsu); % Use Otsu's method to find threshold
        I_bw = im2bw(I_otsu, bound); % Use threshold to make binary image
        final = final | I_bw; % Combine previous binary image with additional imagespace
        I_otsu(final) = 0;
    end
    % Show results of Otsu
    figure('name', 'Otsu','numbertitle','off');
    subplot(2,2,1);
    imshow(final);
    title('Otsu mask')
    subplot(2,2,2);
    I_bbound = Ienh;
    I_bbound(~final) = false;
    %imshow(I_bw)
    imshow(I_bbound, [min(I_bbound(:)), max(I_bbound(:))])
    title('Otsu mask on I enhanced')

    subplot(2,2,3);
    I_edge = edge(final, 'canny');
    imshow(I_edge)
    title('Edge detection');
end


% Show image used for calculation
figure('name',image_string,'numbertitle','off')
subplot(splot_rows,splot_cols,1);
imshow(Ienh, [min(Ienh(:)) max(Ienh(:))]),hold on;
if arg_resize
    title('Resized')
else
    title('Original')
end

% Draw pectoral line
plot(slope_x,slope_y, 'LineWidth', 1,'Color','blue');

% Show image mask is set
if ~strcmp(arg_mask, '')
    subplot(splot_rows,splot_cols,2);
    imshow(Imenh, [min(Imenh(:)) max(Imenh(:))]), hold on
    
    % Draw pectoral line on mask
    plot(slope_x, slope_y, 'LineWidth', 2,'Color','blue');
    if arg_resize
        title('Resized mask')
    else
        title('Original mask')
    end
end

% Show 'all' graphics
if arg_showall
    subplot(splot_rows,splot_cols,4);
    imshow(img, [min(img(:)) max(img(:))]), hold on
    o_mx = x_scale * mx;
    o_my = y_scale * my;
    
    plot(o_mx,o_my, 's','color', 'green');
    patch([crop(1,1) crop(1,1) o_mx o_mx], [crop(1,2) o_my o_my crop(1,2)], 'red', 'FaceColor', 'none');
    plot(o_slope_x, o_slope_y, 'LineWidth', 1, 'Color', 'blue');
    title('Original')
    if ~strcmp(arg_mask, '')
        subplot(splot_rows,splot_cols,5);
        imshow(Im, [min(Im(:)) max(Im(:))]), hold on
        plot(o_slope_x, o_slope_y, 'LineWidth', 1, 'Color', 'blue');
        title('Original mask')
    end
    
    % Show results of hough transformation
    subplot(splot_rows,splot_cols,3);
    imshow(imadjust(mat2gray(H)),[],'XData',theta,'YData',rho,...
        'InitialMagnification','fit');
    xlabel('theta (degrees)'), ylabel('rho');
    axis on, axis normal, hold on;
    %colormap(hot)
    
    plot(x,y,'s','color','blue');
    title('Hough transformation')
end

end