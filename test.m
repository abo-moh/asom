% Meget basal Matlab - starten fra:
% http://initora.wordpress.com/article/learn-matlab-image-processing-edge-232y3pcqwxodx-29/
% Læs et billede ind:

I = imread('images/WH 0038 RAX 99.tif', 'tif');
%I = imread('images/WH 0038 RCC 01.tif', 'tif');
%I = imread('images/coins.gif');

myfilter = fspecial('gaussian', [3 3], 0.5);
Ienhanced = imfilter(I, myfilter, 'replicate');

level = graythresh(I); %Otsu's metode til at thresholde
thresed = im2bw(Ienhanced, level); % bit-billede med thresholden i level.

% Vis billedet. Tror der mangler information i selve billedet/matlab ikke
% kan finde den så man skal selv sætte grænserne for hvad man vil definere
% som højeste og laveste værdi [low high]. Har bare valgt nogle værdi der
% gav et flot billede.

Ienhanced(~thresed) = 0; % Bit billedet "thresed" benyttet som maske
imshow(I, [22 4095]) % Vis det originale billede med low 22 og high 4095
figure, imshow(Ienhanced, [0 4095]) % Vis det enhanced billede.

% Kør en kant detektion (sobel algoritme?)
% Vi skal have den til at se bort fra baggrundsstøjen
BW1 = edge(Ienhanced,'sobel');

% Kør en kant detektion (canny algoritme?)
% Vi skal have den til at se bort fra baggrundsstøjen
BW2 = edge(Ienhanced,'canny');

% Vis resultatet af kant detektionerne
figure, imshow(BW1) % Vis resultatet af sobel kantdetektor
figure, imshow(BW2) % Vis resultater af canny kantdetektor

[H,theta,rho] = hough(BW1);
figure, imshow(imadjust(mat2gray(H)),[],'XData',theta,'YData',rho,...
    'InitialMagnification','fit');
xlabel('theta (degrees)'), ylabel('rho');
axis on, axis normal, hold on;
colormap(hot)

P=houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x = theta(P(:,2));
y = rho(P(:,1));
plot(x,y,'s','color','black');

lines = houghlines(BW1,theta,rho,P,'FillGap',5,'MinLength',7);

figure, imshow(I, [0 4095]), hold on
max_len = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

    % Plot beginnings and ends of lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

    % Determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end

% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');