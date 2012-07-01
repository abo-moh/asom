% Remove all previous calculations and figures
clear all;
close all;

% Get screen size
scrsz = get(0,'ScreenSize');

% Set directory for mammographies
directory = '../../nijmegen-fuckup/no_skinair/';

% List content of directory
images = dir(directory);

for i = 3:length(directory)
    if ~images(i).isdir
        disp(images(i).name);
        im = imread([directory images(i).name]);
        figure(); imshow(im, [min(im(:)) max(im(:))]);
        karssemeijer(directory, images(i).name, 'resize', true, 'showall', true);
        karssemeijer(directory, images(i).name, 'resize', true, 'showall', true, 'canny', true);
    end
end