% Remove all previous calculations and figures
clear all;
close all;

% Get screen size
scrsz = get(0,'ScreenSize');

% Set directory for mammographies
directory = '../../nijmegen/mammographies/';

% List content of directory
images = dir(directory);

for i = 1:3
    if ~images(i).isdir
        element = randi(size(images,1));
        karssemeijer(directory, images(605).name, 'resize', true, 'showall', true, 'canny', true);
        %karssemeijer(directory, images(i).name, 'resize', true, 'showall', true, 'canny', true);
    end
end