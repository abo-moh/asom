%Set up directory to perform karssemeijer on
directory = '../nijmegen/mammographies/';
mask_dir = '../nijmegen/masks/';
%directory = 'images/';
img = dir(directory);
elmCount = numel(img);
i = 0;
%Pick 10 random images and their corresponding masks files
count = 10;
elements = randi(elmCount,1,count);
results = zeros(count,2);
i = 1;
for k = elements
    [pathstr, name, ext] = fileparts([directory img(k).name]);
    % Check if file is a .ini file
    if ~img(k).isdir && ~strcmp(ext,'.ini') && name(1) == 'l'
        % Define mask filename
        mask_name = [mask_dir 'mask' img(k).name(2:end)];
        if exist(mask_name, 'file') == 0
            % Perform Karssemeijer on the image
            karssemeijer([directory img(k).name], 'mask', NaN);
        else
            results = karssemeijer([directory img(k).name], 'mask', mask_name, 'showall', true, 'resize', true, 'otsu', true);
            mask_results = karssemeijer(mask_name, 'showall', true);
            results(i,:) = results ./ mask_results;
            if results(i,1) == Inf || results(i,1) == -Inf
                results(i,1) = 1;
            end
            if results(i,2) == Inf || results(i,2) == -Inf
                results(i,2) = 1;
            end
        end
    end
    i = i+1;
end
disp('Afvigelse i hhv. a og b (ax + b forskrift)');
disp(results);