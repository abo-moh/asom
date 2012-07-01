function [ mx my ] = massemidtpunkt ( img )
    type = 1;
    if type == 1
        %type 1 her. Virker vist ikke helt 100%
        [rows, cols] = size(img);
        x = ones(rows, 1) * [1 : cols];
        y = [1 : rows]' * ones(1, cols);

        area = sum(sum(img));
        mx = sum(sum(double(img).* x)) / area;
        my = sum(sum(double(img).* y)) / area;
    elseif type == 2
        %type 2 her
        level = graythresh(img);
        Ibin = im2bw(img, level);
        props = regionprops(Ibin, 'Centroid')
        mx = props.Centroid(1);
        my = props.Centroid(2);
    end
end