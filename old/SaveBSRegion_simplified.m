function [BS_mask] = SaveBSRegion_simplified(im,iter)
%SaveBSRegion: This function finds Breast regions (air-skin)
%
% Naga P K R Dandu
% Nordicbioscience Imaging, 2008
%
% Modified by
% Andreas Eilschou, 2010

statBShapeLoc = 'StatSahpes';
diameter = 10;
resize_height = 1000;

% Resize image for faster computation
[N,M] = size(im);
if N > resize_height
    im = imresize(im,[resize_height NaN]);
    resized = true;
else
    resized = false;
end

% Find air-skin segmentation
segImage = adaptiveThershold(im,statBShapeLoc,diameter,iter);

% Smooth to remove abrupt edge in skin-air
H=fspecial('gaussian');
im1=filter2(H,single(segImage),'same');
BS_mask=(im1 > 0.5);

% Final brutforce to remove any unneccary connected components other
% than big
masks = bwlabel(BS_mask);
Nmb_Obj = max(max(masks));
sizes=zeros(Nmb_Obj,1);
for k = 1:Nmb_Obj
    sizes(k) = numel(find(masks == k));
end
[dummy,Indices] = max(sizes);
BS_mask = (masks  == Indices);

% Resize air-skin mask to the original size
if resized
    BS_mask = imresize(BS_mask, [N M],'nearest');
end

% Morphology to improve edges
BS_mask = imclose(BS_mask,strel('disk',10));
BS_mask = imopen(BS_mask,strel('disk',10));

end