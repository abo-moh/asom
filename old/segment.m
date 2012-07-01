function [segImage,areaSeg,theImage] = segment(origImage,tune)
%segment(origImage,tune):This function segements brest using otsu
%multilevel thersolding algorithm
% origImag - original image
% tune - level
% segImage - segmented Image (only large connected component)
% areaSeg - number of pixels in segemented region
% theImages - otsus multi level image
%
% Naga P K R Dandu
% Nordicbioscience Imaging, 2008
%
% Modified by
% Andreas Eilschou, 2010

theImage = otsu_mod(origImage,tune);
%theImage = imopen(theImage,strel('disk',diameter));

%find labels
masks = bwlabel(theImage);
%find number of objects
Nmb_Obj = max(max(masks));
% find largest connect component
% MODIFIED to include the intensity. Make sure the brightest connected
% component is chosen
intensityVolumes=zeros(Nmb_Obj,1);
for k = 1:Nmb_Obj
    masks_k = masks == k;
    size_k = sum(sum(masks_k));
    intensityVolumes(k) = size_k * sum(sum(origImage(masks_k)));
end
[maxvalue,Indices] = max(intensityVolumes);
segImage = (masks  == Indices);
areaSeg=sum(sum(segImage));

% sizes=zeros(Nmb_Obj,1);
% for k = 1:Nmb_Obj
%     sizes(k) = numel(find(masks == k));
% end
% [maxvalue,Indices] = max(sizes);
% segImage = (masks  == Indices);
% areaSeg=sum(sum(segImage));
