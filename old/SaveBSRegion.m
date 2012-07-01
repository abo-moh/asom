function SaveBSRegion(inDir,outdir,statBShapeLoc,diameter)
%SaveBSRegion(inDir):This function saves Breast regions (air-skin)
%including pectoral musle in BRSSegmented directry which is created in a direcory
%where image directory (inDir) presents.
% inDir - image ldirectory
% outdir - localtion to save these brest skin masks 
% diameter - diameter for morphological operation (By defauly it is 10)
%
% Naga P K R Dandu
% Nordicbioscience Imaging, 2008
%
% Modified by
% Andreas Eilschou, 2010

%Deault data location is relative data directory (default convention)
if nargin < 1
    inDir = '../NijmegenRAW/good';
end
if nargin < 2
    outdir = '../NijmegenRAW/good';
end
if nargin < 3
    statBShapeLoc = 'StatSahpes';
end
if nargin < 4
    diameter = 10;
end

% %if input from mosix scripts (strict checking)
% if (ischar(diameter))
%     diameter=str2num(diameter);
% end
%     
%look for all tif images
dirc=dir([inDir '/*.tif']);
dirc=dirc(1:end);
 
% if ~exist(outdir,'dir')
%     mkdir(outdir);
% end

%%% MODIFIED START %%%
cutoff = 80;
%%% MODIFIED END %%%

for i=1:numel(dirc)
    %brute
    if strcmp(dirc(i).name(end-3:end),'.tif')
        im=imread([inDir '/' dirc(i).name]);
        dirc(i).name
               
        %%% MODIFIED START %%%
        % Make sure the body is on the left side of the mammogram.
        numberPixels = numel(im);
        flip = (sum(im(1:round(numberPixels/2))) < ...
                sum(im(round(numberPixels/2):numberPixels)));
        if flip
            im = fliplr(im);
            display('Nijmegen flipped')
            %im = im(:,1:765); % cut off the right border
        end
        
        cutoff_im = im(cutoff+1:end-cutoff,:);
        %%% MODIFIED END %%%
        
        segImage = adaptiveThershold(cutoff_im,statBShapeLoc,diameter);
        
        %smooth to remove abrupt edge in skin-air
        H=fspecial('gaussian');
        im1=filter2(H,single(segImage),'same');
        segImage=(im1 > 0.5);
        
        %final brutforce to remove any unneccary connected components other
        %than big
        masks = bwlabel(segImage);

        Nmb_Obj = max(max(masks));
        sizes=zeros(Nmb_Obj,1);
        for k = 1:Nmb_Obj
            sizes(k) = numel(find(masks == k));
        end
        [~,Indices] = max(sizes);
        segRImage = (masks  == Indices);
        
        %morphology to improve edges
        segRImage = imclose(segRImage,strel('disk',10));
        segRImage = imopen(segRImage,strel('disk',10));
        
        %%% MODIFIED START %%%
        % Convert cutoff segmentation to full image.
        zeropadding = zeros(cutoff,size(segRImage,2));
        segRImage = [zeropadding; segRImage; zeropadding];
        %%% MODIFIED END %%%
        segRImage = double(segRImage).*double(im);
        save([outdir '/' dirc(i).name(1:end-4) '.mat'], 'im','segRImage','flip');

    end
end
