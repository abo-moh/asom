self = dir('../../nijmegen-fuckup/');
for i = 3:length(self)
  filename = ['../../nijmegen-fuckup/' self(i).name];
  disp(self(i).name);
  im = imread(filename);
  figure(); imshow(im, [min(im(:)) max(im(:))]);
  karssemeijer('../../nijmegen-fuckup/', self(i).name, 'resize', true, 'false', true, 'canny', true);
end
