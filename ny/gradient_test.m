i = imread('rund_test.png');
i2 = imread('rund_test_stor.png');

i = rgb2gray(i);
i2 = rgb2gray(i2);

[g, d]  = gradients_calc(i);
[g2,d2] = gradients_calc(i2);
subdivisions = 8;
figure;
for i = 1:subdivisions
    min = (i-1)/ceil(subdivisions/2);
    max = i/ceil(subdivisions/2);
    [gf,df] = gradients_filter(g2,d2, 'min', min*pi, 'max', max*pi);
    subplot(ceil(sqrt(subdivisions)), ceil(sqrt(subdivisions)),i);
    imshow(gf);
    title(['Min: ', num2str(min), '. Max: ', num2str(max)]);
end