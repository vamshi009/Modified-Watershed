%same as marker_watershed.m; thresholding is done before morphological
%processing, n-dim smoothing, %%watershed on 3 different channels

close all;
clear all;
clc;

%Read in the Color Image and Convert it to Grayscale
% rgb = imread('traffic.jpg');
% I = rgb2gray(rgb);
I = imread('coins.png');
imshow(I);

%determine adaptive threshold
I = double(I); m= size(I,1); n= size(I,2); 
maxIn = max(max(I));
minIn = min(min(I));
I1 = 1/(maxIn-minIn)*(I - (minIn.*(ones(size(I,1), size(I,2)))));
T1 = graythresh(I1);

I2 = zeros(m,n);
for i = 1:m
    for j=1:n
        if I1(i,j)>T1
            I2(i,j)= I1(i,j);
        end
    end
end

% maxIn2 = max(max(I2));
% minIn2 = min(min(I2));
% I3 = 1/(maxIn2-minIn2)*(I2 - (minIn.*(ones(size(I2,1), size(I2,2)))));
T2 = graythresh(I2);            

%smoothing
Ismooth = imgaussfilt(I1, 1);
figure, imshow(uint8(Ismooth));

%adaptive masking
T1_big = (maxIn-minIn)*T1 + minIn;
range = maxIn-minIn;
Ibin1_new = zeros(m,n);
for i = 1:m
    for j=1:n
        if Ismooth(i,j)>T1
            Ibin1_new(i,j)= Ismooth(i,j);
        end
    end
end
Ibin1 = imbinarize(Ismooth,T1);
Ibin1_new = range.*Ibin1_new + minIn.*ones(m,n);
figure
imshow(uint8(Ibin1_new));

T2_big = (maxIn-minIn)*T2 + minIn;
Ibin2 = zeros(m,n);
for i = 1:m
    for j=1:n
        if Ismooth(i,j)>T2
            Ibin2(i,j)= Ismooth(i,j);
        end
    end
end
Ibin2 = imbinarize(Ismooth,T2);
% Ibin2_new = range.*Ibin2 + minIn.*ones(m,n);
% figure
% imshow(uint8(Ibin2_new));


%impose minima
Ibwater = imimposemin(Ibin1_new, Ibin1); 
figure
imshow(Ibwater);

L = watershed(Ibwater);
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure
imshow(Lrgb)

figure
imshow(uint8(I))
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title('Lrgb superimposed transparently on original image')



