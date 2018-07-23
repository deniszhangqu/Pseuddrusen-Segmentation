close all;
clear all;
%Quelle:image processing of OCT using Matlab%
% load('D:\GoogleDrive\Masterarbeit-Matlab\Dataset\Duke\269AMD\Farsiu_Ophthalmology_2013_AMD_Subject_1004.mat');
% N=50; %the number of N-th sectional Image
% img=images(:,:,N);
img=imread('RPD1.png'); img=rgb2gray(img);

[M,N]=size(img);
%%preprocessing ¡¾A1¡¿
img=double(img)/255;
img_med=medfilt2(img,[3 3]);
img_med=mat2gray(img_med);
% figure,imshow(img);
% figure,imshow(img_med);

%%  RPE segmentation using 0,9 max. intensity in each colum ¡¾B1¡¿
im=img_med;tf=0.9;
[im_bin,y_rpe]=RPE_colummax(im,tf);
figure,imshow(img); hold on,
plot(y_rpe,'r'); hold on
%% ÄâºÏRPEÇúÏß
t=1:1:N;
P = polyfit(t,y_rpe,3); y_rpes = round(polyval(P,t));
plot(y_rpes,'g');
%% ÕÒ³ödrusen
t=1:1:N;
xfill=[t fliplr(t)];
yfill=[y_rpe,fliplr(y_rpes)];
fill(xfill,yfill,'r'); hold off;
