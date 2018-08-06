close all;
clear all;
%Quelle:image processing of OCT using Matlab%
% load('D:\GoogleDrive\Masterarbeit-Matlab\Dataset\Duke\269AMD\Farsiu_Ophthalmology_2013_AMD_Subject_1004.mat');
% N=50; %the number of N-th sectional Image
% img=images(:,:,N);

img=imread('RPD5.jpg'); img=rgb2gray(img);

[M,N]=size(img);
%%preprocessing ¡¾A1¡¿
img=double(img)/255;
img_med=medfilt2(img,[5 5]);
img_med=mat2gray(img_med);

%% RPE segmentation using regiongrow
Ts=max(img(:));
[g,NR,SI,TI]=regiongrow(img,Ts,0.4);

f=zeros(M,N);
for m=1:1:M
    for n=1:1:N
        if g(m,n)==1
            f(m,n)=1;
        end
    end
end

se=strel('disk',3);
f1=imclose(f,se);

%% find boundary and drusendetecction and classification
[y_rpe,y_ez,y_ch]=findboundary(f1,1); 
[peaks,p_drusen,p_normaldrusen,p_rpd]=finddrusen(y_ez,y_ch,7,2);

%% show of boundary and drusen 
figure,imshow(img); title('origin image');

figure,imshow(img); hold on,
plot(y_rpe,'r'); hold on,
plot(y_ez,'g'); hold on,
plot(y_ch,'b'); hold on
plot(p_drusen,y_ez(p_drusen),'o');hold on,
if ~isempty (p_normaldrusen)
    plot(p_normaldrusen,y_ez(p_normaldrusen),'r*'); hold on,
    legend('y\_rpe','y\_ez','y\_ch','kandidate drusen','normal drusen');
    if ~isempty(p_rpd)
        plot(p_rpd,y_ez(p_rpd),'rs');
        legend('y\_rpe','y\_ez','y\_ch','kandidate drusen','normal drusen','reticular pseudedrusen'); hold off
    end
else
    if ~isempty(p_rpd)
        plot(p_rpd,y_ez(p_rpd),'rs');
        legend('y\_rpe','y\_ez','y\_ch','kandidate drusen','reticular pseudedrusen');
    end
end
hold off
title('results after segmentation');