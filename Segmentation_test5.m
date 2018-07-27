%%first calculate the flattend image based on RPE poly linie
%%then calsulate the  average intensity of each row
%%
close all;
clear all;
clc;
% load('D:\GoogleDrive\Masterarbeit-Matlab\Dataset\Duke\269AMD\Farsiu_Ophthalmology_2013_AMD_Subject_1004.mat');
% B=50; %the number of N-th sectional Image
% img=images(:,:,B);
% img(1:10,:)=img(11:20,:);

img=imread('RPD4.jpg');
img=rgb2gray(img); %3 channels to 1 channel
%%
[M,N]=size(img);
img=double(img)/255;
%img=mat2gray(img); %[0 1.0]
re=fspecial('gaussian'); 
img=imfilter(img,re);
figure,imshow(img); title('after gaussian filtering');

%% max. intensity points als rpe and fit two order poly
N1=1;
[img_shift,y_rpe,int_shift]=img_rpe_shift(img,N1);
k=median(y_rpe);
figure,imshow(img_shift); title('fattend image '); hold on,
plot(y_rpe+(k-y_rpe),'r-'); title('fisrt time flattend image and flattend RPE'); hold off

N2=2; img2=img_shift;
[img_shift2,y_rpe2,int_shift2]=img_rpe_shift(img2,N2);
k2=median(y_rpe2);
figure,imshow(img_shift2); title('fattend image '); hold on,
plot(y_rpe2+(k2-y_rpe2),'r-'); title('fisrt time flattend image and flattend RPE'); hold off

%% 找到较大波峰中间的极小值
mean_row=sum(img_shift2,2)';
figure,plot(mean_row);hold on,
[~,locs]=findpeaks(mean_row,'MINPEAKHEIGH',0.5*max(mean_row)); 
plot(locs,mean_row(locs),'ro');hold on

[~,loc_min]=min(mean_row(locs(1):locs(end)));
stem(loc_min+locs(1)-1,mean_row(loc_min+locs(1)-1),'b'); hold off
%% 该极小值左边和右边最大波峰分别对应NFL 和 RPE

[~,locs]=findpeaks(






