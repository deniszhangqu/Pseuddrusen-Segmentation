close all;
clear all;

load('D:\GoogleDrive\Masterarbeit-Matlab\Dataset\Duke\269AMD\Farsiu_Ophthalmology_2013_AMD_Subject_1002.mat');
N=50; %the number of N-th sectional Image
im=images(:,:,N);

%img=imread('RPD1.png'); im=rgb2gray(img); 
im=double(im)/255;
im=medfilt2(im,[6 6]);
se=strel('square',3);
im_d=imdilate(im,se);
im_e=imerode(im,se);
im_new=imsubtract(im_d,im_e);
% figure,imshow(img);
figure,imshow(im);
% figure,imshow(im_d);
% figure,imshow(im_e);
figure,imshow(im_new);

%% NLL segmenttation
[M,N]=size(im_new);
[~,max_index]=max(im_new);
NFL=max_index;
figure,imshow(im); hold on,
plot(NFL,'r'); hold off
%%  RPE segmentation using 0,9 max. intensity in each colum
max_colum=max(im);
T=max_colum*0.95;
img_bin=zeros(M,N);
for i=1:1:N
    v=im(:,i);
    w=v>T(i);
    img_bin(:,i)=w(:);
end
img_3=0.5*img_bin+0.5*im; % fuse the binay image and origin image
figure,imshow(img_bin);
%%
for n=1:1:N
    for m=1:1:M
        if img_bin(m,n)==1
            img_w(m,n)=m;
        else
            img_w(m,n)=0;
        end
    end
    sum_c(n)=sum(img_w(:,n));
    sum_bin(n)=sum(img_bin(:,n));
    y_rpe(n)=sum_c(n)/sum_bin(n);
end
figure,imshow(im); hold on,
plot(y_rpe,'r');
