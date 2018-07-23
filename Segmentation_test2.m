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

%% RPE segmentation using regiongrow
[g,NR,SI,TI]=regiongrow(img,1.0,0.4);
figure,imshow(g,[]);title('after Regiongrow');
f=zeros(M,N);
for m=1:1:M
    for n=1:1:N
        if g(m,n)==1
            f(m,n)=1;
        end
    end
end
figure,imshow(f,[]);title('after region_choise');
se=strel('disk',3);
f1=imclose(f,se);
figure,imshow(f1);title('after close operation');
for n=1:1:N
    for m=1:1:M
        if f1(m,n)==1
            img_w(m,n)=m;
        else
            img_w(m,n)=0;
        end
    end
    sum_c(n)=sum(img_w(:,n));
    sum_bin(n)=sum(f1(:,n));
    y_rpe(n)=sum_c(n)/sum_bin(n);
end

figure,imshow(img); hold on,
plot(y_rpe,'r'); hold on,
%%

for i=1:1:N
    if isnan(y_rpe(i))
        y_ez(i)=NaN;
        y_ch(i)=NaN;
    else
        y_ez(i)=find(img_w(:,i)',1,'first');
        y_ch(i)=find(img_w(:,i)',1,'last');
    end
end
plot(y_ez,'g'); hold on,
plot(y_ch,'b'); hold off
%% ÄâºÏRPEÇúÏß
t=1:1:N;
P = polyfit(t,y_rpe,3); y_rpes = round(polyval(P,t));
plot(y_rpes,'g');
%% ÕÒ³ödrusen
t=1:1:N;
xfill=[t fliplr(t)];
yfill=[y_rpe,fliplr(y_rpes)];
fill(xfill,yfill,'r'); hold off;


