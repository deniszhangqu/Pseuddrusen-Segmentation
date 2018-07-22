close all;
clear all;
%Quelle:image processing of OCT using Matlab%
% load('D:\GoogleDrive\Masterarbeit-Matlab\Dataset\Duke\269AMD\Farsiu_Ophthalmology_2013_AMD_Subject_1004.mat');
% N=50; %the number of N-th sectional Image
% img=images(:,:,N);
img=imread('RPD2.jpg'); img=rgb2gray(img);

[M,N]=size(img);
%%preprocessing 【A1】
img=double(img)/255;
img_med=medfilt2(img,[3 3]);
img_med=mat2gray(img_med);
% figure,imshow(img);
% figure,imshow(img_med);

%%  RPE segmentation using 0,9 max. intensity in each colum 【B1】
max_colum=max(img_med);
T=max_colum*0.90;
img_bin=zeros(M,N);
for i=1:1:N
    v=img_med(:,i);
    w=v>T(i);
    img_bin(:,i)=w(:);
end
img_3=0.5*img_bin+0.5*img_med; % fuse the binay image and origin image
figure,imshow(img_bin);
%%在segemnt RPE的时候收到NFL的干扰
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
    y_rpe(n)=sum_c(n)/sum_bin(n);   %%传统Drusen取中线作为RPE可以，但是 RPD不行
end
figure,imshow(img); hold on,
plot(y_rpe,'r');
%% 通过区域生长消除NFL的干扰，几乎成功 【B2】(先区域生长除去NFL影响，在用B1法提取RPE边缘）
[g,NR,SI,TI]=regiongrow(img,1.0,0.4);
figure,imshow(g,[]);title('after Regiongrow');
f=zeros(M,N);
for m=1:1:M
    for n=1:1:N
        if g(m,n)==1 %有可能不是1，而是2.
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
%% 拟合RPE曲线
t=1:1:N;
P = polyfit(t,y_rpe,3); y_rpes = round(polyval(P,t)); 
plot(y_rpes,'g'); 
%% 找出drusen
t=1:1:N;
xfill=[t fliplr(t)];
yfill=[y_rpe,fliplr(y_rpes)];
fill(xfill,yfill,'r'); hold off;


