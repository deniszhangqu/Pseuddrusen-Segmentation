close all;
clear all;
% load('D:\GoogleDrive\Masterarbeit-Matlab\Dataset\Duke\269AMD\Farsiu_Ophthalmology_2013_AMD_Subject_1004.mat');
% N=50; %the number of N-th sectional Image
% img=images(:,:,N);

img=imread('RPD4.jpg');
img=rgb2gray(img); %3 channels to 1 channel

img=double(img)/255;
img=mat2gray(img); %[0 1.0]
figure,imshow(img);
% img(1:8,:)=img(17:32,:);
re=fspecial('gaussian'); title('original image');
img=imfilter(img,re);
figure,imshow(img); title('after gaussian filtering');
%% max. intensity points als rpe and fit two order poly
[M,N]=size(img); X=zeros(N); Y=zeros(N);
max_column=max(img);
for i=1:1:N
    img_max(:,i)=img(:,i)==max_column(i);
end
figure,imshow(img);hold on,
[X,Y]=find(img_max);
P=polyfit(Y,X,2);
y_rpe=polyval(P,Y);
plot(y_rpe,'r-'); title('Image with fitted RPE'); hold off;
%%  circshift the OCT image to make the RPE geradeaus
k=median(y_rpe);
shift_int=int8(k-y_rpe);
img_shift=zeros(M,N);
for i=1:1:N
    img_shift(:,i)=circshift(img(:,i),shift_int(i));
end
figure,imshow(img_shift); title('fattend image '); hold on,
plot(y_rpe+(k-y_rpe),'r-');
%%
Y(1:N)=0;
gradient_mask_ld=[1;-1]; %like RPE-Choriod
gradient_mask_dl=[-1;1];  %like NFL,OS-RPE
img_gradient=imfilter(img_shift,gradient_mask_dl,'replicate');
img_gradient=img_gradient-min(img_gradient(:))./max(img_gradient(:)-min(img_gradient(:))); %Normierung

%%
im_g=img_gradient(30:80,150:200);
[M,N1]=size(im_g);
N2=N1+2; %左右添加一列
g_min=0.0001;
n=1;
%% calculate the weights bitween nodes
for k=1:1:M*N2
    if rem(k,N2)~=1 && rem(k,N2)~=N2-1 && rem(k,N2)~=0 % 处理2：N2-2列
        if fix(k/N2)+1==1 % 处理第一行
            s(n:n+1)=k;
            t(n:n+1)=[k+1,k+1+N2];
            g_k=im_g(fix(k/N2)+1,rem(k,N2)-1);
            %g_k1=im_g(fix((k+1-N2)./N2)+1,rem(k+1-N2)-1); %最上面一行不能再向上走
            g_k2=im_g(fix((k+1)/N2)+1,rem(k+1,N2)-1);
            g_k3=im_g(fix((k+1+N2)/N2)+1,rem(k+1+N2,N2)-1);
            weights(n:n+1)=[2-(g_k+g_k2)+g_min,2-(g_k+g_k3)+g_min];
            n=n+2;
        elseif fix(k/N2)+1==N1 %处理最后一行
            s(n:n+1)=k;
            t(n:n+1)=[k+1-N2,k+1];
            g_k=im_g(fix((k)./N2)+1,rem(k,N2)-1);
            g_k1=im_g(fix((k+1-N2)./N2)+1,rem(k+1-N2,N2)-1);
            g_k2=im_g(fix((k+1)./N2)+1,rem(k+1,N2)-1);
            %g_k3=im_g(fix((k+1+N2)./N2)+1,rem(k+1+N2,N2)-1); %最下面一行不能再向下
            weights(n:n+1)=[2-(g_k+g_k1)+g_min,2-(g_k+g_k2)+g_min];
            n=n+2;
        else
            s(n:n+2)=k;
            t(n:n+2)=[k+1-N2,k+1,k+1+N2];
            g_k=im_g(fix((k)./N2)+1,rem(k,N2)-1);
            g_k1=im_g(fix((k+1-N2)./N2)+1,rem(k+1-N2,N2)-1);
            g_k2=im_g(fix((k+1)./N2)+1,rem(k+1,N2)-1);
            g_k3=im_g(fix((k+1+N2)./N2)+1,rem(k+1+N2,N2)-1);
            weights(n:n+2)=[2-(g_k+g_k1)+g_min,2-(g_k+g_k2)+g_min,2-(g_k+g_k3)+g_min];
            n=n+3;
        end
        
    elseif  rem(k,N2)==1
        if fix(k/N2)+1==1
            s(n:n+1)=k;
            t(n:n+1)=[k+1,k+1+N2];
            weights(n:n+1)=[g_min,g_min];
            n=n+2;
        elseif fix(k/N2)+1==N1
            s(n:n+1)=k;
            t(n:n+1)=[k+1-N2,k+1];
            weights(n:n+1)=[g_min,g_min];
            n=n+2;
        else
            s(n:n+2)=k;
            t(n:n+2)=[k+1-N2,k+1,k+1+N2];
            weights(n:n+2)=[g_min,g_min,g_min];
            n=n+3;
        end
        
    elseif rem(k,N2)==N2-1
        
        if fix(k/N2)+1==1
            s(n:n+1)=k;
            t(n:n+1)=[k+1,k+1+N2];
            weights(n:n+1)=[g_min,g_min];
            n=n+2;
        elseif fix(k/N2)+1==N2
            s(n:n+1)=k;
            t(n:n+2)=[k+1-N2,k+1];
            weights(n:n+1)=[g_min,g_min];
            n=n+2;
        else
            s(n:n+2)=k;
            t(n:n+2)=[k+1-N2,k+1,k+1+N2];
            weights(n:n+2)=[g_min,g_min,g_min];
            n=n+3;
        end
    else
        continue;
    end
end
%%

G = digraph(s',t',weights);
path=shortestpath(G,1,2703);
%%
j=1;
for i=1:1:N2
    if rem(path(i),N2)~=1 && rem(path(i),N2)~=0
        n(j)=fix((path(i)-1)/N2);
        m(j)=rem(path(i),N2)-1;
        j=j+1;
    else
        continue;
    end
end
%figure,plot(G,'Layout','force')
%%
figure,imshow(im_g); hold on,
plot(m,n,'r*');hold off

