close all;
clear all;
load('D:\GoogleDrive\Masterarbeit-Matlab\Dataset\Duke\269AMD\Farsiu_Ophthalmology_2013_AMD_Subject_1004.mat');
N=50; %the number of N-th sectional Image
img=images(:,:,N);

% img=imread('RPD3.png');
% img=rgb2gray(img); %3 channels to 1 channel

img=double(img)/255;
img=mat2gray(img); %[0 1.0]
figure,imshow(img);
img(1:8,:)=img(9:16,:);
re=fspecial('gaussian'); title('original image');
img=imfilter(img,re);
figure,imshow(img); title('after gaussian filtering');
%% max. intensity points als rpe and fit two order poly
[M,N]=size(img);
max_column=max(img);
for i=1:1:N
    img_max(:,i)=img(:,i)==max_column(i);
end
figure,imshow(img);hold on,
[X,Y]=find(img_max);
P=polyfit(Y,X,2);
y_rpe=polyval(P,Y);
plot(y_rpe,'r-'); title('Image with fitted RPE');
%%
k=median(y_rpe);
shift_int=int8(k-y_rpe);
img_shift=zeros(M,N);
for i=1:1:N
    img_shift(:,i)=circshift(img(:,i),shift_int(i));
end
figure,imshow(img_shift); title('fattend image ');
%%
gradient_mask_ld=[1;-1]; %like RPE-Choriod
gradient_mask_dl=[-1;1];  %like NFL,OS-RPE
img_gradient=imfilter(img_shift,gradient_mask_dl,'replicate');
img_gradient=img_gradient./max(img_gradient(:)); %Normierung
%% calculate the weight
% k=0; g_min=0.0001; l=zeros(8);weig=[];
% for i=31:1:140
%     for j=1:1:1000
%         l0=1000*(i-1)+j;
%         l(1)=l0+1; l(2)=l0+1001; l(3)=l0+1000; l(4)=l0+999; l(5)=l0-1; l(6)=l0-1001; l(7)=l0-1000; l(8)=l0-999;
%         for h=1:1:8
%             k=k+1;
%             weig(k)=2-(img_gradient(l0)+img_gradient(l(h)))+g_min;
%             X(k)=l0;
%             Y(k)=l(h);
%         end
%     end
% end
%%
 g_min=0.0001;
row1=100; row2=150; cost_matrix=[];
for i=1:1:1002*(row2-row1+1);
    if rem(i,1002)==1 || rem(i,1002)==0
        if 0<i+1<=1002*(row2-row1+1)
            cost_matrix(i,i+1)=g_min;
        end
        if 0<i+1001<=1002*(row2-row1+1)
            cost_matrix(i,i+1001)=g_min;
        end
        if 0<i+1002<=1002*(row2-row1+1)
            cost_matrix(i,i+1002)=g_min;
        end
        if 0<i+1003<=1002*(row2-row1+1)
            cost_matrix(i,i+1003)=g_min;
        end
        if 0<i-1 && i-1<=1002*(row2-row1+1)
            cost_matrix(i,i-1)=g_min;
        end
        if 0<i-1001&& i-1001<=1002*(row2-row1+1)
            cost_matrix(i,i-1001)=g_min;
        end
        if 0<i-1002 && i-1002<=1002*(row2-row1+1)
            cost_matrix(i,i-1002)=g_min;
        end
        if 0<i-1003 && i-1003<=1002*(row2-row1+1)
            cost_matrix(i,i-1003)=g_min;
        end
%     else
%         if 0<i+1<=1002*(row2-row1+1)
%             cost_matrix(i,i+1)=;
%         end
%         if 0<i+1001<=1002*(row2-row1+1)
%             cost_matrix(i,i+1001)=g_min;
%         end
%         if 0<i+1002<=1002*(row2-row1+1)
%             cost_matrix(i,i+1002)=g_min;
%         end
%         if 0<i+1003<=1002*(row2-row1+1)
%             cost_matrix(i,i+1003)=g_min;
%         end
%         if 0<i-1<=1002*(row2-row1+1)
%             cost_matrix(i,i-1)=g_min;
%         end
%         if 0<i-1001<=1002*(row2-row1+1)
%             cost_matrix(i,i-1001)=g_min;
%         end
%         if 0<i-1002<=1002*(row2-row1+1)
%             cost_matrix(i,i-1002)=g_min;
%         end
%         if 0<i-1003<=1002*(row2-row1+1)
%             cost_matrix(i,i-1003)=g_min;
%         end
    end
end
%%
for i=1:1:109
    weig(8000*(i-1)+5)=g_min;
    weig(8000*i-7)=g_min;
end
% X=X-28999;Y=Y-28999;
% DG=sparse(X,Y,weig);
% h = view(biograph(DG,[],'ShowWeights','on'))
% Biograph object with 110000 nodes and 880000 edges.



