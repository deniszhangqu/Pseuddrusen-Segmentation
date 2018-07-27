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
img_gradient=imfilter(img,gradient_mask_dl,'replicate');
img_gradient=img_gradient./max(img_gradient(:)); %Normierung


%%
% X=X-28999;Y=Y-28999;
% DG=sparse(X,Y,weig);
% h = view(biograph(DG,[],'ShowWeights','on'))
% Biograph object with 110000 nodes and 880000 edges.



