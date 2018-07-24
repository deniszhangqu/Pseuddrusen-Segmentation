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
img_med=medfilt2(img,[5 5]);
img_med=mat2gray(img_med);
figure,imshow(img);
% figure,imshow(img_med);

%% RPE segmentation using regiongrow
Ts=max(img(:));
[g,NR,SI,TI]=regiongrow(img,Ts,0.4);
figure,imshow(g,[]);title('after Regiongrow');
f=zeros(M,N);
for m=1:1:M
    for n=1:1:N
        if g(m,n)==1
            f(m,n)=1;
        end
    end
end
% figure,imshow(f,[]);title('after region_choise');
se=strel('disk',3);
f1=imclose(f,se);
figure,imshow(f1);title('after close operation');
%% calculate the boundaries
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

for i=1:1:N
    if isnan(y_rpe(i))
        y_ez(i)=NaN;
        y_ch(i)=NaN;
    else
        y_ez(i)=find(img_w(:,i)',1,'first');   % curve of EZ
        y_ch(i)=find(img_w(:,i)',1,'last');    % curve of Choroide
    end
end
figure,imshow(img); hold on,
plot(y_rpe,'r'); hold on,
plot(y_ez,'g'); hold on,
plot(y_ch,'b'); hold off
legend('y\_rpe','y\_ez','y\_ch')

y_rpe=smooth(y_rpe);
y_ez=smooth(y_ez);
y_ch=smooth(y_ch);

figure,imshow(img); hold on,
plot(y_rpe,'r'); hold on,
plot(y_ez,'g'); hold on,
plot(y_ch,'b'); hold off
legend('y\_rpe','y\_ez','y\_ch')
%% detection of Drusen
[pks,locs] = findpeaks(-y_ez);

j=1; p_drusen=[];
y_rpe=uint8(y_rpe);
y_ez=uint8(y_ez);
y_ch=uint8(y_ch);

for i=1:1:size(locs)
    k1=[];
    k2=[];
    m=locs(i);
    if m<7||m>N-6 % the punkt by edge willn't be considered
        continue;
    end
    k1=y_ez(m:m+5) < y_ez(m+1:m+6);
    if sum(k1)>=3
        k2=y_ez(m-5:m) < y_ez(m-6:m-1);
        if sum(k2)>=3
            p_drusen(j)=locs(i);
            j=j+1;
        else
            continue;
        end
    else
        continue;
    end
end
figure,imshow(img); hold on,
plot(y_rpe,'r'); hold on,
plot(y_ez,'g'); hold on,
plot(y_ch,'b'); hold on
plot(p_drusen,y_ez(p_drusen),'o');hold on,
%% detection of RPD
j=1;k=1;p_rpd=[];p_normaldrusen=[];
for i=1:1:size(p_drusen,2)
    k1=[];
    k2=[];
    m=p_drusen(i);
    k1=y_ch(m:m+5) < y_ch(m+1:m+6);
    k2=y_ch(m-5:m) < y_ch(m-6:m-1);
    if (sum(k1)<3) && (sum(k2)<3)
        p_rpd(j)=p_drusen(i);
        j=j+1;
    else
        p_normaldrusen(k)=p_drusen(i);
        k=k+1;
    end
end
%%
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


