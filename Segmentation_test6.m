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

img=imread('RPD2.jpg');
img=rgb2gray(img); %3 channels to 1 channel
%
[M,N]=size(img);
img=double(img)/255;
%img=mat2gray(img); %[0 1.0]
re=fspecial('gaussian'); 
img=imfilter(img,re);
% figure,imshow(img); title('after gaussian filtering');

%% 找到较大波峰中间的极小值
for i=1:N
each_row=img(:,i);
locs=[];
[~,locs]=findpeaks(each_row,'MINPEAKHEIGH',0.6*max(each_row)); 
% figure,plot(each_row),hold on
% plot(locs,each_row(locs),'ro');hold on

[~,loc_min]=min(each_row(locs(1):locs(end)));
point_low=loc_min+locs(1)-1;
% stem(loc_min+locs(1)-1,each_row(point_low),'b'); hold off
% pause(0.5)
%%该极小值左边和右边最大波峰分别对应NFL 和 RPE
[peaks_left,loc_left]=max(each_row(1:point_low));
[peaks_right,loc_right]=max(each_row(point_low:end));
loc_right=loc_right+point_low-1;
rpe(i)=loc_right;
nfl(i)=loc_left;
end
%%
figure, imshow(img);hold on,
plot(rpe,'r'); hold on
plot(nfl,'b');hold off



%%
i=340;
each_row=img(:,i);
[~,locs]=findpeaks(each_row,'MINPEAKHEIGH',0.6*max(each_row)); 
figure,plot(each_row),hold on
plot(locs,each_row(locs),'ro');hold on

[~,loc_min]=min(each_row(locs(1):locs(end)));
point_low=loc_min+locs(1)-1;
stem(loc_min+locs(1)-1,each_row(point_low),'b'); hold off
pause(0.5)
%%该极小值左边和右边最大波峰分别对应NFL 和 RPE
[peaks_left,loc_left]=max(each_row(1:point_low));
[peaks_right,loc_right]=max(each_row(point_low:end));
d=diff(each_row);
figure,plot(d)
%%
for i=1:1:N
    max1(i)=find(img(:,i)==max(img(:,i)),1,'first');
    max2(i)=find(img(:,i)==max(img(:,i)),1,'last');
end
figure, plot(max1),hold on,
plot(max2,'r');

