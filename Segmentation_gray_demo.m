close all;
clear all;
%Quelle:image processing of OCT using Matlab%
%load('D:\GoogleDrive\Masterarbeit-Matlab\Dataset\Duke\269AMD\Farsiu_Ophthalmology_2013_AMD_Subject_1001.mat');
%N=50; %the number of N-th sectional Image
%Img=images(:,:,N);
Im=imread('RPD1.png');
Img=rgb2gray(Im);
Img_gray=Img(1:327,:);
Img_gray=double(Img_gray)/255; Img_org=Img_gray;
%%
Img_med=medfilt2(Img_org,[5 5]);
Img_med=mat2gray(Img_med);
figure;
imshow(Img_med);
n=(1:size(Img_med,2))';
m=(1:size(Img_med,1))';
y_rpe=[];
Img_bin_pre=zeros(size(Img_med));
for ik=1:size(Img_med,2)
    xx_best=[];
    Img_labp=bwlabel(Img_med(:,ik)>(max(Img_med(:,ik))*0.9));
    Img_bin_pre(:,ik)=Img_labp;
    for tt=1:max(Img_labp)
        xxl=m(Img_labp==tt);
        xx_best=[xx_best;mean(xxl)] ;
    end
    if ~isempty(xx_best)
        y_rpe(ik)=max(xx_best);
    else
        y_rpe(ik)=0;
    end
end
figure; imshow(mat2gray(Img_bin_pre*0.5+Img_med));hold on;
plot(y_rpe,'r*-')
figure; imshow(Img_org);hold on;
plot(y_rpe,'r*-')
%%
yg=gradient(y_rpe);
ygg=ones([1 length(y_rpe)]); ygg(abs(yg)>20)=0;
ygl=bwlabel(ygg);
figure; imshow(mat2gray(Img_bin_pre*0.5+Img_med));hold on;
palett=jet(max(ygl));
for iiih=1:max(ygl(:))
    plot(n(ygl==iiih),y_rpe(ygl==iiih),'Color',palett(iiih,:),'LineWidth',4);
end
pam_dl=[];
figure; imshow(mat2gray(Img_bin_pre*0.5+Img_med)); hold on
for iiik=1:max(ygl(:))
    for iiikk=iiik:max(ygl(:))
        if iiik<=iiikk
            ygk=[y_rpe(ygl==iiik),y_rpe(ygl==iiikk)];
            xgk=[n(ygl==iiik);n(ygl==iiikk)];
        else
            ygk=[y_rpe(ygl==iiikk),y_rpe(ygl==iiik)];
            xgk=[n(ygl==iiikk);n(ygl==iiik)];
        end
        if length(ygk)>10
            P = polyfit(xgk',ygk,2); yrpes =round(polyval(P,n));
            plot(yrpes,'g*-')
            pam_dl=[pam_dl;[iiik iiikk sum(abs(y_rpe-yrpes')<20)]];
        end
    end
end
%%
pam_s=sortrows(pam_dl,-3);
if size(pam_s,1)==1
    ygk=[y_rpe(ygl==pam_s(1,1))];
    xgk=[n(ygl==pam_s(1,1))];
else
    ygk=[y_rpe(ygl==pam_s(1,1)),y_rpe(ygl==pam_s(1,2))];
    xgk=[n(ygl==pam_s(1,1));n(ygl==pam_s(1,2))];
end
P = polyfit(xgk',ygk,2); yrpes = round(polyval(P,n));
plot(n,yrpes,'w*-');
y_rpe=y_rpe(:);
plot(n,y_rpe,'m*-');
%%
dx=n; dx(abs(y_rpe-yrpes)>20)=[];
y_rpe(abs(y_rpe-yrpes)>20)=[];
dxl=bwlabel(diff(dx)<125);
pdxl=[];
for qw=1:max(dxl)
    pdxl=[pdxl;[qw, sum(dxl==qw)]];
end
pdxl(pdxl(:,2)<50,:)=[];
dxx=[]; dyy=[];
for wq=1:size(pdxl,1)
    dxx=[dxx; dx(dxl==pdxl(wq,1))];
    dyy=[dyy; y_rpe(dxl==pdxl(wq,1))];
end
dx=dxx; y_rpe=dyy;
plot(dx,y_rpe,'c*-');
figure
imshow(Img_gray); hold on
plot(dx,y_rpe,'c*-');
%% NFL segmentation

L1=rand([201 200]);
xx=-1:0.01:1;
y=gauss(xx+0.5,0.2)+0.5*gauss(xx-0.1,0.05);
Lmed=y'*ones([1 200]);
Lmed=mat2gray(Lmed);
Lmed(:,50:100)=Lmed(:,50:100)*.2;
Lmed = imnoise(Lmed,'gaussian',0.02);
Lmed=medfilt2(Lmed,[3 3]);
figure; imshow(Lmed); hold on
xyinfy=[];
xyinfdl=[];
for ik=1:size(Lmed,2)
    grL1=Lmed(:,ik)>(max(Lmed(:,ik))*0.1);
    lgrL1=bwlabel(grL1);
    for jju=1:max(lgrL1)
        xyinfdl(jju,ik)=sum(lgrL1==jju);
        cuu=1:length(lgrL1);
        cuu(lgrL1~=jju)=[];
        xyinfy(jju,ik)=cuu(1);
        plot(ik,cuu(1),'b*')
    end
end