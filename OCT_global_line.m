function [x,yrpes,dxx,dyy]=OCT_global_line(Lmed)
x=(1:size(Lmed,2))';
yyy=(1:size(Lmed,1))';
yrpe=[];
Lbinrpe=zeros(size(Lmed));
for ik=1:size(Lmed,2)
    xx_best=[];
    Llabp=bwlabel(Lmed(:,ik)>(max(Lmed(:,ik))*0.9));
    Lbinrpe(:,ik)=Llabp;
    for tt=1:max(Llabp)
        xxl=yyy(Llabp==tt);
        xx_best=[xx_best;mean(xxl)] ;
    end
    if ~isempty(xx_best)
        yrpe(ik)=max(xx_best);
    else
        yrpe(ik)=0;
    end
end
figure; imshow(mat2gray(Lbinrpe*0.5+Lmed));hold on;
plot(yrpe,'r*-')
yg=gradient(yrpe);
ygg=ones([1 length(yrpe)]); ygg(abs(yg)>20)=0;
ygl=bwlabel(ygg);
figure; imshow(mat2gray(Lbinrpe*0.5+Lmed));hold on;
palett=jet(max(ygl));
for iiih=1:max(ygl(:))
    plot(x(ygl==iiih),yrpe(ygl==iiih),'Color',palett(iiih,:),'LineWidth',4);
    
end
pam_dl=[];
figure; imshow(mat2gray(Lbinrpe*0.5+Lmed)); hold on
for iiik=1:max(ygl(:))
    for iiikk=iiik:max(ygl(:))
        if iiik<=iiikk
            ygk=[yrpe(ygl==iiik),yrpe(ygl==iiikk)];
            xgk=[x(ygl==iiik);x(ygl==iiikk)];
        else
            ygk=[yrpe(ygl==iiikk),yrpe(ygl==iiik)];
            xgk=[x(ygl==iiikk);x(ygl==iiik)];
        end
        if length(ygk)>10
            P = POLYFIT(xgk',ygk,2); yrpes = round(polyval(P,x));
           
            plot(yrpes,'g*-')
            pam_dl=[pam_dl;[iiik iiikk sum(abs(yrpeyrpes')<20)]];
        end
    end
end
pam_s=sortrows(pam_dl,-3);
if size(pam_s,1)==1
    ygk=[yrpe(ygl==pam_s(1,1))];
    xgk=[x(ygl==pam_s(1,1))];
else
    ygk=[yrpe(ygl==pam_s(1,1)),yrpe(ygl==pam_s(1,2))];
    xgk=[x(ygl==pam_s(1,1));x(ygl==pam_s(1,2))];
end
P = POLYFIT(xgk',ygk,2); yrpes = round(polyval(P,x));
plot(x,yrpes,'w*-');
yrpe=yrpe(:);
plot(x,yrpe,'m*-');
dx=x; dx(abs(yrpe-yrpes)>20)=[];
yrpe(abs(yrpe-yrpes)>20)=[];
dxl=bwlabel(diff(dx)<125);
pdxl=[];
for qw=1:max(dxl)
    pdxl=[pdxl;[qw, sum(dxl==qw)]];
end
pdxl(pdxl(:,2)<50,:)=[];
dxx=[]; dyy=[];
for wq=1:size(pdxl,1)
    dxx=[dxx; dx(dxl==pdxl(wq,1))];
    dyy=[dyy; yrpe(dxl==pdxl(wq,1))];
end
end