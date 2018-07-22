function [xNFL,yNFL,xyinfdl,xyinfy,ggtxnn,ggtynn,ggdlnn,xyinfdl_old,xyinfy_old]=OCT_NFL_line(Lmed,grad_y_punkt)
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
xyinfdl_old=xyinfdl;
xyinfy_old=xyinfy;
ggtxnn=[];
ggtynn=[];
ggdlnn=[];
while sum(sum(xyinfy(:,1:(end-1))))~=0
    ggtx=[];
    ggty=[];
    for hvi=1:(size(xyinfy,2)-1)
        if sum(xyinfy(:,hvi))~=0
            break
        end
    end
    for hv=hvi:(size(xyinfy,2)-1)
        if (min( abs(xyinfy(1,hv)-xyinfy(:,hv+1)))<grad_y_punkt)&(xyinfy(1,hv)~=0)
            vff=1:size(xyinfy,1); vff(abs(xyinfy(1,hv)-xyinfy(:,hv+1))>=grad_y_punkt)=[];
            vff=vff(1);
            xypam=xyinfy(1,hv);
            vff__=1:size(xyinfy,1); vff__(vff)=[];
            xyinfy(1:end,hv+1) = [xyinfy(vff,hv+1);
                xyinfy(vff__,hv+1)];
            xyinfdl(1:end,hv+1)= [xyinfdl(vff,hv+1);
                xyinfdl(vff__,hv+1)];
            xyinfy(1:end,hv)=[xyinfy(2:end,hv);0];
            xyinfdl(1:end,hv)=[xyinfdl(2:end,hv);0];
            ggtx=[ggtx,hv];
            ggty=[ggty,xypam];
        else
            xyinfy(1:end,hv)=[xyinfy(2:end,hv);0];
            xyinfdl(1:end,hv)=[xyinfdl(2:end,hv);0];
            break
        end
    end
    if length(ggty)>10
        ggtxnn(size(ggtxnn,1)+1,1:length(ggtx))=ggtx;
        ggtynn(size(ggtynn,1)+1,1:length(ggty))=ggty;
        ggdlnn=[ggdlnn;[length(ggty) min(ggty)]];
    end
end
ggdlnn_leng=ggdlnn(:,1);
ggdlnn=[(1:size(ggdlnn,1))',ggdlnn];
ggdlnn(:,2)=ggdlnn(:,2)-min(ggdlnn(:,2));
ggdlnn(:,2)=ggdlnn(:,2)./max(ggdlnn(:,2));
ggdlnn_leng(ggdlnn(:,2)<(0.2),:)=[];
ggdlnn(ggdlnn(:,2)<(0.2),:)=[];
for bniewazne=1:(size(ggdlnn,1).^2)
    if size(ggdlnn,1)>=2
        usun_=zeros([1 size(ggdlnn,1)]);
        nr1=ggdlnn(1,1);
        x11=ggtxnn(nr1,:);
        y11=ggtynn(nr1,:);
        x11(y11==0)=[];
        y11(y11==0)=[];
        for nr_=2:size(ggdlnn,1)
            nr2=ggdlnn(nr_,1);
            x22=ggtxnn(nr2,:);
            y22=ggtynn(nr2,:);
            x22(y22==0)=[];
            y22(y22==0)=[];
            for iy=1:length(x11)
                xbn=1:length(x22);
                xbni=xbn(x22==x11(iy));
                if ~isempty(xbni)
                    if y11(iy)<y22(xbni(1))
                        usun_(nr_)=usun_(nr_)+1;
                    end
                end
            end
        end
        if sum(usun_)~=0
            ggdlnn(usun_>(ggdlnn_leng'*0.2),:)=[];
            ggdlnn_leng(usun_>(ggdlnn_leng'*0.2))=[];
            ggdlnn=[ggdlnn(2:end,:);ggdlnn(1,:)];
            ggdlnn_leng=[ggdlnn_leng(2:end);ggdlnn_leng(1,:)];
        else
            ggdlnn=[ggdlnn(2:end,:);ggdlnn(1,:)];
            ggdlnn_leng=[ggdlnn_leng(2:end);ggdlnn_leng(1,:)];
        end
    end
end
ggdlnn_s=sortrows(ggdlnn,-2);
if size(ggdlnn_s,1)==2
    xNFL1=ggtxnn(ggdlnn_s(1,1),:);
    yNFL1=ggtynn(ggdlnn_s(1,1),:);
    xNFL2=ggtxnn(ggdlnn_s(2,1),:);
    yNFL2=ggtynn(ggdlnn_s(2,1),:);
    xNFL1(xNFL1==0)=[];
    yNFL1(yNFL1==0)=[];
    xNFL2(xNFL2==0)=[];
    yNFL2(yNFL2==0)=[];
    yNFL1_poczg=yNFL1(1)+std(yNFL1);
    yNFL1_poczd=yNFL1(1)-std(yNFL1);
    yNFL2_poczg=yNFL2(1)+std(yNFL2);
    yNFL2_poczd=yNFL2(1)-std(yNFL2);
    if min(xNFL1)<min(xNFL2)
        if (abs(yNFL1(end)-yNFL2_poczd)<std(yNFL1))|(abs(yNFL1(end)-yNFL2_poczg)<std(yNFL1))
            xNFL=[xNFL1 xNFL2];
        else
            if length(yNFL1)>length(yNFL2)
                xNFL=[xNFL1];
            else
                xNFL=[xNFL2];
            end
        end
    else
        if (abs(yNFL2(end)- yNFL1_poczd)<std(yNFL2))|(abs(yNFL2(end)- yNFL1_poczg)<std(yNFL2))
           
           
            xNFL=[xNFL2 xNFL1];
        else
            if length(yNFL1)>length(yNFL2)
                xNFL=[xNFL1];
            else
                xNFL=[xNFL2];
            end
        end
    end
else
    xNFL=ggtxnn(ggdlnn_s(1,1),:);
    xNFL(xNFL==0)=[];
end
filtr_med=50;
[xNFL,yNFL]=OCT_NFL_line_end(xNFL,xyinfdl_old,xyinfy_old,grad_y_punkt,filtr_med);
przyci_po_obu_x_proc=0.2;
y_dd=abs(diff(yNFL));
y_dd_lab=bwlabel(y_dd<(grad_y_punkt)/2);
num_1=y_dd_lab(round(length(y_dd_lab)*przyci_po_obu_x_proc));
num_end=y_dd_lab(round(length(y_dd_lab)*(1-przyci_po_obu_x_proc)));
x_sek=1:length(y_dd_lab);
x_sek_1=x_sek(y_dd_lab==num_1); x_sek_1=x_sek_1(1);
x_sek_end=x_sek(y_dd_lab==num_end);
x_sek_end=x_sek_end(end);
xNFL=xNFL(x_sek_1:x_sek_end);
yNFL=yNFL(x_sek_1:x_sek_end);
end