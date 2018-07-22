function [xNFL,yNFL]=OCT_NFL_line_end(xNFL_old,xyinfdl,xyinfy,grad_y_punkt,filtr_med)
x_start=xNFL_old(round(end/2));
xNFL=[];
yNFL=[];
xyinfy(1,:)=medfilt2(xyinfy(1,:),[1 filtr_med]);
for hv=x_start:(size(xyinfy,2)-1)
    if (min( abs(xyinfy(1,hv)-xyinfy(:,hv+1)))<grad_y_punkt)&(xyinfy(1,hv)~=0)
       
        vff=1:size(xyinfy,1); vff(abs(xyinfy(1,hv)-xyinfy(:,hv+1))>=grad_y_punkt)=[]; vff=vff(1);
        
        xypam=xyinfy(1,hv);
        vff__=1:size(xyinfy,1); vff__(vff)=[];
        xyinfy(1:end,hv+1) = [xyinfy(vff,hv+1);
            xyinfy(vff__,hv+1)];
        xyinfdl(1:end,hv+1)= [xyinfdl(vff,hv+1);
            xyinfdl(vff__,hv+1)];
        xyinfy(1:end,hv)=[xyinfy(2:end,hv);0];
        xyinfdl(1:end,hv)=[xyinfdl(2:end,hv);0];
        xNFL=[xNFL;hv];
        yNFL=[yNFL;xypam];
    else
        xyinfy(1:end,hv)=[xyinfy(2:end,hv);0];
        xyinfdl(1:end,hv)=[xyinfdl(2:end,hv);0];
        break
    end
end
for hv=(x_start-1):-1:2
    if (min( abs(xyinfy(1,hv)-xyinfy(:,hv-1)))<grad_y_punkt)&(xyinfy(1,hv)~=0)
        
        vff=1:size(xyinfy,1); vff(abs(xyinfy(1,hv)-xyinfy(:,hv-1))>=grad_y_punkt)=[]; vff=vff(1);
        
        xypam=xyinfy(1,hv);
        vff__=1:size(xyinfy,1); vff__(vff)=[];
        xyinfy(1:end,hv-1) = [xyinfy(vff,hv-1);
            xyinfy(vff__,hv-1)];
        xyinfdl(1:end,hv-1)= [xyinfdl(vff,hv-1);
            xyinfdl(vff__,hv-1)];
        xyinfy(1:end,hv)=[xyinfy(2:end,hv);0];
        xyinfdl(1:end,hv)=[xyinfdl(2:end,hv);0];
        xNFL=[hv;xNFL];
        yNFL=[xypam;yNFL];
    else
        xyinfy(1:end,hv)=[xyinfy(2:end,hv);0];
        xyinfdl(1:end,hv)=[xyinfdl(2:end,hv);0];
        break
    end
end
xNFL=round(xNFL);
yNFL=round(yNFL);
end