function [y_rpe,y_ez,y_ch]=findboundary(img_retina,T)

%%[y_rpe,y_ez,y_ch]=findboundary(img_retina,T)
%this is a function to detection the bounddary of choroid, RPE and EZ in outer retina 
% %input:
% %img_retina: is a binar image, in which the background is 0, the outer
%retina is 0
% %T: is the region number. sometimes the result of regiongrowing have more
%than 1 region. suggest to set T=1
% %outputs:
% %y_rpe:the boundary between RPE and EZ
% %y_ch:the boundary between Choroid and RPE
% %y_ez:the boundary between EZ and ELM
%  ----Qu Zhang, 2018/07/23,  TUI----
[M,N]=size(img_retina);
for n=1:1:N
    for m=1:1:M
        if img_retina(m,n)==T
            img_w(m,n)=m;
        else
            img_w(m,n)=0;
        end
    end
    sum_c(n)=sum(img_w(:,n));
    sum_bin(n)=sum(img_retina(:,n));
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
y_rpe=smooth(y_rpe);
y_ez=smooth(y_ez);
y_ch=smooth(y_ch);
end
