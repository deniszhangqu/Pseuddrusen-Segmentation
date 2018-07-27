function [img_shift,y_rpe,shift_int]=img_rpe_shift(img,K)

[M,N]=size(img); X=zeros(N); Y=zeros(N);
max_column=max(img);
for i=1:1:N
    img_max(:,i)=img(:,i)>=0.9*max_column(i);
end

[X,Y]=find(img_max);
P=polyfit(Y,X,K);
y_rpe=polyval(P,Y);

%%circshift the OCT image to make the RPE geradeaus
k=median(y_rpe);
shift_int=int8(k-y_rpe);
img_shift=zeros(M,N);
for i=1:1:N
    img_shift(:,i)=circshift(img(:,i),shift_int(i));
end

end