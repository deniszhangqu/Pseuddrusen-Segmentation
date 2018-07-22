function [im_bin,y_rpe]=RPE_colummax(im,tf)

%%this function bases on the global maximum intensity of each colum in a
%OCT B-sccan image.
%im: input image, which is a M*N gray image;
%tf: the threshodingsfactor for the RPE. like tf=0.9 means the pixels with
%intensity in a colum more than 0.9*maximum intensity of this colum are
%choise and consist on the RPE
%im_bin: output binar-image, the region with 1 means the RPE
%y_rpe:the rpe boundary will be als the center linie of the segmented RPE
%*author: Qu Zhang   2018/07/21  Tu Ilmenau*

[M,N]=size(im);
max_colum=max(im);
T=max_colum*tf;
im_bin=zeros(M,N);
for i=1:1:N
    v=im(:,i);
    w=v>T(i);
    im_bin(:,i)=w(:);
end
figure,imshow(im_bin);

for n=1:1:N
    for m=1:1:M
        if im_bin(m,n)==1
            im_w(m,n)=m;
        else
            im_w(m,n)=0;
        end
    end
    sum_c(n)=sum(im_w(:,n));
    sum_bin(n)=sum(im_bin(:,n));
    y_rpe(n)=sum_c(n)/sum_bin(n);
end
end