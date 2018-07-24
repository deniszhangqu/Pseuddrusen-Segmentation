function [locs,p_drusen,p_normaldrusen,p_rpd]=finddrusen(y_ez,y_ch,M,T)
  %%this function is to detection the drusen and classify them to RPD or normal drusen
  % % input:
  %   y_ez: the boundary between EZ and ELM
  %   y_ch: the boundary between RPE and Choriod
  %   M: how many punkts on the left or on the right side of one local peak
  %      on the y_ez to culculate their difference.
  %   T: the number of  acsends and decends at least on the left and right 
  %     side of this peaks if this peak  is a drusen
  % %outputs:
  %   locs:             the x-axel vector of peaks;
  %   p_drusen:         the x-axel vector of drusen;
  %   p_normaldrusen:   the x-axel vector of normal drusen;
  %   p_rpd:            the x-axel vector of RPD;
  %
  % -------Qu Zhang,   2018/07/23,   TUI------
     
[~,locs] = findpeaks(-y_ez);
N=numel(y_ez);
j=1; locs_drusen=[];
y_ez=uint8(y_ez);
y_ch=uint8(y_ch);

for i=1:1:size(locs)
    k1=[];
    k2=[];
    m=locs(i);
    if m<=M||m>=N-M % the punkt by edge willn't be considered
        continue;
    end
    k1=y_ez(m:m+M-1) < y_ez(m+1:m+M);
    if sum(k1)>=T
        k2=y_ez(m-M+1:m) < y_ez(m-M:m-1);
        if sum(k2)>=T
            locs_drusen(j)=locs(i);
            j=j+1;
        else
            continue;
        end
    else
        continue;
    end
end
j=1;k=1;locs_rpd=[];locs_normal_drusen=[];
for i=1:1:size(locs_drusen,2)
    k1=[];
    k2=[];
    m=locs_drusen(i);
    k1=y_ch(m:m+M-1) < y_ch(m+1:m+M);
    k2=y_ch(m-M+1:m) < y_ch(m-M:m-1);
    if (sum(k1)<T) && (sum(k2)<T)
        locs_rpd(j)=locs_drusen(i);
        j=j+1;
    else
        locs_normal_drusen(k)=locs_drusen(i);
        k=k+1;
    end
end
p_drusen=locs_drusen;
p_normaldrusen=locs_normal_drusen;
p_rpd=locs_rpd;
end
