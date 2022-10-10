function [SE]=sph(r)
% This function gives back a circle image with radius equal to r
% We are using this function for anlaytical validation of permeabilities
s=ceil(2*r+1); SE=zeros(s,s);
for I=1:s
    for J=1:s
        D=sqrt((I-(r+1))^2+(J-(r+1))^2);
        if D<=r+.5; SE(I,J)=1; end;
    end
end
t=zeros(size(SE)+4); t(3:end-2,3:end-2)=SE; SE=t;
end