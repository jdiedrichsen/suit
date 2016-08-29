function [x,t]=minORmax(v)
% returns the minimum or maximum of a vector, whatever is larger in absolute value
%function [x,t]=minORmax(v)
% x is the min/max of the vector
% t is the time of occurance of this value (when it occurs multiple time, the first
% occurance)
% if v is a matrix, the operation is performed on each column of v
if(isempty(v))
    x=NaN;
    t=NaN;
    return;
end;
[vmin,tmin]=min(v);[vmax,tmax]=max(v);
indmin=find(vmin<-vmax);
indmax=find(vmin>=-vmax);
x(indmin)=vmin(indmin);
x(indmax)=vmax(indmax);
x(find(isnan(vmin) | isnan(vmax)))=NaN;
t(indmin)=tmin(indmin);
t(indmax)=tmax(indmax);
t(find(isnan(tmin) | isnan(tmax)))=NaN;
