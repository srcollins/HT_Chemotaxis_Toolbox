function [r,p] = myNanCorrcoef(x,y)
%computes the r value for correlation of x to y. It can also
%return a p-value for the correlation: [r,p]=ergcorrcoef(x,y)

x=x(:); y=y(:);
ind = (~isnan(x))&(~isnan(y));
x=x(ind);
x=x(:);
y=y(ind);
y=y(:);
if (length(x)>1) & (max(x)~=min(x)) & (max(y)~=min(y))  
    [temp, ptemp]=corrcoef(x,y);
    r=temp(1,2);
    p=ptemp(1,2);
else
    r=NaN;
    p=NaN;
end