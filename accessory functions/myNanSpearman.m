function [r,p]=myNanSpearman(x,y)
%computes the Spearman rank correlation ignoring NaN values

x=x(:); y=y(:);
ind=~isnan(x+y);
x=x(ind);
y=y(ind);

if length(x)>1
    [r,p]=corr(x,y,'type','Spearman');
else
    r=NaN;
    p=NaN;
end
