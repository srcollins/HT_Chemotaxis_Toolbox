function r = myNanCosine(x,y)
%computes the cosine correlation of x to y.

x=x(:); y=y(:);
ind = (~isnan(x))&(~isnan(y));
x=x(ind);
x=x(:);
y=y(ind);
y=y(:);
if (length(x)>1)% & (max(x)~=min(x)) & (max(y)~=min(y))  
    r=(x'*y)/norm(x)/norm(y);
else
    r=NaN;
end