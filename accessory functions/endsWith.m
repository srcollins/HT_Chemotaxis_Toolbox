function bool = endsWith(str1,str2)
%returns 1 if str1 ends with str2, returns 0 otherwise. Can take input
%arguments of type string or type cell. Syntax is bool =
%endswith(str1,str2)
%
%written by Sean Collins (2006) as part of the EMAP toolbox

if ischar(str1)
    bool=singleEndsWith(str1,str2);
else
    len=length(str1);
    bool=zeros(len,1);
    for i=1:len
        bool(i)=singleEndsWith(str1(i),str2);
    end
end

end

%---------------------------------------------------------------------
function bool=singleEndsWith(str1,str2)
str1=char(str1);
str2=char(str2);

len1=length(str1);
len2=length(str2);
bool=0;
if len2<=len1
    if strcmp(str1((len1-len2+1):len1),str2)
        bool=1;
    end
end
if len2==0
    bool=1;
end

end