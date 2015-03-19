function tokens=getTokens(str,expr)
%this function assumes a single match for each element of a cell array. If
%there are multiple, it returns the first one.

tok1=regexp(str, expr, 'tokens');

if ischar(str) %a single string
    %tokens=tok1{1};
    tokens=cellfun(@(x) x{1},tok1,'UniformOutput',false);
else  %a cell array of strings
    %sizes=cellfun(@(x) length(x(:)),tok1);
    tokens=cellfun(@helper, tok1);
end

end

%-------------------------------------------------------------------

function y=helper(x)

if ~isempty(x)
    y=x{1};
else
    y={[]};
end

end