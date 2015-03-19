function bool=boolRegExp(cell_arr,pattern)
%Returns a boolean vector or matrix indicating whether or not the elements
%of cell_arr contain pattern

if ischar(cell_arr)
    bool=length(regexp(cell_arr,pattern))>0;
else  % a cell array was given
    bool=cellfun(@(x) length(x),regexp(cell_arr,pattern))>0;
end