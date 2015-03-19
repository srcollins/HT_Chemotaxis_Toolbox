function tokens=split(str,pattern)
%splits a string into pieces divided by the provided pattern

tokens=regexp(str,pattern,'split');

