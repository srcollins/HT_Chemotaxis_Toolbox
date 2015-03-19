function lines=getNextChunkOfLines(fid,num)
%Read in a chunk of lines from a file. If possible, it will read in num
%lines. If there aren't that many in the file, it will return as many as it
%can.

if nargin<2
    num=10000;
end

lines=cell(num,1);
i=0;
while i<num && feof(fid)==0
    i=i+1;
    lines{i}=fgetl(fid);
end

lines=lines(1:i);
