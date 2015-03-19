function lines=readAllLines(fid)

lines={};
while feof(fid)==0
    lines=[lines; getNextChunkOfLines(fid)];
end
