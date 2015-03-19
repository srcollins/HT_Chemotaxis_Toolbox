function rowlabels=readChemotaxisCoordMap(path)

if endsWith(path,'.csv') || endsWith(path,'.txt')
    fname=path;
else
    fname=[path filesep 'cmap.txt'];
end
fid=fopen(fname,'r');
if fid>0
    lines=readAllLines(fid);
    fclose(fid);
    if endsWith(fname,'.csv')
        symbol=getTokens(lines(2:end),'[A-Z],[0-9]+,([^,]+),[^\n,]+');
        treatment=getTokens(lines(2:end),'[A-Z],[0-9]+,[^,]+,([^\n,]+)');
        rowlabels=cellfun(@(x,y) {[x ' (' y ')']},symbol,treatment);
    else
        rowlabels=getTokens(lines(2:end),'\t([^\n]+)');
    end
else
    rowlabels={};
end