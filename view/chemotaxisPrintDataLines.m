function chemotaxisPrintDataLines(dat,genePat,cols)

if nargin<3
    cols=1:size(dat.s1,2);
end
ind=find(boolRegExp(upper(dat.rowlabels),upper(genePat)));
sDflag=ismember('sD',fieldnames(dat));

if ~isempty(ind)
    if sDflag
        fprintf('%-30s %6s %6s %6s %6s %6s\n','Gene','s0','s1','sD','c1','a1');
    else
        fprintf('%-30s %6s %6s %6s %6s\n','Gene','s0','s1','c1','a1');
    end
    for j=ind(:)'
        if sDflag
            fprintf('<a href="matlab:firefoxNCBIgenePage(''%s'')">%-30s</a> %6.2f %6.2f %6.2f %6.2f %6.2f\n',regexprep(dat.rowlabels{j},' \([A-Za-z\# 0-9]+\)',''),dat.rowlabels{j},nanmean(dat.s0(j,cols),2), ...
                nanmean(dat.s1(j,cols),2),nanmean(dat.sD(j,cols),2),nanmean(dat.c1(j,cols),2),nanmean(dat.a1(j,cols),2));
        else
            fprintf('<a href="matlab:firefoxNCBIgenePage(''%s'')">%-30s</a> %6.2f %6.2f %6.2f %6.2f\n',regexprep(dat.rowlabels{j},' \([A-Za-z\# 0-9]+\)',''),dat.rowlabels{j},nanmean(dat.s0(j,cols),2), ...
                nanmean(dat.s1(j,cols),2),nanmean(dat.c1(j,cols),2),nanmean(dat.a1(j,cols),2));
        end
    end
else
    fprintf('\nNo genes match the given pattern.\n');
end