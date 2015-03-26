function plotCDFs(data)
%overlays cumulative distribution plots for the provided data

figure; hold on;
if iscell(data)
    data=data(:);
    for i=1:length(data)
        d=data{i}; d=d(:);
        cdfplot(d);
        f=(i-1)/(length(data)-1);
        set(findobj('Type','line','Color',[0 0 1]),'Color',[min(1.5*f,1) max(0,f-0.33) 0]);
    end
end
if isnumeric(data)
    [r c]=size(data);
    if r==1
        data=data';
        [r c]=size(data);
    end
    for i=1:c
        d=data(:,i);
        cdfplot(d);
        f=(i-1)/(c-1);
        set(findobj('Type','line','Color',[0 0 1]),'Color',[min(1.5*f,1) max(0,f-0.33) 0]);
    end
end
hold off;
