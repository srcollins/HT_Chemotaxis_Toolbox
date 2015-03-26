function plotFrequencyDistributions(data, bins, colors,figFlag)
%computes histograms for the provided data over the bins, normalizes the
%histograms, and plots the result. figFlag = 1 (default) indicates that a
%new figure should be generated

if nargin<4
    figFlag=1;
end

h=zeros(length(data),length(bins));
f=zeros(length(data),length(bins));
for i=1:length(data)
    d=data{i}; d=d(:);
    h(i,:)=hist(d,bins);
    f(i,:)=h(i,:)/sum(h(i,:));
end

if figFlag>0
    figure;
end
if length(bins)>4
    if nargin<3
        plot(bins,f');
    else
        hold on;
        for i=1:length(data)
            plot(bins,f(i,:),'Color',colors{i},'LineWidth',2);
        end
        hold off;
    end
else
    bar(bins,f');
end
