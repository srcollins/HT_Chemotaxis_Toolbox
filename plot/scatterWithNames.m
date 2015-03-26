function figure1=scatterWithNames(v1,v2,point_labels,varargin)
% This function generates a scatter plot of the data in row1 of
% struct1.data against the data in row2 of struct2.data

%% check first three arguments
if nargin<3, 
    error('need at least three arguments'); 
end

if ~isa(v1,'numeric') || ~isa(v2,'numeric') 
    error('v1 and v2 need to be numeric'); end

if numel(v1) ~= numel(v2)
    error('v1 and v2 must be same length!'); 
end

if numel(v1) ~= numel(point_labels)
    error('data and label vectors must be same length!'); 
end
v1=v1(:);
v2=v2(:);

%% define defaults values for options
opt.xlabel='';
opt.ylabel='';
opt.title='';
opt.corrs={'Pearson','Cosind','Spearman'};
opt.axisColor=[0.25 0.25 0.25];
opt.diag=1;
opt.robustTrend=0;
opt.xGrayRegion=0;
opt.yGrayRegion=0;
opt.grayRegionColor=[0.9 0.9 0.9];
opt.markerColor=[0 0 1];

%% update input arguments based on varargin
if mod(length(varargin),2)~=0
    error('options should always be in pairs')
end
for i=1:2:length(varargin)
    o=varargin{i}; % o is what optio to change
    v=varargin{i+1}; % v is what is its value
    
    % o must be a char
    if ~ischar(o)
        error('Option must be a char!');
    end
    
    % check if option is an allowed one
    if ~ismember(o,fieldnames(opt))
        error('%s is not a valid option!',o)
    end
      
    opt.(o)=v;
end

%% Generate the figure
if sum(~isnan(v1+v2))>0     % Check to make sure there is at least one data point to plot
    figure1 = figure;
%    figure1=figure('PaperSize',[4 6]);
    
    % Choose the point size for the scatter plot
%     if sum(~isnan(v1+v2))>1600
%         pointsize=1;
%     else
%         pointsize=3;
%     end
    pointsize=11;
    
    % Scatter the data
    minVal=min([v1(:); v2(:)])-0.1;
    maxVal=max([v1(:); v2(:)])+0.1;
    ax=gca; hold on;
    if opt.yGrayRegion>0
        patch([minVal minVal maxVal maxVal],[-1*opt.yGrayRegion opt.yGrayRegion opt.yGrayRegion -1*opt.yGrayRegion],opt.grayRegionColor,'LineStyle','none'); 
    end
    if opt.xGrayRegion>0
        patch([-1*opt.xGrayRegion opt.xGrayRegion opt.xGrayRegion -1*opt.xGrayRegion],[minVal minVal maxVal maxVal],opt.grayRegionColor,'LineStyle','none'); 
    end
    if opt.yGrayRegion>0 && opt.xGrayRegion>0
        patch([-1*opt.xGrayRegion opt.xGrayRegion opt.xGrayRegion -1*opt.xGrayRegion],[opt.yGrayRegion opt.yGrayRegion -1*opt.yGrayRegion -1*opt.yGrayRegion],1-2*(1-opt.grayRegionColor),'LineStyle','none'); 
    end
    if opt.robustTrend>0
        p=robustfit(v1,v2);
        plot([minVal maxVal],p(1)+p(2)*[minVal maxVal],'Color',opt.axisColor,'HitTest','off');
    end
    plot(ax,[minVal maxVal],[0 0],'Color',opt.axisColor,'HitTest','off','LineWidth',1);
    plot(ax,[0 0],[minVal maxVal],'Color',opt.axisColor,'HitTest','off','LineWidth',1);
    if opt.diag>0
        plot(ax,[minVal maxVal],[minVal maxVal],'Color',opt.axisColor,'HitTest','off');
    end
    scatter(ax,v1,v2,'Marker','o','MarkerEdgeColor',0.6*opt.markerColor,'MarkerFaceColor',opt.markerColor,'SizeData',pointsize);
    xlim([minVal maxVal]);
    ylim([minVal maxVal]);
    hold off;
    
    % Set up the labeling of the data points in the plot
    dcm_obj = datacursormode(figure1);
    set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on');
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,point_labels})

    % Label the axes
    xlabel(opt.xlabel,'FontSize',10);
    ylabel(opt.ylabel,'FontSize',10);
    
    % Add the title
    title(opt.title);
    
    % Print info to the main window
    if ismember('Pearson',opt.corrs)
        fprintf('Pearson corr = %.3f\n',myNanCorrcoef(v1,v2));
    end
    if ismember('Cosine',opt.corrs)
        fprintf('Cosine corr = %.3f\n',myNanCosine(v1,v2));
    end
    if ismember('Spearman',opt.corrs)
        fprintf('Spearman corr = %.3f\n',myNanSpearman(v1,v2));
    end
    
else
    fprintf('No data to compare.\n');
end

end % of scatterWithNames function

%--------------------------------------------------------------------------
function txt = myupdatefcn(empt,event_obj,labels)

% pos = get(event_obj,'Position');  % you could use the x-y coordinates of
                                    % the data point if it was useful
a = get(event_obj,'DataIndex');
txt = labels(a);
return

end % of myupdatefcn
