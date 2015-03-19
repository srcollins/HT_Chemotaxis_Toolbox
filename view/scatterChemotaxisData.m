function scatterChemotaxisData(m,field,col1,col2,rows,varargin)

if nargin>=4 && isnumeric(col1) && isnumeric(col2) && ~isempty(col1) && ~isempty(col2) ...
        && (ischar(field) || (iscell(field) && length(field)==1))
    % Scatter data for one field from different data columns
    if min(col1)>=1 && min(col2)>=1 && max(col1)<=length(m.collabels) && max(col2)<=length(m.collabels)
        field=char(field);
        if ismember(field,fieldnames(m))
            if nargin<5 || isempty(rows)
                rows=1:length(m.rowlabels);
            end
            dat1=nanmean(m.(field)(rows,col1),2);
            dat2=nanmean(m.(field)(rows,col2),2);
            scatterWithNames(dat1,dat2,m.rowlabels(rows),'corrs',{'Pearson','Spearman','Cosine'},varargin{:});
        else
            fprintf('The specified data field (%s) does not exist in the input variable.\n',field);
        end
    else
        fprintf('Column indices exceed bounds of data matrices.\n');
    end
elseif nargin>=3 && (iscell(field) && length(field)==2)
    % Scatter data for two different data fields
    if min(col1)>=1 && max(col1)<=length(m.collabels)
        if ismember(field{1},fieldnames(m)) && ismember(field{2},fieldnames(m))
            if nargin<4 || isempty(rows)
                rows=1:length(m.rowlabels);
            end
            dat1=nanmean(m.(field{1})(rows,col1),2);
            dat2=nanmean(m.(field{2})(rows,col1),2);
            scatterWithNames(dat1,dat2,m.rowlabels(rows),'corrs',{'Pearson','Spearman','Cosine'},varargin{:});
        else
            if ~ismember(field{1},fieldnames(m))
                field=field{1};
            else
                field=field{2};
            end
            fprintf('The specified data field (%s) does not exist in the input variable.\n',field);
        end
    else
        fprintf('Column indices exceed bounds of data matrices.\n');
    end
else
    fprintf('Must input columns or two data fields to scatter.\n');
end
