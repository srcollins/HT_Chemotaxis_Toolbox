function viewChemotaxisDataInPlateFormat(m,field,cols,bounds,marker)

if nargin<3
    cols=1;
end
if nargin<5
    marker=1;
end
if ismember(field,fieldnames(m))
    datMat=m.(field);
    if min(cols)>=1 && max(cols)<=size(datMat,2)
        dat=nanmean(datMat(:,cols,marker),2);
        switch length(dat)
            case 96
                plateSize=[12 8];
            case 384
                plateSize=[24 8];
        end
        if nargin<4 | isempty(bounds)
            mag2=max(dat);
            mag1=min(dat);
            if mag1<0 & mag2>0
                mag2=max(abs(dat));
                mag1=-1*mag2;
            end
            bounds=[mag1 mag2];
        end
        figure;
        [dat,cmap,bounds]=makeColormapAndAdjustData(dat,bounds);
        imagesc(reshape(dat,plateSize)',bounds);
%        colormap(blue_yellow);
        colormap(cmap);
        colorbar;
        set(gca,'box','off');
        rowlabels=cell(1,plateSize(2));
        for i=1:plateSize(2); rowlabels{i}=char(64+i); end
        set(gca,'YTickLabel',rowlabels);
    else
        fprintf('Column indices exceed bounds in input variable.\n');
    end
else
    fprintf('The specified field (%s) does not exist for the input variable.\n',field);
end

end

%-------------------------------------------------------
function [dat,cmap,bounds]=makeColormapAndAdjustData(dat,bounds)

nrows=100;
cval=zeros(nrows,1);
cmap=zeros(nrows,3);
cval(2:nrows)=bounds(1) + (0:(nrows-2))*(bounds(2)-bounds(1))/(nrows-2);
cval(1)=cval(2) - (cval(3)-cval(2));
cmap(2:nrows,:)=blue_yellow(nrows-1);
cmap(1,:)=[1 1 1]*0.4;

dat(dat<bounds(1))=bounds(1);
dat(dat>bounds(2))=bounds(2);
dat(isnan(dat))=cval(1);

bounds=[cval(1) cval(end)];

end