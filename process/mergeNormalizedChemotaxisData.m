function nOut=mergeNormalizedChemotaxisData(n1,n2)

n1.rowlabels=n1.rowlabels(:);
n2.rowlabels=n2.rowlabels(:);
if length(n1.rowlabels)==length(n2.rowlabels) & sum(strcmp(n1.rowlabels,n2.rowlabels))==length(n1.rowlabels)
    nOut.rowlabels=n1.rowlabels;
else
    nOut.rowlabels=union(n1.rowlabels,n2.rowlabels);
end

nOut.collabels=[n1.collabels(:); n2.collabels(:)];

[~,i1,i2]=intersect(n1.rowlabels,nOut.rowlabels);
fn=fieldnames(n1);
fn=setdiff(fn,{'rowlabels','collabels','set'})
for i=1:length(fn)
    nOut.(fn{i})=nan(length(nOut.rowlabels),length(nOut.collabels),size(n1.(fn{i}),3));
    if size(n1.(fn{i}),3)==1
        nOut.(fn{i})(i2,1:length(n1.collabels))=n1.(fn{i})(i1,:);
    else
        nOut.(fn{i})(i2,1:length(n1.collabels),:)=n1.(fn{i})(i1,:,:);
    end
end

[~,i1,i2]=intersect(n2.rowlabels,nOut.rowlabels);
for i=1:length(fn)
    if size(n1.(fn{i}),3)==1
        nOut.(fn{i})(i2,(length(n1.collabels)+1):(length(n1.collabels)+length(n2.collabels)))=n2.(fn{i})(i1,:);
    else
        nOut.(fn{i})(i2,(length(n1.collabels)+1):(length(n1.collabels)+length(n2.collabels)),:)=n2.(fn{i})(i1,:,:);
    end
end

