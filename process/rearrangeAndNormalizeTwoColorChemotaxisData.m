function [n,mer2]=rearrangeAndNormalizeTwoColorChemotaxisData(mer)
%Assumes the first half of the columns are the controls and the second have
%are the experimental samples

fn=fieldnames(mer);
fn=setdiff(fn,{'rowlabels','collabels','set'});
num=round(length(mer.collabels)/2);

mer2.rowlabels=mer.rowlabels;
mer2.collabels=mer.collabels(1:num);
mer2.set=mer.set;

for i=1:length(fn)
    mer2.(fn{i})(:,:,1)=mer.(fn{i})(:,(num+1):length(mer.collabels));
    mer2.(fn{i})(:,:,2)=mer.(fn{i})(:,1:num);
end

n=normalizeTwoColorChemotaxisData(mer2);
n.a0=n.a0(:,:,1);
n.c0=n.c0(:,:,1);
