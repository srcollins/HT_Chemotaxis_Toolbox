function n=normalizeTwoColorChemotaxisData(mer)
%Normalizes HT chemotaxis data. For the data matrices, column 1 should be
%the experimental and column 2 should be the control samples.

% Remove data with low cell counts
cellCountThresh=20;
ind=mer.num0<cellCountThresh;
fn={'s0','s1','c0','c1','a0','a1'};
for i=1:length(fn)
    mer.(fn{i})(ind)=nan;
end

% Create normalized variable
n=mer;
n.a0=n.s0(:,:,1);
n.c0=n.c0(:,:,1);
fn={'s0','s1','c1','a1'};
for i=1:length(fn)
    n.(fn{i})=nan(length(n.rowlabels),length(n.collabels));
    for j=1:length(mer.collabels)
        c=myNanSpearman(mer.(fn{i})(:,j,2),mer.(fn{i})(:,j,1));
        p=robustfit(mer.(fn{i})(:,j,2),mer.(fn{i})(:,j,1));                     % Compute trendline for experimental data as a function of control data. This trend line is what we will use to normalize the experimental data.
        if c>0.1 && p(2)>0                                                      % Positive correlation -- In well control data is predictive of experimental data
            ref=p(1)+p(2)*mer.(fn{i})(:,j,2);                                   % Compute the expected experimental values, given the control values
            n.(fn{i})(:,j)=mer.(fn{i})(:,j,1)./ref;                             % Normalize by the expected value
            n.(fn{i})(:,j)=n.(fn{i})(:,j)/nanmedian(n.(fn{i})(:,j))-1;          % Further normalize such that the median is 1
        else
            n.(fn{i})(:,j)=mer.(fn{i})(:,j,1)/nanmedian(mer.(fn{i})(:,j,1))-1;  % If the controls had no predictive value, then just normalize by the median experimental value. Using the controls would likely just introduce noise.
        end
        p=robustfit(n.s0(:,j),n.s1(:,j));
        n.sD(:,j)=n.s1(:,j)-(p(1) + p(2)*n.s0(:,j));                            % Compute the chemokinesis parameter from a trendline of s0 vs s1 normalized speed phenotypes
    end
end
