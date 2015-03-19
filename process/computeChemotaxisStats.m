function m=computeChemotaxisStats(t,frames0,framesPerStep,neighborDistThresh,angMinDistMoved)
% Computes chemotaxis statistics from TrackingExperiment objects.
% frames0 is the number of frames before uncaging
% t is expected to be a cell array of TrackingExperiment objects, where the
% number of rows is the number of high-throughput experiments, and there
% can be two columns for control and experimental samples from the same
% wells. In this case, column 1 should be experimental samples, and column 2
% should be control samples.

%% Set the row (e.g. gene name or well name) labels
fn=fieldnames(t{1,1});
if ismember('wellLabel',fn) && ~isempty(t{1,1}.wellLabel)
    m.rowlabels=t{1,1}.wellLabel;
else
    m.rowlabels=t{1,1}.wellName;
end

%% Set the column labels using the folder names
for i=1:size(t,1)
    temp=split(t{i,1}.folderPath,filesep);
    folderName=temp(end);
    m.collabels(i,1)=folderName;
end

%% Load default parameters, and use if necessary
defaultParameters = defaultParametersForComputingStatistics;
if nargin<2
    frames0=defaultParameters.frames0;
end
if nargin<3
    framesPerStep=defaultParameters.framesPerStep;
end
if nargin<4
    neighborDistThresh=defaultParameters.neighborDistThresh;
end
if nargin<5
    angMinDistMoved=defaultParameters.angMinDistMoved;
end
minDistFromCenter=defaultParameters.minDistFromCenter;
maxDistFromCenter=defaultParameters.maxDistFromCenter;

%% Compute Statistics
numWells=length(t{1,1}.wellName);
wells=1:numWells;
numExp=size(t,1);   % The number of experiments being processed
numPop=size(t,2);   % The number of independent populations (e.g., siRNA and Control) with data stored for each experiment

% Initialize fields
m.s0=nan(numWells,numExp,numPop);
m.s1=nan(numWells,numExp,numPop);
m.c0=nan(numWells,numExp,numPop);
m.c1=nan(numWells,numExp,numPop);
m.a0=nan(numWells,numExp,numPop);
m.a1=nan(numWells,numExp,numPop);
m.num0=nan(numWells,numExp,numPop);
m.num1=nan(numWells,numExp,numPop);

% Options sets for filtering step data
% For before gradient steps
optA=t{1,1}.parseStepFilteringInputs('func','nanmean','wells',wells,'minFrame',1,'maxFrame',frames0,'minDistMoved',0,'minDistFromCenter',0, ...
    'framesPerStep',framesPerStep,'maxDistFromCenter',Inf,'neighborDistThresh',neighborDistThresh);
optR1=t{1,1}.parseStepFilteringInputs('func','nanmean','wells',wells,'minFrame',1,'maxFrame',frames0,'minDistMoved',0,'minDistFromCenter',minDistFromCenter, ...
    'framesPerStep',framesPerStep,'maxDistFromCenter',maxDistFromCenter,'neighborDistThresh',neighborDistThresh);
optR2=t{1,1}.parseStepFilteringInputs('func','nanmean','wells',wells,'minFrame',1,'maxFrame',frames0,'minDistMoved',angMinDistMoved,'minDistFromCenter',minDistFromCenter, ...
    'framesPerStep',framesPerStep,'maxDistFromCenter',maxDistFromCenter,'neighborDistThresh',neighborDistThresh);
% For after gradient steps
opt1=t{1,1}.parseStepFilteringInputs('func','nanmean','wells',wells,'minFrame',frames0+1,'maxFrame',t{1,1}.numFrames,'minDistMoved',0,'minDistFromCenter',minDistFromCenter, ...
    'framesPerStep',framesPerStep,'maxDistFromCenter',maxDistFromCenter,'neighborDistThresh',neighborDistThresh);
opt2=t{1,1}.parseStepFilteringInputs('func','nanmean','wells',wells,'minFrame',frames0+1,'maxFrame',t{1,1}.numFrames,'minDistMoved',angMinDistMoved,'minDistFromCenter',minDistFromCenter, ...
    'framesPerStep',framesPerStep,'maxDistFromCenter',maxDistFromCenter,'neighborDistThresh',neighborDistThresh);

% Compute statistics
for i=1:numExp
    for k=1:numPop
        %Before Gradient Steps
        steps0=t{i,k}.stepSizesComputeOnly(t{i,k}.parseStepFilteringInputs('func','nanmean','wells',wells,'minFrame',1,'maxFrame',frames0,'minDistMoved',0,'minDistFromCenter',0, ...
            'framesPerStep',framesPerStep,'maxDistFromCenter',Inf,'neighborDistThresh',neighborDistThresh));
        
        %Find wells with no data
        emptyInd=cellfun(@(x) isempty(x),steps0);
        
        %Compute number of cells per well before gradient
        m.num0(emptyInd,i,k)=0;
        m.num0(~emptyInd,i,k)=cellfun(@(x) size(x(x(:,5)>=0 & x(:,5)<=Inf,:),1),steps0(~emptyInd))/div(frames0-1,framesPerStep);
        
        %Compute before gradient statistics
        for j=1:numWells
            if emptyInd(j)
                m.s0(j,i,k)=nan;
                m.c0(j,i,k)=nan;
                m.a0(j,i,k)=nan;
            else
                ind=t{i,k}.getFilteredStepInd(steps0{j},1,optA);
                m.s0(j,i,k)=feval(optA.func,steps0{j}(ind,1))*t{i,k}.getStepValColInfo(1,'factor',optA.wells(j),optA);
                ind=t{i,k}.getFilteredStepInd(steps0{j},1,optR1);
                m.c0(j,i,k)=feval(optA.func,steps0{j}(ind,7))*t{i,k}.getStepValColInfo(7,'factor',optA.wells(j),optR1);
                ind=t{i,k}.getFilteredStepInd(steps0{j},1,optR2);
                m.a0(j,i,k)=90-feval(optA.func,steps0{j}(ind,2))*t{i,k}.getStepValColInfo(2,'factor',optA.wells(j),optR2);
            end
        end
        
        %After Gradient Steps
        steps1=t{i,k}.stepSizesComputeOnly(t{i,k}.parseStepFilteringInputs('func','nanmean','wells',wells,'minFrame',frames0+1,'maxFrame',t{i,k}.numFrames,'minDistMoved',0,'minDistFromCenter',minDistFromCenter, ...
            'framesPerStep',framesPerStep,'maxDistFromCenter',maxDistFromCenter,'neighborDistThresh',neighborDistThresh));
        
        %Find wells with no data
        emptyInd=cellfun(@(x) isempty(x),steps1);
        
        %Compute number of cells per well before gradient
        m.num1(emptyInd,i,k)=0;
        m.num1(~emptyInd,i,k)=cellfun(@(x) size(x(x(:,5)>=minDistFromCenter & x(:,5)<=maxDistFromCenter,:),1),steps1(~emptyInd))/div(t{i,k}.numFrames-frames0-2,framesPerStep);
        
        %Compute after gradient statistics
        for j=1:numWells
            if emptyInd(j)
                m.s1(j,i,k)=nan;
                m.c1(j,i,k)=nan;
                m.a1(j,i,k)=nan;
            else
                ind=t{i,k}.getFilteredStepInd(steps1{j},1,opt1);
                m.s1(j,i,k)=feval(opt1.func,steps1{j}(ind,1))*t{i,k}.getStepValColInfo(1,'factor',opt1.wells(j),opt1);
                m.c1(j,i,k)=feval(opt1.func,steps1{j}(ind,7))*t{i,k}.getStepValColInfo(7,'factor',opt1.wells(j),opt1);
                ind=t{i,k}.getFilteredStepInd(steps1{j},1,opt2);
                m.a1(j,i,k)=90-feval(opt1.func,steps1{j}(ind,2))*t{i,k}.getStepValColInfo(2,'factor',opt1.wells(j),opt2);
            end
        end
    end
end

%% In case it useful, determine which wells were imaged simultaneously as a group.
for i=1:numExp
    startTimes=zeros(size(t{i,1}.wellName));
    for j=1:length(t{i,1}.wellName)
        [~,datenums,~]=TrackingExperiment.getSortedImageNames([t{i,1}.folderPath filesep t{i,1}.wellName{j}],t{i,1}.nucleusColor);
        startTimes(j)=datenums(1);
        if j==1
            gapSize=datenums(2)-datenums(1);
        end
    end
    rem=1:length(t{i,1}.wellName);
    j=0;
    while ~isempty(rem)
        j=j+1;
        sets{j}=find(abs(startTimes - min(startTimes(rem)))<=1.5*gapSize);
        rem=setdiff(rem,sets{j});
    end
    m.set(:,i)=sets(:);
end


