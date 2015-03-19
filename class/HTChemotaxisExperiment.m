classdef HTChemotaxisExperiment
    properties
        folderPath = '';
        
        %Nucleus detection and measurement properties
        nucleusColor = '';          %Will be entered by user
        cellMinArea = nan;          %Will be overwritten by defaultParametersForTracking
        cellMaxArea = nan;          %Will be overwritten by defaultParametersForTracking
        largeFilt=nan;              %Will be overwritten by defaultParametersForTracking
        thresh0 = nan;              %Will be overwritten by defaultParametersForTracking
        bgBlockSize = nan;          %Will be overwritten by defaultParametersForTracking
        bgPercentile = nan;         %Will be overwritten by defaultParametersForTracking
        imageMicronsPerPixel = nan; %Will be overwritten by defaultParametersForTracking
        localizationError = nan;    %Will be overwritten by defaultParametersForTracking
                                    %Average error in estimate of centroid coordinates (in pixel units)
        gradCenterXY='imageCenter'; %Will be overwritten by defaultParametersForTracking
        
        %Tracking properties
        maxDisp = 15;       %Will be overwritten by defaultParametersForTracking
        
        %Data
        imNumRows;
        imNumCols;
        wellName;           %Standard well name (A01, A02, etc.)
        wellLabel;          %Descriptive name (e.g. siRNA condition)
        numWells;
        numFrames;
        channel;
        wellThresh;
        coord;
        coordCollabels;
        traj;
        trajCollabels;
        filteredTraj;
        steps;
        stepCollabels;
        time;
        numGroups;
        groupName;
        groupInd;
        
    end
    methods
        %Constructor
        function t=HTChemotaxisExperiment(path)
            if strcmp(path,'load') | endsWith(path,'.mat')
                if endsWith(path,'.mat')
                    fullpath=path;
                else
                    [filename,path]=uigetfile('*.mat','Choose HTChemotaxisExperiment data file:');
                    fullpath=[path filename];
                end
                fprintf('Loading data from:\n%s\n\n',fullpath);
                dat=load(fullpath);
                fn=fieldnames(dat);
                fprintf('Found data variable: %s\n',fn{1});
                dat=dat.(fn{1});
                fn=fieldnames(dat);
                fn=intersect(fn,fieldnames(t));
                if ismember('nucleusColor',fn)
                    %load data
                    for i=1:length(fn)
                        t.(fn{i})=dat.(fn{i});
                    end
                else
                    beep;
                    fprintf('Data does not appear to be of the right format.\n');
                    t=[];
                end
            else
                
                %Populate basic data from the folder structure
                t.folderPath=path;
                t.wellName=getSubdirectories(path);
                if isempty(t.wellName)
                    error('No wells detected. Check data folder');
                end
                t.numWells=length(t.wellName);
                t.channel=htChemotaxisGetImageChannels([path filesep t.wellName{1}]);
                if isempty(t.channel)
                    error('No channels detected. Check data folder');
                end
                t.nucleusColor=t.userSelectColor('nucleus color');
                t.numFrames=max(cellfun(@(x) length(htChemotaxisGetSortedImageFilenames([path filesep x],t.nucleusColor)),t.wellName));
                t.wellLabel=readChemotaxisCoordMap(path);

                %Load default parameter values
                defaultParameters=defaultParametersForTracking;
                fn=fieldnames(defaultParameters);
                for i=1:length(fn)
                    t.(fn{i})=defaultParameters.(fn{i});
                end
                t.wellThresh=ones(t.numWells,t.numFrames)*t.thresh0;
            end
            t.gradCenterXY=t.interpretGradCenter;
        end
        
        %Adjust parameters
        function gradCenter=interpretGradCenter(obj)
            temp=obj.gradCenterXY;
            if ischar(temp)
                if strcmp(temp,'imageCenter')    %Read in one image and calculate the coordinates of the image center
                    try
                        im=htChemotaxisReadImage(obj.folderPath,obj.wellName{1},obj.nucleusColor,1);
                        sz=size(im);
                        gradCenter=[(sz(2)+1)/2 (sz(1)+1)/2];
                    catch
                        fprintf(' ----- Could not read corresponding images files.\n');
                        fprintf(' Using default values assuming a 1024 x 1280 image dimension.\n');
                        gradCenter=[640.5 512.5]
                    end
                end
            else
                gradCenter=temp;
            end
        end
        function color=userSelectColor(obj,desc)
            color='';
            message=['Enter ' desc ': '];
            fprintf('\nChannels detected are:\n');
            for i=1:(length(obj.channel)-1)
                fprintf('%s, ',obj.channel{i});
            end
            fprintf('%s\n',obj.channel{end});
            while ~ismember(color,obj.channel)
                color=input(message,'s');
            end
        end
        function obj=adjustThreshForNucleusDetection(obj)
            t0=obj.thresh0;
            fprintf('Current starting threshold value = %.0f\n',obj.thresh0);
            obj.thresh0=input('Enter new value: ');
            ratio=obj.thresh0/t0;
            obj.wellThresh=obj.wellThresh*ratio;
        end
        function obj=userAssignReplicateSets(obj)
            num=input('Enter number experimental groups: ');
            obj.numGroups=num;
            for i=1:num
                message=['Enter name for group ' num2str(i) ': '];
                obj.groupName{i}=input(message,'s');
                indStr=input('Enter indices: ','s');
                tok=getTokens(indStr,'([0-9]+)');
                obj.groupInd{i}=str2double(tok);
            end
        end
        
        %Image display functions
        function showOverlayWithDetectedNucleiSimple(obj,well,frame,bounds)
            %Uses the simple threshold only method for identifying nuclei
            warning('off','Images:initSize:adjustingMag');
            if nargin<4
                bounds=[];
            end
            im=htChemotaxisReadImage(obj.folderPath,obj.wellName{well},obj.nucleusColor,frame);
            coors=obj.getNucleiPositionsSimple(well,frame,im);
            
            figure;
            imshow(im,bounds); hold on;
            scatter(coors(:,1),coors(:,2),coors(:,4));  %use the measured area as the circle area
            hold off;
            colorbar;
            fprintf('Identified %i objects.\n',size(coors,1));
            fprintf('Median area = %.1f pixels.\n',median(coors(:,4)));
            fprintf('Median max intensity = %.0f.\n',10^median(coors(:,6)));
        end
        function showOverlayWithDetectedNucleiWatershed(obj,well,frame,bounds)
            %Uses the slower watershed-based method for identifying nuclei
            warning('off','Images:initSize:adjustingMag');
            if nargin<4
                bounds=[];
            end
            im=htChemotaxisReadImage(obj.folderPath,obj.wellName{well},obj.nucleusColor,frame);
            coors=obj.getNucleiPositionsWatershed(well,frame,im);
            
            figure;
            imshow(imageSubtractBackground(im,obj.bgPercentile,obj.bgBlockSize),bounds); hold on;
            scatter(coors(:,1),coors(:,2),coors(:,4),'MarkerEdgeColor',[0 1 0],'LineWidth',1);  %use the measured area as the circle area
            hold off;
            colorbar;
            fprintf('Identified %i objects.\n',size(coors,1));
            fprintf('Median area = %.1f pixels\n',median(coors(:,4)));
        end
        function plotXYTraj(obj,well,startFrame,endFrame,minDistMoved,minNumFrames)
            if nargin<3
                startFrame=1;
            end
            if nargin<4
                endFrame=Inf;
            end
            if nargin<5
                minDistMoved=0;
            end
            if nargin<6
                minNumFrames=0;
            end
            %Other parameters
            pointsize=1;
            frameCol=size(obj.traj{1}{1},2);
            
            % determine the number of trajectory sets in the input
            n=length(well);
            c=ceil(sqrt(n));
            r=ceil(n/c);
            
            % Set up the figure
            scrsz=get(0,'ScreenSize');
            figure('Position',[20 50 scrsz(3)-50 scrsz(4)-150]); subplot(r,c,1);
            
            % Make the plots
            for j=1:r
                for k=1:c
                    m=(j-1)*c+k;
                    if m<=length(well)
                        subplot(r,c,m);
                        hold on;
                        xyTraj=obj.traj{well(m)};
                        xyTraj=obj.filteredTraj{well(m)};
                        %find trajectories to be included
                        ind=cellfun(@(x) max(x(:,frameCol)),xyTraj) >= startFrame & ...
                            cellfun(@(x) min(x(:,frameCol)),xyTraj) <= endFrame & ...
                            cellfun(@(x) sqrt((x(1,1)-x(end,1))^2 + (x(1,2)-x(end,2))^2)>minDistMoved,xyTraj) & ...
                            cellfun(@(x) max(x(:,frameCol)),xyTraj) - cellfun(@(x) min(x(:,frameCol)),xyTraj) >= minNumFrames;
                        xyTraj=xyTraj(ind);
                        
                        %scatterplot the starting points
                        xval=cellfun(@(x) x(find(x(:,frameCol)>=startFrame,1),1),xyTraj);
                        yval=cellfun(@(x) x(find(x(:,frameCol)<=endFrame,1),2),xyTraj);
                        col=[0 0 0];
                        scatter(xval(:),yval(:),pointsize,col,'filled');
                        
                        %plot the trajectories
                        for i=1:length(xyTraj)
                            frameIDs=xyTraj{i}(:,frameCol);
                            ind=frameIDs >= startFrame & frameIDs <= endFrame;
                            if sum(ind)>1
                                dx=xyTraj{i}(find(ind,1,'last'),1)-xyTraj{i}(find(ind,1),1);
                                dy=xyTraj{i}(find(ind,1,'last'),2)-xyTraj{i}(find(ind,1),2);
                                len=sqrt(dx*dx + dy*dy + 1e-3);
                                dx=dx/len;
                                dy=dy/len;
                                centX=xyTraj{i}(find(ind,1),1)-obj.gradCenterXY(1);
                                centY=xyTraj{i}(find(ind,1),2)-obj.gradCenterXY(2);
                                len=sqrt(centX*centX + centY*centY + 1e-3);
                                centX=centX/len;
                                centY=centY/len;
                                angle=real(acos(dx*centX + dy*centY));
                                if isnan(angle); angle=pi; end
                                rgbVal=max([min(1,2*(pi-angle)/pi) min(1,2*angle/pi) 0],0);
                                %rgbVal=[0 1 0];
                                plot(xyTraj{i}(ind,1),xyTraj{i}(ind,2),'Color',rgbVal,'LineWidth',1);
                            end
                        end
                        xlim([0 obj.imNumCols]);
                        ylim([0 obj.imNumRows]);
                        set(gca,'YDir','reverse');
                        set(gca,'Color',[0 0 0]);
                        set(gca,'XTick',[]);
                        set(gca,'YTick',[]);
                        axis image
                        hold off;
                    end
                end
            end
        end
        function plotTrajPiece(obj,well,startFrame,endFrame,bounds)
            if nargin<3
                startFrame=1;
            end
            if nargin<4
                endFrame=Inf;
            end
            if nargin<5
                bounds=[0 1000];
            end
            frameCol=size(obj.traj{1}{1},2);
            warning('off','Images:initSize:adjustingMag');
            im1=htChemotaxisReadImage(obj.folderPath,obj.wellName{well},obj.nucleusColor,startFrame);
            coors1=obj.getNucleiPositionsWatershed(well,startFrame,im1);
            im2=htChemotaxisReadImage(obj.folderPath,obj.wellName{well},obj.nucleusColor,endFrame);
            coors2=obj.getNucleiPositionsWatershed(well,endFrame,im2);
            imM(:,:,3)=imageSubtractBackground(im1,obj.bgPercentile,obj.bgBlockSize);
            imM(:,:,1)=imageSubtractBackground(im2,obj.bgPercentile,obj.bgBlockSize);
            imM(:,:,2)=imM(:,:,1);
            imM=double(imM);
            imM=(imM-bounds(1))/(bounds(2)-bounds(1));
            imM=min(1,imM);
            
            figure;
            imshow(imM,bounds); hold on;
            xyTraj=obj.traj{well};
            %find trajectories to be included
            ind=cellfun(@(x) max(x(:,frameCol)),xyTraj) >= startFrame & ...
                cellfun(@(x) min(x(:,frameCol)),xyTraj) <= endFrame;
            xyTraj=xyTraj(ind);            
            %plot the trajectories
            for i=1:length(xyTraj)
                frameIDs=xyTraj{i}(:,frameCol);
                ind=frameIDs >= startFrame & frameIDs <= endFrame;
                if sum(ind)>1
                    plot(xyTraj{i}(ind,1),xyTraj{i}(ind,2),'Color',[0 1 0],'LineWidth',1);
                end
            end
            axis image
        end
        
        function plotAlignedTrajectories(obj,wells,startFrame,endFrame,minDistFromCenter,maxDistFromCenter)
            if isnumeric(wells)
                temp=wells;
                wells=cell(size(wells));
                for i=1:length(wells)
                    wells{i}=temp(i);
                end
            end
            if nargin<6
                maxDistFromCenter=650;
            end
            if nargin<5
                minDistFromCenter=100;
            end
            if nargin<3
                startFrame=1;
            end
            if nargin<4
                endFrame=obj.numFrames;
            end
            gradCenter=obj.gradCenterXY;
            frameCol=size(obj.filteredTraj{1,1}{1},2);
            
            figure;
            subplot(1,length(wells),1);
            xMin=0; xMax=0;% yMin=0; yMax=0;
            for i=1:length(wells)
                subplot(1,length(wells),i);
                trajSet={};
                for j=1:length(wells{i})
                    checkTraj=obj.traj{wells{i}(j)};
                    checkTraj=checkTraj(cellfun(@(x) min(x(:,frameCol))<=startFrame & max(x(:,frameCol))>=startFrame,checkTraj));  %startFrame should be included in the trajectory
                    checkTraj=checkTraj(cellfun(@(x) min(x(:,frameCol))<=endFrame & max(x(:,frameCol))>=endFrame,checkTraj));  %startFrame should be included in the trajectory
                    checkTraj=cellfun(@(x) {x(x(:,frameCol)>=startFrame & x(:,frameCol)<=endFrame,:)},checkTraj);
                    distFromCenter=cellfun(@(x) sqrt((gradCenter(1)-x(1,1))^2 + (gradCenter(2)-x(1,2))^2),checkTraj);
                    checkTraj=checkTraj(distFromCenter>minDistFromCenter & distFromCenter<maxDistFromCenter);
                    trajSet=[trajSet; checkTraj];
                end
                hold on
                plotX=[]; plotY=[];
                for j=1:length(trajSet)
                    %Align starting positions and orientations
                    xvals=trajSet{j}(:,1)-trajSet{j}(1,1);
                    yvals=trajSet{j}(:,2)-trajSet{j}(1,2);
                    orientV=[gradCenter(1)-trajSet{j}(1,1) gradCenter(2)-trajSet{j}(1,2)];
                    orientV=orientV/sqrt(sum(orientV.*orientV));
                    if orientV(2)<0
                        yvals=-1*yvals;
                    end
                    orientA=(pi/2)-real(acos(orientV(1)));
                    temp=[cos(orientA) -1*sin(orientA); sin(orientA) cos(orientA)]*[xvals'; yvals'];
                    xvals=temp(1,:)'; yvals=temp(2,:)';
                    
                    plot(xvals,yvals);
                    plotX=[plotX; xvals(end)]; plotY=[plotY; yvals(end)];
                end
                scatter(plotX,plotY,3);
                hold off
                xMin=min([xMin get(gca,'xlim') get(gca,'ylim')]);
                xMax=max([xMax get(gca,'xlim') get(gca,'ylim')]);
            end
            for i=1:length(wells)
                subplot(1,length(wells),i);
                axis image
                xlim([xMin xMax]);
                ylim([xMin xMax]);
            end
        end
        
        function vals=barStepVals(obj,col,varargin)
            if ~isempty(obj.traj)
                opt=obj.parseStepFilteringInputs(varargin{:});

                vals=zeros(length(opt.wells),1);
                tic
                tempSteps=obj.stepSizesComputeOnly(opt);
                toc
                numpts=zeros(length(opt.wells),size(obj.filteredTraj,2));
                for i=1:length(opt.wells)
                    for k=1:size(obj.filteredTraj,2)
                        x=tempSteps{opt.wells(i),k};
                        ind=obj.getFilteredStepInd(x,col,opt);
                        numpts(i,k)=sum(ind);
                        if numpts(i,k)>0
                            vals(i,k)=feval(opt.func,x(ind,col))*obj.getStepValColInfo(col,'factor',opt.wells(i),opt);
                        else
                            vals(i,k)=nan;
                        end
                    end
                end
                toc
                labels=HTChemotaxisExperiment.makeBarGraphLabels(obj.wellName(opt.wells),numpts);
                HTChemotaxisExperiment.generateBarGraph(vals,labels);
            else
                fprintf('Cells must be tracked first.\n');
            end
            if nargout==0
                vals=[];
            end
        end
  
        function vals=barStepPairMeans(obj,stat,varargin)
            if ~isempty(obj.steps)
                opt=obj.parseStepFilteringInputs(varargin{:});
                vals=zeros(size(obj.filteredTraj,2));
                numpts=zeros(length(opt.wells),size(obj.filteredTraj,2));
                s2=obj.stepPairsComputeOnly(opt);
                for i=1:length(opt.wells)
                    for k=1:size(obj.filteredTraj,2)
                        x=s2{opt.wells(i),k};
                        if ~isempty(x)
                            ind=~isnan(x(:,1)) & x(:,5)>=opt.minDistMoved & x(:,6)>=opt.minDistMoved;
                            switch stat
                                case 'pSpd'
                                    vals(i,k)=nanmean(min(x(ind,5),x(ind,6))./max(x(ind,5),x(ind,6)));
                                    ylab='Fraction of faster step';
                                    tlab='Persistence of Speed';
                                case 'pDot'
                                    vals(i,k)=nanmean(x(ind,1).*x(ind,3) + x(ind,2).*x(ind,4));
                                    factor=obj.getStepValColInfo(1,'factor',opt.wells(i),opt)^2;
                                    vals(i,k)=vals(i,k)/factor;
                                    ylab='Dot Product (um/min)^2';
                                    tlab='Persistence Dot Product';
                                case 'pCos'
                                    vals(i,k)=nanmean((x(ind,1).*x(ind,3) + x(ind,2).*x(ind,4))./x(ind,5)./x(ind,6));
                                    ylab='Cosine';
                                    tlab='Persistence of Direction - Cosine';
                                case 'pAng'
                                    vals(i,k)=nanmean((x(ind,1).*x(ind,3) + x(ind,2).*x(ind,4))./x(ind,5)./x(ind,6));
                                    vals(i,k)=real(acos(vals(i,k))/pi*180);
                                    ylab='Angle';
                                    tlab='Persistence of Direction - Angle';
                            end
                            numpts(i,k)=sum(ind);
                        else
                            vals(i,k)=0;
                            numpts(i,k)=0;
                        end
                    end
                end
                labels=HTChemotaxisExperiment.makeBarGraphLabels(obj.wellName(opt.wells),numpts);
                HTChemotaxisExperiment.generateBarGraph(vals,labels);
                title(tlab,'FontSize',14);
                ylabel(ylab,'FontSize',11);
            else
                fprintf('Cells must be tracked first.\n');
            end
        end
        function plotPooledGroupStepValsOverTime(obj,col,groupNums,varargin)
            if ~isempty(obj.traj) && ~isempty(obj.groupInd)
                %Make sure that the requested group indices have indeed been defined
                groupNums0=groupNums;
                groupNums=intersect(groupNums0,1:length(obj.groupInd));
                if length(groupNums)<length(groupNums0)
                    fprintf('Ignoring groups that have not yet been defined.\n');
                end
                
                %Handle input parameters and compute step values
                opt=obj.parseStepFilteringInputs(varargin{:});
                opt=obj.parseStepFilteringInputs(varargin{:});
                if isempty(opt.omitFrame)
                    XVals=1:opt.framesPerStep:(obj.numFrames-opt.framesPerStep);
                else
                    XVals=[];
                    for i=1:(length(opt.omitFrame)+1)
                        XVals=[XVals max([1 max(XVals)+opt.framesPerStep opt.omitFrame(1:(i-1))+1]):opt.framesPerStep:min([obj.numFrames-opt.framesPerStep opt.omitFrame(i:end)-opt.framesPerStep])]; 
                    end
                end
                vals=zeros(length(opt.wells),1);
                tempSteps=obj.stepSizesComputeOnly(opt);
                
                %Generate the figure
                figure; hold on;
                numpts=zeros(length(groupNums)*size(obj.filteredTraj,2),1);
                count=0;
                str=vect(obj.groupName(groupNums))';
                for i=1:length(groupNums)
                    for k=1:size(obj.filteredTraj,2)
                        count=count+1;
                        numpts1=zeros(length(XVals),1);
                        for j=1:length(XVals)                               % time points
                            j1=XVals(j);                                    % this time point
                            x=cell2mat(tempSteps(obj.groupInd{groupNums(i)},k));    % Pull out all steps for wells corresponding to this group
                            x=x(x(:,6)==j1,:);                              % pull out only steps for the right time point
                            if ~isempty(x)
                                ind=obj.getFilteredStepInd(x,col,opt);
                                vals(i,j,k)=feval(opt.func,x(ind,col));
                                numpts1(j)=sum(ind);
                            else
                                vals(i,k)=0;
                                numpts1(j)=0;
                            end
                        end
                        numpts(count)=round(mean(numpts1));
                        % Convert units from pixels to microns, if needed
                        if boolRegExp(obj.getStepValColInfo(col,'units'),'Micron')
                            vals=vals*obj.imageMicronsPerPixel;
                        end
                        % Divide by time intervals, if appropriate
                        timeMat=cell2mat(obj.time(obj.groupInd{i}));
                        timeSteps=nanmean(timeMat(XVals+opt.framesPerStep,:)-timeMat(XVals,:),2)/60; %nanmean(diff(timeMat),2)/60;
                        timeXVals=mean(timeMat(XVals,:)+repmat(timeSteps,[1 size(timeMat,2)])/2,2)/60;
                        if boolRegExp(obj.getStepValColInfo(col,'units'),'Min')
                            vals(i,:,k)=vals(i,:,k)./(vect(timeSteps)');
                        end
                        if k==1
                            plot(timeXVals,vals(i,:,k),'LineWidth',2,'Color',[mod((i+1)/3,1) mod(i/5,1) mod(i/8,1)]);
                        else
                            plot(timeXVals,vals(i,:,k),'LineWidth',2,'LineStyle','- -');
                        end
                    end
                end
                title(obj.getStepValColInfo(col,'title'),'FontSize',14);
                xlabel('Time (min)','FontSize',11);
                ylabel(obj.getStepValColInfo(col,'ylabel'),'FontSize',11);
                labels=HTChemotaxisExperiment.makeBarGraphLabels(str(:),numpts);
                legend(labels);
                
            else
                fprintf('Trajectories must be computed and groups must be defined first.\n');
            end
        end
        function plotStepValsOverTime(obj,col,varargin)
            if ~isempty(obj.traj)
                %Handle input parameters and compute step values
                opt=obj.parseStepFilteringInputs(varargin{:});
                if isempty(opt.omitFrame)
                    XVals=1:opt.framesPerStep:(obj.numFrames-opt.framesPerStep);
                else
                    XVals=[];
                    for i=1:(length(opt.omitFrame)+1)
                        XVals=[XVals max([1 max(XVals)+opt.framesPerStep opt.omitFrame(1:(i-1))+1]):opt.framesPerStep:min([obj.numFrames-opt.framesPerStep opt.omitFrame(i:end)-opt.framesPerStep])]; 
                    end
                end
                vals=zeros(length(opt.wells),1);
                tempSteps=obj.stepSizesComputeOnly(opt);
                
                %Generate the figure
                figure; hold on;
                numpts=zeros(length(opt.wells)*size(obj.filteredTraj,2),1);
                count=0;
                if isnumeric(opt.wells)
                    str=vect(obj.wellName(opt.wells))';
                else   %opt.wells should be a cell array
                    str=cellfun(@(x) {['Group (' obj.wellName{x(1)} ')']},opt.wells);
                end
                for i=1:length(opt.wells)
                    for k=1:size(obj.filteredTraj,2)
                        count=count+1;
                        numpts1=zeros(length(XVals),1);
                        for j=1:length(XVals)
                            j1=XVals(j);
                            if isnumeric(opt.wells)
                                x=tempSteps{opt.wells(i),k};
                            else
                                x=cell2mat(tempSteps(opt.wells{i},k));
                            end
                            x=x(x(:,6)==j1,:);
                            if ~isempty(x)
                                ind=obj.getFilteredStepInd(x,col,opt);
                                vals(i,j,k)=feval(opt.func,x(ind,col));
                                numpts1(j)=sum(ind);
                            else
                                vals(i,k)=0;
                                numpts1(j)=0;
                            end
                        end
                        numpts(count)=round(mean(numpts1));
                        % Convert units from pixels to microns, if needed
                        if boolRegExp(obj.getStepValColInfo(col,'units'),'Micron')
                            vals=vals*obj.imageMicronsPerPixel;
                        end
                        % Divide by time intervals, if appropriate
                        if isnumeric(opt.wells)
                            timeMat=obj.time{i};
                        else
                            timeMat=cell2mat(obj.time(opt.wells{i}));
                        end
                        timeSteps=nanmean(timeMat(XVals+opt.framesPerStep,:)-timeMat(XVals,:),2)/60; %nanmean(diff(timeMat),2)/60;
                        timeXVals=mean(timeMat(XVals,:)+repmat(timeSteps,[1 size(timeMat,2)])/2,2)/60;
                        if boolRegExp(obj.getStepValColInfo(col,'units'),'Min')
                            vals(i,:,k)=vals(i,:,k)./(vect(timeSteps)');
                        end
                        if k==1
                            plot(timeXVals,vals(i,:,k),'LineWidth',2,'Color',[mod((i+1)/3,1) mod(i/5,1) mod(i/8,1)]);
                        else
                            plot(timeXVals,vals(i,:,k),'LineWidth',2,'LineStyle','- -');
                        end
                    end
                end
                title(obj.getStepValColInfo(col,'title'),'FontSize',14);
                xlabel('Time (min)','FontSize',11);
                ylabel(obj.getStepValColInfo(col,'ylabel'),'FontSize',11);
                labels=HTChemotaxisExperiment.makeBarGraphLabels(str(:),numpts);
                legend(labels);
            else
                fprintf('Trajectories must be computed first.\n');
            end
        end
        function cdfStepVals(obj,col,varargin)
            if ~isempty(obj.steps)
                opt=obj.parseStepFilteringInputs(varargin{:});
                
                stp=obj.stepSizesComputeOnly(opt);
                vals=cell(length(opt.wells),size(stp,2));
                numpts=zeros(length(opt.wells),size(obj.filteredTraj,2));
                if iscell(opt.wells)
                    str=vect(obj.wellName(cellfun(@(x) x(1),opt.wells)));
                    str=cellfun(@(x) {['Group ' x]},str);
                else
                    str=vect(obj.wellName(opt.wells));
                end
                for i=1:length(opt.wells)
                    for k=1:size(stp,2)
                        if iscell(opt.wells)
                            x=[];
                            for m=1:length(opt.wells{i})
                                x=[x; stp{opt.wells{i}(m),k}];
                            end
                        else
                            x=stp{opt.wells(i),k};
                        end
                        if ~isempty(x)
                            ind=obj.getFilteredStepInd(x,col,opt);
                            vals(i,k)={x(ind,col)};
                            numpts(i,k)=sum(ind);
                        else
                            vals(i,k)={[]};
                            numpts(i,k)=0;
                        end
                    end
                end
                plotCDFs(vals(:));
                labels=HTChemotaxisExperiment.makeBarGraphLabels(str(:),numpts(:));
                legend(labels);
            else
                fprintf('Cells must be tracked first.\n');
            end
        end
        function [vals,bins]=pdfStepVals(obj,col,bins,varargin)
            if ~isempty(obj.steps)
                opt=obj.parseStepFilteringInputs(varargin{:});

                stp=obj.stepSizesComputeOnly(opt);
                vals=cell(length(opt.wells),size(stp,2));
                numpts=zeros(length(opt.wells),size(obj.filteredTraj,2));
                if iscell(opt.wells)
                    str=vect(obj.wellName(cellfun(@(x) x(1),opt.wells)));
                    str=cellfun(@(x) {['Group ' x]},str);
                else
                    str=vect(obj.wellName(opt.wells));
                end
                for i=1:length(opt.wells)
                    for k=1:size(stp,2)
                        if iscell(opt.wells)
                            x=[];
                            for m=1:length(opt.wells{i})
                                x=[x; stp{opt.wells{i}(m),k}];
                            end
                        else
                            x=stp{opt.wells(i),k};
                        end
                        if ~isempty(x)
                            ind=obj.getFilteredStepInd(x,col,opt);
                            vals(i,k)={x(ind,col)};
                            numpts(i,k)=sum(ind);
                        else
                            vals(i,k)={[]};
                            numpts(i,k)=0;
                        end
                    end
                end
                labels=HTChemotaxisExperiment.makeBarGraphLabels(str(:),numpts(:));
                plotFrequencyDistributions(vals,bins);
                legend(labels);
            else
                fprintf('Step sizes must be computed first.\n');
            end
        end
        function scatterStepVals(obj,col1,col2,varargin)
            if ~isempty(obj.steps)
                opt=obj.parseStepFilteringInputs(varargin{:});
                
                stp=obj.stepSizesComputeOnly(opt);
                vals=cell(length(opt.wells),size(stp,2));
                numpts=zeros(length(opt.wells),size(obj.filteredTraj,2));
                str=vect(obj.wellName(opt.wells));
                for k=1     %     I will have to decide how to handle two columns later
                    x=[];
                    for m=1:length(opt.wells)
                        x=[x; stp{opt.wells(m),k}];
                    end
                    ind=obj.getFilteredStepInd(x,col1,opt);
                    vals1=x(ind,col1);
                    vals2=x(ind,col2);
                    numpts=sum(ind);
                end
                %labels=HTChemotaxisExperiment.makeBarGraphLabels(str(:),numpts(:));
                figure; keyboard;
                dscatter(vals1,vals2);
                %legend(labels);
            else
                fprintf('Step sizes must be computed first.\n');
            end
        end
          
        function barGroupStepVals(obj,col,groupNums,varargin)
            if ~isempty(obj.traj)
                opt=obj.parseStepFilteringInputs(varargin{:});

                stp=obj.stepSizesComputeOnly(opt);
                vals=zeros(length(groupNums),1);
                errvals=vals;
                labels=cell(length(groupNums),1);
                for i=1:length(groupNums)
                    labels{i}=obj.groupName{groupNums(i)};
                    for k=1:size(stp,2)
                        valList=zeros(length(obj.groupInd{groupNums(i)}),1);
                        for m=1:length(obj.groupInd{groupNums(i)})
                            x=stp{obj.groupInd{groupNums(i)}(m),k};
                            if ~isempty(x)
                                ind=obj.getFilteredStepInd(x,col,opt);
                                valList(m,1)=feval(opt.func,x(ind,col))*obj.getStepValColInfo(col,'factor',obj.groupInd{i}(m),opt);
                            else
                                valList(m,1)=0;
                            end
                        end
                        vals(i,k)=mean(valList);
                        errvals(i,k)=std(valList)/sqrt(length(valList));
                    end
                end
                HTChemotaxisExperiment.generateBarGraph(vals,labels,errvals);
                title(obj.getStepValColInfo(col,'title'),'FontSize',14);
                ylabel(obj.getStepValColInfo(col,'ylabel'),'FontSize',11);
            else
                fprintf('Step sizes must be computed first.\n');
            end
        end
        
        %Image processing and calculation
        function [coors,pix]=getNucleiPositionsSimple(obj,well,frame,im)
            %Either well number and frame number can be specified, or an
            %image can be passed in as the argument
            if nargin<4
                im=htChemotaxisReadImage(obj.folderPath,obj.wellName{well},obj.nucleusColor,frame);
            end
            im=imageSubtractBackground(im,obj.bgPercentile,obj.bgBlockSize);
            mask = im > obj.wellThresh(well,frame);
            labeled = logical(bwareaopen(mask,obj.cellMinArea,4));
            dl=regionprops(labeled,'Centroid','PixelIdxList','Area');%get centroid & pixelidxlist
            im=double(im);
            num=size(dl,1);
            coors=nan(num,6);
            if nargout>1
                pix=cell(num,1);
            end
            if isempty(obj.coordCollabels)
                obj=obj.addCoordCollabels;
            end
            for l=1:num
                loc=dl(l).Centroid;
                coors(l,1:2)=loc;
                coors(l,3)=median(log10(double(im(dl(l).PixelIdxList))));
                coors(l,4)=dl(l).Area;
                coors(l,5)=log10(sum(vect(double(im(dl(l).PixelIdxList)))));
                coors(l,6)=log10(max(double(im(dl(l).PixelIdxList))));
                if nargout>1
                    pix{l}=dl(l).PixelIdxList;
                end
            end
            
            coors(coors(:,4)>obj.cellMaxArea,:)=[];
            coors(coors(:,4)<obj.cellMinArea,:)=[];     %remove objects that are too small
%             p1=coors(:,1:2);
%             [~,d]=annsearch(p1',p1',2);
%             coors(:,7)=d(2,:)';
        end
        function [coors,pix]=getNucleiPositionsWatershed(obj,well,frame,im)
            %Either well number and frame number can be specified, or an
            %image can be passed in as the argument
            if nargin<4
                im=htChemotaxisReadImage(obj.folderPath,obj.wellName{well},obj.nucleusColor,frame);
            end
            im=imageSubtractBackground(im,obj.bgPercentile,obj.bgBlockSize);
            mask = im > obj.wellThresh(well,frame);
            if obj.largeFilt>0
                mask=mask-bwareaopen(mask,obj.largeFilt,4);
            end
            D1=bwdist(mask);
            D = -1*double(imopen(im,strel('disk',1))-obj.wellThresh(well,frame));
            D(~mask) = -Inf;
            D(D1>0 & D1<1.2)=0;
            labeled = watershed(D,4);
            dl=regionprops(labeled,'Centroid','PixelIdxList','Area');%get centroid & pixelidxlist
            im=double(im);
            num=size(dl,1);
            coors=nan(num,6);
            if nargout>1
                pix=cell(num,1);
            end
            if isempty(obj.coordCollabels)
                obj=obj.addCoordCollabels;
            end
            for l=1:num
                loc=dl(l).Centroid;
                coors(l,1:2)=loc;
                coors(l,3)=median(log10(double(im(dl(l).PixelIdxList))));
                coors(l,4)=dl(l).Area;
                coors(l,5)=log10(sum(vect(double(im(dl(l).PixelIdxList)))));
                coors(l,6)=log10(max(double(im(dl(l).PixelIdxList))));
                if nargout>1
                    pix{l}=dl(l).PixelIdxList;
                end
            end     %Fill in the data for each object
            
            coors(coors(:,4)>obj.cellMaxArea,:)=[];
            coors(coors(:,4)<obj.cellMinArea,:)=[];     %remove objects that are too small
            p1=coors(:,1:2);
            if size(coors,1)>=3
                [~,d]=annsearch(p1',p1',2);
                coors(:,7)=d(2,:)';
            else
                coors(:,7)=999;
            end
        end
        function obj=addDist2NearestNeighbor(obj)
            %This can be used if nuclei positions have already been
            %identified without this field
            if ismember('dist2nearest',obj.coordCollabels)
                nearestInd=find(strcmp('dist2nearest',obj.coordCollabels));
            else
                nearestInd=size(obj.coord{1,1},2)+1;
                obj.coordCollabels{nearestInd}='dist2nearest';
            end
            for i=1:size(obj.coord,1)
                for j=1:size(obj.coord,2)
                    p1=obj.coord{i,j}(:,1:2);
                    if size(p1,1)>2
                        [~,d]=annsearch(p1',p1',2);
                        obj.coord{i,j}(:,nearestInd)=d(2,:)';
                    end
                    if size(p1,1)==1 | size(p1,1)==2
                        obj.coord{i,j}(:,nearestInd)=999;
                    end
                end
            end
        end
        function obj=findNuclei(obj)
            sourcePath=obj.folderPath;
            nucColor=obj.nucleusColor;
            matlabpool(4);
            obj.coord=cell(obj.numWells,obj.numFrames);
            for well=1:obj.numWells
                tic
                fprintf('Detecting nuclei - Well %s -- ',obj.wellName{well});
                thisWellName=obj.wellName{well};
                [filenames,datenums]=htChemotaxisGetSortedImageFilenames([sourcePath filesep thisWellName],nucColor);
                temp=cell(1,obj.numFrames);
                parfor frame=1:obj.numFrames
                    if frame<=length(filenames)
                        im=htChemotaxisReadImage(sourcePath,thisWellName,nucColor,frame);
                        temp{frame}=obj.getNucleiPositionsWatershed(well,frame,im);
                    else
                        temp{frame}=[];
                    end
                end
                times=(datenums-min(datenums))*24*3600;
                obj.coord(well,:)=temp;
                obj.time{well}=times(:);
                toc
            end
            matlabpool close;
            %Adjust the image size properties
            im=htChemotaxisReadImage(obj.folderPath,obj.wellName{1},obj.nucleusColor,1);
            [obj.imNumRows, obj.imNumCols]=size(im);
        end
        function obj=trackCells(obj)
            obj.traj=cell(obj.numWells,1);
            fprintf('Generating cell tracks\n');
            for i=1:obj.numWells
                fprintf('%s ',obj.wellName{i});
                nums=cellfun(@(x) size(x,1),obj.coord(i,:));
                if nanmedian(nums)>5 %there are at least a few cells to track

                    % Build the set of trajectories
                    trajFound=obj.getTrajectoriesForSingleWell(i);
                    obj.traj{i}=trajFound;
                else
                    obj.traj{i}=[];
                end
            end
            fprintf('\n\n');
            obj.filteredTraj=obj.traj;
            if isempty(obj.coordCollabels)
                obj=obj.addCoordCollabels;
            end
            obj=obj.computeStepSizes;
        end
        function traj=getTrajectoriesForSingleWell(obj,wellNum)
            coord=squeeze(obj.coord(wellNum,:));  % Cell coordinates for tracking
            
            % Prepare a matrix to hold tracked cell IDs across all frames
            trajIndMat=nan(sum(cellfun(@(x) size(x,1),coord(1:end-1))),length(coord));
            trajIndMat(1:size(coord{1},1),1)=vect(1:size(coord{1},1));
            
            numTraj=size(coord{1},1);  % A counter to keep track of the number of trajectories built so far
            for i=1:(length(coord)-1)
                if ~isempty(coord{i}) && ~isempty(coord{i+1})
                    % Forward match
                    [ix,d]=annsearch(coord{i+1}(:,1:2)',coord{i}(:,1:2)',1);
                    M1=[vect(1:length(ix)) ix(:)];
                    M1(d>obj.maxDisp,:)=[];
                    
                    % Backwards match
                    [ix,d]=annsearch(coord{i}(:,1:2)',coord{i+1}(:,1:2)',1);
                    M2=[vect(1:length(ix)) ix(:)];
                    M2(d>obj.maxDisp,:)=[];
                    
                    % Filter for reciprocal matches
                    M=intersect(M1(:,1:2),M2(:,[2 1]),'rows');
                    if isempty(M)
                        M=[0 0];   % If no matches are found, [0 0] is used as a placeholder to avoid an error from intersect
                    end
                else
                    M=[0 0];   % There can't be any matches if there are no cells in one of the frames. [0 0] is used as a placeholder to avoid an error from intersect
                end
                
                % Apply matches
                [junk,i1,i2]=intersect(trajIndMat(:,i),M(:,1));
                trajIndMat(i1,i+1)=M(i2,2);   %Continuing trajectories
                
                if i<(length(coord)-1)        %Only need to add new trajectories if it is not the last frame
                    idVect=vect(1:size(coord{i+1},1));
                    [junk,i1]=setdiff(idVect,M(:,2));
                    trajIndMat((numTraj+1):(numTraj+length(i1)),i+1)=idVect(i1);  %Potential new trajectories
                    numTraj=numTraj+length(i1);
                end
            end
            
            % Clean up trajIndMat
            trajIndMat=trajIndMat(1:numTraj,:);   %Remove empty rows
            hasAtLeastOneLink=sum(~isnan(trajIndMat),2)>1;
            trajIndMat=trajIndMat(hasAtLeastOneLink,:);
            
            % Generate a cell array of trajectories
            numTraj=size(trajIndMat,1);
            numCols=size(coord{1},2)+2;
            traj=cell(numTraj,1);
            for i=1:numTraj
                framesCovered=find(~isnan(trajIndMat(i,:)));
                traj{i}=nan(length(framesCovered),numCols);
                for j=1:length(framesCovered)
                    frameNum=framesCovered(j);
                    traj{i}(j,:)=[coord{frameNum}(trajIndMat(i,frameNum),:) i frameNum];
                end
            end
        end
        function obj=removeDrift(obj)
            trajDistThresh=5;
            if isempty(obj.stepCollabels)
                obj=obj.addStepCollabels;
            end
            trajIDcol=find(strcmp(obj.stepCollabels,'trajID'));
            startFrameCol=find(strcmp(obj.stepCollabels,'startFrame'));
            deltaXcol=find(strcmp(obj.stepCollabels,'deltaX'));
            deltaYcol=find(strcmp(obj.stepCollabels,'deltaY'));
            for i=1:obj.numWells
                steps1=obj.steps{i};  %This should be run before filtering trajectories
                trajlen=cellfun(@(x) size(x,1),obj.traj{i});
                trajdist=cellfun(@(x) sqrt((x(1,1)-x(end,1))^2 + (x(1,2)-x(end,2))^2),obj.traj{i});
                ind=find(trajlen>=0.5*obj.numFrames & trajdist<trajDistThresh);
                fprintf('Well = %s  Using %i trajectories, %i full time course.\n',obj.wellName{i},length(ind),sum(trajlen==obj.numFrames & trajdist<trajDistThresh));
                steps1=steps1(ismember(steps1(:,trajIDcol),ind),:);  %Only use information from the nonmoving cells
                cDrift=zeros(obj.numFrames,2);
                for j=2:obj.numFrames
                    d1=nanmedian(steps1(steps1(:,startFrameCol)==j-1,deltaXcol));
                    d2=nanmedian(steps1(steps1(:,startFrameCol)==j-1,deltaYcol));
                    cDrift(j,1)=cDrift(j-1,1)+d1;
                    cDrift(j,2)=cDrift(j-1,2)+d2;
                    obj.coord{i,j}(:,1)=obj.coord{i,j}(:,1)-cDrift(j,1);
                    obj.coord{i,j}(:,2)=obj.coord{i,j}(:,2)-cDrift(j,2);
                end
            end
            obj=obj.trackCells;
            obj=obj.computeStepSizes;
        end
        function obj=filterTrajectories(obj,minDistMoved)
            for i=1:obj.numWells
                for k=1:size(obj.traj,2)
                    if size(obj.traj{i,k},1)>0
                        ind=cellfun(@(x) sqrt((x(1,1)-x(end,1))^2 + (x(1,2)-x(end,2))^2)>minDistMoved,obj.traj{i,k});
                        obj.filteredTraj{i,k}=obj.traj{i,k}(ind);
                    end
                end
            end
        end
        
        function opt=parseStepFilteringInputs(obj,varargin)
            % defaults
            opt=defaultParametersForComputingStatistics;
            opt=rmfield(opt,'angMinDistMoved');  % Only minDistMoved will be used from this options set
            opt=rmfield(opt,'frames0');          % frames0 will be defined elsewhere, and will not be a part of this options set
            opt.wells=1:obj.numWells;
            opt.minFrame=1;
            opt.maxFrame=obj.numFrames;
            opt.minDistMoved=0;
            opt.frameGap=0;  %Only important for functions using step pairs
            opt.func='mean';
            opt.omitFrame=[];       % Only implemented for plotStepValsOverTime and plotPooledGroupStepValsOverTime
            fn=fieldnames(opt);
            
            % parse inputs
            if mod(length(varargin),2)==1
                fprintf(' ** Number of input parameters and values should match.\nUsing default values.\n');
            else
                for i=1:2:length(varargin)
                    if ismember(varargin{i},fn)
                        opt.(varargin{i})=varargin{i+1};
                    else
                        fprintf(' ** Unrecognized parameter -- %s\n',varargin{i});
                    end
                end
            end
        end
        function steps=stepSizesComputeOnly(obj,opt)
            if nargin<2
                opt=obj.parseStepFilteringInputs;
            end
            if iscell(opt.wells)
                opt.wells=unique(cell2mat(opt.wells));
            end
%            gradCenter=[(1+obj.imNumCols) (1+obj.imNumRows)]/2;
%            gradCenter=[1200 880]/2;
            gradCenter=obj.gradCenterXY;
            steps=cell(size(obj.filteredTraj));
            neighborCol=obj.getCoordColIndex('dist2nearest');
            frameCol=size(obj.filteredTraj{1,1}{1},2);
            trajIDCol=frameCol-1;
            for i=1:obj.numWells
                if ismember(i,opt.wells)
                    for k=1:size(obj.filteredTraj,2)
                        trajMat=cell2mat(obj.filteredTraj{i,k});   % this concatenates across the cell array. A 5x1 cell array of
                        % 10x8 matrices will be converted into a 50x8 matrix. A similar 5x2 cell array would be converted into a 50x16 matrix. xyTraj should be a column vector.
                        steps{i,k}=[];  %in case there are no tracked cells
                        if size(trajMat,1)>1
                            indRem=trajMat(:,frameCol)<opt.minFrame | trajMat(:,frameCol)>opt.maxFrame;
                            trajMat(indRem,:)=[];
                            if size(trajMat,1)>1
                                distbool=trajMat(1:(end-1),neighborCol)'>opt.neighborDistThresh | trajMat(2:end,neighborCol)'>opt.neighborDistThresh;
                                len=size(trajMat,1);
                                
                                %Find the set of lines that would be allowable step starting points. This says the nearest neighbor distance constraint must be satisfied for all lines encompassed in a step starting at this line
                                okInd=true(1,len); okInd((len-opt.framesPerStep+1):len)=false;
                                for mm=1:(opt.framesPerStep)
                                    indval=1:(len-opt.framesPerStep);
                                    okInd(indval)=okInd(indval) & distbool(indval+mm-1);
                                end
                                %Exclude cases where the lines t.plcorresponding to a putative step span from one trajectory into the next
                                okInd(1:(len-opt.framesPerStep))=okInd(1:(len-opt.framesPerStep)) & (trajMat(1:(len-opt.framesPerStep),trajIDCol)==trajMat((1+opt.framesPerStep):len,trajIDCol))';
                                
                                %Find a set of lines for step starting points that avoids duplicate use of any data points
                                indStart=false(1,len);
                                for m=opt.minFrame:(opt.maxFrame-opt.framesPerStep)                 %Iterate through all frames in range
                                    indStart=(indStart+(okInd & trajMat(:,frameCol)'==m))>0;        %If a line with this frame is ok, make it a starting point
                                    okInd=(okInd-indStart)>0;                                       %Exclude these lines for the next iteration
                                    for mm=1:(opt.framesPerStep-1)                                  %Also exclude any lines encompassed in these steps
                                        okInd((1+mm):(len-opt.framesPerStep+mm))=(okInd((1+mm):(len-opt.framesPerStep+mm))-indStart(1:(len-opt.framesPerStep)))>0;
                                    end
                                end
                                
                                %Define the corresponding step end lines
                                indEnd=false(1,len);
                                indEnd((1+opt.framesPerStep):len)=indStart(1:(len-opt.framesPerStep));
                                indKeep=(indStart+indEnd)>0;
                                trajMat=trajMat(indKeep,:); indStart1=indStart(indKeep);
%                                sameTraj=indStart1(1:end-1)' & diff(trajMat(:,trajIDCol))==0;
%                                steps{i,k}=obj.trajMat2steps(trajMat,gradCenter,sameTraj);
                                steps{i,k}=obj.trajMat2steps(trajMat,gradCenter,indStart1(1:end-1)',obj.time{i},opt.framesPerStep); %Start indices where the trajectory ID doesn't match in the following line have already been filtered out
                            end
                        end
                    end
                end
            end
        end
        function s2=stepPairsComputeOnly(obj,opt)
            %This function needs to be finished - it has not been fully
            %corrected for using opt and for filtering based on close
            %neighboring cells
            s2=cell(size(obj.filteredTraj));
            neighborCol=obj.getCoordColIndex('dist2nearest');
            frameCol=size(obj.filteredTraj{1,1}{1},2);
            trajIDCol=frameCol-1;
            for i=1:obj.numWells
                if ismember(i,opt.wells)
                    for k=1:size(obj.filteredTraj,2)
                        trajMat=cell2mat(obj.filteredTraj{i,k});   % this concatenates across the cell array. A 5x1 cell array of
                        % 10x8 matrices will be converted into a 50x8 matrix. A similar 5x2 cell array would be converted into a 50x16 matrix. xyTraj should be a column vector.
                        indRem=trajMat(:,frameCol)<opt.minFrame | trajMat(:,frameCol)>opt.maxFrame;
                        trajMat(indRem,:)=[];
                        
                        sameTraj=diff(trajMat(:,trajIDCol))==0;
                        stepOKbool=sameTraj' & (trajMat(1:(end-1),neighborCol)'>opt.neighborDistThresh | trajMat(2:end,neighborCol)'>opt.neighborDistThresh);
                        len=size(trajMat,1);
                        stepGroupLen=2*opt.framesPerStep+opt.frameGap;
                        groupOK=true(1,len); groupOK((len-stepGroupLen+1):len)=false;
                        for mm=1:stepGroupLen
                            indval=1:(len-stepGroupLen);
                            groupOK(indval)=groupOK(indval) & stepOKbool(indval+mm-1);
                        end
                        indStartS1=false(1,len);
                        for m=opt.minFrame:(opt.maxFrame-stepGroupLen)
                            indStartS1=(indStartS1+(groupOK & trajMat(:,frameCol)'==m))>0;
                            groupOK=(groupOK-indStartS1)>0;
                            for mm=1:(stepGroupLen-1)
                                groupOK((1+mm):(len-stepGroupLen+mm))=(groupOK((1+mm):(len-stepGroupLen+mm))-indStartS1(1:(len-stepGroupLen)))>0;
                            end
                        end
                        indEndS1=false(1,len);
                        indEndS1((1+opt.framesPerStep):len)=indStartS1(1:(len-opt.framesPerStep));
                        indStartS2=false(1,len);
                        indStartS2((1+opt.framesPerStep+opt.frameGap):len)=indStartS1(1:(len-opt.framesPerStep-opt.frameGap));
                        indEndS2=false(1,len);
                        indEndS2((1+2*opt.framesPerStep+opt.frameGap):len)=indStartS1(1:(len-2*opt.framesPerStep-opt.frameGap));

                        s1=nan(sum(indStartS1),6);
                        s1(:,1)=trajMat(indEndS1,1)-trajMat(indStartS1,1);
                        s1(:,2)=trajMat(indEndS1,2)-trajMat(indStartS1,2);
                        s1(:,3)=trajMat(indEndS2,1)-trajMat(indStartS2,1);
                        s1(:,4)=trajMat(indEndS2,2)-trajMat(indStartS2,2);
                        s1(:,5)=sqrt(s1(:,1).*s1(:,1)+s1(:,2).*s1(:,2));    %Dist step 1
                        s1(:,6)=sqrt(s1(:,3).*s1(:,3)+s1(:,4).*s1(:,4));    %Dist step 2
                        s2{i,k}=s1;
                    end
                end
            end
        end
        function obj=computeStepSizes(obj)
            obj.steps=obj.stepSizesComputeOnly;
        end
        function stp=trajMat2steps(obj,trajMat,gradCenter,sameTraj,timeVals,framesPerStep)
            stepMat=diff(trajMat);
            timeSteps=timeVals((1+framesPerStep):obj.numFrames)-timeVals(1:(obj.numFrames-framesPerStep));
            if size(stepMat,1)>0
                if nargin<4
                    sameTraj=stepMat(:,end-1)==0;   % The second to last column should be the trajectory identifier
                end  %otherwise, sameTraj has been provided
                stepMat(~sameTraj,:)=NaN;
                
                if isempty(obj.stepCollabels)
                    obj=obj.addStepCollabels;
                end
                stp(:,1)=sqrt(stepMat(:,1).^2 + stepMat(:,2).^2);
                stp(:,3)=stepMat(:,1);
                stp(:,4)=stepMat(:,2);
                stp(:,6)=trajMat(1:(end-1),end);                                % starting frame
                
                m1=[stp(:,3) stp(:,4)]./[stp(:,1) stp(:,1)];
                m2=[gradCenter(1)-trajMat(1:(end-1),1) gradCenter(2)-trajMat(1:(end-1),2)];
%                m3=[gradCenter(1)-trajMat(2:end,1) gradCenter(2)-trajMat(2:end,2)];
                stp(:,5)=sqrt(m2(:,1).^2 + m2(:,2).^2);
                m2=m2./repmat(stp(:,5),[1 2]);
                stp(:,2)=real(acos(sum(m1.*m2,2))/pi*180);
%                stp(:,7)=stp(:,5) - sqrt(m3(:,1).^2 + m3(:,2).^2);
                stp(:,7)=stp(:,3).*m2(:,1) + stp(:,4).*m2(:,2);
                stp(:,8)=trajMat(1:(end-1),end-1);
                if ismember('dist2nearest',obj.coordCollabels)
                    stp(:,9)=min(trajMat(1:(end-1),7),trajMat(2:end,7));
                    %stp(:,9)=stepMat(:,5);
                    stp(:,10)=max(trajMat(1:(end-1),7),trajMat(2:end,7));
                end
                stp(:,11)=sqrt(max(0,stp(:,1).^2 - 2*obj.localizationError^2));   % Distance moved, corrected for localization errors
                stp(~sameTraj,:)=[];
                stp(:,12)=arrayfun(@(x) timeSteps(x),round(stp(:,6)));            % Add the time interval for steps
            else
                stp=[];
            end
        end
        function trajStatMat=computeTrajProperties(obj,well,group,minFrame,maxFrame)
            %computes statistics for trajectories between the specified
            %frame limits
%             gradCenter=[(1+obj.imNumCols) (1+obj.imNumRows)]/2;
            gradCenter=obj.gradCenterXY;
            xyTraj=obj.traj{well,group};
            stp=obj.steps{well,group};
            trajIDList=unique(stp(stp(:,6)>=minFrame & stp(:,6)<maxFrame,8));
            trajIDTotal=cellfun(@(x) max(x(:,end-1)),xyTraj);
            [~,trajInd]=intersect(trajIDTotal,trajIDList);
            xyTraj=xyTraj(trajInd);
            for i=1:length(xyTraj)
                val=xyTraj{i};
                xyTraj{i}=val(val(:,end)>=minFrame & val(:,end)<=maxFrame,:);
            end
            stp=stp(ismember(stp(:,8),trajIDList),:);
            trajStatMat=nan(length(xyTraj),6);
            
            % 1 = trajectory length
            % 2 = net distance traveled
            % 3 = path length
            % 4 = starting distance from center
            % 5 = ending distance from center
            % 6 = net distance towards the center
            % 7 = #6 / #2
            % 8 = #2 / #3
            
            trajStatMat(:,1)=cellfun(@(x) size(x,1),xyTraj);
            trajStatMat(:,2)=cellfun(@(x) sqrt((x(end,1)-x(1,1))^2 + (x(end,2)-x(1,2))^2),xyTraj);
            for i=1:length(xyTraj)
                trajStatMat(i,3)=sum(stp(stp(:,8)==trajIDList(i),1));
            end
            trajStatMat(:,4)=cellfun(@(x) sqrt((gradCenter(1)-x(1,1))^2 + (gradCenter(2)-x(1,2))^2),xyTraj);
            trajStatMat(:,5)=cellfun(@(x) sqrt((gradCenter(1)-x(end,1))^2 + (gradCenter(2)-x(end,2))^2),xyTraj);
            trajStatMat(:,6)=trajStatMat(:,4)-trajStatMat(:,5);
            trajStatMat(:,7)=trajStatMat(:,6)./(trajStatMat(:,2)+0.1);
            trajStatMat(:,8)=trajStatMat(:,2)./(trajStatMat(:,3)+0.1);
        end
        function obj=addCoordCollabels(obj)
            obj.coordCollabels{1}='centX';
            obj.coordCollabels{2}='centY';
            obj.coordCollabels{3}='medianInt';
            obj.coordCollabels{4}='area';
            obj.coordCollabels{5}='integratedIntLog';
            obj.coordCollabels{6}='maxIntLog';
            obj.coordCollabels{7}='dist2nearest';
        end
        function obj=addStepCollabels(obj)
            obj.stepCollabels{1}='dist';
            obj.stepCollabels{2}='angle';
            obj.stepCollabels{3}='deltaX';
            obj.stepCollabels{4}='deltaY';
            obj.stepCollabels{5}='distFromCent';
            obj.stepCollabels{6}='startFrame';
            obj.stepCollabels{7}='chemDotProd';
            obj.stepCollabels{8}='trajID';
            obj.stepCollabels{9}='minDist2Neighbor';
            obj.stepCollabels{10}='maxDist2Neighbor';
            obj.stepCollabels{11}='distCorrected';
            obj.stepCollabels{12}='deltaTime';
        end
        function ind=getCoordColIndex(obj,label)
            ind=find(strcmp(label,obj.coordCollabels));
        end
        function ind=getStepColIndex(obj,label)
            ind=find(strcmp(label,obj.stepCollabels));
        end
        function out=getStepValColInfo(obj,col,prop,well,opt)
            if strcmp(prop,'title')
                for i=1:1; switch col
                        case 1
                            out='Speed';
                        case 2
                            out='Chemotaxis (Angle)';
                        case 3
                            out='deltaX';
                        case 4
                            out='deltaY';
                        case 5
                            out='Distance From Center';
                        case 6
                            out='Frame # at Start';
%                         case 7
%                             out='Chemotaxis (Progress Towards Center)';
                        case 7
                            out='Chemotaxis (Velocity dot Gradient)';
                        case 8
                            out='Traj ID';
                        case 9
                            out='Nearest Neighbor Dist (min)';
                        case 9
                            out='Nearest Neighbor Dist (max)';
                        case 11
                            out='Speed (corrected)';
                        case 12
                            out='Time Interval (sec)';
                    end;
                end
            end
            if strcmp(prop,'ylabel')
                for i=1:1; switch col
                        case {1,7,11}
                            out='Speed (Microns/Min)';
                        case 2
                            out='Angle (Degrees)';
                        case {3,4,5}
                            out='Distance (Microns)';
                        case 6
                            out='Frame #';
                        case 8
                            out='Traj ID';
                    end;
                end
            end
            if strcmp(prop,'units')
                for i=1:1; switch col
                        case {1,7,11}
                            out='Microns/Min';
                        case {3,4,5}
                            out='Microns';
                        case 2
                            out='Degrees';
                        case 6
                            out='Frames';
                        case 8
                            out='Trajectories';
                    end;
                end
            end
            if strcmp(prop,'factor')
                for i=1:1; switch col
                        case {1,7,11}
                            timeVal=opt.framesPerStep*(obj.time{well}(min(length(obj.time{well}),opt.maxFrame))-obj.time{well}(opt.minFrame))/(min(length(obj.time{well}),opt.maxFrame)-opt.minFrame);
                            out=obj.imageMicronsPerPixel/timeVal*60;
                        case {3,4,5}
                            out=obj.imageMicronsPerPixel;
                        case {2,6,8}
                            out=1;
                    end;
                end;
            end
        end
        function ind=getFilteredStepInd(obj,x,col,opt)
            if size(x,1)>0
                ind=~isnan(x(:,col)) & x(:,6)>=opt.minFrame & x(:,6)<=opt.maxFrame-opt.framesPerStep & x(:,1)>=opt.minDistMoved & x(:,5)>=opt.minDistFromCenter & x(:,5)<=opt.maxDistFromCenter;
            else
                ind=[];
            end
        end
        
        
        %Other functions
        function printCoordCollabels(obj)
            fprintf('\nColumn labels for Coords:\n');
            for i=1:length(obj.coordCollabels)
                fprintf('%2i = %s\n',i,obj.coordCollabels{i});
            end
            fprintf('\n');
        end
        function printStepCollabels(obj)
           fprintf('\nColumn labels for Steps:\n');
           for i=1:length(obj.stepCollabels)
                fprintf('%2i = %s\n',i,obj.stepCollabels{i});
           end
           fprintf('\n');
        end
        function path=getWellPath(obj,well)
            path=[obj.folderPath filesep obj.wellName{well}];
        end
        function im=getImage(obj,well,frame,filter)
            if nargin<3
                filter=obj.nucleusColor;
            end
            im=htChemotaxisReadImage(obj.folderPath,obj.wellName{well},filter,frame);
        end
        function saveData(obj,filename)
            trackData=struct(obj);
            if nargin<1
                filename=input('Output filename: ','s');
            end
            if ~endsWith(filename,'.mat')
                filename=[filename '.mat'];
            end
            save([obj.folderPath filesep filename],'trackData');
        end
    end
    methods (Static)
        function labels=makeBarGraphLabels(str,numpts)
            labels=cell(size(str));
            for i=1:length(str)
                labels{i}=[str{i} ' ('];
                for j=1:size(numpts,2)
                    if j>1
                        labels{i}=[labels{i} ', '];
                    end
                    labels{i}=[labels{i} num2str(numpts(i,j))];
                end
                labels{i}=[labels{i} ')'];
            end
        end
        function generateBarGraph(vals,labels,errVals)
            figure;
            bar(vals,'FaceColor',[0.6 0.6 1]);
            %            [vals(:) errVals(:)]
            if nargin>2
                hold on;
                [s1 s2]=size(vals);
                if s2==1
                    errorbar(vals,errVals,'Marker','none','Color',[0 0 0],'LineWidth',2,'LineStyle','none');
                else
                    errorbar((1:s1)-0.15,vals(:,1),errVals(:,1),'+');
                    errorbar((1:s1)+0.15,vals(:,2),errVals(:,2),'+');
                end
            end
            set(gca,'XTick',1:length(vals));
            xlim([0.4 length(vals)+0.6]);
            set(gca,'XTickLabel',labels);
            set(gca,'Box','off');
        end
    end
end