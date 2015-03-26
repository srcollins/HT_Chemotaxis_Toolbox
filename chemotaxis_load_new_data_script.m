%% Script for loading and processing new ht chemotaxis experiments
% This script is an example of how to load data from HT Chemotaxis
% experiments. It is set up to load data stored in a master folder root,
% in which each subfolder represents a different 96-well plate experiment.

root = 'c:\ht-chemotaxis-data\project1';
folders=getSubdirectories(root);    % Each of these subfolders is expected to contain data for one 96-well plate experiment
for i=1:length(folders)
    folders{i}=[root folders{i}];   % Build full paths to each of the experiment data folders
end
folders(:)

for i=1:length(folders)
    fprintf('-----\nCont cells: %s\n', folders{i});     % Indicate that we are loading the data for the control cells
    t{i,2}=HTChemotaxisExperiment(folders{i});              % Load the control cell data - the user will be prompted to enter the correct imaging channel (e.g. YFP)
end
for i=1:length(folders)
    fprintf('-----\nExp cells: %s\n', folders{i});     % Indicate that we are loading the data for the experimental cells
    t{i,1}=HTChemotaxisExperiment(folders{i});             % Load the control cell data - the user will be prompted to enter the correct imaging channel (e.g. RFP)
end
%%  Test the nucleus detection threshold
i=4; well=36; frame=1;
t{i,1}.showOverlayWithDetectedNucleiWatershed(well,frame,[0 500]);  % This will show an image with detected nuclei. The displayed image will be a background subtracted image, with an image intensity scale determined by the third input.
%%  Adjust the nucleus detection threshold based on results from above
for i=4
    t{i,1}=t{i,1}.adjustThreshForNucleusDetection;      % If visible cells are being missed, the threshold should be lowered. If imaging noise is erroneously being detected as nuclei, the threshold should be raised.
end
%%  Process the data to detect nuclei, track cells, remove trajectories for cells that never move, and finally save the processed data to disk
for i=1:size(t,1)                                                   % Iterate over the different experiments
    for k=1:size(t,2)                                               % Iterate to process both the control and experimental cells
        fprintf('%s  %s\n',t{i,k}.nucleusColor,t{i,k}.folderPath);  
        t{i,k}=t{i,k}.findNuclei;
        t{i,k}=t{i,k}.trackCells;
        t{i,k}=t{i,k}.filterTrajectories(5);
        t{i,k}.saveData(['tracks_' t{i,k}.nucleusColor '.mat']);
    end
end
%%  Compute statistics of cell movement
frames0=17;                                                         % This will be the number of frames before the initial gradient generation
mer=computeChemotaxisStats(t,frames0);                              % Note that many important parameters, including properties of the imaging system, are set in defaultParametersForComputingStatistics. Alternately, a number of them can be specified as inputs to this function.
mer.rowlabels=readChemotaxisCoordMap([root filesep 'example coordmap.csv']);    % This will load labels for the wells from a coordinate map file
[n,mer2]=rearrangeAndNormalizeTwoColorChemotaxisData(mer);
%% ------ View some of the data -----

%% Plot trajectories for one well
wellNum=1;
t{1,1}.plotplotXYTraj(wellNum);

%% Look at CDFs of basal speed
plotCDFs(mer.s0);

%% Scatter plot raw speed data for before and after attractant stimulation
scatterChemotaxisData(mer,{'s0','s1'},1);  % This plot will have clickable data points, which will show you the corresponding rowlabels

%% Scatter plot normalized angular bias for two different experiments
scatterChemotaxisData(n,'a1',1,2);
