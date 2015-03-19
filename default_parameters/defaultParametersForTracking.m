function defaultParameters = defaultParametersForTracking
% This function loads default parameters for cell tracking which are
% hard-coded in this file

% Parameters with units of microns/pixel
defaultParameters.imageMicronsPerPixel = 2.5;

% The following parameters are defined in pixel units. These parameters
% worked well for imaging conditions (4X) where 1 pixel in the image corresponds to 2.5 microns in the sample.
defaultParameters.cellMinArea = 6;      % The minimum number of pixels for an object to be considered a cell
defaultParameters.cellMaxArea = 75;     % Objects larger than this will be filtered out as unlikely to be single cells (this filtering is after object detection)
defaultParameters.largeFilt = 200;      % Objects larger than this size will be filtered out prior to object detection (which includes watershed separation of potentially touching cells). This filter is intended to remove large debris which may occasionally make it into a well.
defaultParameters.bgBlockSize = 32;     % This parameter is for local background subtraction. Square blocks of N x N pixels are used to compute local background intensity.
defaultParameters.maxDisp=15;           % The maximum displacement a cell could move between consecutive frames. This is a search limit for the algorithm matching cells between frames.
defaultParameters.localizationError = 0.5;                % The average error in estimate of centroid coordinates (in pixel units). This will be used to apply a correction to measured distances for localization error

% Parameters with units of percent
defaultParameters.bgPercentile = 80;    % This is the pixel intensity percentile which will be used to determine the local background. The corresponding pixel intensity will be subtracted from the local region of the image. Resulting negative values will then be set to zero. A value greater than 50 is typically used to suppress spurious noise.

% Parameters with units of intensity
defaultParameters.thresh0=180;          % This is an intensity threshold used to distinguish cell from non-cell. It will be highly dependent on your imaging conditions, but should be fairly consistent from day to day if the same conditions are used. I typically adjust this value slightly using the adjustThreshForNucleusDetection function for every experiment.

% Parameter for defining the gradient direction
defaultParameters.gradCenterXY='imageCenter';   % This parameter defines the gradient direction. 'imageCenter' will result in automatic calculation of the image center, and this value will be used as the center of the gradient. Alternately, numerical coordinates in pixels may be defined.

