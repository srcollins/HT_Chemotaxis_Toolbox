function defaultParameters = defaultParametersForComputingStatistics
% This function loads default parameters for computing chemotaxis
% statistics from cell tracks

% The following parameters are defined in pixel units. These parameters
% worked well for imaging conditions (4X) where 1 pixel in the image corresponds to 2.5 microns in the sample.
defaultParameters.neighborDistThresh = 16;      % This is a filter to remove from analysis pairs of cells that are too close to each other to be reliably tracked and distinguished in the following frame.
defaultParameters.angMinDistMoved = 4;          % The minimum allowed distance moved by a cell in order to use that movement for the calculation of the angle of movement. This is to avoid noise from non-moving cells.
defaultParameters.minDistFromCenter = 100;      % As gradient steepness approaches zero in the center of the gradient, this filter removes cells near the center of the image which experience little gradient.
defaultParameters.maxDistFromCenter = 650;      % As gradient intensity falls off far from the gradient center, this filter removes cells near the edge of the image that experience a lower chemoattractant concentration.

% Parameters with units of frames
defaultParameters.frames0 = 17;         % This is the number of frames prior to gradient generation by uncaging.
defaultParameters.framesPerStep = 1;    % This defines the time interval used to compute statistics such as speed and angle of movement. A value of 1 indicates that the movement between consecutive frames will be used. Larger values are often used to minimize the effect of cell localization noise. I typically use a value of 2 when imaging with a 30 second interval between frames.

% A microns/pixel parameter is defined in defaultParametersForTracking