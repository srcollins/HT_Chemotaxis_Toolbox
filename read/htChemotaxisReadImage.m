function im=htChemotaxisReadImage(folderPath,well,nucleusColor,frameNumber)
% htChemotaxisReadImage is a function to read a single image from the disk
% and return the result. This version is written with a particular file
% organization and naming strategy in mind. You should either copy this
% strategy, or edit this function (and/or the functions called by this
% function) to make sure that your images are read in properly.
%
% folderPath -- is the source folder for all of the images in the experiment.
%       It is expected to have one subfolder for each well.
% well -- should indicate the well for the desired image (e.g. 'A02')
% nucleusColor -- should indicate the color channel for the image (e.g. 'RFP')
% frameNumber -- should indicate the desired frame (e.g. 7)

wellFolder = [folderPath filesep well];
imageFilenames = htChemotaxisGetSortedImageFilenames(wellFolder,nucleusColor);

im=imread([wellFolder filesep imageFilenames{frameNumber}]);
