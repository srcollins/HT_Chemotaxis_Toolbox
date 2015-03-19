function channel=htChemotaxisGetImageChannels(folder)
% htChemotaxisGetImageChannels returns a cell array of the detected image
% channels present in a folder. It assumes a particular file naming
% convention which it uses to detect the image channels. You should either
% use that same convention, or adjust this function so that the image
% channels are detected properly.
% folder should be the path to a folder that contains a representative set
% of images from the experiment (including all imaging channels used).

filenames=htChemotaxisGetSortedImageFilenames(folder);              % Calling this function without specifying nucleusColor will return all filenames from the folder
channel=getTokens(filenames,'([A-Za-z0-9 ]+)_[A-Za-z0-9]+\.tif');   % Looks for an expected pattern at the end of the filename (e.g. xxx_EYFP_20110107T141024.tif)
channel=unique(channel);

