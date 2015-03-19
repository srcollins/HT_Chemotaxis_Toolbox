function [filenames,datenums]=htChemotaxisGetSortedImageFilenames(wellFolder,nucleusColor)
% htChemotaxisGetImageFilenames returns a list of the image filenames from
% the specified folder, filtering to return only those images taken in the
% color channel specified by nucleusColor. This function is written with a
% particular file naming convention in mind. You should either copy this
% naming convention, or adjust this function so that your images can be
% read in properly.
%
% datenums is an optional second output which should return the timestamp
% from the images. This function uses the date modified returned from
% Matlab's dir function to get the timestamp. The timestamp could
% alternately be extracted from either the image filename or from a
% metadata file.

% Get all the filenames
fileStruct=dir(wellFolder);
fileStruct=fileStruct([fileStruct.isdir]==0);
fileStruct=fileStruct(~ismember({fileStruct.name},{'Thumbs.db','thumbs.db'}));

% Filter out only ones with the correct color channel
if nargin>1
    fileStruct=fileStruct(boolRegExp({fileStruct.name},nucleusColor));
end

% Sort the images according to frame order
[sortedDateNums,sortOrder]=sort([fileStruct.datenum]);
fileStruct=fileStruct(sortOrder);

% Extract the filenames
filenames={fileStruct.name};
filenames=filenames(:);
datenums=vect([fileStruct.datenum]);
