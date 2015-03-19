function filenames=getFilenames(path,regExp)
% returns a list of the filenames in a specified directory
% the second parameter can be used to specify a regular expression required
% in the filenames

d0=dir(path);
filenames={d0([d0.isdir]==0).name};
filenames=filenames(:);
filenames=setdiff(filenames,{'Thumbs.db','thumbs.db'});

if nargin>1
    filenames=filenames(boolRegExp(filenames,regExp));
end
