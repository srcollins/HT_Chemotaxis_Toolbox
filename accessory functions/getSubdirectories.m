function subdir=getSubdirectories(folder,regExp)
%Returns the names of the subfolders in folder

d0=dir(folder);
subdir={d0([d0.isdir]>0).name};
subdir=setdiff(subdir,{'.','..'});

if nargin>1
    subdir=subdir(boolRegExp(subdir,regExp));
end
