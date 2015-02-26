function path = GetFilePath(type, varargin) 
%path GetFilePath(type,...)  returns path for datafiles, preferences, etc
%So that future changes in server disk names are less disruptive.
%Also allows users to customize these, by replacing this function
%GetFilePath('data') returns prefix for data files ('/b/data')
%GetFilePath('bgcmatlab') returns prefix for bruces matlab direcotyr ('/b/bgc/matlab')
%GetFilePath('preferences') returns prefix for config/layout files ('/b/group/matlab/prefernces')
%GetFilePath('binoclean') returns prefix for binoclean source  ('/b/binoclean')
%GetFilePath('perl') returns prefix for perl scripts  ('/b/bgc/perl')
%GetFilePath('anal') returns prefix for perl scripts  ('/b/bgc/anal')
%GetFilePath('group') returns prefix for perl scripts  ('/b/group')

path = [];
prefix = '/b';

if strcmp(type,'data')
    path = [prefix '/data'];
end
if strcmp(type,'preferences')
    path = [prefix '/group/matlab/preferences'];
end
if strcmp(type,'binoclean')
    path = [prefix '/binoclean'];
end
if strcmp(type,'group')
    path = [prefix '/group'];
end
if strcmp(type,'perl')
    path = [prefix '/bgc/perl'];
end
if strcmp(type,'anal')
    path = [prefix '/bgc/anal'];
end
if strcmp(type,'bgcmatlab')
    path = [prefix '/bgc/matlab'];
end

path = CheckNameBug(path);