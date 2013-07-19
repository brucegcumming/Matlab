function compdirectory(dir1, dir2)

%COMPDIRECTORY   GUI tool for directory comparison.
%   COMPDIRECTORY opens a GUI for comparing the contents of 2 directories.
%   It performs recursive comparison for all subdirectories. It checks for
%   file/directory names, file sizes, and modification dates.
%
%   In the GUI, if an entry does not exist in one of the directories, then
%   it marks an 'x' in front of it. If an entry has a different file size,
%   then it marks an 's' in front of it. If an entry has the same file size
%   but a different modification date, then it marks a 'd' in front of it.
%
%   COMPDIRECTORY(DIR1, DIR2) starts the GUI with DIR1 and DIR2
%

% Copyright 2006-2011 The MathWorks, Inc.

% Version:
%   v1.0 - original version
%   v1.1 - cosmetic changes to the GUI (December 2, 2006)
%   v1.2 - fixed bug that emerges when a folder name shows up above . or ..
%          (March 30, 2011)
%
% Jiro Doke
% November 2006-2011

versionNum = '1.2';

delete(findall(0, 'tag', 'compdirFig'));

if nargin == 0
  dir1 = '';
  dir2 = '';
end

if ~ischar(dir1) || ~isdir(dir1)
  dir1 = '';
end
if ~ischar(dir2) || ~isdir(dir2)
  dir2 = '';
end

directories.dir1  = dir1;
directories.dir2  = dir2;
str               = '';
fontsize          = 12;
bgcolor           = [.8 .8 .8];

figH = figure(...
  'name'                            , 'Directory Comparison', ...
  'numbertitle'                     , 'off', ...
  'color'                           , bgcolor, ...
  'tag'                             , 'compdirFig', ...
  'units'                           , 'normalized', ...
  'position'                        , [.1 .1 .8 .8], ...
  'resize'                          , 'on', ...
  'resizefcn'                       , @setUIPos, ...
  'keypressfcn'                     , @keypressFcn, ...
  'menubar'                         , 'none', ...
  'toolbar'                         , 'none', ...
  'handlevisibility'                , 'callback', ...
  'visible'                         , 'off', ...
  'defaultuicontrolunits'           , 'pixels', ...
  'defaultuicontrolfontunits'       , 'pixels', ...
  'defaultuicontrolfontsize'        , fontsize, ...
  'defaultuicontrolfontname'        , 'Verdana', ...
  'defaultuicontrolbackgroundcolor' , 'white');

uh(1) = uicontrol(...
  'style'               , 'pushbutton', ...
  'string'              , 'Compare Directories', ...
  'fontweight'          , 'bold', ...
  'enable'              , 'off', ...
  'backgroundcolor'     , bgcolor, ...
  'parent'              , figH, ...
  'callback'            , @compareDirectories);
uh(2) = uicontrol(...
  'style'               , 'popupmenu', ...
  'string'              , [dir1;{'Select a directory...'}], ...
  'value'               , 1, ...
  'tag'                 , 'dir1', ...
  'parent'              , figH, ...
  'callback'            , @selectDirectory);
uh(3) = uicontrol(...
  'style'               , 'popupmenu', ...
  'string'              , [dir2;{'Select a directory...'}], ...
  'value'               , 1, ...
  'tag'                 , 'dir2', ...
  'parent'              , figH, ...
  'callback'            , @selectDirectory);
uh(4) = uicontrol(...
  'style'               , 'listbox', ...
  'max'                 , 3, ...
  'min'                 , 1, ...
  'fontname'            , 'FixedWidth', ...
  'parent'              , figH, ...
  'callback'            , @selectEntryCallback);
uh(5) = uicontrol(...
  'style'               , 'listbox', ...
  'max'                 , 3, ...
  'min'                 , 1, ...
  'fontname'            , 'FixedWidth', ...
  'parent'              , figH, ...
  'callback'            , @selectEntryCallback);
uh(6) = uicontrol(...
  'style'               , 'listbox', ...
  'max'                 , 3, ...
  'min'                 , 1, ...
  'enable'              , 'inactive', ...
  'horizontalalignment' , 'left', ...
  'fontname'            , 'FixedWidth', ...
  'parent'              , figH);
uh(7) = uicontrol(...
  'style'               , 'text', ...
  'fontweight'          , 'bold', ...
  'horizontalalignment' , 'left', ...
  'backgroundcolor'     , bgcolor, ...
  'string'              , 'Comparison Details:', ...
  'parent'              , figH);
st = {'Compare the contents of 2 directories, including subdirectories.', ...
  '', 'Differences are denoted by the following labels:', ...
  'x:  file or directory does not exist', ...
  's:  different file size', ...
  'd:  different modification date'};
uh(8) = uicontrol(...
  'style'               , 'text', ...
  'fontname'            , 'FixedWidth', ...
  'horizontalalignment' , 'left', ...
  'backgroundcolor'     , bgcolor, ...
  'string'              , sprintf('%s\n',st{:}), ...
  'parent'              , figH);

uh(9) = uipanel(...
  'units'               , 'pixels', ...
  'backgroundcolor'     , bgcolor, ...
  'bordertype'          , 'beveledin', ...
  'parent'              , figH);
uh(10) = uicontrol(...
  'style'               , 'text', ...
  'horizontalalignment' , 'center', ...
  'backgroundcolor'     , bgcolor, ...
  'fontweight'          , 'bold', ...
  'string'              , sprintf('Ver %s', versionNum), ...
  'parent'              , uh(9));
uh(11) = uipanel(...
  'units'               , 'pixels', ...
  'backgroundcolor'     , bgcolor, ...
  'bordertype'          , 'beveledin', ...
  'parent'              , figH);
uh(12) = uicontrol(...
  'style'               , 'text', ...
  'horizontalalignment' , 'right', ...
  'backgroundcolor'     , bgcolor, ...
  'fontweight'          , 'bold', ...
  'string'              , '', ...
  'parent'              , uh(11));

setUIPos;

compareDirectories;

set(figH, 'visible', 'on');

%--------------------------------------------------------------------------
% Nested Subfunctions
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% setUIPos
%   Set the position of the uicontrols (for resizing)
%--------------------------------------------------------------------------
  function p = setUIPos(varargin)

    set(figH, 'units', 'pixels');
    figPos = get(figH, 'position');
        
    % figure can't be too small or off the screen
    if figPos(3) < 750 || figPos(4) < 500
      figPos(3) = max([750 figPos(3)]);
      figPos(4) = max([500 figPos(4)]);
      screenSize = get(0, 'screensize');
      if figPos(1)+figPos(3) > screenSize(3)
        figPos(1) = screenSize(3) - figPos(3) - 50;
      end
      if figPos(2)+figPos(4) > screenSize(4)
        figPos(2) = screenSize(4) - figPos(4) - 50;
      end
      
      set(figH, 'position', figPos);
      
    end
    
    p = [...
      5, figPos(4)-5-fontsize*2, figPos(3)-10, fontsize*2; ...
      5, figPos(4)-5-fontsize*4, (figPos(3)-15)/2, fontsize*1.5; ...
      (figPos(3)+5)/2, figPos(4)-5-fontsize*4, (figPos(3)-15)/2, fontsize*1.5; ...
      5, fontsize*17.5, (figPos(3)-15)/2, figPos(4)-10-fontsize*21.5; ...
      (figPos(3)+5)/2, fontsize*17.5, (figPos(3)-15)/2, figPos(4)-10-fontsize*21.5; ...
      5, 30, figPos(3)-300, fontsize*13.5; ...
      5, 30+fontsize*13.5, figPos(3)-10, fontsize*1.5; ...
      figPos(3)-290, 30, 285, fontsize*13.5; ...
      2, 1, 100, 25; ...
      2, 2, 96, 20;
      102, 1, figPos(3)-102, 25; ...
      2, 2, figPos(3)-106, 20];
    
    set(uh, 'units', 'pixels');
    set(uh, {'position'}, mat2cell(p, ones(1, size(p, 1)), 4));
    
  end


%--------------------------------------------------------------------------
% keypressFcn
%   Increase or decrease the font size in the GUI
%--------------------------------------------------------------------------
  function keypressFcn(varargin)
    k = get(figH, 'currentkey');
    
    switch k
      case 'equal'
        fontsize = min([fontsize + 2, 24]);
      case 'hyphen'
        fontsize = max([10, fontsize - 2]);
      otherwise
        return;
    end
    
    set(uh(1:end-1), 'fontsize', fontsize);
    setUIPos;
    
  end


%--------------------------------------------------------------------------
% compareDirectories
%   Compare the contents of 2 directories
%--------------------------------------------------------------------------
  function compareDirectories(varargin)

    str = '';

    if ~isempty(directories.dir1) && ~isempty(directories.dir2)
      
      set(uh(12), 'string', 'Comparing ...');drawnow;
      
      [l1, l2, ll1, ll2] = compareDirectoriesEngine(directories.dir1, directories.dir2);

      l1 = cellstr([char(ll1),num2str((1:length(l1))'),repmat(': ', length(l1), 1), char(l1)]);
      l2 = cellstr([char(ll2),num2str((1:length(l2))'),repmat(': ', length(l2), 1), char(l2)]);

      set(uh(4) , 'string', l1 , 'value', 1);
      set(uh(5) , 'string', l2 , 'value', 1);
      set(uh(6) , 'string', str, 'value', 1);
      set(uh(12), 'string', '');
    end

    set(uh(1), 'backgroundcolor', bgcolor, 'enable', 'off');


    %----------------------------------------------------------------------
    % compareDirectoriesEngine
    %   The engine for comparing 2 directories
    %----------------------------------------------------------------------
    function [list1, list2, diff1, diff2] = compareDirectoriesEngine(dirname1, dirname2)

      % LIST1 and LIST2 are list of files/directories in dir1 and dir2,
      % respectively. If an entry does not exist in the directory, an empty
      % string is placed. DIFF1 and DIFF2 are cell arrays that indicate
      % differences. No difference is '  ', missing entry is 'x ', different
      % file size is 's ', and different modification date is 'd '.
      
      d1 = dir(dirname1); d1(ismember({d1.name}, {'.', '..'})) = '';
      d2 = dir(dirname2); d2(ismember({d2.name}, {'.', '..'})) = '';
      
      str = sprintf('%sComparing:|(dir1) %s|(dir2) %s||', ...
        str, dirname1, dirname2);

      % separate into files and directories
      d1dirs  = d1([d1.isdir]);
      d2dirs  = d2([d2.isdir]);
      d1files = d1(~[d1.isdir]);
      d2files = d2(~[d2.isdir]);

      list1 = {};
      list2 = {};
      diff1 = {};
      diff2 = {};

      % traverse through all files in dir1
      for id1 = 1:length(d1files)
        list1 = [list1;{d1files(id1).name}];
        id2 = strcmp(d1files(id1).name, {d2files.name});
        if ~any(id2)
          diff1 = [diff1;{'x  '}];
          diff2 = [diff2;{'x  '}];
          list2 = [list2;{''}];
          str = sprintf('%s    "%s" does not exist in dir2.|', str, d1files(id1).name);
        else
          list2 = [list2;{d2files(id2).name}];
          if d1files(id1).bytes ~= d2files(id2).bytes
            diff1 = [diff1;{' s '}];
            diff2 = [diff2;{' s '}];
            str = sprintf('%s    "%s" has a different file size.|', str, d1files(id1).name);
          elseif ~strcmp(d1files(id1).date, d2files(id2).date)
            diff1 = [diff1;{' d '}];
            diff2 = [diff2;{' d '}];
            str = sprintf('%s    "%s" has a different modification date.|', str, d1files(id1).name);
          else % same file
              diff1 = [diff1;{'   '}];
              diff2 = [diff2;{'   '}];
          end
          % remove from the list
          d2files(id2) = '';
        end
      end

      % traverse through the remainding files in dir2
      % these files do not exist in dir1
      list2 = [list2;{d2files.name}'];
      list1 = [list1;repmat({''},length(d2files),1)];
      diff1 = [diff1;repmat({'x  '},length(d2files),1)];
      diff2 = [diff2;repmat({'x  '},length(d2files),1)];
      for id2 = 1:length(d2files)
        str = sprintf('%s    "%s" does not exist in dir1.|', str, d2files(id2).name);
      end
      str = sprintf('%s|', str);

      % traverse through all directories in dir1
      subDirs = [];
      for id1 = 1:length(d1dirs)
        id2 = strcmp(d1dirs(id1).name, {d2dirs.name});
        if ~any(id2)
          list1 = [list1;{['\', d1dirs(id1).name]}];
          list2 = [list2;{''}];
          diff1 = [diff1;{'x  '}];
          diff2 = [diff2;{'x  '}];
          str = sprintf('%s    directory "%s" does not exist in dir2.|', str, d1dirs(id1).name);
        else
          subDirs = [subDirs; id1, id2];
        end
      end

      % traverse through the remainding directories in dir2
      for id2 = 1:length(d2dirs)
        id1 = strcmp(d2dirs(id2).name, {d1dirs.name});
        if ~any(id1)
          list1 = [list1;{''}];
          list2 = [list2;{['\', d2dirs(id2).name]}];
          diff1 = [diff1;{'x  '}];
          diff2 = [diff2;{'x  '}];
          str = sprintf('%s    directory "%s" does not exist in dir1.|', str, d2dirs(id2).name);
        end
      end

      % recursive comparison of subdirectories
      for ii = 1:size(subDirs, 1)
        list1 = [list1;{['\', d1dirs(subDirs(ii, 1)).name]}];
        list2 = [list2;{['\', d2dirs(subDirs(ii, 2)).name]}];
        if ~strcmp(d1dirs(subDirs(ii, 1)).date, d2dirs(subDirs(ii, 2)).date)
          diff1 = [diff1;{' d '}];
          diff2 = [diff2;{' d '}];
          str = sprintf('%s    directory "%s" has a different modification date.|', str, d1dirs(subDirs(ii, 1)).name);
        else
          diff1 = [diff1;{'   '}];
          diff2 = [diff2;{'   '}];
        end
        str = sprintf('%s|', str);
        [l1, l2, ll1, ll2] = compareDirectoriesEngine(...
          fullfile(dirname1, d1dirs(subDirs(ii, 1)).name), ...
          fullfile(dirname2, d2dirs(subDirs(ii, 2)).name));
        if ~isempty(l1)
          list1 = [list1;cellstr([repmat(' ', length(l1), 2), char(l1)])];
          list2 = [list2;cellstr([repmat(' ', length(l2), 2), char(l2)])];
          diff1 = [diff1;ll1];
          diff2 = [diff2;ll2];
        end
      end
    end
  end


%--------------------------------------------------------------------------
% selectEntryCallback
%   Update selection to show the same entry on both lists
%--------------------------------------------------------------------------
  function selectEntryCallback(varargin)

    val = get(varargin{1}, 'value');
    listtop = get(varargin{1}, 'listboxtop');

    set(uh(4:5), 'value', val, 'listboxtop', listtop);

  end


%--------------------------------------------------------------------------
% selectDirectory
%   Select directory for comparison
%--------------------------------------------------------------------------
  function selectDirectory(varargin)

    dirs = cellstr(get(varargin{1}, 'string'));
    val = get(varargin{1}, 'value');

    tagname = get(varargin{1}, 'tag');

    if strcmp(dirs{val}, 'Select a directory...')  % select new directory
      d = uigetdir(directories.(tagname), 'Select a directory');
      if ischar(d)
        dirs = [dirs(1:end-1);{d};dirs(end)];
        set(varargin{1}, 'string', dirs, 'value', length(dirs)-1);

        directories.(tagname) = d;
        set(uh(4:6), 'string', {}, 'value', 1);
      else % canceled
        if ~isempty(directories.(tagname))
          id = find(strcmp(directories.(tagname), dirs));
        else
          id = length(dirs);
        end
        set(varargin{1}, 'value', id);
        return;
      end
    else % selected existing directory from the list
      directories.(tagname) = dirs{val};
      set(uh(4:6), 'string', {}, 'value', 1);
    end
    
    if ~isempty(directories.dir1) && ~isempty(directories.dir2)
      set(uh(1), 'backgroundcolor', [1 .75 .75], 'enable', 'on');
    else
      set(uh(1), 'backgroundcolor', bgcolor, 'enable', 'off');
    end
    
  end
end




%#ok<*AGROW>