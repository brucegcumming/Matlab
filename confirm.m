function ok = confirm(str, varargin)
% confirm(str) calles questdlg, returns 0 or 1
T = 'Confirm popup';
ystr = 'Yes';
nstr = 'No';
j = 1;

ok = 0;
yn = questdlg(str, T, 'Yes', 'No','Yes');
if strcmp(yn,'Yes')
    ok = 1;
end