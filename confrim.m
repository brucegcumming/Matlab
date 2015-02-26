function ox = confirm(str, varargin)
% confirm(str) calles questdlg, returns 0 or 1
T = 'Confirm popup';
j = 1;

ok = 0;
yn = questdlg(str, T, 'Yes', 'No');
if strcmp(yn,'Yes')
    ok = 1;
end