function DATA = CheckState(DATA, varargin)
%check for incompatible options, and set options defined
%by ocnjuntion of GUI options

if DATA.state.showspkxy && DATA.state.showspikes == 0 && DATA.plot.quickspks == 0
    DATA.state.usexycache = 1;
else
    DATA.state.usexycache = 0;
end
