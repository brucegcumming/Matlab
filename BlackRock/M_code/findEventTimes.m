function [timestamps snippets indices] = findEventTimes(NEV, Electrode, Unit)

% Given a NEV file, an electrode number, and a unit number this function
% will return indices corresponding to those spikes, timestamps
% corresponding to those spikes (optional), and snippets corresponding to
% those spikes (optional). 
%
% Use [timestamps snippets indices] = findEventTimes(NEV, Electrode, Unit)
% 
% INPUTS
%
%   NEV:          This corresponds to the NEV file the desired data is
%                 being extracted from.
%
%   Electrode:    The electrode number the data is being extracted for.
%
%   Unit:         The unit number the data is being extracted for.
%
% OUTPUT
%
%   indices:      An array of all indices that correspond to neural data
%                 for electrode and unit passed.
%
%   timestamps:   An array of all timestamps that correspond to neural data
%                 for electrode and unit passed.
%
%   snippets:     A matrix of all indices that correspond to neural data
%                 for electrode and unit passed.
%
%   IF OUTPUT IS NOT SPECIFIED ONLY only "timestamps" WILL BE PASSED TO THE
%   CALLING FUNCTION.
%
%   Example: 
%   
%   [timestamps snippets indices] = findEventTimes(NEV, 3, 5);
%
%   In the example above, the indices, timestamps, and snippets for
%   electrode #3 and unit #5 will be passed to the calling function.
%
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 1.1.0 - August 27, 2009

% Find indices that correspond to "Electrode"
ElectrodeIndices = find([NEV.Data.Spikes.Electrode] == Electrode);

% Find indices that correspond to "Unit" within "Electrode" indices
UnitIndices = find([NEV.Data.Spikes.Unit(ElectrodeIndices)] == Unit);

% Updating the indices so they correspond to the original NEV indices
indices = ElectrodeIndices(UnitIndices);

% Finding the timestamps corresponding to the indices
timestamps = NEV.Data.Spikes.Timestamps(indices);

% Finding the snippets corresponding to the indices
if isfield(NEV.Data.Spikes, 'Snippets')
    snippets = NEV.Data.Spikes.Snippets(indices,:);
elseif nargout == 3
    display('Snippet data was not retrieved because the NEV does not contain any snippets.');
    snippets = [];
end
