function varargout = AllVPcs(varargin)
%AllVPcs(V, ...)  Extract and Classify Spikes from continuous voltage record
%Takes an MxN matrix (electrode by voltage) of continuous
%voltages, extracts segments triggered off one row, and plots PCS.
%  V can also be a filename of a FullV file.  Will only load if new. 
%
%AllVPcs(V,'tchan',c,..    uses Channel c to find trigger points.
%AllVPcs(V,'tchan',c,'reclassify')  uses saved clustering and exaclty
%recaptiulates
%AllVPcs(V,'tchan',c,'reapply')  uses saved clustering parameters, but
%appies them to new data (e.g. changes in trigger, new probe).
%AllVPcs(V,'tchan',c,'reapply', Clusters{p}) uses cluster given 
%AllVPcs(V,'tchan',c,'usecluster') Doesn't recalculate space. Just usesAuto
%                    times and classifications from ClusterTimes
%AllVPcs(V,'tchan',c,'usecluster', force) uses a list of time indices
%        in 'force' as the event list. If force is a strcutre, then
%        force.tid gives time indices
%        force.clst gives classification
%
%AllVPcs(name, 'tchan', c, 'refcut')    Applies the cluster defined in RefClusters.mat, if no cluster is yet defined
%AllVPcs(name, 'tchan', c, 'refclusterforce') forces application of cluster defined in RefClusters.mat;
%
%AllVPcs(V, ..., 'spkrate', R) sets the threshold (what evet type) to
%produce a mean of R events per second
%AllVPcs(V, ..., 'th+') triggers off maxima
%AllVPcs(V, ..., 'th-') triggers off minima (default)
%
%AllVPcs(V, 'tchan',c,'trigdt') or AllV.AllVPcs(V, 'tchan',c,'dtthr')  trigggers off dVdt
%AllVPcs(V, 'tchan',c,'dtthr2') Triggers off second temporal derivative
%
%AllVPcs(V, 'tchan',c,'dtthr3') Triggers off energy.  N.B. This is dagnerous, because 
%  of an attempt to align spikes. The trigger point is moved to the nearest max (th+) or min (th-)
%  in the voltage record.  The means that the histogram of trigger values is altered (measures energy at the new
%  trigger point, whihc is different, so no longer get a hard edge to the histogram.
%   ....'dtthr3','thboth') disables the voltage alignment so that the tirgger values are correct
%  but this can produce artifactual clustering...
%
%
%AllVPcs(V, 'tchan',c, 'triggertemplate', Clusters{p}) convolves trgger
%    channel with template in cluster, triggers on peaks
%
%
%AllVPcs(V, 'tchan',c,'pcchan',P)  uses the channels defined in vector P
%  to build PCA scores, and summed template scores. Make the second element
%  of P your trigger channel - plots by probe are P(2) vs P(1) and P(2) vs
%  P(end)
%
%AllVPcs(V, 'tchan',c,'smooth', sigma)  Smooths the trigger channel with a
%Gaussain before triggering. NB this smoothing is NOT applied to the data
%then used after triggering. 
%
%By default, AllV.AllVPcs only includes events that happened in a trial
%        (+pre/postperiod). To include more events:
%...,'usealltrials') includes trials terminated by badfix etc
%...,'allspikes')  includes all spikes, ignoring teh expt
%       (But be sure teh FullV has everything , ProcessGridFullV(...,'nochopfile'

varargout{:} = AllV.AllVPcs(varargin{:});
