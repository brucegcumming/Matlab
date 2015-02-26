function RebuildImages(Expt, varargin)
%RebuildImages(Expt) Builds images for an ORBW Expt

or = [];
filterargs = {};
px  = 0.0188; %was size on Berlioz for V1 expts, and on Ravel.
sf = Expt.Stimvals.sf;
wi = Expt.Stimvals.wi;
outname = [];
if isfield(Expt.Header,'loadname')
    outname = regexprep(Expt.Header.loadname,'ORBW.*','ORBW.Images.mat');
    if strcmp(outname,Expt.Header.loadname)
        outname = [];
    end
end

if isfield(Expt.Stimvals,'impx') & Expt.Stimvals.impx < 1
    px = Expt.Stimvals.impx;
end


if isempty(or)
or = Expt.Stimvals.or;
ors = unique([Expt.Trials.or]);
or = min(ors);
end
obs = unique([Expt.Trials.ob]);

nf = 0;
Images = {};
for j = 1:length(obs)
    if obs(j) > 120
        I = BuildImages(Expt, ors(1), obs(j), px, sf, wi);
        Images = {Images{:} I{:}};
    else
        for o = 1:length(ors)
            I = BuildImages(Expt, ors(o), obs(j), px, sf, wi);
            Images = {Images{:} I{:}};
        end
    end
end


if ~isempty(outname) 
    save(outname,'Images','-v7.3');
end


function Images = BuildImages(Expt, or, ob, px, sf, wi)
zid = find([Expt.Trials.or] == or & [Expt.Trials.ob] == ob & abs([Expt.Trials.RespDir]) ==1);
[seedcounts, details.seedoffsets] = Counts([Expt.Trials(zid).imseed]);

if length(seedcounts) > 1
    fprintf('Seeds for %.0f,%.0f: %s Counts %s\n', or, ob, sprintf(' %d', details.seedoffsets),sprintf(' %d', seedcounts));
end

for s = 1:length(seedcounts)
    nfts = filterim([sf sf/2], [or ob], wi, 'seedoffset', details.seedoffsets(s), 'pix2deg', px,'nseeds', 1000, 'getimages','noplot');
    Images{s}.or = or;
    Images{s}.ob = ob;
    Images{s}.imseed = details.seedoffsets(s);
    Images{s}.px = px;
    Images{s}.Images = nfts;
end
