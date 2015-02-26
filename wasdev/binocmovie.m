function mov = binocmovie(dir,varargin)
% BINOCMOVIE make avi file out of image frames from Binoc
%   mov = BINOCMOVIE(dir) read images from 'dir' and 
%   output frame stack to 'mov'. The filenames should be numbered sequentially,
%   and EM data (if desired) should be read in by Binoc and output
%   as a black cross.
%
%   mov = BINOCMOVIE(dir,'prefix',str) search for image 
%   files that begin with str (current default is 'ORBWmovie.*')
%
%   mov = BINOCMOVIE(dir,'filename',name) write output avi file
%   to 'name', default is 'ORBWmovie.avi' in 'dir' 
%
%   mov = BINOCMOVIE(dir,'play') play movie after writing
%
%   mov = BINOCMOVIE(dir,'emcolor') find black eye position 
%   cross and change to red
%
%   mov = BINOCMOVIE(dir,'framerate',rate) set frame rate 
%   (in Hz) to 'rate'
%
%   mov = BINOCMOVIE(dir,'repeat',n) repeat frame sequence n times
%
%   Adrian Bondy, 2013

emcolor=0;
j=1;
filename='ORBWmovie.avi';
play=0;
nrpts=1;
framerate=30;
verbose = 0;
prefix='ORBWmovieim.*';
Lframestack = [];
cmap = [];

while j<=length(varargin)
    if strncmpi(varargin{j},'emcolor',3)
        emcolor=1;
    elseif strncmpi(varargin{j},'filename',5)
        j=j+1;
        filename=varargin{j};
    elseif strncmpi(varargin{j},'play',4)
        play=1;
    elseif strncmpi(varargin{j},'prefix',5)
        j=j+1;
        prefix=varargin{j};
    elseif strncmpi(varargin{j},'framerate',6)
        j=j+1;
        framerate=varargin{j};
    elseif strncmpi(varargin{j},'repeat',5)
        j=j+1;
        nrpts=varargin{j};
    elseif strncmpi(varargin{j},'verbose',5)
        verbose = 1;
    end
    j=j+1;
    
end

PrintMsg(0,'Finding image files.');
framefiles=TreeFind(dir,'name',prefix);

% get frame number from file names and read frames
PrintMsg(0,'Sorting and reading in %s images.',num2str(length(framefiles)));
nf = 0;
for f=1:length(framefiles)
    im=imread(framefiles{f});
    if regexp(framefiles{f},'[0-9]R.pgm')
        nf = nf+1;
        fnum(nf) = str2num(regexprep(framefiles{f},'.*m([0-9]+)[A-Z].pgm','$1')) + 1;
        framestack(:,:,nf)=flipud(im);
        lfile = regexprep(framefiles{f},'R.pgm','L.pgm');
        if exist(lfile)
            im=imread(lfile);
            Lframestack(:,:,nf)=flipud(im);
        elseif ~isempty(Lframestack)
            fprintf('Missing L Image %s\n',lfile);
        end
    end
end
[~,order] = sort(fnum);
framestack = framestack(:,:,order);  % reorder frame stack    
    
% get frames
PrintMsg(0,'Getting movie frames.');
vidobj = VideoWriter([dir,'/',filename]);
set(vidobj,'FrameRate',framerate,'Quality',100);
open(vidobj);
figure;
if isempty(cmap) && isempty(Lframestack)
    colormap('gray');
end
set(gcf,'Position',[100 100 size(framestack,2) size(framestack,1)],'renderer','openGL');
counter=0;
for n=1:nrpts
for f=1:1:size(framestack,3)
    counter=counter+1;
    currframe=squeeze(framestack(:,:,f));
    if ~isempty(Lframestack) && sum(size(Lframestack) == size(framestack)) == 3 %all 3 dimensions match
        currframe = cat(3, currframe, zeros(size(currframe)), squeeze(Lframestack(:,:,f)));
    end
    a=find(currframe==128);
    currframe(a)=127;    
    if emcolor
        [i j] = find(currframe==min(currframe(:)));
        currframe=cat(3,currframe,currframe,currframe);
        if length(i)<500 & length(i)>40
            for k=1:length(i)
                currframe(i(k),j(k),1)=255;
                currframe(i(k),j(k),[2 3])=0;
            end
        end
    end
    image(currframe);
    if verbose
        fprintf('%d:%s\n',f,framefiles{order(f)});
    end
    set(gca,'XTick',[],'YTick',[]);    
    mov(counter)=getframe;
    writeVideo(vidobj,mov(counter));
end
end
close(vidobj);

if play
    PrintMsg(0,'Playing movie.');
    movie(mov);
end
    