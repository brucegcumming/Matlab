function [expname, exptypename, suff, stimname]  = Expt2Name(Expt, varargin)
% take an expt and return a name identifying the
% exptype, of the form rds.dxXce  
SpkDefs;
suff = [];
if isfield(Expt.Stimvals,'st')
    if ischar(Expt.Stimvals.st)
    stimname = Expt.Stimvals.st;
    else
    stimname = stimnames{Expt.Stimvals.st+1};
    end
else
    stimname = '';
end

    if strmatch(Expt.Stimvals.e2, 'e0')
        exptypename = Expt.Stimvals.et;
        if strmatch(Expt.Stimvals.e3, 'e0')
            exptypename = Expt.Stimvals.et;
        else
            exptypename = [Expt.Stimvals.et 'X' Expt.Stimvals.e3];
        end
    else
        if strmatch('or',Expt.Stimvals.et) & Expt.Stimvals.ei > 100 & strcmp(Expt.Stimvals.e2,'me')
            exptypename = 'dirXme';
        else
        if strmatch(Expt.Stimvals.e3, 'e0')
            exptypename = [Expt.Stimvals.et 'X' Expt.Stimvals.e2];
        else
            exptypename = [Expt.Stimvals.et 'X' Expt.Stimvals.e2 'X' Expt.Stimvals.e3];
        end
        end
    end
    
    if isfield(Expt,'Header')
     
    if isfield(Expt.Header,'psych') && Expt.Header.psych
        exptypename = [exptypename 'P'];
    end
        if isfield(Expt.Header,'rc') && Expt.Header.rc
        suff = 'RC';
    end
    end
    if strcmp(Expt.Stimvals.et,'or') && Expt.Stimvals.ei > 100
        exptypename = strrep(exptypename, 'or','dir');
    end
    
    if strcmp(Expt.Stimvals.et,'tf') & crtrial > nt/2 %Counterphase
        exptypename = ['C' exptypename];
    end
    if isfield(Expt.Stimvals,'rb') & Expt.Stimvals.rb ~= 0
        exptypename = [exptypename 'RB'];
    end

    if strmatch(exptypename,{'dxXId' 'dxXIdP' 'orXId' 'dirXId'},'exact')
        jx = GetEval(Expt,'jx');
        if jx > 0
         exptypename = [exptypename 'D'];
        end
        if strmatch(Expt.Stimvals.Bs,'cylinder')
         exptypename = [exptypename 'B'];
        end
    end
    if strmatch(exptypename,{'OpXdx' 'PpXdx'},'exact')
        Bs = GetEval(Expt,'Bs');
        bh = GetEval(Expt,'bh');
        st = GetEval(Expt,'st');
        sz = GetEval(Expt,'sz');
        if Bs == st && bh >= sz
            exptypename = [exptypename];
        else
            exptypename = [exptypename 'noback'];
        end
            
        
    end
    if strncmp(exptypename,'dxXce',5)
        if Expt.Stimvals.n2 > 2
            exptypename = strrep(exptypename,'dxXce','dxXces');
        elseif Expt.Stimvals.i2 < 2 %not +- 1
            exptypename = strrep(exptypename,'dxXce','dxXces');
        end
    end

    expname = [stimname '.' exptypename];


    
