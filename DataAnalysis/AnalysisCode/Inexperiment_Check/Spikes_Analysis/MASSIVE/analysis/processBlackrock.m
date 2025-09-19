function processBlackrock(chk,par)
dbstop if error
if ~nargin
    chk = 0;
    par = 0;
elseif nargin == 1
    par = 0;
end
%% Step One: Initialise
% If specified, hard-delete old versions of files
if (chk)
    fclose('all');
    if ~isempty(dir('*ns6_dnBR.dat'))
        delete('*ns6_dnBR.dat');
    end
    if ~isempty(dir('*ns6.muBR.dat'))
        delete('*ns6.muBR.dat');
    end    
end
if (par)
    hPool = gcp;
    if isempty(hPool)
        hPool=parpool('local');
    end
end
% Work out what input type we used.
ifile = pwd;
if ~isempty(dir('*.ns6'))
    fileinfo = dir([ifile '\*.ns6']);
    dName = fileinfo.name;
end
try
    header = openNSx('noread',dName);
    nChn = header.MetaTags.ChannelCount;
    FS = header.MetaTags.SamplingFreq;
    % Calculate an appropriate memory size. Bigger is faster.
    t_len = fileinfo.bytes/(nChn * 2 * FS);
catch
    % We're probably running a conjoined dataset - there is no nsx file
    nChn = 96; FS = 30000;
    t_len = 1e9;
end
try
    % This is just to make sure that running over the entire datafile at
    % once (which is theoretically fastest) doesn't exceed the known
    % available system memory.
    mem = memory;
    T = mem.MaxPossibleArrayBytes ./ (2 * nChn * FS);
    if T > t_len
        T = t_len + 1;
    end
    if T > 256
        T = 256;
    end
catch
    T = 256; % If we aren't in Windows, default to a few minutes
end
%% Step Two: Blank the Artefact
disp(['Analysing the ' dName ' file']);
if ~isempty(dir('*_exp_datafile_*')) && isempty(dir('*ns6_dnBR.dat'))
    denoiseBlackrock(ifile,dName,nChn,FS,T,par);
end
%% Step Three: Generate mu, lfp, and spiking files
allExtractBR(dName,nChn,FS,T,par);
%% Step Four: Delete unnecessary files
fclose('all');
if (par)
    delete(hPool);
end
end