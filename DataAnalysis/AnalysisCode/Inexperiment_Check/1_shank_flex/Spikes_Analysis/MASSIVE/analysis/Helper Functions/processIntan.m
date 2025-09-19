function processIntan(chk,par)
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
    if ~isempty(dir('amplifier_dn.dat'))
        delete('amplifier_dn.dat');
    end
    if ~isempty(dir('*.mu.dat'))
        delete('*.mu.dat');
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
if exist('analogin.dat','file')
    dName = 'analogin';
elseif exist('amplifier.dat','file')
    dName = 'amplifier';
end
% Calculate an appropriate memory size. Bigger is faster.
fileinfo = dir([ifile '\' dName '.dat']);
t_len = fileinfo.bytes/(32 * 2 * 30000);
try
    % This is just to make sure that running over the entire datafile at
    % once (which is theoretically fastest) doesn't exceed the known
    % available system memory.
    mem = memory;
    T = mem.MaxPossibleArrayBytes ./ (2 * 32 * 30000);
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
if ~isempty(dir('*_exp_datafile_*')) && isempty(dir([dName '_dn.dat']))
    denoiseIntan(ifile,dName,T,par);
end
%% Step Three: Generate mu, lfp, and spiking files
allExtract(dName,T,par);
%% Step Four: Delete unnecessary files
fclose('all');
if ~isempty(dir('time.dat'))
    delete('time.dat');
end
if ~isempty(dir('stim.dat'))
    delete('stim.dat');
end
if ~isempty(dir('amplifier_lfpraw.dat'))
    delete('amplifier_lfpraw.dat');
end
if (par)
   delete(hPool);
end
end