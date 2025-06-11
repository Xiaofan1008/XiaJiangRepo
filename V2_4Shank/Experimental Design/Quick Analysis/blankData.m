function blankData(nChn,FS)
%% This function attempts to blank the raw data
% You must be in the right directory, and have an amplifier data file, a
% trigger data file, and an experimental data file present
if nargin<2
    FS = 30000;
end
if nargin<1
    nChn = 32;
end
%% Load Data
cleanTrig;
trig = loadTrig(0);
theseTrig = trig;
mS = isMultiStim;
ms_shift = 6000; % Milliseconds post stimulation to reintroduce shift artefact
info = dir([pwd filesep 'amplifier.dat']);
info = info.bytes/2;
nL = (ceil(info / (nChn*FS*double(256)))+1)-1;
vFID = fopen([pwd filesep 'amplifier.dat'],'r');
vdnFID = fopen([pwd filesep 'amplifier_dn.dat'],'W');
N = 1;
CHK = split(pwd,'_');
blank = 1;
for i = 1:size(CHK,1)
    if strcmp(CHK{i},'CSD') || isempty(dir('*exp_datafile*.mat'))
        blank = 0;
    end
end
BREAK = 1;
%% Loop through the Data
while (N) && (BREAK)
    v = fread(vFID,[nChn, (FS * 256)],'int16') .* 0.195;
    if ~size(v,2)
        BREAK = 0;
    end
    if (N ~= 1)
        v = [hv v]; %#ok<*AGROW>
    end
    if (blank)
        for iChn = 1:nChn
            [v(iChn,:),~,~] = simpleBlank(v(iChn,:),N,256,theseTrig,1,mS,ms_shift);
        end
    end
    % Save a little data at the start of each loop
    hv = v(:,end-FS+1:end);
    % Don't save the last second, we'll deal with that next time
    if (BREAK)
        v = v(:,1:end-FS);
    end
    fwrite(vdnFID,v./0.195,'int16');
    N = N + 1;
end
%% Close Down
fclose(vFID);
fclose(vdnFID);
end

