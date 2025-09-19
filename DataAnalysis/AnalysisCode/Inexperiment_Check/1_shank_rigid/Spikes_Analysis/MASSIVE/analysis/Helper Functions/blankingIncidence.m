function chnTotal = blankingIncidence
% This function assesses the incidence rate of various blanking metrics
% across a datafile
%% Variables
datafile = 'X:\Tim\RAT0017\estim_pen1_010_190401_201325\amplifier_dn.dat';
FS = 30000;
nChn = 32;
d = Depth(3);
fFactor = 8; % Flag factor - larger is harder to trip
%% Loading data
% Calculate an appropriate memory size. Bigger is faster.
fileinfo = dir(datafile);
t_len = fileinfo.bytes/(32 * 2 * 30000);
try
    % This is just to make sure that running over the entire datafile at
    % once (which is theoretically fastest) doesn't exceed the known
    % available system memory.
    mem = memory;
    T = cast((mem.MaxPossibleArrayBytes ./ (2 * 32 * 30000)),'int64');
    if T > t_len
        T = t_len + 1;
    end
catch
    T = 124; % If we aren't in Windows, default to a few minutes
end
T = cast(T,'int32');
nT = t_len/T;
trig = loadTrig .* 30;
nTrig = length(trig);
tStart = 1;
%% Output data structures
local_minima = zeros(nChn,nTrig);
%% Logic
for l = 1:nT
    % For all data
    dataFID = fopen(datafile,'r');
    v = fread(dataFID,[nChn, T*FS],'int16') .* 0.195;
    for c = 1:nChn
        % For all channels
        for t = tStart:nTrig
            try
                % For all trigger lines
                thisData = v(d(c),trig(t)-120:trig(t)+600);
                thisDiff = abs(diff(thisData));
                chk(1) = mean(thisDiff(400:end));
                chk(2) = std(thisDiff(400:end));
                pt = chk(1) + fFactor*chk(2);
                if any(thisDiff(1:225) > pt)
                    local_minima(c,t) = 1;
                end
            catch
                tStart = t;
                break;
            end
        end
    end
end
%% Print summary stats
nTotal = sum(sum(local_minima));
chnTotal = sum(local_minima,2);
disp([num2str((100.*nTotal)./(nChn*nTrig)) '% of trials triggered a flag']);
for c = 1:nChn
    disp([num2str((100.*chnTotal(c))./(nTrig)) '% of trials triggered a flag on channel ' num2str(c)]);
end
%% Closing
fclose(dataFID);
end