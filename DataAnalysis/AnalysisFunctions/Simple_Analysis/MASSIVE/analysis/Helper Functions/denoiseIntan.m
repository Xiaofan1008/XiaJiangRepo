function denoiseIntan(filepath,dName,T,par)
dbstop if error
%% This function denoises an Intan datafile
%% Initial Setup
if strcmp(dName,'analogin')
    nChn = 1;
elseif strcmp(dName,'amplifier')
    nChn = 32;
end
FS = 30000;
%% Deal with the digital lines
cleanTrig;
%% Blank stimulation artefact
trig = loadTrig(0);
theseTrig = trig;
if ~isempty(theseTrig)
    info = dir([filepath '\' dName '.dat']);
    info = info.bytes/2;
    nL = (ceil(info / (nChn*FS*double(T)))+1);
    %lvFID = fopen([filepath '\' dName '.dat'],'r');
    vFID = fopen([filepath '\' dName '.dat'],'r');
    vdnFID = fopen([filepath '\' dName '_dn.dat'],'W');
    %lbFID = fopen([filepath '\' dName '_lfpraw.dat'],'w');
    N = 1;
    CHK = split(filepath,'_');
    blank = 1;
    for i = 1:size(CHK,1)
        if strcmp(CHK{i},'CSD') || strcmp(dName,'analogin') || isempty(dir('*exp_datafile*.mat'))
            blank = 0;
        end
    end
    dispstat('','init');
    if (blank)
        dispstat(sprintf('Blanking artefact . . .'),'keepthis','n');
    else
        dispstat(sprintf('Skipping blanking . . .'),'keepthis','n');
    end    
    BREAK = 1;
    while (N) && (BREAK)                
        if strcmp(dName,'analogin')
            v = fread(vFID,[nChn, (FS * T)],'uint16');
            v = (v - 32768) .* 0.0003125;
        elseif strcmp(dName,'amplifier')
            v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195;
        end
        if ~size(v,2)
            BREAK = 0;
        end
        dispstat(sprintf('Progress %03.2f%%',100*((N-1)/nL)),'timestamp');
        % Append saved data to the start of the next loop
        if (N ~= 1)
            v = [hv v]; %#ok<*AGROW>
        end
        if (blank)
            if (par)
                parfor iChn = 1:32                   
                    v(iChn,:) = simpleBlank(v(iChn,:),N,T,theseTrig,1);
                end
            else
                d = Depth;
                for iChn = 1:32
                    if d(iChn) == 13
                        pause(0.01);
                    end
                    v(iChn,:) = simpleBlank(v(iChn,:),N,T,theseTrig,1);
                end
            end
        end
        % Save a little data at the start of each loop
        hv = v(:,end-FS+1:end);
        % Don't save the last second, we'll deal with that next time
        if (BREAK)
            v = v(:,1:end-FS);
        end
        fwrite(vdnFID,v./0.195,'int16');        
        dispstat(sprintf('Progress %03.2f%%',100*((N)/nL)),'timestamp');
        N = N + 1;
    end
end
%% Clean up
fclose(vFID);
fclose(vdnFID);
end