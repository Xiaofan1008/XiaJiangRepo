function denoiseBlackrock(filepath,dName,nChn,FS,T,par)
dbstop if error
%% This function denoises a Blackrock datafile
%% Deal with the digital lines
if isempty(dir('*trigBR.dat'))
    cleanTrigBR;
end
%% Blank stimulation artefact
trig = loadTrigBR(0);
theseTrig = trig;
if ~isempty(theseTrig)
    info = openNSx('noread',dName);
    info = info.MetaTags.DataPoints; % Number of samples in the datafile
    if length(info) > 1
        neuInfo = 0;
        for i = 1:length(info)
            neuInfo = neuInfo + info(i);
        end
        info = neuInfo;
    end
    nL = ceil(info / (FS*double(T)));
    %lvFID = fopen([filepath '\' dName '.dat'],'r');
    if (nPauses == 1)
        vFID = [filepath '\' dName];
    else
        NSx = openNSx(name,['e:1:' num2str(nChn)],'p:double');
        v = [];
        for n = 1:nPauses
            v = [v, NSx.Data{n}];
        end
        total_v = v;
    end
    vdnFID = fopen([filepath '\' dName '_dnBR.dat'],'W');
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
        if ((1+(N-1)*(FS*T))) > info
            BREAK = 0;
            continue;
        end
        if (nPauses == 1)
            NSx = openNSx(vFID,['e:1:' num2str(nChn)],['t:' num2str((1+(N-1)*(FS*T))) ':' num2str(((FS*T)+(N-1)*(FS*T)))],'p:double');
            v = NSx.Data;
        else
            v = total_v(:,(1+(N-1)*(FS*T)):(FS*T)+(N-1)*(FS*T));
        end
        dispstat(sprintf('Progress %03.2f%%',100*((N-1)/nL)),'timestamp');
        % Append saved data to the start of the next loop
        if (N ~= 1)
            v = [hv v]; %#ok<*AGROW>
        end
        if (blank)
            if (par)
                parfor iChn = 1:nChn
                    v(iChn,:) = simpleBlank(v(iChn,:),N,T,theseTrig,1);
                end
            else
                for iChn = 1:nChn
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
%% Clean up
fclose(vdnFID);
end
end