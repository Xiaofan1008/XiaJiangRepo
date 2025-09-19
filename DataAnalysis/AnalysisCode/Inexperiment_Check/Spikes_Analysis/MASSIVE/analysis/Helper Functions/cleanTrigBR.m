%% Initial Setup
here = pwd;
fileinfo = dir('*.nev');
name = fileinfo.name;
TrialParams = loadTrialParams; TrialParams = cell2mat(TrialParams(:,2))';
%% Deal with the digital lines
disp('Cleaning the digital lines');
header = openNEV([here '\' name],'nosave');
FS = header.MetaTags.SampleRes;
time_stamps = header.Data.SerialDigitalIO.TimeStamp;
try
    ee = dir('*exp_datafile*.mat');
    ee = ee.name;
    nStamps = load(ee,'n_Trials');
    nStamps = nStamps.n_Trials;
catch
    fprintf('Please include an Intan experimental datafile in the working directory\n');
    keyboard;
end
time_stamps = time_stamps(1:2:end);
trig = time_stamps;
fprintf('There are %.0f digital lines here.\n',length(trig));
% Correct for movement
fileinfo = dir('*.ns6');
name = fileinfo.name;
loadCHN; AMP = loadAMP;
zT = zeros(1,length(CHN));
for iC = 1:length(CHN)
    zT(1,iC) = 1+(iC-1)*length(AMP);
end
warning off signal:findpeaks:largeMinPeakHeight
try
    if (par)
        parfor t = 1:length(trig)
            sT = trig(t)-FS/10;
            eT = trig(t)+FS/10;
            NSx = openNSx(name,'e:1',['t:' num2str(sT) ':' num2str(eT)],'p:double');
            v = abs(NSx.Data);
            v = v - mean(v);
            if any(zT == TrialParams(t))
                L = 3000;
            else
                [~,L] = findpeaks(thisV,'MinPeakHeight',5000);
                if isempty(L)
                    [~,L] = findpeaks(thisV,'MinPeakHeight',4000);
                    if isempty(L)
                        [~,L] = findpeaks(thisV,'MinPeakHeight',3000);
                        if isempty(L)
                            [~,L] = findpeaks(thisV,'MinPeakHeight',2000);
                            if isempty(L)
                                [~,L] = findpeaks(thisV,'MinPeakHeight',1000);
                            end
                        end
                    end
                end
            end
            %L(L < 12000 | L > 18000) = [];
            if ~isempty(L)
                L = L(1);
                while thisV(L) > 1e3
                    L = L - 1;
                    if (L == 0)
                        L = 3000;
                        break;
                    end
                end
            else
                L = 3000;
            end
            trig(t) = trig(t) - (3000 - L);
        end
    else
        for t = 1:length(trig)
            sT = trig(t)-FS/10;
            eT = trig(t)+FS/10;
            NSx = openNSx(name,'e:1',['t:' num2str(sT) ':' num2str(eT)],'p:double');
            v = abs(NSx.Data);
            v = v - mean(v);
            if any(zT == TrialParams(t))
                L = 3000;
            else
                [~,L] = findpeaks(v,'MinPeakHeight',5000);
                if isempty(L)
                    [~,L] = findpeaks(v,'MinPeakHeight',4000);
                    if isempty(L)
                        [~,L] = findpeaks(v,'MinPeakHeight',3000);
                        if isempty(L)
                            [~,L] = findpeaks(v,'MinPeakHeight',2000);
                            if isempty(L)
                                [~,L] = findpeaks(v,'MinPeakHeight',1000);
                            end
                        end
                    end
                end
            end
            %L(L < 12000 | L > 18000) = [];
            if ~isempty(L)
                L = L(1);
                while thisV(L) > 1e3
                    L = L - 1;
                    if (L == 0)
                        L = 3000;
                        break;
                    end
                end
            else
                L = 3000;
            end
            trig(t) = trig(t) - (3000 - L);
            fprintf('Done with line %.0f of %.0f\n',t,length(trig));
        end
    end
catch
    % Datafile might contain pauses
    NSx = openNSx(name,'e:1','p:double');
    nPauses = length(NSx.Data);
    if nPauses > 1
        v = [];
        for nP = 1:nPauses
            v = [v, NSx.Data{nP}]; %#ok<AGROW>
        end
        v = abs(v);
        if (par)
            parfor t = 1:length(trig)
                sT = trig(t)-FS/10;
                eT = trig(t)+FS/10;
                thisV = v(sT:eT); %#ok<PFBNS>
                thisV = thisV - mean(thisV);
                if any(zT == TrialParams(t))
                    L = 3000;
                else
                    [~,L] = findpeaks(thisV,'MinPeakHeight',5000);
                    if isempty(L)
                        [~,L] = findpeaks(thisV,'MinPeakHeight',4000);
                        if isempty(L)
                            [~,L] = findpeaks(thisV,'MinPeakHeight',3000);
                            if isempty(L)
                                [~,L] = findpeaks(thisV,'MinPeakHeight',2000);
                                if isempty(L)
                                    [~,L] = findpeaks(thisV,'MinPeakHeight',1000);
                                end
                            end
                        end
                    end
                end
                %L(L < 12000 | L > 18000) = [];
                if ~isempty(L)
                    L = L(1);
                    while thisV(L) > 1e3
                        L = L - 1;
                        if (L == 0)
                            L = 3000;
                            break;
                        end
                    end
                else
                    L = 3000;
                end
                trig(t) = trig(t) - (3000 - L);
            end
        else
            for t = 1:length(trig)
                sT = trig(t)-FS/10;
                eT = trig(t)+FS/10;
                thisV = v(sT:eT);
                thisV = thisV - mean(thisV);
                if any(zT == TrialParams(t))
                    L = 3000;
                else
                    [~,L] = findpeaks(thisV,'MinPeakHeight',5000);
                    if isempty(L)
                        [~,L] = findpeaks(thisV,'MinPeakHeight',4000);
                        if isempty(L)
                            [~,L] = findpeaks(thisV,'MinPeakHeight',3000);
                            if isempty(L)
                                [~,L] = findpeaks(thisV,'MinPeakHeight',2000);
                                if isempty(L)
                                    [~,L] = findpeaks(thisV,'MinPeakHeight',1000);
                                end
                            end
                        end
                    end
                end
                %L(L < 12000 | L > 18000) = [];
                if ~isempty(L)
                    L = L(1);
                    while thisV(L) > 1e3
                        L = L - 1;
                        if (L == 0)
                            L = 3000;
                            break;
                        end
                    end
                else
                    L = 3000;
                end
                trig(t) = trig(t) - (3000 - L);
                fprintf('Done with line %.0f of %.0f\n',t,length(trig));
            end
        end
    end
end
trig_fid = fopen([name '.trigBR.dat'],'W');
fwrite(trig_fid,trig,'double');
fclose(trig_fid);