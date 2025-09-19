function online_processIntan(par)
%% This function is written to be as fast as possible. It excludes some features of the offline version to do this, and is not intended as a replacement to the offline function
% Note that this function requires a substantial amount of memory.
dbstop if error
if ~nargin
    par = 0;
end
%% Step One: Initialise
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
    T = cast((mem.MaxPossibleArrayBytes ./ (2 * 32 * 30000)),'int64');
    if T > t_len
        T = t_len + 1;
    end
catch
    T = 124; % If we aren't in Windows, default to a few minutes
end
T = cast(T,'int32');
if strcmp(dName,'analogin')
    nChn = 1;
elseif strcmp(dName,'amplifier')
    nChn = 32;
end
FS = 30000;
threshfac = -7;
sp = cell(1, nChn);
thresh = cell(1, nChn);
muN = 0; time = 0;
NSp = zeros(1,nChn);
ARTMAX = 1.5e3;
%name = pwd;
%name = strsplit(name,'\');
%name = name{end};
%name = name(1:end-14);
artchk = zeros(1,nChn);
munoise = cell(1,nChn);
nRun = 0;
CHANNEL = 12; % The Channel(s) of Interest
%% Generate filters
[Mufilt,Lfpfilt] = generate_Filters;
MuNf = length(Mufilt);
%LfpNf = length(Lfpfilt);
%% Step Two: Do Everything All At Once
disp(['Analysing the ' dName ' file']);
if ~isempty(dir('*_exp_datafile_*'))
    % Blank the artefact
    %% Deal with the digital lines
    cleanTrig;
    %% Blank stimulation artefact
    trig = loadTrig;
    theseTrig = trig;
    info = dir([ifile '\' dName '.dat']);
    info = info.bytes/2;
    nL = (ceil(info / (nChn*FS*double(T)))+1)-1;
    nSam = info/(nChn);
    ntimes = ceil(info / nChn / FS / double(T));
    vFID = fopen([ifile '\' dName '.dat'],'r');
    N = 1;
    CHK = split(ifile,'_');
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
    %lfp_fid = fopen([name '.lfp.dat'],'W');
    %lfp3 = zeros(nChn,1);
    while (N) && (BREAK)
        dispstat(sprintf('Progress %03.2f%%',100*((N-1)/nL)),'timestamp');
        if strcmp(dName,'analogin')
            v = fread(vFID,[nChn, (FS * T)],'uint16');
            v = (v - 32768) .* 0.0003125;
        elseif strcmp(dName,'amplifier')
            v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195;
        end
        if ~size(v,2)
            BREAK = 0;
        end
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
                for iChn = 1:32
                    v(iChn,:) = simpleBlank(v(iChn,:),N,T,theseTrig,1);
                end
            end
        end
        % Save a little data at the start of each loop
        hv = v(:,end-FS+1:end);
        if (BREAK)
            v = v(:,1:end-FS);
        end
        dispstat(sprintf('Progress %03.2f%%',100*((N)/nL)),'timestamp');
        %% Move to extraction
        if sum(artchk) < nChn
            dispstat('','init');
            dispstat(sprintf('Processing thresholds . . .'),'keepthis','n');
            for iChn = 1:nChn
                sp{iChn} = zeros(ntimes * T * 20, FS * 1.6 / 1e3 + 1 + 1);
            end
            tRun = ceil(nSam / FS / 10) + 1; % Number of 10 second chunks of time in the data
            while sum(artchk) < nChn
                nRun = nRun + 1;
                dispstat(sprintf('Progress %03.2f%%',(100*(nRun/tRun))),'timestamp');
                t_v = v(:,1+((nRun-1)*FS*10):1+((nRun)*FS*10));
                if ~isempty(t_v)
                    % Checks for artifact
                    for iChn = 1:nChn
                        if isempty(thresh{iChn})
                            munoise{iChn} = [];
                            mu = conv(t_v(iChn,:),Mufilt);
                            if max(abs(mu(MuNf+1:end-MuNf))) < ARTMAX
                                artchk(iChn) = 1;
                                sd = median(abs(mu(MuNf+1:end-MuNf)))./0.6745;
                                thresh{iChn} = threshfac*sd;
                            else
                                munoise{iChn} = [munoise{iChn} mu(MuNf+1:end-MuNf)];
                            end
                        end
                    end
                else
                    for iChn = 1:nChn
                        if isempty(thresh{iChn}) % If threshold is still 0 - noisy channel
                            artchk(iChn) = 1;
                            sd = median(abs(munoise{iChn}))./0.6745;
                            thresh{iChn} = threshfac*sd;
                        end
                    end
                end
            end
            clear t_v
            disp(['Total recording time: ' num2str(nSam / FS) ' seconds.']);
            disp(['Time analysed per loop: ' num2str(T) ' seconds.']);
        end
        %% Loop through the data
        dispstat('','init');
        dispstat(sprintf('Processing SPMULFP . . .'),'keepthis','n');
        muN = muN + 1;
        dispstat(sprintf('Progress %03.2f%%',100*((muN-1)/ntimes)),'timestamp');
        if (size(v,2))
            Ndata = size(v,2);
            %LfpOut = cell(1,nChn);
            if (par)
                parfor iChn = 1:nChn
                    tmp = conv(fliplr(v(iChn,:)),Mufilt);
                    mu = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
                    %tmp = conv(v(iChn,:),Lfpfilt);
                    %lfp = tmp(1,LfpNf/2:FS/1e3:Ndata+LfpNf/2-1);
                    [Sp_tmp, Spktimes_tmp] = spikeextract(mu, thresh{iChn}, FS);
                    NSp_tmp = length(Spktimes_tmp);
                    if NSp_tmp > 1
                        SpMat = [Spktimes_tmp+double(time) Sp_tmp];
                        NSp_tmp = size(SpMat,1);
                        sp{iChn}(NSp(iChn)+1:NSp(iChn)+NSp_tmp,:) = SpMat;
                        NSp(iChn) = NSp(iChn) + NSp_tmp;
                    end
                    %LfpOut{iChn} = lfp;
                end
            else
                for iChn = 1:nChn
                    tmp = conv(fliplr(v(iChn,:)),Mufilt);
                    mu = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
                    %tmp = conv(v(iChn,:),Lfpfilt);
                    %lfp = tmp(1,LfpNf/2:FS/1e3:Ndata+LfpNf/2-1);
                    [Sp_tmp, Spktimes_tmp] = spikeextract(mu, thresh{iChn}, FS);
                    NSp_tmp = length(Spktimes_tmp);
                    if NSp_tmp > 1
                        SpMat = [Spktimes_tmp+double(time) Sp_tmp];
                        NSp_tmp = size(SpMat,1);
                        sp{iChn}(NSp(iChn)+1:NSp(iChn)+NSp_tmp,:) = SpMat;
                        NSp(iChn) = NSp(iChn) + NSp_tmp;
                    end
                    %LfpOut{iChn} = lfp;
                end
            end
            clear mu tmp
            %lfp2 = zeros(nChn,size(LfpOut{1},2));
            %for iChn = 1:nChn
            %    lfp2(iChn,:) = LfpOut{iChn};
            %end
            %lfp3 = [lfp3, lfp2];
            time = time + T*1e3;
        end
        dispstat(sprintf('Progress %03.2f%%',100*((muN)/ntimes)),'timestamp');
        N = N + 1;
    end
    %lfp3(:,1) = [];
    %% Save the datafiles
    disp('Saving spikes');
    for iChn = 1:nChn
        sp{iChn} = sp{iChn}(1:NSp(iChn),:);
    end
    %save([name '.sp.mat'],'sp','thresh','threshfac','-v7.3');
    %fwrite(lfp_fid,lfp3,'float')
    %fclose(lfp_fid);
    %% Go Into Analysis
    n_REP = []; AMP = [];% Basic initialisation
    %% Variables
    stimSt = 3.2;
    baselineBIN = [-175, -21];
    stimulationBIN = [stimSt, 50];
    suppressionBIN = [51,100];
    dead_channels = [2,9];
    mChn = 1; stimChn = 0; loadStimChn;
    if size(CHANNEL,2) > 1
        mChn = 0;
        setupGraphsBasic;
        collapse = zeros(size(CHANNEL,2),4);
        collapse(:,1) = CHANNEL;
    end
    responseFig = figure('units','normalized','outerposition',[0 0 1 1]);
    if ~(mChn)
        s_axes = tight_subplot(6,6,[.025 .025],[.075 .080],[.05 .01]);
        delete(s_axes([1,6,31,36]));
    end
    name = split(pwd,'\');
    n1 = name{end-1};
    n2 = split(name{end},'_');
    n3 = n2{3};
    n2 = n2{2};
    for thisChn = CHANNEL
        %% Loading
        TrialParams = loadTrialParams; trig = loadTrig; loadNTrials; loadAMP; d = Depth(3);
        %% Initialisation of data structures
        nTr = size(TrialParams,1)/n_REP;
        nT = zeros(n_REP,nTr);
        TrialParams = cell2mat(TrialParams);
        COLOR = cell(1:nTr);COLOR{1} = [0.7 0.7 0.7];COLOR{2} = [0.6 0.6 0.6];COLOR{3} = [0.5 0.5 0.5];COLOR{4} = [0.4 0.4 0.4];
        COLOR{5} = [0.3 0.3 0.3];COLOR{6} = [0.2 0.2 0.2];COLOR{7} = [0.1 0.1 0.1];COLOR{8} = [0 0 0];
        for i = 1:nTr
            nT(:,i) = TrialParams(TrialParams(:,2) == i,1);
        end
        preStim = zeros(n_REP,nTr); postStim = zeros(n_REP,nTr); suppStim = zeros(n_REP,nTr);
        chn = d(thisChn);
        theseSp = cast(sp{chn},'int32');
        theseSp = denoiseSpikes(theseSp);
        SPIKES = cell(1,nTr); spWaves = cell(1,nTr); psthBIN = [-200 200]; rate = zeros(nTr,diff(psthBIN)+1);
        TICKLABELS = cell(1,length(AMP)); LEGEND = cell(1,length(AMP));
        for i = 1:length(AMP)
            TICKLABELS{i} = AMP(i);
            LEGEND{i} = [num2str(AMP(i)) ' uA'];
        end
        TICKLABELS{1} = 'No Stim';
        LEGEND{1} = 'No Stim';
        %% PSTH Methodology for window adjustment
        if (mChn)
            rateFig = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
            % Line at time point 0
            line([0-psthBIN(1) 10*nTr-psthBIN(1)],[0 100*nTr],'Color','b','HandleVisibility','off');
            if (mChn)
                title([n1 ' | ' n2 ' | ' n3 ' | Chn ' num2str(thisChn)]);
            else
                title(['Chn: ' num2str(thisChn)]);
            end
            ylabel('Firing Rate (Sp/s) + Offset 100 S/s/Condition');
            xlabel('Time (ms)');
            xticks(0:10:400);
            xticklabels(-200:10:200);
            rate_x = 1:1:diff(psthBIN)+1+(nTr-1)*10;
            box on; grid on;
            spikeFig = figure('units','normalized','outerposition',[0 0 1 1]);
        end
        for r = 1:nTr
            meanWave = zeros(1,49); count = 0;
            for t = 1:size(TrialParams,1)/nTr
                thisSp = theseSp(theseSp(:,1) >= (trig(nT(t,r))+psthBIN(1)) & theseSp(:,1) <= (trig(nT(t,r))+psthBIN(2)),:);
                SPIKES{t} = thisSp(:,1) - (trig(nT(t,r))+psthBIN(1));
            end
            rate(r,:) = psth(SPIKES(:),psthBIN,1,500);
            baseRateMean = mean(rate(r,21:abs(psthBIN(1))-21));
            stStart = find(rate(r,abs(psthBIN(1)):end)>=baseRateMean*1.2,1) - 1;
            stStop = find(rate(r,abs(psthBIN(1))+stStart:end)<=baseRateMean*1.2,1) + stStart - 1;
            stimulationBIN(r,1:2) = [stimSt, min([max([stStop,7]),10])];
            suppressionBIN(r,1:2) = [stimulationBIN(r,2)+1, 50];
            if r == 1
                stimulationBIN(r,1:2) = [stimSt,8];
                suppressionBIN(r,1:2) = [51,201];
            end
            if (mChn)
                for t = 1:size(TrialParams,1)/nTr
                    thisSp = theseSp(theseSp(:,1) >= (trig(nT(t,r))+stimulationBIN(r,1)) & theseSp(:,1) <= (trig(nT(t,r))+stimulationBIN(r,2)),:);
                    spWaves{t} = thisSp(:,2:end);
                end
                figure(rateFig);
                tmp = rate(r,196:204);
                rate(r,197:203) = NaN;
                plot(rate_x(20+(r-1)*10:diff(psthBIN)+1+(r-1)*10),rate(r,20:end) + (r-1)*100,...
                    'Color',COLOR{r});
                plot(rate_x(20+(r-1)*10+176:20+(r-1)*10+184),tmp + (r-1)*100,...
                    'Color','b');
                figure(spikeFig);
                subplot(1,nTr,r); hold on; grid on; box on;
                for t = 1:size(TrialParams,1)/nTr
                    plot(spWaves{t}');
                    meanWave = meanWave + sum(spWaves{t},1);
                    count = count + size(spWaves{t},1);
                end
                xticks(0:15:49);
                xticklabels(0:0.5:49/30);
                meanWave = meanWave ./ count;
                plot(meanWave,'Color','k','LineWidth',4);
                text(15,200,[num2str(count) ' spikes'],'FontSize',12);
                text(15,180,[num2str(n_REP) ' trials'],'FontSize',12);
                text(15,160,[num2str(diff(stimulationBIN(r,1:2))) ' ms window']);
                title([num2str(TICKLABELS{r}) ' uA']);
                ylabel('Voltage (uV)');
                xlabel('Time (ms)');
                ylim([-300 250]);
                xlim([0 49]);
                set(gca,'FontSize',14);
            end
            sgtitle([n1 ' | ' n2 ' | ' n3 ' | Chn ' num2str(thisChn)]);
        end
        if (mChn)
            figure(rateFig);
            axis([-50-psthBIN(1) 50-psthBIN(1)+(nTr*(10)) -10 max(max(rate))+(nTr*(100))]);
            yticks(0:50:max(max(rate))+(nTr*(100)));
            set(gca,'FontSize',18);
            legend(LEGEND,'Location','northwest');
        end
        %% Logic
        for r = 1:size(TrialParams,1)/nTr
            for t = 1:nTr
                thisTrig = trig(nT(r,t));
                preStim(r,t) = sum(theseSp(:,1) >= (thisTrig + baselineBIN(1)) & theseSp(:,1) <= (thisTrig + baselineBIN(2)));
                postStim(r,t) = sum(theseSp(:,1) >= (thisTrig + stimulationBIN(t,1)) & theseSp(:,1) <= (thisTrig + stimulationBIN(t,2)));
                suppStim(r,t) = sum(theseSp(:,1) >= (thisTrig + suppressionBIN(t,1)) & theseSp(:,1) <= (thisTrig + suppressionBIN(t,2)));
            end
        end
        %% Statistics
        preMean = mean(preStim,1) ./ (diff(baselineBIN)/1e3);
        postMean = mean(postStim,1) ./ (diff(stimulationBIN,1,2)/1e3)';
        suppMean = mean(suppStim,1) ./ (diff(suppressionBIN,1,2)/1e3)';
        preSte = (std(preStim,1) ./ (diff(baselineBIN)/1e3)) ./ sqrt(n_REP);
        postSte = (std(postStim,1) ./ (diff(stimulationBIN,1,2)/1e3)') ./ sqrt(n_REP);
        suppSte = (std(suppStim,1) ./ (diff(suppressionBIN,1,2)/1e3)') ./ sqrt(n_REP);
        [preP,~,stats] = anova1(preStim ./ (diff(baselineBIN)/1e3),[],'off'); %#ok<*ASGLU>
        [st,m] = multcompare(stats,'alpha',0.05,'Display','off');
        [postP,~] = anova1(postStim,[],'off');
        [suppP,~] = anova1(suppStim,[],'off');
        % Confidence intervals
        ci(1) = mean(m(:,1)) + mean(st(:,5));
        ci(2) = mean(m(:,1)) + mean(st(:,3));
        H = zeros(1,nTr);
        P = zeros(1,nTr);
        for n = 1:nTr
            [h,p] = ttest2((preStim(:,n) ./ (diff(baselineBIN)/1e3)),(postStim(:,n) ./ (diff(stimulationBIN(n,1:2)/1e3))));
            if isnan(h)
                h = 0;
            end
            if isnan(p)
                p = 1;
            end
            H(n) = h;
            P(n) = p;
        end
        %% Plotting
        figure(responseFig); hold on;
        if (mChn)
            dotSize = 100;
            ax = gca;
            hold(ax,'on');
        else
            INDEX = find(A' == thisChn);
            ax = s_axes(INDEX);
            hold(ax,'on');
            dotSize = 40;
        end
        if any(dead_channels == thisChn)
            xlim(ax,[0 1]);
            xticks(ax,'');
            ylim(ax,[0 1]);
            yticks(ax,'');
            line(ax,[0 1],[0 1],'Color','k');
            line(ax,[0 1],[1 0],'Color','k');
            beautifyPlot(ax);
            drawnow;
        else
            scatter(ax,[-1,1:nTr-1],preMean,dotSize,[0,0,0]);
            errorbar(ax,[-1,1:nTr-1],preMean,preSte,'Color',[0,0,0],'LineStyle','none','LineWidth',1,'HandleVisibility','off');
            scatter(ax,[-1,1:nTr-1],postMean,dotSize,[0,0,0],'filled','r');
            errorbar(ax,[-1,1:nTr-1],postMean,postSte,'Color',[0,0,0],'LineStyle','none','LineWidth',1,'HandleVisibility','off');
            scatter(ax,[-1,1:nTr-1],suppMean,dotSize,[0,0,0],'filled','b','d');
            errorbar(ax,[-1,1:nTr-1],suppMean,suppSte,'Color',[0,0,0],'LineStyle','none','LineWidth',1,'HandleVisibility','off');
            line(ax,[-2 nTr+2],[ci(1) ci(1)],'Color','r','LineStyle','--','LineWidth',1);
            line(ax,[-2 nTr+2],[ci(2) ci(2)],'Color','r','LineStyle','--','LineWidth',1);
            HX = [-1,1:nTr-1];
            for nT = 1:nTr
                if (H(nT))
                    %text(HX(nT)-0.04,postMean(nT)*1.2,['* ' num2str(P(nT))],'FontSize',20);
                    text(ax,HX(nT)-0.04,postMean(nT)+2*postSte(nT),'*','FontSize',20);
                end
            end
            if (mChn)
                xlabel(ax,'Current (uA)');
                ylabel(ax,'Response (Sp/s)');
            else
                if find(XLABEL == INDEX,1)
                    xlabel(ax,'Current (uA)');
                    xlim(ax,[-2 nTr]);
                    xticks(ax,[-1,1:nTr-1]);
                    xticklabels(ax,TICKLABELS);
                else
                    xlim(ax,[-2 nTr]);
                    xticks(ax,[-1,1:nTr-1]);
                    xticklabels(ax,'');
                end
                if find(YLABEL == INDEX,1)
                    ylabel(ax,'Response (Sp/s)');
                end
            end
            if max(postMean) == 0
                ylim(ax,[0 20]);
                yticks(ax,0:10:20);
            else
                ylim(ax,[0 max(postMean)+3*max(postSte)]);
                y = ceil(max(postMean)+3*max(postSte)) / 3;
                if (y < 100)
                    if (y < 10)
                        y = round(y);
                    else
                        y = round(y,-1);
                    end
                else
                    y = round(y,-2);
                end
                if y == 0
                    y = 1;
                end
                yticks(ax,0:y:max(postMean)+3*max(postSte));
            end
            if (mChn)
                legend({'Baseline','Stimulation-Evoked','Post-Stim Suppression'},'Location','northwest');
            else
                % Assess the channel for a response
                chkChn;
                if collapse(thisChn,2)
                    set(ax,'YColor','b');
                    set(ax,'XColor','b');
                end
                if thisChn == stimChn
                    set(ax,'XColor','r');
                end
            end
            beautifyPlot(ax);
            drawnow;
        end
    end
end
%% Step Three: Close Down
fclose('all');
if (par)
    delete(hPool);
end
end