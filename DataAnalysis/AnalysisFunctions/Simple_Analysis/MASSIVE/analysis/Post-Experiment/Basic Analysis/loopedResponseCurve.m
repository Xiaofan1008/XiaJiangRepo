%% Function takes in a filepath and outputs stimulus response curves
function loopedResponseCurve(file)
if nargin
    cd(file);
end
%% Initialisation
BIN = [-100, 200]; % Extraction window for reading in data, in msec.
spBIN = [4, 9]; % Spike extraction window relative to the pulse, in msec.
bpBIN = [16, 316];
dotSize = 50;
ARTIFACT_CUTOFF = 5e2; % 500 uV cutoff for discriminating between a spike and artifact
threshsign = 1; % 1 for negative thresholds, 2 for positive thresholds
TrialParams = []; AMP = []; DUR = []; sp = []; spExist = []; StimParams = [];% Basic initialisation
setupIntan; % Basic Intan data setup
TrialParams = cell2mat(TrialParams(2:end,:));
if ~isempty(spExist)
    disp('Loading spikes . . .');
    loadSpikes;
end
% Set up the graphs
setupGraphsBasic;
% Load the collapse matrix
if exist('collapse.mat','file')
    load('collapse.mat','collapse');
else
    collapse = zeros(nChn,2);
    for c = 1:nChn
        collapse(c,1) = c;
    end
end
%% Loop logic
nTr = length(time_stamps)/length(TrialParams(TrialParams(:,2) == 1));
nT = zeros(n_REP,nTr);
for i = 1:nTr
    nT(:,i) = TrialParams(TrialParams(:,2) == i,1);
end
% Stim Channel
stimChn = StimParams{2,1};
stimChn = str2double(stimChn(end-1:end))+1;
stimChn = find(depth == stimChn,1);
preMEAN = cell(1,nChn);
preSTE = preMEAN;
postMEAN = preSTE;
postSTE = postMEAN;
KEEP = cell(1,nChn);
figure; hold on;
waveforms = nan(20000,49);
count = 1;
for LOOP = 1:nChn
    thisCHN = LOOP;
    % Loading, filtering, and blanking.
    preStim = zeros(n_REP,nTr);
    postStim = zeros(n_REP,nTr);
    chan = depth(thisCHN);
    if ~isempty(spExist)
        theseSp = sp{chan};
        theseSp = denoiseSpikes(theseSp);
    end
    dispstat('','init');
    dispstat(sprintf('Processing data . . .'),'keepthis','n');
    for r = 1:n_REP
        dispstat(sprintf('Progress %03.2f%%',(100*(r/n_REP))),'timestamp');
        for t = 1:nTr
            %% Calculate the trial in question
            thisTrig = time_stamps(nT(r,t));
            if isempty(spExist)
                OFFSET = cast(nChn*2*(FS/1e3),'int64')*cast((thisTrig+BIN(1)),'int64');
                fseek(v_fid,OFFSET,'bof'); % Find the relevant timestamp
                v = fread(v_fid,[nChn,(FS/1e3)*diff(BIN)],'int16') * 0.195;
            end
            %% Extract spikes from v. OUTPUT will contain n_REP cell arrays
            if isempty(spExist)
                % Blank artefact
                if (t ~= 1)
                    b_v = BlankArte(v(chan,:),BIN,1);
                else
                    b_v = v(chan,:);
                end
                % Filter data
                tmp = conv(b_v,Mufilt);
                filt_v = tmp(1,MuNf/2:nData+MuNf/2-1);
                % Extract spikes
                [WAVES,Sp] = spikeextract(filt_v,thresh(threshsign,LOOP),FS);
                % Sort spikes to remove artefact
                if sum(Sp) ~= 0
                    for l = 1:length(Sp)
                        % Chuck away any spikes that exceed the threshold
                        if min(WAVES(l,:)) < -ARTIFACT_CUTOFF || max(WAVES(l,:)) > ARTIFACT_CUTOFF
                            Sp(l) = NaN;
                        end
                    end                    
                    Sp(isnan(Sp)) = [];
                else
                    Sp(Sp == 0) = [];
                end
                if size(Sp,1) == 1 && isempty(Sp)
                    Sp = Sp';
                end
                Sp = Sp(:,1);
                preStim(r,t) = sum(Sp >= (abs(BIN(1))-spBIN(2)) & Sp <= (abs(BIN(1))-spBIN(1)));
                postStim(r,t) = sum(Sp >= (abs(BIN(1))+spBIN(1)) & Sp <= (abs(BIN(1))+spBIN(2)));
            else
                preSp = theseSp(theseSp(:,1) >= (thisTrig - bpBIN(2)) & theseSp(:,1) <= (thisTrig - bpBIN(1)),:);
                postSp = theseSp(theseSp(:,1) >= (thisTrig + spBIN(1)) & theseSp(:,1) <= (thisTrig + spBIN(2)),:);
                preStim(r,t) = length(preSp(:,1));
                postStim(r,t) = length(postSp(:,1));
                if postStim(r,t) > 0
                    plot(postSp(:,2:end)');
                    waveforms(count:count+size(postSp(:,1),1)-1,:) = postSp(:,2:end);
                    count = count + size(postSp(:,1),1);
                end
            end
        end
    end
    KEEP{LOOP} = preStim;
    preMEAN{LOOP} = mean(preStim,1) ./ (diff(bpBIN)/1e3);
    preSTD = std(preStim,1) ./ (diff(bpBIN)/1e3);
    preSTE{LOOP} = preSTD ./ sqrt(length(preStim(:,1)));
    postMEAN{LOOP} = mean(postStim,1) ./ (diff(spBIN)/1e3);
    postSTD = std(postStim,1) ./ (diff(spBIN)/1e3);
    postSTE{LOOP} = postSTD ./ sqrt(length(postStim(:,1)));
end
plot(nanmean(waveforms),'Color','k','LineWidth',5);
waveforms(isnan(waveforms(:,1)),:) = [];
size(waveforms,1)
figure; hold on;
if length(AMP) > 1
    TICKLABELS = cell(1,length(AMP));
    for i = 1:length(AMP)
        TICKLABELS{i} = AMP(i);
    end
elseif length(DUR) > 1
    TICKLABELS = cell(1,length(DUR));
    for i = 1:length(DUR)
        TICKLABELS{i} = DUR(i);
    end
end
% Update the collapse matrix
collapse(:,2) = 0;
for c = 1:nChn
    chkChn;
end
save('collapse.mat','collapse');
TICKLABELS{1} = 'No Stim';
MAX = zeros(1,nChn);
for c = 1:nChn
    MAX(c) = max([postMEAN{c},preMEAN{c}]);
end
for c = 1:nChn
    INDEX = find(A' == c);
    subplot(ROW,COL,INDEX); hold on;
    set(gca,'FontSize',12)
    scatter([-1,1:nTr-1],preMEAN{c},dotSize,[0,0,0]);
    scatter([-1,1:nTr-1],postMEAN{c},dotSize,[0,0,0],'filled');    
    testStat = reshape(KEEP{c},[size(KEEP{c},1)*size(KEEP{c},2), 1]);
    [~,~,ci,~] = ttest(testStat,mean(testStat));
    threshold = ci(2) ./ (diff(bpBIN)/1e3);
    line([-2 nTr],[threshold threshold],'Color','r','LineStyle','--','LineWidth',2);
    threshold = ci(1) ./ (diff(bpBIN)/1e3);
    line([-2 nTr],[threshold threshold],'Color','r','LineStyle','--','LineWidth',2);
    errorbar([-1,1:nTr-1],preMEAN{c},preSTE{c},'Color',[0,0,0],'LineStyle','none','LineWidth',2,'HandleVisibility','off');
    errorbar([-1,1:nTr-1],postMEAN{c},postSTE{c},'Color',[0,0,0],'LineStyle','none','LineWidth',2,'HandleVisibility','off');
    xlim([-2 nTr-1]);
    ylim([0 MAX(c)*1.1]);
    yticks(0:25:MAX(c)*1.1);
    box on;
    grid on;
    if find(XLABEL == INDEX,1)
        xticks([-1,1:nTr-1]);
        xticklabels(TICKLABELS);
        xlabel('Current (uA)');
    else
        xticks('');
    end
    if find(YLABEL == INDEX,1)
        ylabel('Response (S/s)');
    else
        ylabel('');
    end
    if collapse(c,2) == 1
        set(gca,'YColor','b');
        set(gca,'XColor','b');
        title(['\color{blue}Electrode ' num2str(c)]);
    else
        title(['Electrode ' num2str(c)]);
    end
    if c == stimChn
        set(gca,'YColor','r');
        set(gca,'XColor','r');
        if collapse(c,2) == 1
            title(['\color{blue}Electrode ' num2str(c) ' | Stim']);
        else
            title(['\color{red}Electrode ' num2str(c) ' | Stim']);
        end
    end
end
sgtitle('300 \mum penetration 2 | RAT0011')
end