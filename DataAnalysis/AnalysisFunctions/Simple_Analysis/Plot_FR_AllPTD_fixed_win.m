%% ==================== USER INPUTS ====================
seq_folders = {
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_8ms_251125_113541'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_10ms_251125_114346'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_15ms_251125_115256'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_20ms_251125_120014'
    '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_25ms_251125_120717'
};

Electrode_Type = 1;
target_channel = 18;
plot_amps      = 10;

FS = 30000;
win_ms = 5;     % analysis window after each pulse

%% ==================== LOAD DATA ====================
Seq = struct;
for f = 1:numel(seq_folders)
    folder = seq_folders{f};
    cd(folder);

    % ---- Load spikes ----
    if ~isempty(dir('*sp_xia_FirstPulse.mat'))
        load(dir('*sp_xia_FirstPulse.mat').name,'sp_seq');
        Seq(f).sp = sp_seq;
    else
        load(dir('*sp_xia.mat').name,'sp_clipped');
        Seq(f).sp = sp_clipped;
    end

    % ---- Trigger ----
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
    Seq(f).trig = loadTrig(0);

    % ---- StimParams ----
    S = load(dir('*_exp_datafile_*.mat').name, ...
             'StimParams','simultaneous_stim','E_MAP','n_Trials');

    Stim = S.StimParams;
    simN = S.simultaneous_stim;

    Seq(f).nTrials   = S.n_Trials;
    Seq(f).simN      = simN;
    Seq(f).E_MAP     = S.E_MAP;

    % Amp (col 16)
    amps_all = cell2mat(Stim(2:end,16));
    Seq(f).trialAmps = amps_all(1:simN:end);

    % PTD (col 6, µs → ms)
    PTD_all = cell2mat(Stim(2:end,6));      
    Seq(f).PTD = PTD_all(2:simN:end) / 1000;

    % ----- Ordered stimulation sets -----
    stimNames = Stim(2:end,1);
    [~, idx]  = ismember(stimNames, Seq(f).E_MAP(2:end));

    stimChPerTrial = cell(Seq(f).nTrials,1);
    for t = 1:Seq(f).nTrials
        rr = (t-1)*simN + (1:simN);
        v = idx(rr); 
        v = v(v>0);
        stimChPerTrial{t} = v(:).';  % keep order
    end
    Seq(f).stimChPerTrial = stimChPerTrial;

    % ----- Unique stimulation sets -----
    comb = zeros(Seq(f).nTrials, simN);
    for t=1:Seq(f).nTrials
        v = stimChPerTrial{t};
        comb(t,1:numel(v)) = v;
    end
    [Seq(f).uniqueComb,~,Seq(f).combClass] = unique(comb,'rows','stable');
    Seq(f).nSets = size(Seq(f).uniqueComb,1);
end

%% ==================== DEPTH MAP ====================
d = Depth_s(Electrode_Type);
recCh = d(target_channel);

%% ==================== FIRING RATE (GLOBAL UNIQUE SETS) ====================
fprintf('\n===== FIRING RATE RESULTS (Corrected stim-set grouping) =====\n');

nF = numel(seq_folders);

% --- Step 1: Collect all unique stimulation sets across all PTDs ---
allSets = [];
for f = 1:nF
    tmp = Seq(f).uniqueComb;
    allSets = [allSets; tmp];
end
% remove zero columns, remove duplicates
allSets(allSets==0) = nan;
[globalSets,~,~] = unique(allSets,'rows','stable'); 
nGlobal = size(globalSets,1);

% --- Allocate ---
FR_total  = nan(nF, nGlobal);
SEM_total = nan(nF, nGlobal);
PTDvals   = zeros(nF,1);

for f = 1:nF
    PTD = Seq(f).PTD(1);
    PTDvals(f) = PTD;
    sp_ch = Seq(f).sp{recCh};
    trig  = Seq(f).trig;

    % for each global stimulation pattern
    for g = 1:nGlobal
        thisSet = globalSets(g, :);
        thisSet = thisSet(~isnan(thisSet));     % remove NaNs

        % find matching set index within THIS folder
        setIdx = [];
        for s = 1:Seq(f).nSets
            localSet = Seq(f).uniqueComb(s,:);
            localSet = localSet(localSet>0);
            if isequal(thisSet, localSet)
                setIdx = s;
                break;
            end
        end

        if isempty(setIdx)
            % this folder does not contain this stim set
            continue;
        end

        % extract trials for this stimulation set AND amplitude
        trial_ids = find(Seq(f).combClass == setIdx & Seq(f).trialAmps == plot_amps);
        if isempty(trial_ids)
            continue;
        end

        nT = numel(trial_ids);
        fr_all = zeros(nT,1);

        for ii = 1:nT
            tr = trial_ids(ii);
            t0 = trig(tr)/FS*1000;

            mask1 = sp_ch(:,1) >= t0      & sp_ch(:,1) < t0 + win_ms;
            mask2 = sp_ch(:,1) >= t0+PTD  & sp_ch(:,1) < t0 + PTD + win_ms;

            n_spk = sum(mask1) + sum(mask2);
            fr_all(ii) = n_spk / ((2*win_ms)/1000);
        end

        FR_total(f,g)  = mean(fr_all);
        SEM_total(f,g) = std(fr_all) / sqrt(nT);

        fprintf("PTD=%2g ms | Set [%s] | FR %.2f ± %.2f\n", ...
            PTD, strjoin(string(thisSet),'→'), FR_total(f,g), SEM_total(f,g));
    end
end


%% ==================== PLOT (CORRECTED LEGEND + LINES) ====================
figure('Color','w','Position',[200 200 900 550]); hold on;

colors = lines(nGlobal);
hSet = gobjects(nGlobal,1);

for g = 1:nGlobal
    y  = FR_total(:,g);
    se = SEM_total(:,g);

    if all(isnan(y)), continue; end

    hLine = plot(PTDvals, y, '-o', ...
        'Color', colors(g,:), 'LineWidth', 2.5, ...
        'MarkerSize', 7, 'MarkerFaceColor', colors(g,:), ...
        'MarkerEdgeColor', 'k');

    hSet(g) = hLine;

    errorbar(PTDvals, y, se, ...
        'LineStyle','none', ...
        'Color', colors(g,:), 'LineWidth',1.3);
end

xlabel('Inter Stimulus Interval (ISI) /ms');
ylabel('Firing Rate (spikes/s)');
title(sprintf('Firing Rate — Ch %d (win = %d ms)', target_channel, win_ms));
box off;

% --- Correct legend ---
legendStrings = cell(nGlobal,1);
for g = 1:nGlobal
    setG = globalSets(g,~isnan(globalSets(g,:)));
    legendStrings{g} = sprintf('Set [%s]', strjoin(string(setG),'→'));
end

legend(hSet, legendStrings, 'Box','off', 'Location','northwest');