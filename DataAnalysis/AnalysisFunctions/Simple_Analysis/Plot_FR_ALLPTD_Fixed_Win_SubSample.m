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
target_channel = 19;
plot_amps      = 10;

FS = 30000;
win_ms = 5;     % analysis window after each pulse

subsample_N = 25;     % (set [] to disable subsampling)

%% ==================== LOAD DATA ====================
Seq = struct;
for f = 1:numel(seq_folders)
    folder = seq_folders{f};
    cd(folder);

    if ~isempty(dir('*sp_xia_FirstPulse.mat'))
        load(dir('*sp_xia_FirstPulse.mat').name,'sp_seq');
        Seq(f).sp = sp_seq;
    else
        load(dir('*sp_xia.mat').name,'sp_clipped');
        Seq(f).sp = sp_clipped;
    end

    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
    Seq(f).trig = loadTrig(0);

    S = load(dir('*_exp_datafile_*.mat').name, ...
             'StimParams','simultaneous_stim','E_MAP','n_Trials');
    Stim = S.StimParams;
    simN = S.simultaneous_stim;

    Seq(f).nTrials   = S.n_Trials;
    Seq(f).simN      = simN;
    Seq(f).E_MAP     = S.E_MAP;

    amps_all = cell2mat(Stim(2:end,16));
    Seq(f).trialAmps = amps_all(1:simN:end);

    PTD_all = cell2mat(Stim(2:end,6));      
    Seq(f).PTD = PTD_all(2:simN:end) / 1000;

    stimNames = Stim(2:end,1);
    [~, idx]  = ismember(stimNames, Seq(f).E_MAP(2:end));

    stimChPerTrial = cell(Seq(f).nTrials,1);
    for t = 1:Seq(f).nTrials
        rr = (t-1)*simN + (1:simN);
        v = idx(rr); 
        v = v(v>0);
        stimChPerTrial{t} = v(:).';
    end
    Seq(f).stimChPerTrial = stimChPerTrial;

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

%% ==================== FIRING RATE ====================
fprintf('\n===== FIRING RATE RESULTS (Corrected stim-set grouping + Subsampling) =====\n');

nF = numel(seq_folders);

allSets = [];
for f = 1:nF
    tmp = Seq(f).uniqueComb;
    allSets = [allSets; tmp];
end

allSets(allSets==0) = nan;
[globalSets,~,~] = unique(allSets,'rows','stable'); 
nGlobal = size(globalSets,1);

FR_total  = nan(nF, nGlobal);
SEM_total = nan(nF, nGlobal);
PTDvals   = zeros(nF,1);

for f = 1:nF
    PTD = Seq(f).PTD(1);
    PTDvals(f) = PTD;

    sp_ch = Seq(f).sp{recCh};
    trig  = Seq(f).trig;

    for g = 1:nGlobal
        thisSet = globalSets(g, :);
        thisSet = thisSet(~isnan(thisSet));

        % find matching local stimulation set
        setIdx = [];
        for s = 1:Seq(f).nSets
            localSet = Seq(f).uniqueComb(s,:);
            localSet = localSet(localSet>0);
            if isequal(thisSet, localSet)
                setIdx = s;
                break;
            end
        end
        if isempty(setIdx), continue; end

        % original trial selection
        trial_ids = find(Seq(f).combClass == setIdx & Seq(f).trialAmps == plot_amps);

        % ========= SUBSAMPLE OPTION (only new code) =========
        if ~isempty(subsample_N) && numel(trial_ids) > subsample_N
            trial_ids = trial_ids(randperm(numel(trial_ids), subsample_N));
        end
        % ====================================================

        if isempty(trial_ids), continue; end

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

        fprintf("PTD=%2g ms | Set [%s] | FR %.2f ± %.2f | n=%d\n", ...
            PTD, strjoin(string(thisSet),'→'), FR_total(f,g), SEM_total(f,g), nT);
    end
end

%% ==================== PLOT ====================
figure('Color','w','Position',[200 200 900 550]); 
hold on;

colors = lines(nGlobal);
hSet = gobjects(nGlobal,1);

for g = 1:nGlobal
    y  = FR_total(:,g);
    se = SEM_total(:,g);

    % skip stim sets that are missing for all PTDs
    if all(isnan(y)), continue; end

    % ---- Main curve (line + markers) ----
    hLine = plot(PTDvals, y, '-o', ...
        'Color', colors(g,:), ...
        'LineWidth', 2.5, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors(g,:), ...
        'MarkerEdgeColor', 'k');

    hSet(g) = hLine;

    % ---- Error bars ----
    errorbar(PTDvals, y, se, ...
        'LineStyle', 'none', ...
        'Color', colors(g,:), ...
        'LineWidth', 1.2);
end

% ---- Labels ----
xlabel('Inter-Stimulus Interval (ISI) / ms');
ylabel('Firing Rate (spikes/s)');
title(sprintf('Firing Rate — Ch %d (window = %d ms)', target_channel, win_ms));
box off;

% ---- Legend: use the global stimulation-set labeling ----
legendStrings = {};
legendHandles = [];

for g = 1:nGlobal
    y = FR_total(:,g);
    if all(isnan(y)), continue; end   % skip unused sets

    setG = globalSets(g, ~isnan(globalSets(g,:)));  % remove padding NaN
    labelStr = sprintf('Set [%s]', strjoin(string(setG),'→'));

    legendStrings{end+1} = labelStr;
    legendHandles(end+1) = hSet(g);
end

legend(legendHandles, legendStrings, ...
    'Box','off', 'Location','northwest');