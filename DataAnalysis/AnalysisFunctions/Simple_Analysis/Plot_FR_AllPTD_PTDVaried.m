%% ---------------- USER INPUTS ----------------
data_folder   = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim4_251128_150648';
Electrode_Type = 2;
target_channel = 31;
plot_amps      = 6;

FS     = 30000;
win_ms = 4;           % analysis window after each pulse
subsample_N = [];     % set [] to disable subsampling

%% ---------------- LOAD DATA ----------------
cd(data_folder);

% --- spikes (FirstPulse preferred) ---
if ~isempty(dir('*sp_xia_FirstPulse.mat'))
    load(dir('*sp_xia_FirstPulse.mat').name, 'sp_seq');
    sp = sp_seq;
else
    load(dir('*sp_xia.mat').name, 'sp_clipped');
    sp = sp_clipped;
end

% --- triggers ---
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);

% --- StimParams ---
S = load(dir('*_exp_datafile_*.mat').name, ...
         'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams;
simN = S.simultaneous_stim;
nTrials = S.n_Trials;
E_MAP   = S.E_MAP;

% --- amplitude per trial (col 16) ---
amps_all  = cell2mat(Stim(2:end,16));
trialAmps = amps_all(1:simN:end);

% --- PTD per trial (col 6, µs → ms) ---
PTD_all_us = cell2mat(Stim(2:end,6));   % one row per pulse
% row 2 of each trial block is the delayed pulse
PTD_trial_ms = PTD_all_us(2:simN:end) / 1000;   % size: [nTrials x 1]
PTDvals = unique(PTD_trial_ms);                % all PTD values present (ms)

% --- ordered stimulation sets per trial ---
stimNames = Stim(2:end,1);
[~, idx]  = ismember(stimNames, E_MAP(2:end));

stimChPerTrial = cell(nTrials,1);
for t = 1:nTrials
    rr = (t-1)*simN + (1:simN);
    v  = idx(rr);
    v  = v(v>0);
    stimChPerTrial{t} = v(:).';       % keep order
end

comb = zeros(nTrials, simN);
for t = 1:nTrials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end

[uniqueComb,~,combClass] = unique(comb,'rows','stable');
nSets = size(uniqueComb,1);

%% ---------------- DEPTH MAP ----------------
d = Depth_s(Electrode_Type);
recCh = d(target_channel);

%% ---------------- FIRING RATE ----------------
fprintf('\n===== FIRING RATE RESULTS (one folder, multi-PTD, subsampled) =====\n');

nPTD = numel(PTDvals);

FR_total  = nan(nPTD, nSets);
SEM_total = nan(nPTD, nSets);

sp_ch = sp{recCh};

for g = 1:nSets
    thisSet = uniqueComb(g,:);
    thisSet = thisSet(thisSet>0);

    for p = 1:nPTD
        PTD = PTDvals(p);

        % trials with: this set, chosen amplitude, this PTD
        trial_ids = find( combClass == g & ...
                          trialAmps == plot_amps & ...
                          abs(PTD_trial_ms - PTD) < 1e-6 );

        if isempty(trial_ids)
            continue;
        end

        % ---- optional subsampling per PTD & set ----
        if ~isempty(subsample_N) && numel(trial_ids) > subsample_N
            trial_ids = trial_ids(randperm(numel(trial_ids), subsample_N));
        end

        nT = numel(trial_ids);
        fr_all = zeros(nT,1);

        for ii = 1:nT
            tr = trial_ids(ii);
            t0 = trig(tr)/FS*1000;

            % 1st pulse window: [0, win_ms]
            mask1 = sp_ch(:,1) >= t0          & sp_ch(:,1) < t0 + win_ms;

            % 2nd pulse window: [PTD, PTD+win_ms]
            mask2 = sp_ch(:,1) >= t0 + PTD    & sp_ch(:,1) < t0 + PTD + win_ms;

            n_spk = sum(mask1) + sum(mask2);
            fr_all(ii) = n_spk / ((2*win_ms)/1000);   % spikes / second
        end

        FR_total(p,g)  = mean(fr_all);
        SEM_total(p,g) = std(fr_all) / sqrt(nT);

        fprintf('Set [%s] | PTD = %2g ms | FR = %.2f ± %.2f sp/s | n=%d\n', ...
            strjoin(string(thisSet),'→'), PTD, FR_total(p,g), SEM_total(p,g), nT);
    end
end

%% ---------------- PLOT ----------------
figure('Color','w','Position',[200 200 900 550]); 
hold on;

colors = lines(nSets);
hSet   = gobjects(nSets,1);

for g = 1:nSets
    y  = FR_total(:,g);
    se = SEM_total(:,g);

    if all(isnan(y)), continue; end

    hLine = plot(PTDvals, y, '-o', ...
        'Color', colors(g,:), ...
        'LineWidth', 2.5, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors(g,:), ...
        'MarkerEdgeColor', 'k');

    hSet(g) = hLine;

    errorbar(PTDvals, y, se, ...
        'LineStyle','none', ...
        'Color', colors(g,:), ...
        'LineWidth',1.2);
end

xlabel('Post-Trigger Delay (PTD, ms)');
ylabel('Firing Rate (spikes/s)');
title(sprintf('Firing Rate — Ch %d | Amp = %.0f µA | win = %d ms', ...
    target_channel, plot_amps, win_ms));
box off;

% legend
legendStrings = cell(nSets,1);
for g = 1:nSets
    setG = uniqueComb(g, uniqueComb(g,:)>0);
    legendStrings{g} = sprintf('Set [%s]', strjoin(string(setG),'→'));
end
legend(hSet, legendStrings, 'Box','off', 'Location','northwest');