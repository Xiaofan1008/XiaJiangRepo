%% Plot rasters for single vs sequential (before & after spike insertion)
clear all;

%% ===== User parameters =====
single_folder     = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Single1_251203_121428';
sequential_folder = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Seq_Sim1_251203_121823';

rec_ch_plot       = 58;          % recording channel index in sp{ch}
amps_to_plot      = [];          % [] = all amplitudes
ptd_to_plot_ms    = [];          % [] = all PTDs (ms), non-zero only
max_trials        = 40;          % max trials per condition to show in raster
ras_win_ms        = [-10 50];    % raster window around first pulse (ms)

FS = 30000;

%% ===== Load single-pulse data =====
cd(single_folder);
Ssp = load(dir('*sp_xia.mat').name, 'sp_clipped');
sp_single   = Ssp.sp_clipped;
nChn_single = numel(sp_single);
trig_single = loadTrig(0);

S1 = load(dir('*_exp_datafile_*.mat').name, ...
          'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_single  = S1.StimParams;
E_MAP              = S1.E_MAP;
n_Trials_single    = S1.n_Trials;
simultaneous_stim1 = S1.simultaneous_stim;

% Trial amplitudes for single
trialAmps_single = cell2mat(StimParams_single(2:end,16));
trialAmps_single = trialAmps_single(1:simultaneous_stim1:end);

% Stim channels per single trial (usually 1)
stimNames_single  = StimParams_single(2:end,1);
[~, idx_all_single] = ismember(stimNames_single, E_MAP(2:end));
stimChPerTrial_single = cell(n_Trials_single,1);
for t = 1:n_Trials_single
    rr = (t-1)*simultaneous_stim1 + (1:simultaneous_stim1);
    v  = idx_all_single(rr);
    v  = v(v>0);
    stimChPerTrial_single{t} = v(:).';
end

%% ===== Load sequential data: BEFORE insertion =====
cd(sequential_folder);
Sseq_before = load(dir('*sp_xia.mat').name, 'sp_clipped');
sp_seq_before = Sseq_before.sp_clipped;
trig_seq      = loadTrig(0);

S2 = load(dir('*_exp_datafile_*.mat').name, ...
          'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_seq    = S2.StimParams;
n_Trials_seq      = S2.n_Trials;
simultaneous_stim2 = S2.simultaneous_stim;

% Trial amplitudes in sequential data
trialAmps_seq = cell2mat(StimParams_seq(2:end,16));
trialAmps_seq = trialAmps_seq(1:simultaneous_stim2:end);

% Stim channels per sequential trial (order-sensitive)
stimNames_seq = StimParams_seq(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP(2:end));
stimChPerTrial_seq = cell(n_Trials_seq,1);
for t = 1:n_Trials_seq
    rr = (t-1)*simultaneous_stim2 + (1:simultaneous_stim2);
    v  = idx_all_seq(rr);   % keep order
    v  = v(v>0);
    stimChPerTrial_seq{t} = v(:).';
end

% Order-sensitive stimulation sets
comb_seq = zeros(n_Trials_seq, simultaneous_stim2);
for t = 1:n_Trials_seq
    v = stimChPerTrial_seq{t};
    comb_seq(t,1:numel(v)) = v;
end
[uniqueComb_seq, ~, combClass_seq] = unique(comb_seq, 'rows', 'stable');
nSeqSets = size(uniqueComb_seq,1);

% PTD per sequential trial (µs → ms)
PTD_us_all = cell2mat(StimParams_seq(2:end,6));
PTD_us     = PTD_us_all(2:simultaneous_stim2:end);   % second row of each trial block
PTD_ms     = PTD_us / 1000;
isSimultaneous = (PTD_ms == 0);  % PTD = 0 → simultaneous, no injection

%% ===== Load sequential data: AFTER insertion =====
fsp = dir('*sp_xia_FirstPulse.mat');
assert(~isempty(fsp), 'No *_FirstPulse spike file found in sequential folder.');
Sseq_after = load(fsp(1).name, 'sp_seq');
sp_seq_after = Sseq_after.sp_seq;

%% ===== Selection of amps and PTDs =====
all_amps = unique(trialAmps_seq);
if isempty(amps_to_plot)
    amps_sel = all_amps;
else
    amps_sel = intersect(all_amps, amps_to_plot);
end

all_PTD_ms = unique(PTD_ms(~isSimultaneous));  % non-zero only
if isempty(ptd_to_plot_ms)
    PTD_sel = all_PTD_ms;
else
    PTD_sel = intersect(all_PTD_ms, ptd_to_plot_ms);
end

% fprintf('Amplitudes to plot: '); disp(amps_sel');
% fprintf('PTDs (ms) to plot:  '); disp(PTD_sel');

%% ===== Helper for raster plotting =====
edges = ras_win_ms(1):1:ras_win_ms(2);     % 1 ms bins
ctrs  = edges(1:end-1) + 0.5;

make_raster = @(ax, sp_cell, trig_vec, trials, color) ...
    raster_channel(ax, sp_cell{rec_ch_plot}, trig_vec, trials, ...
                   ras_win_ms, color);

%% ===== Main plotting loops =====
for a = 1:numel(amps_sel)
    amp = amps_sel(a);

    trials_amp = find(trialAmps_seq == amp);
    if isempty(trials_amp), continue; end

    for ptd_val = PTD_sel(:)'
        trials_ptd = trials_amp(PTD_ms(trials_amp) == ptd_val & ...
                                ~isSimultaneous(trials_amp));
        if isempty(trials_ptd), continue; end

        % Stimulation sets present at this amp+PTD
        sets_here = unique(combClass_seq(trials_ptd));

        for ss = sets_here(:)'
            stimVec = uniqueComb_seq(ss,:);
            stimVec = stimVec(stimVec>0);
            if isempty(stimVec), continue; end

            stim_label = sprintf('[%s]', strjoin(string(stimVec),'→'));
            ch_first   = stimVec(1);

            % Sequential trials in this condition
            cond_tr_seq = find(trialAmps_seq==amp & ...
                               PTD_ms==ptd_val & ...
                               combClass_seq==ss & ...
                               ~isSimultaneous);
            if isempty(cond_tr_seq), continue; end
            if numel(cond_tr_seq) > max_trials
                cond_tr_seq = cond_tr_seq(1:max_trials);
            end

            % Find matching single trials: same amp, stim channel contains ch_first
            match_single = find( ...
                cellfun(@(x) ismember(ch_first,x), stimChPerTrial_single) & ...
                (trialAmps_single == amp) );
            if isempty(match_single)
                warning('No matching single trial for amp %d, first stim ch %d', amp, ch_first);
                continue;
            end
            if numel(match_single) > max_trials
                cond_tr_single = match_single(1:max_trials);
            else
                cond_tr_single = match_single;
            end

            %% ----- Make figure for this (amp, PTD, set) -----
            figName = sprintf('Ch%d | Amp %d µA | PTD %g ms | Set %s', ...
                              rec_ch_plot, amp, ptd_val, stim_label);
            figure('Name',figName, 'Color','w', 'Position',[100 100 900 800]);
            tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
            title(tl, sprintf('Channel %d | Amp %d µA | PTD %g ms | Set %s', ...
                              rec_ch_plot, amp, ptd_val, stim_label), ...
                  'Interpreter','none');

            % ----- 1) Single stimulation raster -----
            ax1 = nexttile(tl); hold(ax1,'on');
            make_raster(ax1, sp_single, trig_single, cond_tr_single, [0 0 0]);
            xline(ax1, 0, 'r--');
            xlim(ax1, ras_win_ms);
            ylabel(ax1, 'Trial');
            title(ax1, sprintf('Single stim | stim ch %d', ch_first));
            box(ax1,'off');

            % ----- 2) Sequential BEFORE insertion -----
            ax2 = nexttile(tl); hold(ax2,'on');
            make_raster(ax2, sp_seq_before, trig_seq, cond_tr_seq, [0 0 1]);
            xline(ax2, 0, 'r--');
            xlim(ax2, ras_win_ms);
            ylabel(ax2, 'Trial');
            title(ax2, sprintf('Sequential BEFORE | set %s', stim_label));
            box(ax2,'off');

            % ----- 3) Sequential AFTER insertion -----
            ax3 = nexttile(tl); hold(ax3,'on');
            make_raster(ax3, sp_seq_after, trig_seq, cond_tr_seq, [0 0.6 0]);
            xline(ax3, 0, 'r--');
            xlim(ax3, ras_win_ms);
            xlabel(ax3, 'Time (ms)');
            ylabel(ax3, 'Trial');
            title(ax3, sprintf('Sequential AFTER | set %s', stim_label));
            box(ax3,'off');
        end
    end
end

fprintf('Finished plotting.\n');

%% ===== Local function: raster for one channel =====
function raster_channel(ax, sp_ch, trig_vec, trials, win_ms, color)
    if isempty(sp_ch), return; end
    FS = 30000;  % sampling rate in Hz (fixed here)
    nTr = numel(trials);
    for k = 1:nTr
        tr = trials(k);
        t0 = trig_vec(tr)/FS*1000;
        tt = sp_ch(:,1);
        rel_t = tt - t0;
        rel_t = rel_t(rel_t >= win_ms(1) & rel_t <= win_ms(2));
        if isempty(rel_t), continue; end
        y = k * ones(size(rel_t));
        plot(ax, rel_t, y, '.', 'Color', color, 'MarkerSize', 6);
    end
    ylim(ax, [0 nTr+1]);
end