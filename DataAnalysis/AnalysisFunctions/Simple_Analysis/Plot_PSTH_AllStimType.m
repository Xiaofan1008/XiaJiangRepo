%% === PSTH Comparison per Channel, Stim Set and Amplitude === %%
clear; clc;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% --- USER INPUTS --- %%
single_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single1';
simul_folder  = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim1';
seq_folder    = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq1';

target_channels = [28];  
target_amps     = [5];  % µA
ras_win   = [-20 100];    % ms
bin_ms    = 2;            % ms
smooth_ms = 3;            % ms
FS = 30000;               % Hz
cmap = lines(numel(target_amps));

%% === Helper: Load spike + stim info === %%
function D = loadStimData(folder, useFirstPulse)
    cd(folder);
    if useFirstPulse
        load(dir('*sp_xia_FirstPulse.mat').name,'sp_seq');
        D.sp = sp_seq;
    else
        load(dir('*sp_xia.mat').name,'sp_clipped');
        D.sp = sp_clipped;
    end
    D.trig = loadTrig(0);
    S = load(dir('*_exp_datafile_*.mat').name, ...
        'StimParams','simultaneous_stim','E_MAP','n_Trials');
    StimParams = S.StimParams;
    D.simultaneous_stim = S.simultaneous_stim;
    D.n_Trials = S.n_Trials;
    D.E_MAP = S.E_MAP;

    trialAmps_all = cell2mat(StimParams(2:end,16));
    D.trialAmps = trialAmps_all(1:D.simultaneous_stim:end);
    stimNames = StimParams(2:end,1);
    [~, idx_all] = ismember(stimNames, D.E_MAP(2:end));
    stimChPerTrial = cell(D.n_Trials,1);
    for t = 1:D.n_Trials
        rr = (t-1)*D.simultaneous_stim + (1:D.simultaneous_stim);
        v = unique(idx_all(rr)); v = v(v>0).';
        stimChPerTrial{t} = v;
    end
    comb = zeros(D.n_Trials, D.simultaneous_stim);
    for i = 1:D.n_Trials
        v = stimChPerTrial{i};
        comb(i,1:numel(v)) = v;
    end
    [D.uniqueComb,~,D.combClass] = unique(comb,'rows');
    D.combClass_win = D.combClass;
end

%% === Load data === %%
fprintf('\nLoading datasets...\n');
D_single = loadStimData(single_folder,false);
D_sim    = loadStimData(simul_folder,false);
D_seq    = loadStimData(seq_folder,true);
datasets = {D_single, D_sim, D_seq};

%% === PSTH setup === %%
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms / 1000;
g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
g = g / sum(g);

d = Depth_s(1); % channel depth mapping

%% === Loop over channels === %%
for ch = target_channels
    n_single_sets = size(D_single.uniqueComb,1);
    total_subplots = n_single_sets + 2;

    %% --- Pass 1: find global max firing rate across all subplots ---
    global_maxRate = 0;
    all_D = {D_single, D_sim, D_seq};
    all_set_counts = [size(D_single.uniqueComb,1), 1, 1]; % single has multiple sets

    for group_id = 1:3
        D = all_D{group_id};
        nSets = all_set_counts(group_id);
        for s = 1:nSets
            for a = 1:numel(target_amps)
                amp_val = target_amps(a);
                amp_trials = find(D.trialAmps == amp_val & D.combClass_win == s);
                if isempty(amp_trials), continue; end
                counts = zeros(1,numel(edges)-1);
                for t = 1:numel(amp_trials)
                    t0 = D.trig(amp_trials(t))/FS*1000;
                    S_ch = D.sp{d(ch)};
                    if isempty(S_ch), continue; end
                    tt = S_ch(:,1);
                    tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                    counts = counts + histcounts(tt, edges);
                end
                rate = counts / (numel(amp_trials)*bin_s);
                rate_s = filter(g,1,rate);
                global_maxRate = max(global_maxRate, max(rate_s));
            end
        end
    end
    global_ylim = [0 ceil(global_maxRate*1.1/10)*10]; % nice rounded upper limit

    %% --- Pass 2: actual plotting ---
    figure('Color','w','Position',[100 100 950 750]);
    sgtitle(sprintf('PSTH Channel %d', ch));

    % --- SINGLE STIMULATION (each set separately) ---
    for s = 1:n_single_sets
        subplot(total_subplots,1,s); hold on; box off;
        D = D_single;
        stimIdx = D.uniqueComb(s,:);
        stimIdx = stimIdx(stimIdx>0);
        setLabel = strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimIdx,'UniformOutput',false),' + ');

        for a = 1:numel(target_amps)
            amp_val = target_amps(a);
            amp_trials = find(D.trialAmps==amp_val & D.combClass_win==s);
            if isempty(amp_trials), continue; end
            counts = zeros(1,numel(edges)-1);
            total_spikes = 0;
            for t = 1:numel(amp_trials)
                t0 = D.trig(amp_trials(t))/FS*1000;
                S_ch = D.sp{d(ch)};
                if isempty(S_ch), continue; end
                tt = S_ch(:,1);
                tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                counts = counts + histcounts(tt, edges);
                total_spikes = total_spikes + numel(tt);
            end
            rate = counts / (numel(amp_trials)*bin_s);
            rate_s = filter(g,1,rate);
            plot(ctrs, rate_s, 'Color', cmap(a,:), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%.0f µA (%d spikes)', amp_val, total_spikes));
        end
        xline(0,'r--','LineWidth',2,'HandleVisibility','off');
        xlim(ras_win); ylim(global_ylim);
        ylabel('Rate (sp/s)');
        title(sprintf('Separate stim | Set: %s', setLabel));
        legend('show','Location','northeast','Box','off');
    end

    % --- SIMULTANEOUS STIMULATION ---
    subplot(total_subplots,1,n_single_sets+1); hold on; box off;
    D = D_sim;
    for set_id = 1:size(D.uniqueComb,1)
        stimIdx = D.uniqueComb(set_id,:);
        stimIdx = stimIdx(stimIdx>0);
        setLabel = strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimIdx,'UniformOutput',false),' + ');
        for a = 1:numel(target_amps)
            amp_val = target_amps(a);
            amp_trials = find(D.trialAmps==amp_val & D.combClass_win==set_id);
            if isempty(amp_trials), continue; end
            counts = zeros(1,numel(edges)-1);
            total_spikes = 0;
            for t = 1:numel(amp_trials)
                t0 = D.trig(amp_trials(t))/FS*1000;
                S_ch = D.sp{d(ch)};
                if isempty(S_ch), continue; end
                tt = S_ch(:,1);
                tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                counts = counts + histcounts(tt, edges);
                total_spikes = total_spikes + numel(tt);
            end
            rate = counts / (numel(amp_trials)*bin_s);
            rate_s = filter(g,1,rate);
            plot(ctrs, rate_s, 'Color', cmap(a,:), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%.0f µA | %s (%d spikes)', amp_val, setLabel, total_spikes));
        end
    end
    xline(0,'r--','LineWidth',2,'HandleVisibility','off');
    xlim(ras_win); ylim(global_ylim);
    ylabel('Rate (sp/s)');
    title('Simultaneous stimulation');
    legend('show','Location','northeast','Box','off');

    % --- SEQUENTIAL STIMULATION ---
    subplot(total_subplots,1,n_single_sets+2); hold on; box off;
    D = D_seq;
    for set_id = 1:size(D.uniqueComb,1)
        stimIdx = D.uniqueComb(set_id,:);
        stimIdx = stimIdx(stimIdx>0);
        setLabel = strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimIdx,'UniformOutput',false),' + ');
        for a = 1:numel(target_amps)
            amp_val = target_amps(a);
            amp_trials = find(D.trialAmps==amp_val & D.combClass_win==set_id);
            if isempty(amp_trials), continue; end
            counts = zeros(1,numel(edges)-1);
            total_spikes = 0;
            for t = 1:numel(amp_trials)
                t0 = D.trig(amp_trials(t))/FS*1000;
                S_ch = D.sp{d(ch)};
                if isempty(S_ch), continue; end
                tt = S_ch(:,1);
                tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                counts = counts + histcounts(tt, edges);
                total_spikes = total_spikes + numel(tt);
            end
            rate = counts / (numel(amp_trials)*bin_s);
            rate_s = filter(g,1,rate);
            plot(ctrs, rate_s, 'Color', cmap(a,:), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%.0f µA | %s (%d spikes)', amp_val, setLabel, total_spikes));
        end
    end
    xline(0,'r--','LineWidth',2,'HandleVisibility','off');
    xlim(ras_win); ylim(global_ylim);
    xlabel('Time (ms)');
    ylabel('Rate (sp/s)');
    title('Sequential stimulation');
    legend('show','Location','northeast','Box','off');
end

fprintf('\nGenerated PSTHs for channels: %s\n', num2str(target_channels));