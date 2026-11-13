%% === PSTH Comparison per Channel, Stim Set and Amplitude === %%
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% --- USER INPUTS --- %%
single_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single5';
simul_folder  = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim5';
seq_folder    = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq5_5ms';

target_channels = [30];  
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

   %% === PASS 2: PSTH OVERLAY PLOT WITH CHANNEL LABELS IN LEGEND ===
        for ch = target_channels
            figure('Color','w','Position',[100 100 950 550]); hold on;
            title(sprintf('PSTH Comparison — Recording Ch %d', ch), 'FontSize', 14);
            xlabel('Time (ms)'); ylabel('Rate (sp/s)');
            box off;
        
            legend_entries = {};
            legend_handles = [];
        
            %% === Color settings ===
            color_single = lines(size(D_single.uniqueComb,1));  % each single set different color
            color_sum    = [0 0 0];
            color_sim    = [0.1 0.7 0.2];
            color_seq    = [0.5 0.3 0.5];
        
            linestyle_single = '-';
            linestyle_sum    = ':';
            linestyle_sim    = '-';
            linestyle_seq    = '-';
        
            single_psth_store = {};
        
            % 1. SINGLE STIMULATION SETS (each set has its own color)

            for s = 1:size(D_single.uniqueComb,1)
                D = D_single;
        
                stimIdx = D.uniqueComb(s, D.uniqueComb(s,:)>0);
                stimChLabel = strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimIdx,'UniformOutput',false),' + ');
        
                col_this = color_single(s,:);
        
                for a = 1:numel(target_amps)
                    amp_val = target_amps(a);
        
                    amp_trials = find(D.trialAmps==amp_val & D.combClass_win==s);
                    if isempty(amp_trials), continue; end
        
                    counts = zeros(1,numel(edges)-1);
                    for t = 1:numel(amp_trials)
                        t0 = D.trig(amp_trials(t))/FS*1000;
                        S_ch = D.sp{d(ch)};
                        tt = S_ch(:,1);
                        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                        counts = counts + histcounts(tt, edges);
                    end
        
                    rate = counts / (numel(amp_trials)*bin_s);
                    rate_s = filter(g,1,rate);
                    single_psth_store{s} = rate_s;
        
                    h = plot(ctrs, rate_s, ...
                        'Color', col_this, ...
                        'LineStyle', linestyle_single, ...
                        'LineWidth', 1.7);
        
                    legend_handles(end+1) = h;
                    legend_entries{end+1} = sprintf( ...
                        'Single Set %d — %s — %.0f µA', ...
                        s, stimChLabel, amp_val);
                end
            end
        
            % % 2. SUM OF TWO SINGLE SETS (dashed black)
            % if size(D_single.uniqueComb,1) == 2
            %     stimA = D_single.uniqueComb(1, D_single.uniqueComb(1,:)>0);
            %     stimB = D_single.uniqueComb(2, D_single.uniqueComb(2,:)>0);
            %     stimSumLabel = [ ...
            %         strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimA,'UniformOutput',false),' + ') ...
            %         '  +  ' ...
            %         strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimB,'UniformOutput',false),' + ') ...
            %     ];
            % 
            %     psth_sum = single_psth_store{1} + single_psth_store{2};
            % 
            %     h = plot(ctrs, psth_sum, ...
            %         'Color', color_sum, ...
            %         'LineStyle', linestyle_sum, ...
            %         'LineWidth', 1.5);
            % 
            %     legend_handles(end+1) = h;
            %     legend_entries{end+1} = sprintf('Sum of Singles — %s', stimSumLabel);
            % end
        
            % 3. SIMULTANEOUS STIM (solid green)
            D = D_sim;
        
            for set_id = 1:size(D.uniqueComb,1)
                stimIdx = D.uniqueComb(set_id,D.uniqueComb(set_id,:)>0);
                stimChLabel = strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimIdx,'UniformOutput',false),' + ');
        
                for a = 1:numel(target_amps)
                    amp_val = target_amps(a);
        
                    amp_trials = find(D.trialAmps==amp_val & D.combClass_win==set_id);
                    if isempty(amp_trials), continue; end
        
                    counts = zeros(1,numel(edges)-1);
                    for t = 1:numel(amp_trials)
                        t0 = D.trig(amp_trials(t))/FS*1000;
                        tt = D.sp{d(ch)}(:,1);
                        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                        counts = counts + histcounts(tt, edges);
                    end
        
                    rate_s = filter(g,1,counts/(numel(amp_trials)*bin_s));
        
                    h = plot(ctrs, rate_s, ...
                        'Color', color_sim, ...
                        'LineStyle', linestyle_sim, ...
                        'LineWidth', 2);
        
                    legend_handles(end+1) = h;
                    legend_entries{end+1} = sprintf( ...
                        'Simultaneous — %s — %.0f µA', stimChLabel, amp_val);
                end
            end
        
            % 4. SEQUENTIAL STIM (red dotted)
            D = D_seq;
        
            for set_id = 1:size(D.uniqueComb,1)
                stimIdx = D.uniqueComb(set_id,D.uniqueComb(set_id,:)>0);
                stimChLabel = strjoin(arrayfun(@(x)sprintf('Ch%d',x),stimIdx,'UniformOutput',false),' + ');
        
                for a = 1:numel(target_amps)
                    amp_val = target_amps(a);
        
                    amp_trials = find(D.trialAmps==amp_val & D.combClass_win==set_id);
                    if isempty(amp_trials), continue; end
        
                    counts = zeros(1,numel(edges)-1);
                    for t = 1:numel(amp_trials)
                        t0 = D.trig(amp_trials(t))/FS*1000;
                        tt = D.sp{d(ch)}(:,1);
                        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                        counts = counts + histcounts(tt, edges);
                    end
        
                    rate_s = filter(g,1,counts/(numel(amp_trials)*bin_s));
        
                    h = plot(ctrs, rate_s, ...
                        'Color', color_seq, ...
                        'LineStyle', linestyle_seq, ...
                        'LineWidth', 2);
        
                    legend_handles(end+1) = h;
                    legend_entries{end+1} = sprintf( ...
                        'Sequential — %s — %.0f µA', stimChLabel, amp_val);
                end
            end
        
            %% FINAL FORMATTING
            xline(0,'r--','LineWidth',1.3,'HandleVisibility','off');
            xlim(ras_win); ylim(global_ylim);
            box off;
        
            legend(legend_handles, legend_entries, 'Location','northeast','Box','off');
        end
end

fprintf('\nGenerated PSTHs for channels: %s\n', num2str(target_channels));