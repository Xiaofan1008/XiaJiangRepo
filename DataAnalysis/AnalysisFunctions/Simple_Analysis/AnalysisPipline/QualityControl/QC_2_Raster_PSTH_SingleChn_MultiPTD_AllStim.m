%% ============================================================
%   Unified Raster + PSTH Comparison
%   (Combined Dataset: Distinct Colors per PTD)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ===================== USER INPUTS ===================== %%
folder_single = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Single2';
folder_seq    = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim2'; 

Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% --- Analysis Settings ---
target_channels = [5:10]; 
plot_amps      = [10];       % Amplitudes to plot
target_PTDs    = [0 5 10];  % Specific Seq PTDs to plot. Leave [] for ALL > 0.

ras_win        = [-10 50];  % ms
bin_ms         = 1;
smooth_ms      = 5;
FS             = 30000;

%% Colors
pastel_colors = [
    0.4 0.6 0.85;   % pastel blue
    0.9 0.6 0.5;    % pastel salmon
    0.6 0.8 0.6;    % pastel green
    0.85 0.65 0.85; % pastel purple
    0.95 0.8 0.5;   % pastel orange
];
col_sim = [0.3 0.65 0.3];    % soft green (Simultaneous)

% [MODIFIED] PTD Color Map (Distinct colors for up to 7 PTDs)
ptd_base_colors = lines(7); 

%% ============================================================
% LOAD DATASETS
% ============================================================

% ---------- Load Single ----------
cd(folder_single)
f = dir('*sp_xia_SSD.mat');
if isempty(f), error('No SSD file in Single folder'); end
tmp = load(f(1).name);
if isfield(tmp, 'sp_corr'), D_single.sp = tmp.sp_corr; 
elseif isfield(tmp, 'sp_SSD'), D_single.sp = tmp.sp_SSD; 
else, D_single.sp = tmp.sp_in; end
if isempty(dir('*.trig.dat')); cleanTrig_sabquick; end
D_single.trig = loadTrig(0);
S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams; simN = S.simultaneous_stim; D_single.E_MAP = S.E_MAP; D_single.nTrials = S.n_Trials;
amps_all = cell2mat(Stim(2:end,16)); D_single.trialAmps = amps_all(1:simN:end);
stimNames = Stim(2:end,1); [~, idx]  = ismember(stimNames, D_single.E_MAP(2:end));
comb = zeros(D_single.nTrials, simN);
for t = 1:D_single.nTrials, rr = (t-1)*simN + (1:simN); v = unique(idx(rr)); v = v(v>0).'; comb(t,1:numel(v)) = v; end
[D_single.uniqueComb,~,D_single.combClass] = unique(comb,'rows');

% ---------- Load Combined (Sim + Seq) ----------
cd(folder_seq)
f = dir('*sp_xia_FirstPulse.mat');
if isempty(f), error('No SSD/FirstPulse file in Combined folder'); end
tmp = load(f(1).name);
if isfield(tmp, 'sp_seq'), D_seq.sp = tmp.sp_seq; 
elseif isfield(tmp, 'sp_SSD'), D_seq.sp = tmp.sp_SSD; 
else, D_seq.sp = tmp.sp_in; end
if isempty(dir('*.trig.dat')); cleanTrig_sabquick; end
D_seq.trig = loadTrig(0);
S_seq = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim_seq = S_seq.StimParams; simN_seq = S_seq.simultaneous_stim; D_seq.E_MAP = S_seq.E_MAP; D_seq.nTrials = S_seq.n_Trials;
amps_all = cell2mat(Stim_seq(2:end,16)); D_seq.trialAmps = amps_all(1:simN_seq:end);
stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, D_seq.E_MAP(2:end));
comb_seq = zeros(D_seq.nTrials, simN_seq);
for t = 1:D_seq.nTrials, rr = (t-1)*simN_seq + (1:simN_seq); v  = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[D_seq.uniqueComb,~,D_seq.combClass] = unique(comb_seq,'rows','stable');
D_seq.nSets = size(D_seq.uniqueComb,1);

% Extract PTDs (0 = Sim, >0 = Seq)
PTD_all_us = cell2mat(Stim_seq(2:end,6)); 
D_seq.trialPTD = PTD_all_us(2:simN_seq:end) / 1000; 

% Determine PTDs to plot
available_PTDs = unique(D_seq.trialPTD);
available_PTDs(available_PTDs == 0) = []; % Remove Sim
if isempty(target_PTDs)
    use_PTDs = available_PTDs;
else
    use_PTDs = intersect(available_PTDs, target_PTDs);
end
fprintf('Plotting Sequential PTDs: %s ms\n', num2str(use_PTDs'));

%% PSTH KERNEL & DEPTH
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;
g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g/sum(g);
d = Depth_s(Electrode_Type);

%% ===================== MAIN PLOTTING LOOPS =====================
for ch_idx = 1:length(target_channels)
    target_channel = target_channels(ch_idx);
    for amp_val = plot_amps
        figure('Color','w','Position',[200 100 800 1400]);
        
        nSingleSets = size(D_single.uniqueComb,1);
        nSeqSets    = D_seq.nSets;
        nPTD        = length(use_PTDs);
        
        % Tiles: Single + Sim + (SeqSets * PTDs) + PSTHs
        total_tiles = nSingleSets + 1 + (nSeqSets * nPTD) + 2; 
        
        tl = tiledlayout(total_tiles, 1, 'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('Channel %d — %.0f µA',target_channel,amp_val),'FontSize',16);
        max_psth_rate = 0;
        
        % ============ 1. SINGLE RASTERS =============
        for s = 1:nSingleSets
            ax = nexttile(tl); set(ax,'NextPlot','add'); box(ax,'off');
            stimCh = D_single.uniqueComb(s, D_single.uniqueComb(s,:)>0);
            S_ch = D_single.sp{d(target_channel)};
            trial_ids = find(D_single.combClass==s & D_single.trialAmps==amp_val);
            y = 0;
            for tr = trial_ids'
                t0 = D_single.trig(tr)/FS*1000;
                tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                for k=1:numel(tt), plot(ax,[tt(k) tt(k)],[y y+0.8], 'Color', pastel_colors(mod(s-1,5)+1,:), 'LineWidth',1.5); end
                y=y+1;
            end
            title(ax,sprintf('Single Set %d — Ch%s', s, num2str(stimCh)),'FontSize',11);
            xline(ax,0,'r--','HandleVisibility','off');
            xlim(ax,ras_win); ylim(ax,[0 max(1,y)]);
        end
        
        % ============ 2. SIMULTANEOUS RASTER (PTD = 0) =============
        ax = nexttile(tl); set(ax,'NextPlot','add'); box(ax,'off');
        S_ch = D_seq.sp{d(target_channel)};
        trial_ids = find(D_seq.trialAmps==amp_val & D_seq.trialPTD == 0);
        
        if ~isempty(trial_ids)
            sim_set_id = D_seq.combClass(trial_ids(1));
            stimCh = D_seq.uniqueComb(sim_set_id, D_seq.uniqueComb(sim_set_id,:)>0);
        else
            stimCh = [NaN];
        end
        
        y = 0;
        for tr = trial_ids'
            t0 = D_seq.trig(tr)/FS*1000;
            tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
            for k=1:numel(tt), plot(ax,[tt(k) tt(k)],[y y+0.8], 'Color', col_sim, 'LineWidth',1.5); end
            y=y+1;
        end
        title(ax, sprintf('Simultaneous (PTD 0) — Ch%s', num2str(stimCh)), 'FontSize',11);
        xline(ax,0,'r--','HandleVisibility','off');
        xlim(ax,ras_win); ylim(ax,[0 max(1,y)]);
        
        % ============ 3. SEQUENTIAL RASTERS (Color by PTD) =============
        for s = 1:nSeqSets
            stimCh = D_seq.uniqueComb(s, D_seq.uniqueComb(s,:)>0);
            
            for p_idx = 1:nPTD
                ptd_val = use_PTDs(p_idx);
                
                % [MODIFIED] Use PTD-Specific Color
                this_ptd_col = ptd_base_colors(mod(p_idx-1, size(ptd_base_colors,1))+1, :);
                
                ax = nexttile(tl); set(ax,'NextPlot','add'); box(ax,'off');
                S_ch = D_seq.sp{d(target_channel)};
                
                trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val & abs(D_seq.trialPTD - ptd_val) < 0.01);
                
                y=0;
                for tr = trial_ids'
                    t0 = D_seq.trig(tr)/FS*1000;
                    tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                    for k=1:numel(tt), plot(ax,[tt(k) tt(k)],[y y+0.8], 'Color', this_ptd_col, 'LineWidth',1.5); end
                    y=y+1;
                end
                
                title(ax, sprintf('Seq Set %d (%.0f ms) — Order: %s', s, ptd_val, num2str(stimCh)), 'FontSize',11);
                xline(ax,0,'r--','HandleVisibility','off');
                xline(ax,ptd_val,'k:','HandleVisibility','off'); 
                xlim(ax,ras_win); ylim(ax,[0 max(1,y)]);
            end
        end
        
        %% ============ 4. PSTH (Single Sets) ============
        ax_single_psth = nexttile(tl); box(ax_single_psth,'off'); hold(ax_single_psth,'on');
        title(ax_single_psth,'PSTH — Single Sets');    
        for s = 1:nSingleSets
            S_ch = D_single.sp{d(target_channel)};
            trial_ids = find(D_single.combClass==s & D_single.trialAmps==amp_val);
            if isempty(trial_ids), continue; end
            counts = zeros(1,length(edges)-1);
            for tr = trial_ids'
                t0 = D_single.trig(tr)/FS*1000;
                tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                counts = counts + histcounts(tt,edges);
            end 
            rate = counts / (numel(trial_ids) * bin_s);
            rate_s = filter(g,1,rate);
            plot(ax_single_psth,ctrs,rate_s,'Color',pastel_colors(mod(s-1,5)+1,:), 'LineWidth',2, 'DisplayName',sprintf('Single %d',s));
            max_psth_rate = max(max_psth_rate, max(rate_s));
        end    
        xline(ax_single_psth,0,'r--','HandleVisibility','off');
        xlim(ax_single_psth,ras_win); ylabel(ax_single_psth,'Rate (sp/s)');
        legend(ax_single_psth,'Box','off','Location','northeast');
        
        %% ============ 5. PSTH (Sim + All Seq Sets) ============
        ax_simseq_psth = nexttile(tl); box(ax_simseq_psth,'off'); hold(ax_simseq_psth,'on');
        title(ax_simseq_psth,'PSTH — Simultaneous vs Sequential (Color by PTD)');
        
        % --- Sim (PTD = 0) ---
        S_ch = D_seq.sp{d(target_channel)};
        trial_ids = find(D_seq.trialAmps==amp_val & D_seq.trialPTD == 0);
        if ~isempty(trial_ids)
            counts = zeros(1,length(edges)-1);
            for tr = trial_ids'
                t0 = D_seq.trig(tr)/FS*1000;
                tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                counts = counts + histcounts(tt,edges);
            end
            rate_sim = filter(g,1, counts/(numel(trial_ids)*bin_s));
            plot(ax_simseq_psth,ctrs,rate_sim,'Color',col_sim,'LineWidth',3,'DisplayName','Simultaneous');
            max_psth_rate = max(max_psth_rate, max(rate_sim));
        end
        
        % --- Seq (Loop PTD > 0) ---
        for s = 1:nSeqSets
            for p_idx = 1:nPTD
                ptd_val = use_PTDs(p_idx);
                
                % [MODIFIED] Use PTD-Specific Color
                this_ptd_col = ptd_base_colors(mod(p_idx-1, size(ptd_base_colors,1))+1, :);
                
                S_ch = D_seq.sp{d(target_channel)};
                trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val & abs(D_seq.trialPTD - ptd_val) < 0.01);
                
                if isempty(trial_ids), continue; end
                
                counts = zeros(1,length(edges)-1);
                for tr = trial_ids'
                    t0 = D_seq.trig(tr)/FS*1000;
                    tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                    counts = counts + histcounts(tt,edges);
                end
                rate_seq = filter(g,1, counts/(numel(trial_ids)*bin_s));
                
                % Solid lines for all, distinguished by Color
                plot(ax_simseq_psth,ctrs,rate_seq,'Color',this_ptd_col,'LineStyle','-', 'LineWidth',2, ...
                    'DisplayName',sprintf('Seq Set %d (%.0f ms)', s, ptd_val));
                
                max_psth_rate = max(max_psth_rate, max(rate_seq));
            end
        end
        
        xline(ax_simseq_psth,0,'r--','HandleVisibility','off');
        xlim(ax_simseq_psth,ras_win); ylabel(ax_simseq_psth,'Rate (sp/s)'); xlabel(ax_simseq_psth,'Time (ms)');
        legend(ax_simseq_psth,'Box','off','Location','northeast');
        if max_psth_rate > 0
            ymax = max_psth_rate * 1.1;
            ylim(ax_single_psth, [0 ymax]);
            ylim(ax_simseq_psth, [0 ymax]);
        end
    end
end