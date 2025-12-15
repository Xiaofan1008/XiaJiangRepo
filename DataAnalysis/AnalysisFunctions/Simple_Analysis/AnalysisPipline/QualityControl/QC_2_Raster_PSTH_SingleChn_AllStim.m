%% ============================================================
%   Unified Raster + PSTH Comparison
%   (SSD Loading + Multi-Channel + Separate Sequential Orders)
%   + Clearer Seq Colors + Unified Y-Limits
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ===================== USER INPUTS ===================== %%
folder_single = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Single6_251125_180744';
folder_sim    = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim6_251125_181554';
folder_seq    = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq6_5ms_251125_182437';

Electrode_Type = 1; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% --- List of channels to plot ---
target_channels = [1:32]; 

plot_amps      = [10];      % amplitudes → one figure per amp
ras_win        = [-10 30];  % ms
bin_ms         = 1;
smooth_ms      = 3;
FS             = 30000;

%% Colors
pastel_colors = [
    0.4 0.6 0.85;   % pastel blue
    0.9 0.6 0.5;    % pastel salmon
    0.6 0.8 0.6;    % pastel green
    0.85 0.65 0.85; % pastel purple
    0.95 0.8 0.5;   % pastel orange
];
col_sim = [0.3 0.65 0.3];    % soft green

% MODIFIED: Brighter/clearer sequential colors
seq_colors = [
    0.85 0.10 0.10;  % Bright red
    0.90 0.60 0.10;  % Bright orange/gold
    0.60 0.20 0.20;  % Dark red (fallback)
];

%% ============================================================
% LOAD DATASETS
% ============================================================
% ... (Loading sections for Single, Sim, Seq remain unchanged) ...
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

% ---------- Load Simultaneous ----------
cd(folder_sim)
f = dir('*sp_xia_SSD.mat');
if isempty(f), error('No SSD file in Sim folder'); end
tmp = load(f(1).name);
if isfield(tmp, 'sp_corr'), D_sim.sp = tmp.sp_corr; 
elseif isfield(tmp, 'sp_SSD'), D_sim.sp = tmp.sp_SSD; 
else, D_sim.sp = tmp.sp_in; end
if isempty(dir('*.trig.dat')); cleanTrig_sabquick; end
D_sim.trig = loadTrig(0);
S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams; simN = S.simultaneous_stim; D_sim.E_MAP = S.E_MAP; D_sim.nTrials = S.n_Trials;
amps_all = cell2mat(Stim(2:end,16)); D_sim.trialAmps = amps_all(1:simN:end);
stimNames = Stim(2:end,1); [~, idx]  = ismember(stimNames, D_sim.E_MAP(2:end));
comb = zeros(D_sim.nTrials, simN);
for t = 1:D_sim.nTrials, rr = (t-1)*simN + (1:simN); v = unique(idx(rr)); v = v(v>0).'; comb(t,1:numel(v)) = v; end
[D_sim.uniqueComb,~,D_sim.combClass] = unique(comb,'rows');

% ---------- Load Sequential ----------
cd(folder_seq)
f = dir('*sp_xia_SSD.mat');
if isempty(f), error('No SSD file in Seq folder'); end
tmp = load(f(1).name);
if isfield(tmp, 'sp_corr'), D_seq.sp = tmp.sp_corr; 
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
        total_tiles = nSingleSets + 1 + nSeqSets + 2; 
        
        tl = tiledlayout(total_tiles, 1, 'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('Channel %d — %.0f µA',target_channel,amp_val),'FontSize',16);

        % Variable to track maximum firing rate across all PSTHs
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
        
        % ============ 2. SIMULTANEOUS RASTER =============
        ax = nexttile(tl); set(ax,'NextPlot','add'); box(ax,'off');
        stimCh = D_sim.uniqueComb(1, D_sim.uniqueComb(1,:)>0);
        S_ch = D_sim.sp{d(target_channel)};
        trial_ids = find(D_sim.trialAmps==amp_val);
        y = 0;
        for tr = trial_ids'
            t0 = D_sim.trig(tr)/FS*1000;
            tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
            for k=1:numel(tt), plot(ax,[tt(k) tt(k)],[y y+0.8], 'Color', col_sim, 'LineWidth',1.5); end
            y=y+1;
        end
        title(ax, sprintf('Simultaneous — Ch%s', num2str(stimCh)), 'FontSize',11);
        xline(ax,0,'r--','HandleVisibility','off');
        xlim(ax,ras_win); ylim(ax,[0 max(1,y)]);
        
        % ============ 3. SEQUENTIAL RASTERS =============
        for s = 1:nSeqSets
            ax = nexttile(tl); set(ax,'NextPlot','add'); box(ax,'off');
            stimCh = D_seq.uniqueComb(s, D_seq.uniqueComb(s,:)>0);
            S_ch = D_seq.sp{d(target_channel)};
            trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val);
            this_col = seq_colors(mod(s-1,size(seq_colors,1))+1, :);
            y=0;
            for tr = trial_ids'
                t0 = D_seq.trig(tr)/FS*1000;
                tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                for k=1:numel(tt), plot(ax,[tt(k) tt(k)],[y y+0.8], 'Color', this_col, 'LineWidth',1.5); end
                y=y+1;
            end
            title(ax, sprintf('Sequential Set %d — Order: %s', s, num2str(stimCh)), 'FontSize',11);
            xline(ax,0,'r--','HandleVisibility','off');
            xlim(ax,ras_win); ylim(ax,[0 max(1,y)]);
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
            % Update max rate
            max_psth_rate = max(max_psth_rate, max(rate_s));
        end    
        xline(ax_single_psth,0,'r--','HandleVisibility','off');
        xlim(ax_single_psth,ras_win); ylabel(ax_single_psth,'Rate (sp/s)');
        legend(ax_single_psth,'Box','off','Location','northeast');
        
        %% ============ 5. PSTH (Sim + All Seq Sets) ============
        ax_simseq_psth = nexttile(tl); box(ax_simseq_psth,'off'); hold(ax_simseq_psth,'on');
        title(ax_simseq_psth,'PSTH — Simultaneous vs Sequential');
        
        % --- Sim ---
        S_ch = D_sim.sp{d(target_channel)};
        trial_ids = find(D_sim.trialAmps==amp_val);
        if ~isempty(trial_ids)
            counts = zeros(1,length(edges)-1);
            for tr = trial_ids'
                t0 = D_sim.trig(tr)/FS*1000;
                tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                counts = counts + histcounts(tt,edges);
            end
            rate_sim = filter(g,1, counts/(numel(trial_ids)*bin_s));
            plot(ax_simseq_psth,ctrs,rate_sim,'Color',col_sim,'LineWidth',2,'DisplayName','Simultaneous');
            % Update max rate
            max_psth_rate = max(max_psth_rate, max(rate_sim));
        end
        
        % --- Seq (Loop over sets) ---
        for s = 1:nSeqSets
            S_ch = D_seq.sp{d(target_channel)};
            trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val);
            if isempty(trial_ids), continue; end
            counts = zeros(1,length(edges)-1);
            for tr = trial_ids'
                t0 = D_seq.trig(tr)/FS*1000;
                tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                counts = counts + histcounts(tt,edges);
            end
            rate_seq = filter(g,1, counts/(numel(trial_ids)*bin_s));
            this_col = seq_colors(mod(s-1,size(seq_colors,1))+1, :);
            plot(ax_simseq_psth,ctrs,rate_seq,'Color',this_col,'LineWidth',2, 'DisplayName',sprintf('Seq Set %d', s));
            % Update max rate
            max_psth_rate = max(max_psth_rate, max(rate_seq));
        end
        
        xline(ax_simseq_psth,0,'r--','HandleVisibility','off');
        xlim(ax_simseq_psth,ras_win); ylabel(ax_simseq_psth,'Rate (sp/s)'); xlabel(ax_simseq_psth,'Time (ms)');
        legend(ax_simseq_psth,'Box','off','Location','northeast');

        if max_psth_rate > 0
            ymax = max_psth_rate * 1.1; % Add 10% buffer
            ylim(ax_single_psth, [0 ymax]);
            ylim(ax_simseq_psth, [0 ymax]);
        end

    end  % END amplitude loop
end % END channel loop