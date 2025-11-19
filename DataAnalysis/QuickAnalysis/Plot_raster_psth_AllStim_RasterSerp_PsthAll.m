%% ============================================================
%   Unified Raster + PSTH Comparison (Single / Simultaneous / Sequential)
%   No local functions – fully expanded for script execution
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ===================== USER INPUTS ===================== %%
folder_single = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single1';
folder_sim    = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim1';
folder_seq    = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq1';

Electrode_Type = 1; % 0:single shank rigid; 1:single shank flex; 2:four shank flex
target_channel = 19;         % channel to plot
plot_amps      = [5];      % amplitudes → one figure per amp

ras_win        = [-10 50];  % ms
bin_ms         = 2;
smooth_ms      = 3;
FS             = 30000;

%% Soft pastel colors
pastel_colors = [
    0.4 0.6 0.85;   % pastel blue
    0.9 0.6 0.5;    % pastel salmon
    0.6 0.8 0.6;    % pastel green
    0.85 0.65 0.85; % pastel purple
    0.95 0.8 0.5;   % pastel orange
];

col_sim = [0.3 0.65 0.3];    % soft green
col_seq = [0.85 0.35 0.35];  % soft red

%% ============================================================
% LOAD DATASETS (EXPANDED – NO FUNCTIONS)
% ============================================================

% ---------- Load Single ----------
cd(folder_single)
load(dir('*sp_xia.mat').name,'sp_clipped'); D_single.sp = sp_clipped;
if isempty(dir('*.trig.dat')); cleanTrig_sabquick; end
D_single.trig = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams;
simN = S.simultaneous_stim;
D_single.E_MAP = S.E_MAP;
D_single.nTrials = S.n_Trials;
D_single.simN = simN;

amps_all = cell2mat(Stim(2:end,16));
D_single.trialAmps = amps_all(1:simN:end);

stimNames = Stim(2:end,1);
[~, idx]  = ismember(stimNames, D_single.E_MAP(2:end));
stimChPerTrial = cell(D_single.nTrials,1);
for t = 1:D_single.nTrials
    rr = (t-1)*simN + (1:simN);
    v = unique(idx(rr)); v = v(v>0).';
    stimChPerTrial{t} = v;
end
comb = zeros(D_single.nTrials, simN);
for t = 1:D_single.nTrials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end
[D_single.uniqueComb,~,D_single.combClass] = unique(comb,'rows');

pp_all = cell2mat(Stim(2:end,9));
D_single.pulseIdx = pp_all(1:simN:end);
D_single.PulsePeriods = unique(D_single.pulseIdx);

% ---------- Load Simultaneous ----------
cd(folder_sim)
load(dir('*sp_xia.mat').name,'sp_clipped'); D_sim.sp = sp_clipped;
if isempty(dir('*.trig.dat')); cleanTrig_sabquick; end
D_sim.trig = loadTrig(0);
S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams; simN = S.simultaneous_stim;
D_sim.E_MAP = S.E_MAP; D_sim.nTrials = S.n_Trials; D_sim.simN = simN;

amps_all = cell2mat(Stim(2:end,16));
D_sim.trialAmps = amps_all(1:simN:end);

stimNames = Stim(2:end,1);
[~, idx]  = ismember(stimNames, D_sim.E_MAP(2:end));
stimChPerTrial = cell(D_sim.nTrials,1);
for t = 1:D_sim.nTrials
    rr = (t-1)*simN + (1:simN);
    v = unique(idx(rr)); v = v(v>0).';
    stimChPerTrial{t} = v;
end
comb = zeros(D_sim.nTrials, simN);
for t = 1:D_sim.nTrials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end
[D_sim.uniqueComb,~,D_sim.combClass] = unique(comb,'rows');

pp_all = cell2mat(Stim(2:end,9));
D_sim.pulseIdx = pp_all(1:simN:end);
D_sim.PulsePeriods = unique(D_sim.pulseIdx);

% ---------- Load Sequential (FirstPulse) ----------
cd(folder_seq)
load(dir('*sp_xia_FirstPulse.mat').name,'sp_seq'); D_seq.sp = sp_seq;
if isempty(dir('*.trig.dat')); cleanTrig_sabquick; end
D_seq.trig = loadTrig(0);
S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams; simN = S.simultaneous_stim;
D_seq.E_MAP = S.E_MAP; D_seq.nTrials = S.n_Trials; D_seq.simN = simN;

amps_all = cell2mat(Stim(2:end,16));
D_seq.trialAmps = amps_all(1:simN:end);

stimNames = Stim(2:end,1);
[~, idx]  = ismember(stimNames, D_seq.E_MAP(2:end));
stimChPerTrial = cell(D_seq.nTrials,1);
for t = 1:D_seq.nTrials
    rr = (t-1)*simN + (1:simN);
    v = unique(idx(rr)); v = v(v>0).';
    stimChPerTrial{t} = v;
end
comb = zeros(D_seq.nTrials, simN);
for t = 1:D_seq.nTrials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end
[D_seq.uniqueComb,~,D_seq.combClass] = unique(comb,'rows');

pp_all = cell2mat(Stim(2:end,9));
D_seq.pulseIdx = pp_all(1:simN:end);
D_seq.PulsePeriods = unique(D_seq.pulseIdx);

%% PSTH KERNEL

edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g/sum(g);

d = Depth_s(Electrode_Type);

%% MAIN LOOP: ONE FIGURE PER AMPLITUDE

for amp_val = plot_amps

    figure('Color','w','Position',[250 0 750 1300]);
    tl = tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('Channel %d — %.0f µA',target_channel,amp_val),'FontSize',16);

    % ============ Helper inline raster (expanded) =============
    % RASTER for single sets
    for s = 1:size(D_single.uniqueComb,1)
    
        ax = nexttile(tl);
        set(ax,'NextPlot','add'); box(ax,'off');
    
        % Stim channels in this set
        stimCh = D_single.uniqueComb(s, D_single.uniqueComb(s,:)>0);
        labelStr = sprintf('Separate Ch%s', num2str(stimCh));
    
        S_ch = D_single.sp{d(target_channel)};
        trial_ids = find(D_single.combClass==s & D_single.trialAmps==amp_val);
    
        y = 0;
        for tr = trial_ids'
            t0 = D_single.trig(tr)/FS*1000;
            tt = S_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
    
            for k=1:numel(tt)
                plot(ax,[tt(k) tt(k)],[y y+1], ...
                     'Color', pastel_colors(s,:), 'LineWidth',1.5);
            end
            y=y+1;
        end
    
        title(ax,labelStr,'FontSize',11,'Interpreter','none');
        xline(ax,0,'r--','HandleVisibility','off');
        xlim(ax,ras_win); ylim(ax,[0 y]);
    end

    % ============ SIMULTANEOUS RASTER =============
    ax = nexttile(tl); set(ax,'NextPlot','add'); box(ax,'off');
    stimCh = D_sim.uniqueComb(1, D_sim.uniqueComb(1,:)>0);
    title(ax, sprintf('Simultaneous — Ch%s', num2str(stimCh)), 'FontSize',11);
    
    S_ch = D_sim.sp{d(target_channel)};
    trial_ids = find(D_sim.trialAmps==amp_val);
    
    y = 0;
    for tr = trial_ids'
        t0 = D_sim.trig(tr)/FS*1000;
        tt = S_ch(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
    
        for k=1:numel(tt)
            plot(ax,[tt(k) tt(k)],[y y+1], ...
                 'Color', col_sim, 'LineWidth',1.5);
        end
        y=y+1;
    end
    
    xline(ax,0,'r--','HandleVisibility','off');
    xlim(ax,ras_win); ylim(ax,[0 y]);

    % ============ SEQUENTIAL RASTER =============
    ax = nexttile(tl); set(ax,'NextPlot','add'); box(ax,'off');
    stimCh = D_seq.uniqueComb(1, D_seq.uniqueComb(1,:)>0);
    title(ax, sprintf('Sequential — Ch%s', num2str(stimCh)), 'FontSize',11);
    
    S_ch = D_seq.sp{d(target_channel)};
    trial_ids = find(D_seq.trialAmps==amp_val);
    
    y=0;
    for tr = trial_ids'
        t0 = D_seq.trig(tr)/FS*1000;
        tt = S_ch(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
    
        for k=1:numel(tt)
            plot(ax,[tt(k) tt(k)],[y y+1], ...
                 'Color', col_seq, 'LineWidth',1.5);
        end
        y=y+1;
    end
    
    xline(ax,0,'r--','HandleVisibility','off');
    xlim(ax,ras_win); ylim(ax,[0 y]);

    %% ====== Combined PSTH: Single (dashed), Sim (solid), Seq (solid) ======
    ax = nexttile(tl); hold(ax,'on'); box(ax,'off');
    title(ax,'PSTH','FontSize',12);
    
    % ----- PSTH for Single Sets (dashed) -----
    for s = 1:size(D_single.uniqueComb,1)
        
        % spike data
        S_ch = D_single.sp{d(target_channel)};
        trial_ids = find(D_single.combClass==s & D_single.trialAmps==amp_val);
    
        % accumulate counts
        counts = zeros(1,length(edges)-1);
        for tr = trial_ids'
            t0 = D_single.trig(tr)/FS*1000;
            tt = S_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
            counts = counts + histcounts(tt,edges);
        end 
        % convert to firing rate
        bin_s = bin_ms/1000;
        rate = counts / (numel(trial_ids) * bin_s);
        rate_s = filter(g,1,rate);
    
        plot(ax,ctrs,rate_s, ...
            'Color', pastel_colors(s,:), ...
            'LineWidth', 1.8, ...
            'LineStyle','--', ...
            'DisplayName', sprintf('Separate Ch%s', num2str(D_single.uniqueComb(s))));
    end
    
    
    % ----- PSTH: Simultaneous (solid) -----
    S_ch = D_sim.sp{d(target_channel)};
    trial_ids = find(D_sim.trialAmps==amp_val);
    
    counts = zeros(1,length(edges)-1);
    for tr = trial_ids'
        t0 = D_sim.trig(tr)/FS*1000;
        tt = S_ch(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
        counts = counts + histcounts(tt,edges);
    end
    
    % rate_sim = counts / (max(1,numel(trial_ids)) * bin_s);
    rate_sim = counts / (numel(trial_ids) * bin_s);
    rate_sim = filter(g,1,rate_sim);
    
    plot(ax,ctrs,rate_sim, ...
        'Color', col_sim, ...
        'LineWidth', 2.2, ...
        'LineStyle','-', ...
        'DisplayName','Sim');
    
    %% ----- PSTH: Sequential (solid) -----
    S_ch = D_seq.sp{d(target_channel)};
    trial_ids = find(D_seq.trialAmps==amp_val);
    
    counts = zeros(1,length(edges)-1);
    for tr = trial_ids'
        t0 = D_seq.trig(tr)/FS*1000;
        tt = S_ch(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
        counts = counts + histcounts(tt,edges);
    end
    
    % rate_seq = counts / (max(1,numel(trial_ids)) * bin_s);
    rate_seq = counts / (numel(trial_ids) * bin_s);
    rate_seq = filter(g,1,rate_seq);
    
    plot(ax,ctrs,rate_seq, ...
        'Color', col_seq, ...
        'LineWidth', 2.2, ...
        'LineStyle','-', ...
        'DisplayName','Seq');
    
    
    %% ----- Final formatting -----
    xline(ax,0,'r--','HandleVisibility','off');
    xlim(ax,ras_win);
    ylabel(ax,'Rate (sp/s)');
    xlabel(ax,'Time (ms)');
    legend(ax,'Box','off','Location','northeast');

end  % END amplitude loop