%% ============================================================
%   Unified Raster + PSTH Comparison (Single / Simultaneous / Sequential)
%   REFINED VERSION (Xia 2025-11)
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ===================== USER INPUTS ===================== %%
folder_single = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Single1_251125_110714';
folder_sim    = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq    = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

Electrode_Type = 1;
target_channel = 16;
plot_amps      = [10];

ras_win        = [-10 50];  % ms
bin_ms         = 1;
smooth_ms      = 2;
FS             = 30000;

%% Colors
pastel_colors = [
    0.4 0.6 0.85;
    0.9 0.6 0.5;
    0.6 0.8 0.6;
    0.85 0.65 0.85;
    0.95 0.8 0.5;
];

sim_colors = pastel_colors(3,:);        % SIMULTANEOUS (unique)
seq_colors = pastel_colors(4:5,:);       % SEQUENTIAL (unique)

%% ============================================================
% LOAD DATA
% (I will not comment to save space — logic unchanged)
% ============================================================

%% ---------- SINGLE ----------
cd(folder_single)
load(dir('*sp_xia.mat').name,'sp_clipped');
D_single.sp = sp_clipped;
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
D_single.trig = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams; simN = S.simultaneous_stim;
D_single.E_MAP = S.E_MAP; D_single.nTrials = S.n_Trials; D_single.simN = simN;

amps_all = cell2mat(Stim(2:end,16));
D_single.trialAmps = amps_all(1:simN:end);

stimNames = Stim(2:end,1);
[~, idx] = ismember(stimNames, D_single.E_MAP(2:end));

% Unordered (old ST behaviour)
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
D_single.nSets = size(D_single.uniqueComb,1);

%% ---------- SIMULTANEOUS ----------
cd(folder_sim)
load(dir('*sp_xia.mat').name,'sp_clipped'); 
D_sim.sp = sp_clipped;
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
D_sim.trig = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams; simN = S.simultaneous_stim;
D_sim.E_MAP=S.E_MAP; D_sim.nTrials=S.n_Trials; D_sim.simN=simN;

amps_all = cell2mat(Stim(2:end,16));
D_sim.trialAmps = amps_all(1:simN:end);
stimNames = Stim(2:end,1);
[~, idx] = ismember(stimNames, D_sim.E_MAP(2:end));

% ORDERED (Sim)
stimChPerTrial = cell(D_sim.nTrials,1);
for t = 1:D_sim.nTrials
    rr = (t-1)*simN + (1:simN);
    v = idx(rr); v = v(v>0);
    stimChPerTrial{t} = v(:).';
end

comb = zeros(D_sim.nTrials, simN);
for t = 1:D_sim.nTrials, 
    v=stimChPerTrial{t}; comb(t,1:numel(v))=v; 
end

[D_sim.uniqueComb,~,D_sim.combClass] = unique(comb,'rows','stable');
D_sim.nSets = size(D_sim.uniqueComb,1);

%% ---------- SEQUENTIAL ----------
cd(folder_seq)
load(dir('*sp_xia_FirstPulse.mat').name,'sp_seq'); 
D_seq.sp = sp_seq;
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
D_seq.trig = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams; simN=S.simultaneous_stim;
D_seq.E_MAP=S.E_MAP; D_seq.nTrials=S.n_Trials; D_seq.simN=simN;

amps_all = cell2mat(Stim(2:end,16));
D_seq.trialAmps = amps_all(1:simN:end);
stimNames = Stim(2:end,1);
[~, idx] = ismember(stimNames, D_seq.E_MAP(2:end));

% ORDERED (Seq)
stimChPerTrial = cell(D_seq.nTrials,1);
for t=1:D_seq.nTrials
    rr=(t-1)*simN+(1:simN);
    v=idx(rr); v=v(v>0);
    stimChPerTrial{t}=v(:).';
end

comb=zeros(D_seq.nTrials,simN);
for t=1:D_seq.nTrials
    v=stimChPerTrial{t}; comb(t,1:numel(v))=v;
end

[D_seq.uniqueComb,~,D_seq.combClass] = unique(comb,'rows','stable');
D_seq.nSets = size(D_seq.uniqueComb,1);

%% ===================== PSTH KERNEL =====================
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;
g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2); g=g/sum(g);

%% Depth → channel index
d = Depth_s(Electrode_Type);

%% ===================== MAIN LOOP =====================
for amp_val = plot_amps

    nRows_raster = D_single.nSets + D_sim.nSets + D_seq.nSets;
    nRows_total  = nRows_raster + 2;

    figure('Color','w','Position',[300 50 900 1400]);
    tl = tiledlayout(nRows_total,1,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('Channel %d — %.0f µA',target_channel,amp_val),'FontSize',16);

    %% -----------------------------------------------------------
    %                    R A S T E R   P L O T S
    % -----------------------------------------------------------

    raster_marker = 12;     % larger, clear
    raster_lineH  = 0.7;

    %% -------- SINGLE SETS --------
    S_ch = D_single.sp{d(target_channel)};
    for s = 1:D_single.nSets

        ax = nexttile(tl);
        hold(ax,'on'); box(ax,'off');

        stimCh = D_single.uniqueComb(s, D_single.uniqueComb(s,:)>0);
        title(ax, sprintf('Single — Ch %s', strjoin(string(stimCh),',')));

        trial_ids = find(D_single.combClass==s & D_single.trialAmps==amp_val);

        y=0;
        for tr = trial_ids'
            t0 = D_single.trig(tr)/FS*1000;
            tt = S_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;

            scatter(ax, tt, y*ones(size(tt)), raster_marker, ...
                pastel_colors(min(s,end),:), 'filled');

            y=y+1;
        end

        xline(ax,0,'r--','LineWidth',1);
        ylim(ax,[0 max(1,y)]);
        xlim(ax,ras_win);
        ylabel(ax,'Trials');
    end

    %% -------- SIMULTANEOUS SETS --------
    S_ch = D_sim.sp{d(target_channel)};
    for s = 1:D_sim.nSets

        ax = nexttile(tl);
        hold(ax,'on'); box(ax,'off');

        stimCh = D_sim.uniqueComb(s, D_sim.uniqueComb(s,:)>0);
        title(ax, sprintf('Sim — Ch %s', strjoin(string(stimCh),'→')));

        trial_ids = find(D_sim.combClass==s & D_sim.trialAmps==amp_val);

        y=0;
        for tr = trial_ids'
            t0 = D_sim.trig(tr)/FS*1000;
            tt = S_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;

            scatter(ax, tt, y*ones(size(tt)), raster_marker, ...
                sim_colors(s,:), 'filled');

            y=y+1;
        end

        xline(ax,0,'r--','LineWidth',1);
        ylim(ax,[0 max(1,y)]);
        xlim(ax,ras_win);
        ylabel(ax,'Trials');
    end

    %% -------- SEQUENTIAL SETS --------
    S_ch = D_seq.sp{d(target_channel)};
    for s = 1:D_seq.nSets

        ax = nexttile(tl);
        hold(ax,'on'); box(ax,'off');

        stimCh = D_seq.uniqueComb(s, D_seq.uniqueComb(s,:)>0);
        title(ax, sprintf('Seq — Ch %s', strjoin(string(stimCh),'→')));

        trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val);
        y=0;

        for tr = trial_ids'
            t0 = D_seq.trig(tr)/FS*1000;
            tt = S_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;

            scatter(ax, tt, y*ones(size(tt)), raster_marker, ...
                seq_colors(s,:), 'filled');

            y=y+1;
        end

        xline(ax,0,'r--','LineWidth',1);
        ylim(ax,[0 max(1,y)]);
        xlim(ax,ras_win);
        ylabel(ax,'Trials');
    end

    %% -----------------------------------------------------------
    %                           P S T H
    % -----------------------------------------------------------

    %% ----- PSTH: SINGLE -----
    ax = nexttile(tl); hold(ax,'on'); box(ax,'off');
    title(ax,'PSTH — Single Sets');

    for s = 1:D_single.nSets
        S_ch = D_single.sp{d(target_channel)};
        trial_ids = find(D_single.combClass==s & D_single.trialAmps==amp_val);
        if isempty(trial_ids), continue; end

        counts = zeros(1,length(edges)-1);
        for tr = trial_ids'
            t0 = D_single.trig(tr)/FS*1000;
            tt = S_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
            counts = counts + histcounts(tt,edges);
        end

        rate = filter(g,1, counts/(numel(trial_ids)*bin_s));
        stimCh = D_single.uniqueComb(s, D_single.uniqueComb(s,:)>0);

        plot(ax, ctrs, rate, 'LineWidth',2, ...
            'Color', pastel_colors(min(s,end),:), ...
            'DisplayName', sprintf('Ch %s', strjoin(string(stimCh),',')));
    end

    xline(ax,0,'r--','LineWidth',1,'HandleVisibility','off');
    legend(ax,'Box','off');
    ylabel(ax,'Rate (sp/s)');
    xlim(ax,ras_win);

    %% ----- PSTH: SIM + SEQ -----
    ax = nexttile(tl); hold(ax,'on'); box(ax,'off');
    title(ax,'PSTH — Simultaneous + Sequential');

    % SIM
    for s = 1:D_sim.nSets
        S_ch = D_sim.sp{d(target_channel)};
        trial_ids = find(D_sim.combClass==s & D_sim.trialAmps==amp_val);
        if isempty(trial_ids), continue; end

        counts = zeros(1,length(edges)-1);
        for tr = trial_ids'
            t0 = D_sim.trig(tr)/FS*1000;
            tt = S_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
            counts = counts + histcounts(tt,edges);
        end

        rate = filter(g,1, counts/(numel(trial_ids)*bin_s));
        stimCh = D_sim.uniqueComb(s, D_sim.uniqueComb(s,:)>0);

        plot(ax, ctrs, rate, 'LineWidth',2, 'Color', sim_colors(s,:), ...
            'DisplayName', sprintf('Sim %s', strjoin(string(stimCh),'→')));
    end

    % SEQ
    for s = 1:D_seq.nSets
        S_ch = D_seq.sp{d(target_channel)};
        trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val);
        if isempty(trial_ids), continue; end

        counts = zeros(1,length(edges)-1);
        for tr = trial_ids'
            t0 = D_seq.trig(tr)/FS*1000;
            tt = S_ch(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
            counts = counts + histcounts(tt,edges);
        end

        rate = filter(g,1, counts/(numel(trial_ids)*bin_s));
        stimCh = D_seq.uniqueComb(s, D_seq.uniqueComb(s,:)>0);

        plot(ax, ctrs, rate, 'LineWidth',2, 'Color', seq_colors(s,:), ...
            'DisplayName', sprintf('Seq %s', strjoin(string(stimCh),'→')));
    end

    xline(ax,0,'r--','LineWidth',1,'HandleVisibility','off');
    ylabel(ax,'Rate (sp/s)');
    xlabel(ax,'Time (ms)');
    legend(ax,'Box','off');
    xlim(ax,ras_win);

end