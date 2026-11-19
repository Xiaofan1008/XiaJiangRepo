%% ============================================================
%   SINGLE vs SEQUENTIAL Raster + PSTH 
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ===================== USER INPUTS ===================== %%
folder_single = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Single1_251128_113846';
folder_seq    = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim1_251128_114443';

Electrode_Type = 2;
target_channel = 10;
plot_amps      = [10];

ras_win   = [-10 50];
bin_ms    = 1;
smooth_ms = 2;
FS        = 30000;

cmap_ptd = lines(20);   % Enough colours for many PTDs


%% ============================================================
%                         LOAD SINGLE
% ============================================================
cd(folder_single)
load(dir('*sp_xia.mat').name,'sp_clipped');
D_single.sp = sp_clipped;

if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
D_single.trig = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams;
simN = S.simultaneous_stim;

D_single.E_MAP  = S.E_MAP;
D_single.nTrials = S.n_Trials;

amps_all = cell2mat(Stim(2:end,16));
D_single.trialAmps = amps_all(1:simN:end);

% unordered set decoding for single
stimNames = Stim(2:end,1);
[~, idx] = ismember(stimNames, D_single.E_MAP(2:end));

stimChPerTrial = cell(D_single.nTrials,1);
for t=1:D_single.nTrials
    v = unique(idx((t-1)*simN + (1:simN))); 
    v = v(v>0);
    stimChPerTrial{t} = v;
end

comb = zeros(D_single.nTrials,simN);
for t=1:D_single.nTrials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end

[D_single.uniqueComb,~,D_single.combClass] = unique(comb,'rows');
D_single.nSets = size(D_single.uniqueComb,1);


%% ============================================================
%                         LOAD SEQUENTIAL
% ============================================================
cd(folder_seq)
load(dir('*sp_xia_FirstPulse.mat').name,'sp_seq');
D_seq.sp = sp_seq;

if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
D_seq.trig = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams;
simN = S.simultaneous_stim;

D_seq.E_MAP  = S.E_MAP;
D_seq.nTrials = S.n_Trials;

amps_all = cell2mat(Stim(2:end,16));
D_seq.trialAmps = amps_all(1:simN:end);

stimNames = Stim(2:end,1);
[~, idx] = ismember(stimNames, D_seq.E_MAP(2:end));

stimChPerTrial = cell(D_seq.nTrials,1);
for t=1:D_seq.nTrials
    v = idx((t-1)*simN+(1:simN));
    v = v(v>0);
    stimChPerTrial{t} = v(:).';
end

comb = zeros(D_seq.nTrials,simN);
for t=1:D_seq.nTrials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end

[D_seq.uniqueComb,~,D_seq.combClass] = unique(comb,'rows','stable');
D_seq.nSets = size(D_seq.uniqueComb,1);

% PTDs
ptd_all_us = cell2mat(Stim(2:end,6));
PTD_us  = ptd_all_us(2:simN:end);
PTD_vals = unique(PTD_us);

%% ============================================================
%                     PSTH KERNEL
% ============================================================
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

d = Depth_s(Electrode_Type);


%% ============================================================
%                    MAIN LOOP
% ============================================================

for amp_val = plot_amps
for s = 1:D_seq.nSets

    stimCh = D_seq.uniqueComb(s, D_seq.uniqueComb(s,:)>0);

    % number of subplots:
    nRaster = 1 + numel(PTD_vals);  % single raster + one per PTD
    nPSTH   = 2;                    % single PSTH + seq PSTH
    totalRows = nRaster + nPSTH;

    figure('Color','w','Position',[200 50 950 1400])
    tl = tiledlayout(totalRows,1,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Ch %d  |  %.0f µA  |  Set %s', ...
        target_channel, amp_val, strjoin(string(stimCh),'→')));


    %% ---------------------------------------------------
    %                RASTERS — ALL AT TOP
    % ---------------------------------------------------

    %% --- SINGLE RASTER ---
    ax = nexttile(tl); hold(ax,'on'); box(ax,'off');
    title(ax,'Single — Raster');

    S_ch_single = D_single.sp{d(target_channel)};
    trial_ids = find(D_single.combClass==s & D_single.trialAmps==amp_val);

    y=0;
    for tr = trial_ids'
        t0 = D_single.trig(tr)/FS*1000;
        tt = S_ch_single(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
        scatter(ax,tt,y*ones(size(tt)),12,'k','filled');
        y=y+1;
    end
    xline(ax,0,'r--'); ylim(ax,[0 max(1,y)]); xlim(ax,ras_win);


    %% --- SEQUENTIAL RASTERS FOR EACH PTD ---
    S_ch_seq = D_seq.sp{d(target_channel)};

    for p = 1:numel(PTD_vals)
        ptd = PTD_vals(p);

        ax = nexttile(tl); hold(ax,'on'); box(ax,'off');
        title(ax, sprintf('Seq — Raster | PTD %d µs', ptd));

        trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val & PTD_us==ptd);

        y=0;
        for tr = trial_ids'
            t0 = D_seq.trig(tr)/FS*1000;
            tt = S_ch_seq(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;

            scatter(ax,tt,y*ones(size(tt)),12,cmap_ptd(p,:),'filled');
            y=y+1;
        end

        xline(ax,0,'r--');
        ylim(ax,[0 max(1,y)]);
        xlim(ax,ras_win);
    end


    %% ---------------------------------------------------
    %             PSTH — ALL AT BOTTOM
    % ---------------------------------------------------

    %% --- SINGLE PSTH ---
    ax = nexttile(tl); hold(ax,'on'); box(ax,'off');
    title(ax,'Single — PSTH');

    counts = zeros(1,length(edges)-1);
    for tr = find(D_single.combClass==s & D_single.trialAmps==amp_val)'
        t0 = D_single.trig(tr)/FS*1000;
        tt = S_ch_single(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
        counts = counts + histcounts(tt,edges);
    end

    rate = filter(g,1,counts/(numel(trial_ids)*bin_s));
    plot(ax, ctrs, rate,'k','LineWidth',2);
    xline(ax,0,'r--'); xlim(ax,ras_win); ylabel(ax,'sp/s');


    %% --- SEQ PSTH (overlaid PTDs) ---
    ax = nexttile(tl); hold(ax,'on'); box(ax,'off');
    title(ax,'Sequential — PSTH (all PTDs)');

    for p = 1:numel(PTD_vals)
        ptd = PTD_vals(p);
        trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val & PTD_us==ptd);

        if isempty(trial_ids), continue; end

        counts = zeros(1,length(edges)-1);
        for tr = trial_ids'
            t0 = D_seq.trig(tr)/FS*1000;
            tt = S_ch_seq(:,1);
            tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
            counts = counts + histcounts(tt,edges);
        end

        rate = filter(g,1,counts/(numel(trial_ids)*bin_s));

        plot(ax, ctrs, rate,'LineWidth',2,'Color',cmap_ptd(p,:), ...
            'DisplayName',sprintf('%d µs',ptd));
    end

    xline(ax,0,'r--');
    legend(ax,'Box','off');
    xlim(ax,ras_win); xlabel(ax,'Time (ms)'); ylabel(ax,'sp/s');

end
end