%% =============================================================
%   Raster + PSTH (All Channels) per Set × Amp × PTD
%   - Uses same spike loading style as RespondingChn_identify.m
%   - Highlights responding channels (from *_RespondingChannels.mat)
%   - Marks bad channels (from *.BadChannels.mat) but does NOT remove them
%   - Does NOT remove bad trials; they are only noted in the title
% =============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ========================

data_folder      = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq5_5ms';

Electrode_Type   = 1;          % 0 rigid, 1 single-shank flex, 2 four-shank flex
raster_chn_start = 1;          % Depth_s index
raster_chn_end   = 32;

% ---- plotting windows ----
ras_win       = [-50 80];      % ms, time relative to first pulse
bin_ms_raster = 1;             % ms, PSTH bin size
smooth_ms     =5;             % ms, Gaussian smoothing width

% Which amplitudes to plot (µA). If empty → plot ALL amplitudes.
Plot_Amps = [];     % e.g. [4 6 8];  [] = all

% Fixed figure size (all figures the same)
fig_position = [50 50 1600 900];

%% ===================== INITIAL SETUP =========================

if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end
cd(data_folder);
fprintf('Raster+PSTH plotting in folder:\n%s\n\n', data_folder);

% ---- basename (same rule as your other scripts) ----
parts       = split(data_folder, filesep);
last_folder = parts{end};
u           = strfind(last_folder, '_');
if numel(u) >= 4
    base_name = last_folder(1:u(end-1)-1);   % e.g. 'Xia_Exp1_Seq1'
else
    base_name = last_folder;
end

FS = 30000;

%% ===================== LOAD SPIKES (as in RespondingChn) =====
ssd_file  = [base_name '.sp_xia_SSD.mat'];
base_file = [base_name '.sp_xia.mat'];

if isfile(ssd_file)
    S = load(ssd_file);
    if     isfield(S,'sp_corr'), sp = S.sp_corr;
    elseif isfield(S,'sp_SSD'),  sp = S.sp_SSD;
    else, error('SSD file missing usable variable (sp_corr / sp_SSD).');
    end
elseif isfile(base_file)
    S = load(base_file);
    if     isfield(S,'sp_clipped'), sp = S.sp_clipped;
    elseif isfield(S,'sp'),         sp = S.sp;
    else, error('Base spike file missing usable variable (sp_clipped / sp).');
    end
else
    error('No spike file %s or %s found.', ssd_file, base_file);
end

nCh = numel(sp);
fprintf('  -> loaded %d spike channels.\n', nCh);

%% ===================== LOAD BadChannels & BadTrials ==========

BadCh_perSet = {};
BadTrialsPerCh = {};

badch_file = [base_name '.BadChannels.mat'];
if isfile(badch_file)
    tmp = load(badch_file);
    if isfield(tmp,'BadCh_perSet')
        BadCh_perSet = tmp.BadCh_perSet;   % cell{set_id} = [depth_idx ...]
        fprintf('Loaded BadChannels from %s\n', badch_file);
    else
        warning('BadChannels.mat found but BadCh_perSet missing – ignored.');
    end
else
    fprintf('No BadChannels.mat found – no channels marked as bad.\n');
end

badtr_file = [base_name '.BadTrials.mat'];
if isfile(badtr_file)
    tmp = load(badtr_file);
    if isfield(tmp,'BadTrials')
        BadTrialsPerCh = tmp.BadTrials;    % cell{depth_idx} = [trial list]
        fprintf('Loaded BadTrials from %s\n', badtr_file);
    else
        warning('BadTrials.mat found but BadTrials missing – ignored.');
    end
else
    fprintf('No BadTrials.mat found – no bad trials annotated.\n');
end

%% ===================== LOAD Responding Channels ==============

Resp = [];
resp_file = [base_name '_RespondingChannels.mat'];
hasResp = false;
if isfile(resp_file)
    tmp = load(resp_file);
    if isfield(tmp,'Responding')
        Resp    = tmp.Responding;
        hasResp = true;
        fprintf('Loaded responding-channel data from %s\n', resp_file);
    else
        warning('RespondingChannels.mat found but Responding struct missing.');
    end
else
    fprintf('No RespondingChannels file – channels will not be highlighted as RESP.\n');
end

%% ===================== LOAD TRIGGERS =========================

if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);   % in samples
nTrig = numel(trig);

%% ===================== LOAD StimParams & DECODE ==============

fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat file found.');
S = load(fileDIR(1).name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');

StimParams = S.StimParams;
sim_stim   = S.simultaneous_stim;
E_MAP      = S.E_MAP;
n_Trials   = S.n_Trials;

% ---- Amplitudes (one per trial) ----
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:sim_stim:end);
[Amps,~,ampIdx] = unique(trialAmps);
Amps(Amps == -1) = 0;
nAMP = numel(Amps);

% Reduce Plot_Amps if user left it empty
if isempty(Plot_Amps)
    Plot_Amps = Amps(:).';
end

% ---- PTD decoding ----
if sim_stim > 1
    % column 6, starting at third row for pulse 2
    PTD_all = cell2mat(StimParams(3:sim_stim:end,6));   % µs
else
    PTD_all = zeros(n_Trials,1);    % single pulse: PTD = 0
end
[PTDs,~,ptdIdx] = unique(PTD_all);
nPTD = numel(PTDs);

% ---- Stimulation sets (order-sensitive) ----
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));

stimSeq = zeros(n_Trials, sim_stim);
for t = 1:n_Trials
    rr = (t-1)*sim_stim + (1:sim_stim);
    v  = idx_all(rr);
    v  = v(v>0);
    stimSeq(t,1:numel(v)) = v;
end
[uniqueComb,~,combClass] = unique(stimSeq, 'rows','stable');
nSets = size(uniqueComb,1);

fprintf('\nDetected %d stimulation sets (order-sensitive):\n', nSets);
for k = 1:nSets
    v = uniqueComb(k, uniqueComb(k,:)>0);
    fprintf('  Set %d: %s\n', k, mat2str(v));
end

%% ===================== ELECTRODE MAP =========================

d = Depth_s(Electrode_Type);   % depth index -> Intan channel index
depth_range = raster_chn_start : min(raster_chn_end, numel(d));
nChPlot = numel(depth_range);

%% ===================== PSTH KERNEL ===========================

edges = ras_win(1) : bin_ms_raster : ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;

g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
g = g / sum(g);

%% ===================== MAIN CONDITION LOOPS ==================
for si = 1:nSets
    
    stimIdx = uniqueComb(si, uniqueComb(si,:)>0);
    setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, ...
                  'UniformOutput', false), ' + ');
    
    % Bad channels for THIS set (depth indices)
    if ~isempty(BadCh_perSet) && si <= numel(BadCh_perSet)
        badDepthThisSet = BadCh_perSet{si};
        if isempty(badDepthThisSet), badDepthThisSet = []; end
    else
        badDepthThisSet = [];
    end
    
    for aVal = Plot_Amps
        
        ai = find(abs(Amps - aVal) < 1e-6, 1);
        if isempty(ai)
            fprintf('Set %d | %.1f µA: amplitude not present, skipping.\n', si, aVal);
            continue;
        end
        
        for pi = 1:nPTD
            
            trials_this = find(combClass==si & ampIdx==ai & ptdIdx==pi);
            if isempty(trials_this)
                continue;
            end
            
            PTD_us = PTDs(pi);
            PTD_ms = PTD_us/1000;
            
            % ---- create figure for this (set,amp,PTD) ----
            figTitle = sprintf('Set %d (%s) | Amp %.1f µA | PTD %.1f ms | nTrials=%d | Sequential', ...
                               si, setLabel, aVal, PTD_ms, numel(trials_this));
            figure('Color','w','Name',figTitle,'Position',fig_position);
            tiledlayout('flow','TileSpacing','compact','Padding','compact');
            sgtitle(figTitle, 'FontSize', 14, 'FontWeight','bold');
            
            % ---- channel loop ----
            for idxDepth = 1:nChPlot
                
                ich = depth_range(idxDepth);   % depth index
                ch  = d(ich);                  % Intan index
                
                if ch < 1 || ch > nCh || isempty(sp{ch})
                    nexttile; axis off;
                    continue;
                end
                
                % spikes for this channel
                sp_times = sp{ch}(:,1);   % assumed in ms already
                
                % bad trials for this *channel* (for info only)
                if ~isempty(BadTrialsPerCh) && ich <= numel(BadTrialsPerCh)
                    badTr_ch = BadTrialsPerCh{ich};
                    if isempty(badTr_ch), badTr_ch = []; end
                else
                    badTr_ch = [];
                end
                
                % collect spikes and PSTH
                allTrialSpikes = cell(numel(trials_this),1);
                counts         = zeros(1, numel(edges)-1);
                
                for ti = 1:numel(trials_this)
                    tr = trials_this(ti);
                    t0 = trig(tr)/FS*1000;  % ms
                    tt = sp_times;
                    tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                    allTrialSpikes{ti} = tt;
                    counts = counts + histcounts(tt, edges);
                end
                
                if ~any(counts)
                    rate_s = zeros(size(ctrs));
                else
                    rate   = counts / (numel(trials_this)*bin_s);
                    rate_s = filter(g,1,rate);
                end
                maxRate   = max(rate_s);
                yMaxPSTH  = max(50, ceil(maxRate*1.1/10)*10);
                nTr       = numel(trials_this);
                
                % ================= PLOT SUBPLOT ==================
                ax = nexttile; hold(ax,'on');
                
                % ---- highlight responding / bad channel ----
                respTag = '';
                isResp = false;
                
                if hasResp && si <= numel(Resp.set) && ...
                        ai <= numel(Resp.set(si).amp) && ...
                        pi <= numel(Resp.set(si).amp(ai).ptd) && ...
                        ich <= numel(Resp.set(si).amp(ai).ptd(pi).channel)
                
                    R = Resp.set(si).amp(ai).ptd(pi).channel(ich);
                    if isfield(R,'is_responsive') && R.is_responsive
                        isResp = true;
                        respTag = 'RESP';
                    end
                end
                
                isBadCh = ismember(ich, badDepthThisSet);
                
                % ----- COLOR LOGIC -----
                if isResp && ~isBadCh
                    % responding channel only → light red background
                    set(ax,'Color',[1.0 0.88 0.88]);    % pale pink
                elseif isResp && isBadCh
                    % both responding and bad → orange warning color
                    set(ax,'Color',[1.0 0.90 0.75]);
                elseif isBadCh
                    % bad only → existing light gray
                    set(ax,'Color',[0.95 0.95 0.95]);
                end
                
                % Axis color for responding channels
                if isResp
                    ax.XColor = [0.6 0 0];
                    ax.YColor = [0.6 0 0];
                    ax.LineWidth = 1.4;
                end
                
                % ---- left y-axis: PSTH ----
                yyaxis(ax,'left');
                if any(rate_s)
                    plot(ax, ctrs, rate_s, 'LineWidth', 1.4);
                end
                xlim(ax, ras_win);
                ylim(ax, [0 yMaxPSTH]);
                ylabel(ax,'Rate (sp/s)');
                
                % ---- right y-axis: raster ----
                yyaxis(ax,'right');
                for ti = 1:nTr
                    tt = allTrialSpikes{ti};
                    if isempty(tt), continue; end
                    plot(ax, tt, ti*ones(size(tt)), '.', ...
                         'Color',[0 0 0], 'MarkerSize',4);
                end
                ylim(ax,[0 nTr+1]);
                set(ax,'YTick',[]);
                
                xline(ax,0,'r--','LineWidth',1);
                xlim(ax, ras_win);
                
                % ---- title & labels ----
                chLabel = sprintf('Ch %d', ich);
                if isBadCh
                    chLabel = [chLabel ' BAD'];
                end
                chLabel = [chLabel respTag];
                
                % title(ax, chLabel, 'FontSize',10,'Interpreter','none');
                ht = title(ax, chLabel, 'FontSize',11, 'FontWeight','bold', 'Interpreter','none');
                
                % Color the title for responding channels
                if exist('isResp','var') && isResp
                    set(ht, 'Color', [0.7 0 0]);   % dark red title for RESP channels
                end
                
                if idxDepth > nChPlot - ceil(sqrt(nChPlot))
                    xlabel(ax,'Time (ms)');
                end
                
                % add small text with count of bad trials for this channel
                if ~isempty(badTr_ch)
                    text(ax, ras_win(1)+1, nTr, ...
                        sprintf('badTr: %d', numel(intersect(trials_this,badTr_ch))), ...
                        'FontSize',7,'Color',[0.3 0.3 0.3], ...
                        'HorizontalAlignment','left','VerticalAlignment','top');
                end
            end
            
        end % PTD
    end % amplitude
end % set