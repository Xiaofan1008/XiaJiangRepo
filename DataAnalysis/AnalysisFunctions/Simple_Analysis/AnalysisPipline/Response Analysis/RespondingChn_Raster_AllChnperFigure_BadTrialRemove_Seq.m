%% =============================================================
%   Raster + PSTH (All Channels) per Set × Amp × PTD
%   - EXCLUDES bad trials from Plot and PSTH calculation
% =============================================================
clear; 
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));
%% ====================== USER SETTINGS ========================
data_folder      = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq4_5ms_251125_154235';
Electrode_Type   = 1;          
raster_chn_start = 1;          
raster_chn_end   = 32;
% ---- plotting windows ----
ras_win       = [-20 40];      
bin_ms_raster = 1;             
smooth_ms     = 5;             
Plot_Amps = [];     
fig_position = [50 50 1600 900];
%% ===================== INITIAL SETUP =========================
if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end
cd(data_folder);
fprintf('Raster+PSTH plotting in folder:\n%s\n\n', data_folder);
parts       = split(data_folder, filesep);
last_folder = parts{end};
u           = strfind(last_folder, '_');
if numel(u) >= 4
    base_name = last_folder(1:u(end-1)-1);   
else
    base_name = last_folder;
end
FS = 30000;
%% ===================== LOAD SPIKES  ==========================
ssd_file  = [base_name '.sp_xia_SSD.mat'];
base_file = [base_name '.sp_xia.mat'];
if isfile(ssd_file)
    S = load(ssd_file);
    if     isfield(S,'sp_corr'), sp = S.sp_corr;
    elseif isfield(S,'sp_SSD'),  sp = S.sp_SSD;
    else, error('SSD file missing usable variable.'); end
elseif isfile(base_file)
    S = load(base_file);
    if     isfield(S,'sp_clipped'), sp = S.sp_clipped;
    elseif isfield(S,'sp'),         sp = S.sp;
    else, error('Base spike file missing usable variable.'); end
else
    error('No spike file found.');
end
nCh = numel(sp);
fprintf('loaded %d spike channels.\n', nCh);
%% ===================== LOAD BadChannels & BadTrials ==========
BadCh_perSet = {};
BadTrialsPerCh = {};
badch_file = [base_name '.BadChannels.mat'];
if isfile(badch_file)
    tmp = load(badch_file);
    if isfield(tmp,'BadCh_perSet'), BadCh_perSet = tmp.BadCh_perSet; end
end
badtr_file = [base_name '.BadTrials.mat'];
if isfile(badtr_file)
    tmp = load(badtr_file);
    if isfield(tmp,'BadTrials'), BadTrialsPerCh = tmp.BadTrials; end
end
%% ===================== LOAD Responding Channels ==============
Resp = [];
resp_file = [base_name '_RespondingChannels.mat'];
hasResp = false;
if isfile(resp_file)
    tmp = load(resp_file);
    if isfield(tmp,'Responding'), Resp = tmp.Responding; hasResp = true; end
end
%% ===================== LOAD TRIGGERS =========================
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);   
nTrig = numel(trig);
%% ===================== LOAD StimParams & DECODE ==============
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat file found.');
S = load(fileDIR(1).name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams = S.StimParams;
sim_stim   = S.simultaneous_stim;
E_MAP      = S.E_MAP;
n_Trials   = S.n_Trials;
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:sim_stim:end);
[Amps,~,ampIdx] = unique(trialAmps);
Amps(Amps == -1) = 0;
if isempty(Plot_Amps), Plot_Amps = Amps(:).'; end
if sim_stim > 1
    PTD_all = cell2mat(StimParams(3:sim_stim:end,6)); 
else
    PTD_all = zeros(n_Trials,1);    
end
[PTDs,~,ptdIdx] = unique(PTD_all);
nPTD = numel(PTDs);
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));
stimSeq = zeros(n_Trials, sim_stim);
for t = 1:n_Trials
    rr = (t-1)*sim_stim + (1:sim_stim);
    v  = idx_all(rr); v = v(v>0); stimSeq(t,1:numel(v)) = v;
end
[uniqueComb,~,combClass] = unique(stimSeq, 'rows','stable');
nSets = size(uniqueComb,1);
%% ===================== ELECTRODE MAP =========================
d = Depth_s(Electrode_Type);   
depth_range = raster_chn_start : min(raster_chn_end, numel(d));
nChPlot = numel(depth_range);
%% ===================== PSTH KERNEL ===========================
edges = ras_win(1) : bin_ms_raster : ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;
g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2); g = g / sum(g);
%% ===================== MAIN CONDITION LOOPS ==================
for si = 1:nSets
    stimIdx = uniqueComb(si, uniqueComb(si,:)>0);
    setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, 'UniformOutput', false), ' + ');
    
    if ~isempty(BadCh_perSet) && si <= numel(BadCh_perSet)
        badDepthThisSet = BadCh_perSet{si};
    else
        badDepthThisSet = [];
    end
    
    for aVal = Plot_Amps
        ai = find(abs(Amps - aVal) < 1e-6, 1);
        if isempty(ai), continue; end
        
        for pi = 1:nPTD
            trials_this = find(combClass==si & ampIdx==ai & ptdIdx==pi);
            if isempty(trials_this), continue; end
            
            PTD_ms = PTDs(pi)/1000;
            figTitle = sprintf('Set %d (%s) | Amp %.1f µA | PTD %.1f ms | Sequential', ...
                               si, setLabel, aVal, PTD_ms);
            figure('Color','w','Name',figTitle,'Position',fig_position);
            tiledlayout('flow','TileSpacing','compact','Padding','compact');
            sgtitle(figTitle, 'FontSize', 14, 'FontWeight','bold');
            
            for idxDepth = 1:nChPlot
                ich = depth_range(idxDepth);   
                ch  = d(ich);                  
                
                if ch < 1 || ch > nCh || isempty(sp{ch})
                    nexttile; axis off; continue;
                end
                
                sp_times = sp{ch}(:,1);
                
                % Get bad trials for this channel
                if ~isempty(BadTrialsPerCh) && ich <= numel(BadTrialsPerCh)
                    badTr_ch = BadTrialsPerCh{ich};
                else
                    badTr_ch = [];
                end
                
                % collect spikes (SKIPPING BAD TRIALS)
                allTrialSpikes = cell(numel(trials_this),1);
                counts         = zeros(1, numel(edges)-1);
                nValid = 0;
                
                for ti = 1:numel(trials_this)
                    tr = trials_this(ti);
                    
                    % --- FILTER: Skip if this trial is bad for this channel ---
                    if ismember(tr, badTr_ch)
                        continue; 
                    end
                    
                    nValid = nValid + 1;
                    t0 = trig(tr)/FS*1000; 
                    tt = sp_times;
                    tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                    allTrialSpikes{nValid} = tt; % Store in next valid slot
                    counts = counts + histcounts(tt, edges);
                end
                
                % Trim empty cells
                allTrialSpikes = allTrialSpikes(1:nValid);
                
                if nValid == 0
                    rate_s = zeros(size(ctrs));
                else
                    rate   = counts / (nValid*bin_s); % Normalized by CLEAN trials only
                    rate_s = filter(g,1,rate);
                end
                maxRate   = max(rate_s);
                yMaxPSTH  = max(50, ceil(maxRate*1.1/10)*10);
                
                % ================= PLOT ==================
                ax = nexttile; hold(ax,'on');
                
                respTag = ''; isResp = false;
                if hasResp && si <= numel(Resp.set) && ai <= numel(Resp.set(si).amp) && ...
                   pi <= numel(Resp.set(si).amp(ai).ptd) && ich <= numel(Resp.set(si).amp(ai).ptd(pi).channel)
                    R = Resp.set(si).amp(ai).ptd(pi).channel(ich);
                    if isfield(R,'is_responsive') && R.is_responsive
                        isResp = true; respTag = 'RESP';
                    end
                end
                
                isBadCh = ismember(ich, badDepthThisSet);
                if isResp && ~isBadCh, set(ax,'Color',[1.0 0.88 0.88]);
                elseif isResp && isBadCh, set(ax,'Color',[1.0 0.90 0.75]);
                elseif isBadCh, set(ax,'Color',[0.95 0.95 0.95]); end
                
                if isResp, ax.XColor=[0.6 0 0]; ax.YColor=[0.6 0 0]; ax.LineWidth=1.4; end
                
                yyaxis(ax,'left');
                if any(rate_s), plot(ax, ctrs, rate_s, 'LineWidth', 1.4); end
                xlim(ax, ras_win); ylim(ax, [0 yMaxPSTH]); ylabel(ax,'Rate (sp/s)');
                
                yyaxis(ax,'right');
                for ti = 1:nValid
                    tt = allTrialSpikes{ti};
                    if isempty(tt), continue; end
                    plot(ax, tt, ti*ones(size(tt)), '.', 'Color',[0 0 0], 'MarkerSize',4);
                end
                ylim(ax,[0 nValid+1]); set(ax,'YTick',[]);
                
                xline(ax,0,'r--','LineWidth',1); xlim(ax, ras_win);
                
                chLabel = sprintf('Ch %d', ich);
                if isBadCh, chLabel = [chLabel ' BAD']; end
                chLabel = [chLabel respTag];
                ht = title(ax, chLabel, 'FontSize',11, 'FontWeight','bold', 'Interpreter','none');
                if isResp, set(ht, 'Color', [0.7 0 0]); end
                
                if idxDepth > nChPlot - ceil(sqrt(nChPlot)), xlabel(ax,'Time (ms)'); end
            end
        end 
    end 
end