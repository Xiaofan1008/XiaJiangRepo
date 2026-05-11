%% =============================================================
%   Raster + PSTH (One Figure per ISI per Channel)
%   - For selected Set(s) × one amplitude
%   - Each figure = one ISI / PTD condition
%   - PSTH + raster are overlaid within the same axes
%   - Uses ONE common PSTH y-limit across all selected ISIs for each channel
%   - EXCLUDES bad trials from Plot and PSTH calculation
%   - Can plot all responsive channels OR only user-selected channels
%   - Paper-style export formatting
% =============================================================
clear; 
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ========================
data_folder      = '/Volumes/MACData/Data/Data_Xia/DX021/Xia_ISI_SimSeq1';
Electrode_Type   = 1;          
raster_chn_start = 1;          
raster_chn_end   = 32;

% ---- select one or multiple sets ----
plot_all_sets = false;   % true = plot all sets, false = use Plot_Sets below
Plot_Sets     = [1 2];   % only used when plot_all_sets = false

% ---- select amplitude ----
target_amp  = 10;        

% ---- select channels to plot ----
plot_all_responsive_channels = false;   % true = plot all responsive channels
Plot_Channels = [20];                % only used when plot_all_responsive_channels = false

% ---- plotting windows ----
ras_win       = [-40 40];      
bin_ms_raster = 1; 
smooth_ms     = 6;             

% Which PTDs to plot (ms). If empty -> plot ALL PTDs.
% 0 = Simultaneous. Example: [0 5 9 10 12 15]
Plot_PTDs = [];    

% ---- export settings ----
save_figures   = false;
save_dir       = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/ISI_PSTH_Raster';
tiff_dpi       = 600;

% ---- paper figure formatting ----
fig_width_cm   = 4.4;
fig_height_cm  = 4.4;
font_name      = 'Arial';
axis_fontsize  = 6;
label_fontsize = 6;
axis_linewidth = 1.0;
fig_bg_color   = 'w';

% ---- plotting style settings ----
psth_color   = [0 0 0];          % black PSTH
raster_color = [0.35 0.35 0.35]; % dark gray raster
stim1_color  = [0 0 0];          % black dashed line at 0 ms
stim2_color  = [0 0 0];          % black dotted line at 2nd pulse
raster_msize = 3;
psth_lw      = 1;

%% ===================== INITIAL SETUP =========================
if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end
cd(data_folder);
fprintf('Raster+PSTH plotting in folder:\n%s\n\n', data_folder);

if save_figures && ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

parts       = split(data_folder, filesep);
last_folder = parts{end};
u           = strfind(last_folder, '_');
if numel(u) > 4
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
    else, error('SSD file missing usable variable.'); 
    end
elseif isfile(base_file)
    S = load(base_file);
    if     isfield(S,'sp_clipped'), sp = S.sp_clipped;
    elseif isfield(S,'sp'),         sp = S.sp;
    else, error('Base spike file missing usable variable.'); 
    end
else
    error('No spike file found.');
end

nCh = numel(sp);
fprintf('loaded %d spike channels.\n', nCh);

%% ===================== LOAD BadChannels & BadTrials ==========
BadCh_perSet = {};
BadTrialsPerCh = {};

badch_file = [base_name '.MultiISIsBadChannels.mat'];
if isfile(badch_file)
    tmp = load(badch_file);
    if isfield(tmp,'BadCh_perSet'), BadCh_perSet = tmp.BadCh_perSet; end
end

badtr_file = [base_name '.MultiISIsBadTrials.mat'];
if isfile(badtr_file)
    tmp = load(badtr_file);
    if isfield(tmp,'BadTrials'), BadTrialsPerCh = tmp.BadTrials; end
end

%% ===================== LOAD Responding Channels ==============
Resp = [];
resp_file = [base_name '_MultiISIRespondingChannels.mat'];
hasResp = false;
if isfile(resp_file)
    tmp = load(resp_file);
    if isfield(tmp,'Responding')
        Resp = tmp.Responding; 
        hasResp = true; 
    end
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
    v  = idx_all(rr); 
    v  = v(v>0); 
    stimSeq(t,1:numel(v)) = v;
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

g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2); 
g = g / sum(g);

%% ===================== SET SELECTION =========================
if plot_all_sets
    sets_to_plot = 1:nSets;
else
    sets_to_plot = Plot_Sets;
    sets_to_plot = sets_to_plot(sets_to_plot >= 1 & sets_to_plot <= nSets);
end

if isempty(sets_to_plot)
    error('No valid stimulation sets selected.');
end

ai = find(abs(Amps - target_amp) < 1e-6, 1);
if isempty(ai)
    error('target_amp = %.1f uA not found in this dataset.', target_amp);
end

fprintf('Selected sets: %s\n', num2str(sets_to_plot));
fprintf('Selected amp : %.1f uA\n\n', target_amp);

%% ===================== MAIN LOOP OVER SETS ===================
for target_set = sets_to_plot

    %% ===================== TARGET CONDITION SELECTION ============
    stimIdx = uniqueComb(target_set, uniqueComb(target_set,:)>0);
    setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, 'UniformOutput', false), ' + ');

    if ~isempty(BadCh_perSet) && target_set <= numel(BadCh_perSet)
        badDepthThisSet = BadCh_perSet{target_set};
    else
        badDepthThisSet = [];
    end

    % ----- Select which PTDs / ISIs to plot -----
    selected_pi = [];
    selected_ptd_ms = [];

    for pi = 1:nPTD
        PTD_ms = PTDs(pi) / 1000;
        
        if ~isempty(Plot_PTDs)
            if ~any(abs(Plot_PTDs - PTD_ms) < 1e-4)
                continue;
            end
        end
        
        trials_this = find(combClass==target_set & ampIdx==ai & ptdIdx==pi);
        if isempty(trials_this)
            continue;
        end
        
        selected_pi(end+1) = pi; %#ok<SAGROW>
        selected_ptd_ms(end+1) = PTD_ms; %#ok<SAGROW>
    end

    if isempty(selected_pi)
        fprintf('No valid PTDs found for Set %d and Amp %.1f uA. Skipped.\n', target_set, target_amp);
        continue;
    end

    fprintf('\nSelected Set = %d (%s)\n', target_set, setLabel);
    fprintf('Selected Amp = %.1f uA\n', target_amp);
    fprintf('Selected PTDs (ms) = %s\n\n', num2str(selected_ptd_ms));

    %% ===================== BUILD RESPONSIVE CHANNEL POOL =========
    responsive_mask = false(numel(d), 1);

    if hasResp
        for jj = 1:length(selected_pi)
            pi = selected_pi(jj);
            
            if target_set <= numel(Resp.set) && ...
               ai <= numel(Resp.set(target_set).amp) && ...
               pi <= numel(Resp.set(target_set).amp(ai).ptd)
                
                this_resp = Resp.set(target_set).amp(ai).ptd(pi).channel;
                
                for ich = depth_range
                    if ich <= numel(this_resp)
                        if isfield(this_resp(ich), 'is_responsive') && this_resp(ich).is_responsive
                            responsive_mask(ich) = true;
                        end
                    end
                end
            end
        end
    else
        warning('No RespondingChannels file found. Using all channels in selected range.');
        responsive_mask(depth_range) = true;
    end

    responsive_channels = find(responsive_mask);
    responsive_channels = responsive_channels(responsive_channels >= raster_chn_start & responsive_channels <= raster_chn_end);
    responsive_channels = responsive_channels(:)';

    fprintf('Responsive channels in selected pool: %s\n', num2str(responsive_channels));

    if isempty(responsive_channels)
        fprintf('No responsive channels found for Set %d × Amp %.1f × selected PTDs.\n', target_set, target_amp);
        continue;
    end

    if plot_all_responsive_channels
        channels_to_plot = responsive_channels;
    else
        channels_to_plot = Plot_Channels;
        channels_to_plot = channels_to_plot(channels_to_plot >= raster_chn_start & channels_to_plot <= raster_chn_end);
        channels_to_plot = intersect(channels_to_plot, responsive_channels, 'stable');
    end

    fprintf('Channels that will be plotted: %s\n', num2str(channels_to_plot));

    if isempty(channels_to_plot)
        fprintf('No valid selected channels to plot for Set %d × Amp %.1f.\n', target_set, target_amp);
        continue;
    end

    %% ===================== MAIN LOOP: ONE CHANNEL, ONE FIGURE PER ISI =====
    for ich = channels_to_plot
        
        ch = d(ich);
        if ch < 1 || ch > nCh || isempty(sp{ch})
            continue;
        end
        
        isBadCh = ismember(ich, badDepthThisSet);
        sp_times = sp{ch}(:,1);
        
        if ~isempty(BadTrialsPerCh) && ich <= numel(BadTrialsPerCh)
            badTr_ch = BadTrialsPerCh{ich};
        else
            badTr_ch = [];
        end
        
        % =========================================================
        % FIRST PASS: compute one common PSTH y-limit across all selected ISIs
        % =========================================================
        maxRate_all = 0;
        
        for jj = 1:length(selected_pi)
            pi = selected_pi(jj);
            trials_this = find(combClass==target_set & ampIdx==ai & ptdIdx==pi);
            
            counts = zeros(1, numel(edges)-1);
            nValid = 0;
            
            for ti = 1:numel(trials_this)
                tr = trials_this(ti);
                
                if ismember(tr, badTr_ch)
                    continue;
                end
                
                nValid = nValid + 1;
                t0 = trig(tr)/FS*1000; 
                tt = sp_times;
                tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                counts = counts + histcounts(tt, edges);
            end
            
            if nValid > 0
                rate   = counts / (nValid * bin_s);
                rate_s = filter(g, 1, rate);
                maxRate_all = max(maxRate_all, max(rate_s));
            end
        end
        
        yMaxPSTH = max(50, ceil(maxRate_all*1.1/10)*10);
        
        % =========================================================
        % SECOND PASS: one figure per ISI
        % =========================================================
        for jj = 1:length(selected_pi)
            pi = selected_pi(jj);
            PTD_ms = PTDs(pi) / 1000;
            
            trials_this = find(combClass==target_set & ampIdx==ai & ptdIdx==pi);
            
            allTrialSpikes = cell(numel(trials_this),1);
            counts         = zeros(1, numel(edges)-1);
            nValid = 0;
            
            for ti = 1:numel(trials_this)
                tr = trials_this(ti);
                
                if ismember(tr, badTr_ch)
                    continue;
                end
                
                nValid = nValid + 1;
                t0 = trig(tr)/FS*1000; 
                tt = sp_times;
                tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                
                allTrialSpikes{nValid} = tt;
                counts = counts + histcounts(tt, edges);
            end
            
            allTrialSpikes = allTrialSpikes(1:nValid);
            
            if nValid == 0
                rate_s = zeros(size(ctrs));
            else
                rate   = counts / (nValid * bin_s);
                rate_s = filter(g,1,rate);
            end
            
            % =====================================================
            % CREATE FIGURE FOR THIS ISI
            % =====================================================
            % fig = figure('Color', fig_bg_color, ...
            %              'Units','centimeters', ...
            %              'Position',[5 5 fig_width_cm fig_height_cm]);

            safe_setLabel = regexprep(setLabel, '[^a-zA-Z0-9]+', '_');
            fig_name = sprintf('StimSet%d_Ch%d_Amp%.1fuA_ISI%.1fms', ...
                target_set, ich, target_amp, PTD_ms);
            
            fig = figure('Color', fig_bg_color, ...
                         'Units','centimeters', ...
                         'Position',[5 5 fig_width_cm fig_height_cm], ...
                         'Name', fig_name, ...
                         'NumberTitle', 'off');
            ax = axes(fig);
            hold(ax,'on');
            
            if isBadCh
                set(ax,'Color',[0.95 0.95 0.95]); 
            end
            
            % ---- PSTH axis (left) ----
            yyaxis(ax,'left');
            ax.YColor = [0 0 0];
            if any(rate_s)
                plot(ax, ctrs, rate_s, 'Color', psth_color, 'LineWidth', psth_lw);
            end
            xlim(ax, ras_win); 
            ylim(ax, [0 yMaxPSTH]); 
            ylabel(ax,'Rate (sp/s)', 'FontSize', label_fontsize, 'FontName', font_name);
            
            % ---- Raster axis (right) ----
            yyaxis(ax,'right');
            ax.YColor = [1 1 1];   % make right axis invisible on white background
            for ti = 1:nValid
                tt = allTrialSpikes{ti};
                if isempty(tt), continue; end
                plot(ax, tt, ti*ones(size(tt)), '.', 'Color', raster_color, 'MarkerSize', raster_msize);
            end
            ylim(ax,[0 nValid+1]); 
            set(ax,'YTick',[]);
            ylabel(ax,'');
            
            % ---- Stimulation timing markers ----
            xline(ax, 0, '--', 'Color', stim1_color, 'LineWidth', 1);
            if PTD_ms > 0
                xline(ax, PTD_ms, ':', 'Color', stim2_color, 'LineWidth', 1);
            end
            xlim(ax, ras_win);
            
            % ---- Axis labels / format ----
            xlabel(ax,'Time (ms)', 'FontSize', label_fontsize, 'FontName', font_name);
            set(ax, 'FontSize', axis_fontsize, ...
                    'FontName', font_name, ...
                    'LineWidth', axis_linewidth, ...
                    'TickDir', 'out');
            box(ax,'off');
            axis(ax,'square');
            
            % =====================================================
            % SAVE FIGURE
            % =====================================================
            if save_figures
                safe_setLabel = regexprep(setLabel, '[^a-zA-Z0-9]+', '_');
                filename = sprintf('%s_Set%d_%s_Amp%.1fuA_Ch%d_ISI%.1fms.tiff', ...
                    base_name, target_set, safe_setLabel, target_amp, ich, PTD_ms);
                exportgraphics(fig, fullfile(save_dir, filename), 'Resolution', tiff_dpi);
            end
        end
    end
end

fprintf('\nDone.\n');