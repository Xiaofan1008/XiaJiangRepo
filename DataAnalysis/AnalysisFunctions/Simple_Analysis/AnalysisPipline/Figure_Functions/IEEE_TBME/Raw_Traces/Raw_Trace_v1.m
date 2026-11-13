%% ============================================================
%   Figure 1 Panel: Filtered Trace (standalone)
%   - depth_channel can be a SCALAR (single clean panel) or a VECTOR
%     (stacked montage across channels, to scan for a good example).
%   - Change USER SETTINGS and re-run to browse for a good example.
%   - Based on the loading logic in Raw_trace_Monitor.m
% ============================================================
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
data_folder     = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim3';
Electrode_Type  = 1;              % 0: rigid; 1: single-shank flex; 2: four-shank flex
depth_channel   = 26;             % SCALAR for one channel (final panel), or a VECTOR (e.g. 1:32) to scan several at once

% ---- Which segment to show ----
% false = just grab a fixed window from the recording (simplest -- pick a
%         clean spontaneous stretch away from stimulation artifacts)
% true  = center the window on a specific stimulation trial, selected the
%         same way Raw_trace_Monitor.m does (by set/PTD/amplitude)
use_trigger_window = true;

% -- If use_trigger_window = true --
amp_to_plot    = 5;               % uA
ptd_to_plot    = 0;               % ms (0 = simultaneous)
set_to_plot    = 1;
trial_pick     = 1;               % which matching trial (1st, 2nd, ...) to use
plot_window_ms = [4 600];        % window relative to trigger (ms)

% -- If use_trigger_window = false --
window_start_sec = 30;            % start time in the recording (s) -- pick a quiet stretch
window_length_ms = 5000;           % length of the segment to display (ms)

% ---- Which file to read ----
% 'raw' = amplifier.dat (wideband)
% 'dn'  = amplifier_dn_sab.dat (artifact-blanked)
% 'mu'  = <base_name>.mu_sab.dat (filtered 300-6000 Hz multi-unit band)
trace_type = 'mu';

% mu_sab.dat is written by allExtract_sab_1.m as SCALEFACTOR*mu2 (SCALEFACTOR = 10),
% where mu2 is already in uV -- so the correct read-back for 'mu' is value/10.
% raw/dn files are read as value*0.195 (standard Intan bit->uV scaling).
% This is left explicit (rather than silently assumed) so you can sanity-check
% the resulting amplitude before trusting the panel.
mu_scale_mode = 'divide10';       % 'divide10' (matches allExtract_sab_1.m) or 'times0195' (matches old Raw_trace_Monitor.m)

% ---- Optional: overlay spike-time ticks on the trace ----
% Requires a *.sp_xia_QC.mat / *.sp_xia_SSD.mat / *.sp_xia.mat file in data_folder.
% Set to false to plot the trace alone with no spike-file dependency.
show_spike_ticks     = false;      % set true if you want ticks while exploring; reference panel style has none
amp_reject_threshold = 300;       % drop spikes with any waveform sample beyond +/- this (uV)

% ---- Display (single-channel final panel) ----
show_title              = true;   % set false for the polished/final export
show_amplitude_scalebar = false;  % target style only shows the time bar (amplitude comes from the waveform panel)
auto_ylim               = false;   % tightly fit y-range to this channel's actual amplitude, like the reference panel
auto_ylim_margin        = 1.15;   % headroom multiplier on top of max(abs(trace)) when auto_ylim = true
trace_ylim_uv           = [-250 250];   % used only when auto_ylim = false
trace_scalebar          = [100 100];    % [amplitude_uV, time_ms] -- amplitude part only drawn if show_amplitude_scalebar = true
line_color              = 'k';
line_width              = 1;   % thinner line reads better at this wide/short aspect ratio
panel_size_px           = [900 130];    % [width height] for the single-channel final panel

% ---- Display (multi-channel scan/montage) ----
montage_row_height_px = 90;       % only used when depth_channel has >1 entry

save_fig_path  = '';              % e.g. 'Fig1_trace_example.pdf'; leave '' to skip saving

%% ====================== CHECK FOLDER ======================
if ~isfolder(data_folder), error('Invalid folder: %s', data_folder); end
cd(data_folder);

%% ====================== BASE NAME ======================
parts = split(data_folder, filesep);
lastfld = parts{end};
u = strfind(lastfld,'_');
if numel(u) >= 4, base_name = lastfld(1:u(end-1)-1); else, base_name = lastfld; end

%% ====================== CHOOSE FILE ======================
switch trace_type
    case 'raw', data_file = 'amplifier.dat';
    case 'dn',  data_file = 'amplifier_dn_sab.dat';
    case 'mu',  data_file = [base_name '.mu_sab.dat'];
    otherwise, error('Unknown trace_type');
end
if ~isfile(data_file), error('Trace file not found: %s', data_file); end

%% ====================== HEADER / MAPPING ======================
[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS   = freq_params.amplifier_sample_rate;
nChn = numel(amp_channels);
d    = Depth_s(Electrode_Type);
ch_intan_list = d(depth_channel(:)');     % row vector, one intan channel per requested depth channel

%% ====================== DETERMINE TIME WINDOW ======================
if use_trigger_window
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
    trig = loadTrig(0);

    fDIR = dir('*_exp_datafile_*.mat');
    Sx = load(fDIR(1).name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');
    StimParams = Sx.StimParams; simN = Sx.simultaneous_stim; E_MAP = Sx.E_MAP;

    trialAmps = cell2mat(StimParams(2:end,16)); trialAmps = trialAmps(1:simN:end);
    postTrig  = cell2mat(StimParams(2:end,6));  postTrig  = postTrig(2:simN:end);
    stimNames = StimParams(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP(2:end));
    comb = zeros(Sx.n_Trials, simN);
    for t = 1:Sx.n_Trials
        v = idx_all((t-1)*simN + (1:simN)); v = v(v>0);
        comb(t,1:numel(v)) = v(:)';
    end
    [~,~,combClass] = unique(comb,'rows','stable');

    ptd_val_us = ptd_to_plot * 1000;
    tidx = intersect(intersect(find(combClass==set_to_plot), find(postTrig==ptd_val_us)), ...
                      find(trialAmps==amp_to_plot));
    if isempty(tidx) || trial_pick > numel(tidx)
        error('No matching trial found for the requested set/PTD/amp/trial_pick.');
    end
    tr = tidx(trial_pick);

    samp_win = round(plot_window_ms/1000 * FS);
    Nsamp    = samp_win(2) - samp_win(1) + 1;
    t0_samp  = trig(tr) + samp_win(1);
    time_ms  = (samp_win(1):samp_win(2)) / FS * 1000;
else
    Nsamp   = round(window_length_ms/1000 * FS);
    t0_samp = round(window_start_sec * FS);
    time_ms = (0:Nsamp-1) / FS * 1000;
end
t0_ms = t0_samp / FS * 1000;

%% ====================== READ TRACE (all requested channels at once) ======================
fid = fopen(data_file,'r');
if fid < 0, error('Cannot open %s', data_file); end
byte_pos = t0_samp * nChn * 2;
fseek(fid, byte_pos, 'bof');
switch trace_type
    case {'raw','dn'}
        data_block = fread(fid, [nChn, Nsamp], 'int16');
        trace_block_uv = data_block(ch_intan_list,:) .* 0.195;   % rows = requested channels
    case 'mu'
        data_block = fread(fid, [nChn, Nsamp], 'short');
        switch mu_scale_mode
            case 'divide10',  trace_block_uv = data_block(ch_intan_list,:) ./ 10;
            case 'times0195', trace_block_uv = data_block(ch_intan_list,:) .* 0.195;
            otherwise, error('Unknown mu_scale_mode');
        end
end
fclose(fid);

%% ====================== OPTIONAL SPIKE TICKS (load once, reuse per channel) ======================
sp_use = [];
if show_spike_ticks
    qc_file   = [base_name '.sp_xia_QC.mat'];
    ssd_file  = [base_name '.sp_xia_SSD.mat'];
    base_file = [base_name '.sp_xia.mat'];
    if isfile(qc_file)
        S = load(qc_file);
        if isfield(S,'sp_qc'), sp_use = S.sp_qc;
        elseif isfield(S,'sp_corr'), sp_use = S.sp_corr;
        elseif isfield(S,'sp_clipped'), sp_use = S.sp_clipped;
        end
    elseif isfile(ssd_file)
        S = load(ssd_file);
        if isfield(S,'sp_pca'), sp_use = S.sp_pca;
        elseif isfield(S,'sp_corr'), sp_use = S.sp_corr;
        elseif isfield(S,'sp_SSD'), sp_use = S.sp_SSD;
        elseif isfield(S,'sp_in'), sp_use = S.sp_in;
        end
    elseif isfile(base_file)
        S = load(base_file);
        if isfield(S,'sp_clipped'), sp_use = S.sp_clipped;
        elseif isfield(S,'sp'), sp_use = S.sp;
        end
    else
        fprintf('No spike file found -- plotting traces without ticks.\n');
    end
end

%% ====================== PLOT ======================
nCh_plot = numel(depth_channel);
win_ms_total = Nsamp/FS*1000;

if nCh_plot == 1
    % ---- Single clean panel, styled to match the target reference (e.g. colleague's Fig. panel) ----
    figure('Color','w','Position',[200 200 panel_size_px(1) panel_size_px(2)]);
    ax = axes; hold on;
    trace_uv = trace_block_uv(1,:);
    plot(time_ms, trace_uv, line_color, 'LineWidth', line_width);

    if auto_ylim
        ylim_max = max(abs(trace_uv)) * auto_ylim_margin;
        this_ylim = [-ylim_max ylim_max];
    else
        this_ylim = trace_ylim_uv;
    end

    if show_spike_ticks
        spike_times_rel = get_spike_ticks(sp_use, ch_intan_list(1), t0_ms, win_ms_total, amp_reject_threshold);
        if ~isempty(spike_times_rel)
            tick_y = this_ylim(2) * 0.92;
            plot(spike_times_rel, tick_y*ones(size(spike_times_rel)), 'v', ...
                'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',5);
        end
    end

    ylim(this_ylim); xlim([time_ms(1) time_ms(end)]);
    axis off;

    if show_amplitude_scalebar
        add_scalebar(ax, time_ms(1), this_ylim(1)*0.95, trace_scalebar(2), trace_scalebar(1), ...
            sprintf('%g ms', trace_scalebar(2)), sprintf('%g \\muV', trace_scalebar(1)));
    else
        add_xscalebar(ax, time_ms(1), this_ylim(1)*0.95, trace_scalebar(2), sprintf('%g ms', trace_scalebar(2)));
    end

    if show_title
        title(sprintf('%s | Ch %d (depth) | %s', base_name, depth_channel, trace_type), ...
            'FontWeight','normal', 'Interpreter','none');
    end
else
    % ---- Stacked montage across channels, to scan for a good example ----
    figure('Color','w','Position',[100 50 900 min(montage_row_height_px*nCh_plot, 1400)]);
    tl = tiledlayout(nCh_plot, 1, 'TileSpacing','none', 'Padding','compact');
    for k = 1:nCh_plot
        ax = nexttile; hold on;
        trace_uv = trace_block_uv(k,:);
        plot(time_ms, trace_uv, line_color, 'LineWidth', 1);

        spike_times_rel = get_spike_ticks(sp_use, ch_intan_list(k), t0_ms, win_ms_total, amp_reject_threshold);
        n_sp = numel(spike_times_rel);
        if n_sp > 0
            tick_y = trace_ylim_uv(2) * 0.9;
            plot(spike_times_rel, tick_y*ones(size(spike_times_rel)), 'v', ...
                'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',3);
        end

        ylim(trace_ylim_uv); xlim([time_ms(1) time_ms(end)]);
        box off; set(gca,'YColor','none');
        text(time_ms(end), 0, sprintf('Ch %d (n=%d)', depth_channel(k), n_sp), ...
            'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',8);

        if k < nCh_plot
            set(gca,'XColor','none');
        else
            xlabel('Time (ms)');
            set(gca,'XColor','k');
        end
    end
    sgtitle(sprintf('%s | %s | scanning %d channels', base_name, trace_type, nCh_plot), 'Interpreter','none');
end

if ~isempty(save_fig_path)
    exportgraphics(gcf, save_fig_path, 'Resolution', 300);
    fprintf('Saved to %s\n', save_fig_path);
end

%% ====================== LOCAL FUNCTIONS ======================
function spike_times_rel = get_spike_ticks(sp_use, ch_intan, t0_ms, win_ms, amp_reject_threshold)
    spike_times_rel = [];
    if isempty(sp_use) || ch_intan > numel(sp_use) || isempty(sp_use{ch_intan}), return; end
    sp_times_all = sp_use{ch_intan}(:,1);
    sp_wave_all  = sp_use{ch_intan}(:,2:end);
    valid_idx    = all(abs(sp_wave_all) <= amp_reject_threshold, 2);
    sp_times_all = sp_times_all(valid_idx);
    win_mask = sp_times_all >= t0_ms & sp_times_all < (t0_ms + win_ms);
    spike_times_rel = sp_times_all(win_mask) - t0_ms;
end

function add_scalebar(ax, x0, y0, xlen, ylen, xlabelstr, ylabelstr)
    axes(ax); %#ok<LAXES>
    plot([x0 x0+xlen], [y0 y0], 'k-', 'LineWidth', 1.5);
    plot([x0 x0], [y0 y0+ylen], 'k-', 'LineWidth', 1.5);
    text(x0+xlen/2, y0 - 0.06*range(ylim), xlabelstr, ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
    text(x0 - 0.02*range(xlim), y0+ylen/2, ylabelstr, ...
        'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9);
end

function add_xscalebar(ax, x0, y0, xlen, xlabelstr)
    % Time-only scale bar (no amplitude tick), matching the reference panel style.
    axes(ax); %#ok<LAXES>
    plot([x0 x0+xlen], [y0 y0], 'k-', 'LineWidth', 1.5);
    text(x0+xlen/2, y0 - 0.06*range(ylim), xlabelstr, ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
end