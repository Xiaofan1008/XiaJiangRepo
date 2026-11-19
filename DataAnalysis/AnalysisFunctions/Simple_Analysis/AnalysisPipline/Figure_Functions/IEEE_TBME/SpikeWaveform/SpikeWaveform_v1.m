%% ============================================================
%   Figure 1 Panel: Population Spike Waveforms (standalone)
%   - Overlaid spike waveforms (aligned at trough) + mean, for one channel,
%     pooled across the whole recording (not just one trial window).
%   - Carries the amplitude scale bar for BOTH this panel and the
%     neighboring trace panel (Fig1_FilteredTrace_only.m) -- see the
%     "MUST MATCH" notes below.
%   - Based on the alignment logic in QC_Spikewaveform_MultiPTDs_Seq.m
% ============================================================
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
data_folder    = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim3';
Electrode_Type = 1;               % 0: rigid; 1: single-shank flex; 2: four-shank flex
depth_channel  = 25;              % SCALAR -- use the SAME channel as the trace panel for narrative consistency

amp_reject_threshold = 300;       % drop waveforms with any sample beyond +/- this (uV), same as QC_Spikewaveform_MultiPTDs_Seq.m

% Population overlay: pull spikes from the whole recording for this channel
% (not just the ~500 ms window shown in the trace), so this reads as a
% genuine "population of spike waveforms," not just the handful visible in
% one snippet. Capped and randomly subsampled if there are a lot, so the
% overlay doesn't turn into unreadable gray mush.
max_waveforms_overlay = 200;
rng_seed = 1;                     % fixed seed so re-running gives the same subsample while you tune styling

% ---- Display ----
% *** trace_ylim_uv below MUST MATCH trace_ylim_uv in Fig1_FilteredTrace_only.m ***
% *** panel_height_in below MUST MATCH panel_height_in in Fig1_FilteredTrace_only.m ***
% Both are required for the shared "100 uV" scale bar to mean the same thing
% in both panels (same uV range over the same physical height).
wave_ylim_uv    = [-250 250];
wave_scalebar   = [100 0.5];      % [amplitude_uV, time_ms] -- 0.5 ms = 500 us
line_color_indiv = [0.6 0.6 0.6]; % individual waveforms
line_alpha_indiv = 0.25;
line_width_indiv = 0.5;
line_color_mean  = 'k';
line_width_mean  = 1.5;

panel_width_in  = 2.36;           % *** trace_width + gap + this should sum to the target row width (e.g. 7.16 in) ***
panel_height_in = 1.0;            % *** must match panel_height_in in Fig1_FilteredTrace_only.m ***
export_dpi      = 600;            % lineart minimum per IEEE TBME template (mixed figure meets the stricter rule)

show_title = true;                % set false for the polished/final export
save_fig_path = '';   % '' to skip saving

%% ====================== CHECK FOLDER ======================
if ~isfolder(data_folder), error('Invalid folder: %s', data_folder); end
cd(data_folder);

%% ====================== BASE NAME ======================
parts = split(data_folder, filesep);
lastfld = parts{end};
u = strfind(lastfld,'_');
if numel(u) >= 4, base_name = lastfld(1:u(end-1)-1); else, base_name = lastfld; end

%% ====================== HEADER / MAPPING ======================
[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;
d  = Depth_s(Electrode_Type);
ch_intan = d(depth_channel);

%% ====================== LOAD SPIKES (QC > SSD > base) ======================
qc_file   = [base_name '.sp_xia_QC.mat'];
ssd_file  = [base_name '.sp_xia_SSD.mat'];
base_file = [base_name '.sp_xia.mat'];
sp_use = [];
if isfile(qc_file)
    fprintf('Loading QC spikes: %s\n', qc_file);
    S = load(qc_file);
    if isfield(S,'sp_qc'), sp_use = S.sp_qc;
    elseif isfield(S,'sp_corr'), sp_use = S.sp_corr;
    elseif isfield(S,'sp_clipped'), sp_use = S.sp_clipped;
    else, error('No usable spike variable in %s', qc_file); end
elseif isfile(ssd_file)
    fprintf('Loading SSD spikes: %s\n', ssd_file);
    S = load(ssd_file);
    if isfield(S,'sp_pca'), sp_use = S.sp_pca;
    elseif isfield(S,'sp_corr'), sp_use = S.sp_corr;
    elseif isfield(S,'sp_SSD'), sp_use = S.sp_SSD;
    elseif isfield(S,'sp_in'), sp_use = S.sp_in;
    else, error('No usable spike variable in %s', ssd_file); end
elseif isfile(base_file)
    fprintf('Loading base spikes: %s\n', base_file);
    S = load(base_file);
    if isfield(S,'sp_clipped'), sp_use = S.sp_clipped;
    elseif isfield(S,'sp'), sp_use = S.sp;
    else, error('No usable spike variable in %s', base_file); end
else
    error('No spike file found in %s', data_folder);
end

if ch_intan > numel(sp_use) || isempty(sp_use{ch_intan})
    error('Channel %d (depth idx %d) has no spikes in the loaded spike file.', ch_intan, depth_channel);
end
sp_wave_all = sp_use{ch_intan}(:,2:end);
valid_idx   = all(abs(sp_wave_all) <= amp_reject_threshold, 2);
sp_wave_all = sp_wave_all(valid_idx,:);
n_total     = size(sp_wave_all,1);
wf_len      = size(sp_wave_all,2);
t_wave      = ((0:wf_len-1) / FS * 1000) - (wf_len/2)/FS*1000;   % centered around 0 for display

fprintf('Channel %d: %d spikes pass QC/amplitude screen.\n', depth_channel, n_total);
if n_total == 0
    error('No spikes survive the amplitude-reject screen on this channel.');
end

%% ====================== SUBSAMPLE FOR OVERLAY ======================
rng(rng_seed);
if n_total > max_waveforms_overlay
    keep_idx = randperm(n_total, max_waveforms_overlay);
else
    keep_idx = 1:n_total;
end
sp_wave_sub = sp_wave_all(keep_idx,:);

%% ====================== ALIGN TO TROUGH (same trick as QC script) ======================
aligned = zeros(size(sp_wave_sub));
for k = 1:size(sp_wave_sub,1)
    [~, min_idx] = min(sp_wave_sub(k,:));
    shift = ceil(wf_len/2) - min_idx;
    aligned(k,:) = circshift(sp_wave_sub(k,:), shift, 2);
end
mean_wave = mean(aligned, 1);

%% ====================== PLOT ======================
fig = figure('Color','w','Units','inches','Position',[1 1 panel_width_in panel_height_in], ...
    'PaperUnits','inches','PaperPosition',[0 0 panel_width_in panel_height_in]);
ax = axes; hold on;

plot(t_wave, aligned', 'Color', [line_color_indiv line_alpha_indiv], 'LineWidth', line_width_indiv);
plot(t_wave, mean_wave, line_color_mean, 'LineWidth', line_width_mean);

ylim(wave_ylim_uv); xlim([t_wave(1) t_wave(end)]);
axis off;

add_scalebar(ax, t_wave(1), wave_ylim_uv(1)*0.95, wave_scalebar(2), wave_scalebar(1), ...
    sprintf('%g \\mus', wave_scalebar(2)*1000), sprintf('%g \\muV', wave_scalebar(1)));

if show_title
    title(sprintf('%s | Ch %d (depth) | n=%d (of %d)', base_name, depth_channel, size(aligned,1), n_total), ...
        'FontWeight','normal', 'Interpreter','none');
end

if ~isempty(save_fig_path)
    exportgraphics(fig, save_fig_path, 'Resolution', export_dpi);
    fprintf('Saved to %s\n', save_fig_path);
end

%% ====================== LOCAL FUNCTION ======================
function add_scalebar(ax, x0, y0, xlen, ylen, xlabelstr, ylabelstr)
    axes(ax); %#ok<LAXES>
    plot([x0 x0+xlen], [y0 y0], 'k-', 'LineWidth', 1.5);
    plot([x0 x0], [y0 y0+ylen], 'k-', 'LineWidth', 1.5);
    text(x0+xlen/2, y0 - 0.06*range(ylim), xlabelstr, ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
    text(x0 - 0.02*range(xlim), y0+ylen/2, ylabelstr, ...
        'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9);
end