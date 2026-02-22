%% ============================================================
%   Unified Raster + PSTH (Stacked Subplots, Sequential Only)
%   - Top Plot: Shaded bands per ISI + Raster dots
%   - Bottom Plot: PSTH curves overlaid
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ===================== USER INPUTS ===================== %%
folder_seq    = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_1'; 
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% --- Analysis Settings ---
target_channels = [1:40];     % Choose which channel(s) to plot
plot_amps      = [10];       % Choose which Amplitude(s) to plot
target_PTDs    = [5 8 10 12 15 17 20]; % Only sequential ISIs > 0
num_trials_plot = 20;        % Set the number of random trials to plot per ISI
ras_win        = [-20 60];   % ms (Set your desired time window here)
bin_ms         = 1;
smooth_ms      = 10;
FS             = 30000;

%% Colors
% Distinct colors for up to 10 ISIs
ptd_base_colors = turbo(10); 

%% ============================================================
% LOAD DATASETS
% ============================================================
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
available_PTDs(available_PTDs == 0) = []; % Remove Sim to guarantee Seq only
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
        for s = 1:D_seq.nSets
            stimCh = D_seq.uniqueComb(s, D_seq.uniqueComb(s,:)>0);
            
            % [MODIFIED] Create Figure with Tiled Layout for Top/Bottom split
            fig = figure('Color','w','Position',[200 100 800 800]);
            tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            title(tl, sprintf('Seq Set %d (Ch %s) | Rec: Ch %d | Amp: %.0f µA', s, num2str(stimCh), target_channel, amp_val), 'FontSize', 14, 'FontWeight', 'bold');
            
            % [MODIFIED] Setup Top Axis (Rasters)
            ax_raster = nexttile(1); hold(ax_raster, 'on'); box(ax_raster, 'off');
            ylabel(ax_raster, 'Trial Number', 'FontWeight', 'bold');
            
            % [MODIFIED] Setup Bottom Axis (PSTH)
            ax_psth = nexttile(2); hold(ax_psth, 'on'); box(ax_psth, 'off');
            ylabel(ax_psth, 'Firing Rate (Sp/s)', 'FontWeight', 'bold');
            xlabel(ax_psth, 'Time (ms)', 'FontWeight', 'bold');
            
            y_raster_offset = 0; % Tracks vertical position for stacking rasters
            
            % --- LOOP OVER ISIs ---
            for p_idx = 1:length(use_PTDs)
                ptd_val = use_PTDs(p_idx);
                this_col = ptd_base_colors(mod(p_idx-1, size(ptd_base_colors,1))+1, :);
                
                S_ch = D_seq.sp{d(target_channel)};
                trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val & abs(D_seq.trialPTD - ptd_val) < 0.01);
                
                % Randomly select a specific number of trials
                if length(trial_ids) > num_trials_plot
                    rand_idx = randperm(length(trial_ids), num_trials_plot);
                    trial_ids = trial_ids(rand_idx);
                end
                
                num_trials = length(trial_ids);
                if num_trials == 0, continue; end
                
                % ============ TOP PLOT: SHADING & RASTERS ============
                y_bottom = y_raster_offset;
                y_top = y_raster_offset + num_trials;
                
                % 1. Draw Shaded Background Box on Raster Axis
                patch(ax_raster, [ras_win(1) ras_win(2) ras_win(2) ras_win(1)], ...
                      [y_bottom y_bottom y_top y_top], ...
                      this_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                
                % 2. Draw Raster Dots
                y_curr = y_raster_offset;
                for tr = trial_ids'
                    t0 = D_seq.trig(tr)/FS*1000;
                    tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                    
                    if ~isempty(tt)
                        % [MODIFIED] Plot dots directly onto ax_raster
                        scatter(ax_raster, tt, repmat(y_curr + 0.5, size(tt)), 6.5, 'k', 'filled', 'HandleVisibility', 'off');
                    end
                    y_curr = y_curr + 1;
                end
                y_raster_offset = y_raster_offset + num_trials;
                
                % ============ BOTTOM PLOT: PSTH ============
                counts = zeros(1,length(edges)-1);
                for tr = trial_ids'
                    t0 = D_seq.trig(tr)/FS*1000;
                    tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                    counts = counts + histcounts(tt,edges);
                end
                rate_seq = filter(g,1, counts/(num_trials*bin_s));
                
                % [MODIFIED] Plot smoothed PSTH line directly onto ax_psth
                plot(ax_psth, ctrs, rate_seq, 'LineStyle', '-', 'Marker', 'none', 'Color', this_col, 'LineWidth', 3, 'DisplayName', sprintf('ISI %.0f ms', ptd_val));
                
            end % End PTD loop
            
            % --- FORMATTING & AESTHETICS ---
            % Format Raster Axis (Top)
            xlim(ax_raster, ras_win);
            ylim(ax_raster, [0, max(1, y_raster_offset)]);
            xline(ax_raster, 0, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            set(ax_raster, 'XTickLabel', []); % Remove x-labels from top plot so they don't crowd the bottom plot
            
            % Format PSTH Axis (Bottom)
            xlim(ax_psth, ras_win);
            xline(ax_psth, 0, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            legend(ax_psth, 'Location', 'northeast', 'Box', 'off');
            
        end % End Set loop
    end % End Amp loop
end % End Channel loop