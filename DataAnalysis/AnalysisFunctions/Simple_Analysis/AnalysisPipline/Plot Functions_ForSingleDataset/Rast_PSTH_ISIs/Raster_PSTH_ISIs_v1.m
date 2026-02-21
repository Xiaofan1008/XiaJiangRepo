%% ============================================================
%   Unified Raster + PSTH Overlay (Sequential Only)
%   - Layer 1 (Background): Shaded bands per ISI
%   - Layer 2 (Right Y-Axis): Raster dots stacked by ISI
%   - Layer 3 (Left Y-Axis): PSTH curves overlaid
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ===================== USER INPUTS ===================== %%
% [MODIFIED] Removed Single folder, kept only Seq folder
folder_seq    = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_1'; 
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% --- Analysis Settings ---
target_channels = [1:20];     % Choose which channel(s) to plot
plot_amps      = [10];       % Choose which Amplitude(s) to plot
target_PTDs    = [5 8 10 12 15 17 20 25]; % [MODIFIED] Only sequential ISIs > 0
ras_win        = [-20 60];  % ms (Set your desired time window here)
bin_ms         = 1;
smooth_ms      = 20;
FS             = 30000;

%% Colors
% Distinct colors for up to 7 ISIs
ptd_base_colors = turbo(10); 

%% ============================================================
% LOAD DATASETS
% ============================================================
% [MODIFIED] Removed Single Dataset loading block entirely.

% ---------- Load Combined (Seq) ----------
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
        % [MODIFIED] Loop through Sets. Each Set gets its own layered figure.
        for s = 1:D_seq.nSets
            stimCh = D_seq.uniqueComb(s, D_seq.uniqueComb(s,:)>0);
            
            % Create Figure
            figure('Color','w','Position',[200 100 700 600]);
            ax = gca; hold(ax, 'on');
            title(ax, sprintf('Seq Set %d (Ch %s) | Rec: Ch %d | Amp: %.0f µA', s, num2str(stimCh), target_channel, amp_val), 'FontSize', 14);
            
            y_raster_offset = 0; % Tracks vertical position for stacking rasters
            
            % --- LOOP OVER ISIs ---
            for p_idx = 1:length(use_PTDs)
                ptd_val = use_PTDs(p_idx);
                this_col = ptd_base_colors(mod(p_idx-1, size(ptd_base_colors,1))+1, :);
                
                S_ch = D_seq.sp{d(target_channel)};
                trial_ids = find(D_seq.combClass==s & D_seq.trialAmps==amp_val & abs(D_seq.trialPTD - ptd_val) < 0.01);
                num_trials = length(trial_ids);
                
                if num_trials == 0, continue; end
                
                % ============ LAYER 1 & 2: SHADING & RASTERS (Right Y-Axis) ============
                yyaxis right
                
                % 1. Draw Shaded Background Box
                y_bottom = y_raster_offset;
                y_top = y_raster_offset + num_trials;
                patch([ras_win(1) ras_win(2) ras_win(2) ras_win(1)], ...
                      [y_bottom y_bottom y_top y_top], ...
                      this_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                
                % 2. Draw Raster Dots
                y_curr = y_raster_offset;
                for tr = trial_ids'
                    t0 = D_seq.trig(tr)/FS*1000;
                    tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                    
                    % Plot spikes as tiny black dots for maximum contrast
                    if ~isempty(tt)
                        scatter(ax, tt, repmat(y_curr + 0.5, size(tt)), 5, 'k', 'filled', 'HandleVisibility', 'off');
                    end
                    y_curr = y_curr + 1;
                end
                
                % Update offset for the next ISI block
                y_raster_offset = y_raster_offset + num_trials;
                
                % ============ LAYER 3: PSTH (Left Y-Axis) ============
                yyaxis left
                counts = zeros(1,length(edges)-1);
                for tr = trial_ids'
                    t0 = D_seq.trig(tr)/FS*1000;
                    tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                    counts = counts + histcounts(tt,edges);
                end
                rate_seq = filter(g,1, counts/(num_trials*bin_s));
                
                % Plot smoothed PSTH line
               plot(ax, ctrs, rate_seq, 'LineStyle', '-', 'Marker', 'none', 'Color', this_col, 'LineWidth', 2.5, 'DisplayName', sprintf('ISI %.0f ms', ptd_val));
                
            end % End PTD loop
            
            % --- FORMATTING & AESTHETICS ---
            % X-Axis
            xlim(ax, ras_win);
            xlabel(ax, 'Time (ms)', 'FontWeight', 'bold');
            xline(ax, 0, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Pulse 1 Marker
            
            % Left Y-Axis (PSTH Rate)
            yyaxis left
            ylabel(ax, 'Firing Rate (Sp/s)', 'FontWeight', 'bold');
            ax.YColor = 'k'; % Force axis color to black
            
            % Right Y-Axis (Trial Stacks)
            yyaxis right
            ylabel(ax, 'Trial Number', 'FontWeight', 'bold');
            ylim(ax, [0, max(1, y_raster_offset)]);
            ax.YColor = 'k'; % Force axis color to black
            
            legend(ax, 'Location', 'northeast', 'Box', 'off');
            box(ax, 'off');
            
        end % End Set loop
    end % End Amp loop
end % End Channel loop