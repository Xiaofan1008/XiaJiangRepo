%% ============================================================
%   Spatial Diversity Analysis (Bars + Logistic Regression)
%   - Figure 1: Grouped Bar Chart (Binned Data)
%   - Figure 2: Logistic Regression Curves (Probability Trend)
%   - Output: Saves to 'Result_Spatial_... .mat'
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

% 1 = Single Shank (Linear 32)
% 2 = Four Shank (4x16)
Electrode_Type = 1; 
target_amp = 10; 
target_ptd_seq = 5; 

% Distance Bins (e.g. 0-50, 50-100...)
dist_edges = 0:50:800; 
dist_centers = dist_edges(1:end-1) + 25; 

%% =================== LOAD DATA ====================
cd(folder_sim); 
Rsim=load(dir('*RespondingChannels.mat').name).Responding; 
Ssim=load(dir('*_exp_datafile_*.mat').name); 
amps_all_sim=cell2mat(Ssim.StimParams(2:end,16)); 
trialAmps_sim=amps_all_sim(1:Ssim.simultaneous_stim:end); 
[Amps_sim,~,ampIdx_sim]=unique(trialAmps_sim); 
Amps_sim(Amps_sim==-1)=0;

cd(folder_seq); 
Rseq=load(dir('*RespondingChannels.mat').name).Responding; 
Sseq=load(dir('*_exp_datafile_*.mat').name); 
simN_seq=Sseq.simultaneous_stim; 
amps_all_seq=cell2mat(Sseq.StimParams(2:end,16)); 
trialAmps_seq=amps_all_seq(1:simN_seq:end); 
[Amps_seq,~,ampIdx_seq]=unique(trialAmps_seq); 
Amps_seq(Amps_seq==-1)=0;

PTD_all_us=cell2mat(Sseq.StimParams(3:simN_seq:end,6)); PTDs_ms=unique(PTD_all_us/1000); stimNames=Sseq.StimParams(2:end,1); [~,idx_all]=ismember(stimNames,Sseq.E_MAP(2:end)); if isfield(Sseq,'n_Trials'), nTr_seq=Sseq.n_Trials; else, nTr_seq=(size(Sseq.StimParams,1)-1)/simN_seq; end; comb_seq=zeros(nTr_seq,simN_seq); for t=1:nTr_seq, rr=(t-1)*simN_seq+(1:simN_seq); v=idx_all(rr); v=v(v>0); comb_seq(t,1:numel(v))=v(:).'; end; [uniqueComb_seq,~,combClass_seq]=unique(comb_seq,'rows','stable'); nSets_seq=size(uniqueComb_seq,1);

%% ================= GEOMETRY & TARGETS =================
Coords = get_electrode_coordinates(Electrode_Type);
if Electrode_Type == 1, valid_channels = 1:32; else, valid_channels = [1:16, 33:48, 49:64, 17:32]; end

idx_sim_amp = find(Amps_sim == target_amp, 1);
idx_seq_amp = find(Amps_seq == target_amp, 1);
idx_ptd = find(abs(PTDs_ms - target_ptd_seq) < 0.1, 1);
if isempty(idx_ptd), idx_ptd = 1; end

%% ================= STEP 1: CALCULATE EXACT DISTANCES =================
% A. SIMULTANEOUS
stim_chans_sim = Rsim.set(1).stimChannels;
stim_locs_sim  = Coords(stim_chans_sim, :); 
Sim_Results = []; % [Distance, IsResponsive]
chan_structs = Rsim.set(1).amp(idx_sim_amp).ptd(1).channel;

for ch = valid_channels
    rec_loc = Coords(ch, :);
    dists = sqrt(sum((stim_locs_sim - rec_loc).^2, 2));
    min_dist = round(min(dists));
    
    if min_dist < 1, continue; end 
    
    is_resp = 0;
    if isfield(chan_structs(ch), 'is_responsive') && chan_structs(ch).is_responsive
        is_resp = 1;
    end
    Sim_Results = [Sim_Results; min_dist, is_resp];
end

% B. SEQUENTIAL
Seq_Results = cell(nSets_seq, 1);
for ss = 1:nSets_seq
    stim_chans_seq = uniqueComb_seq(ss, :); stim_chans_seq = stim_chans_seq(stim_chans_seq>0);
    stim_locs_seq  = Coords(stim_chans_seq, :);
    chan_structs = Rseq.set(ss).amp(idx_seq_amp).ptd(idx_ptd).channel;
    
    for ch = valid_channels
        rec_loc = Coords(ch, :);
        dists = sqrt(sum((stim_locs_seq - rec_loc).^2, 2));
        min_dist = round(min(dists));
        
        if min_dist < 1, continue; end
        
        is_resp = 0;
        if isfield(chan_structs(ch), 'is_responsive') && chan_structs(ch).is_responsive
            is_resp = 1;
        end
        Seq_Results{ss} = [Seq_Results{ss}; min_dist, is_resp];
    end
end

%% ================= STEP 2: GROUP BY UNIQUE DISTANCE =================
all_dists = Sim_Results(:, 1);
for ss=1:nSets_seq, all_dists = [all_dists; Seq_Results{ss}(:,1)]; end
Unique_Dists = unique(all_dists);
Unique_Dists = Unique_Dists(Unique_Dists <= 600); % Cap distance

Plot_Data = zeros(length(Unique_Dists), 1 + nSets_seq);

% Fill Sim
for i = 1:length(Unique_Dists)
    d = Unique_Dists(i);
    mask = (Sim_Results(:,1) == d);
    if sum(mask) > 0
        Plot_Data(i, 1) = mean(Sim_Results(mask, 2)) * 100;
    end
end

% Fill Seq
for ss = 1:nSets_seq
    dat = Seq_Results{ss};
    for i = 1:length(Unique_Dists)
        d = Unique_Dists(i);
        mask = (dat(:,1) == d);
        if sum(mask) > 0
            Plot_Data(i, 1+ss) = mean(dat(mask, 2)) * 100;
        end
    end
end

% --- DEFINE PROFESSIONAL COLORS (Muted/Pastel) ---
prof_colors = [ ...
    0.25, 0.45, 0.65;  ... % Blue (Sim)
    0.80, 0.45, 0.30;  ... % Orange (Seq 1)
    0.55, 0.40, 0.65];     % Purple (Seq 2)

%% ================= STEP 3: PLOT GROUPED BARS =================
figure('Color','w', 'Position',[200 600 800 500]); hold on;

% Labels
leg_names = {'Simultaneous'};
for ss = 1:nSets_seq
    stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
    leg_names{end+1} = sprintf('Seq Set %d (Ch:%s)', ss, num2str(stimCh));
end

% Plot
b = bar(Unique_Dists, Plot_Data, 'grouped');

% Apply Colors
for i = 1:length(b)
    if i <= size(prof_colors,1)
        b(i).FaceColor = prof_colors(i,:);
    end
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.85;
end

ylabel('% of Responding Channels', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Distance from Stimulation Site (µm)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Spatial Spread (Bar Chart at %.0f µA)', target_amp), 'FontSize', 14);

ax = gca; ax.Box = 'off';
xticks(Unique_Dists); xtickangle(0);
ylim([0 105]);
legend(leg_names, 'Location', 'northeast', 'Box', 'off');

%% ================= STEP 4: PLOT LOGISTIC REGRESSION =================
% This plot fits a probability curve to the raw data (0s and 1s)
figure('Color','w', 'Position',[200 100 800 500]); hold on;

x_fit = linspace(0, 600, 200)'; % Smooth X-axis for curve plotting

% --- A. Plot Simultaneous ---
% Fit Logistic Regression (Binomial Distribution)
[b_sim, ~, stats_sim] = glmfit(Sim_Results(:,1), Sim_Results(:,2), 'binomial', 'link', 'logit');
[p_sim, ci_lo, ci_hi] = glmval(b_sim, x_fit, 'logit', stats_sim);

% Plot Confidence Interval (Shaded)
fill([x_fit; flipud(x_fit)], [p_sim - ci_lo; flipud(p_sim + ci_hi)], ...
    prof_colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% Plot Curve
plot(x_fit, p_sim, '-', 'Color', prof_colors(1,:), 'LineWidth', 3, 'DisplayName', leg_names{1});

% --- B. Plot Sequential Sets ---
for ss = 1:nSets_seq
    dat = Seq_Results{ss};
    if isempty(dat), continue; end
    
    col_idx = 1 + ss;
    if col_idx > size(prof_colors,1), col_idx = size(prof_colors,1); end
    col = prof_colors(col_idx, :);
    
    try
        [b_seq, ~, stats_seq] = glmfit(dat(:,1), dat(:,2), 'binomial', 'link', 'logit');
        [p_seq, ci_lo, ci_hi] = glmval(b_seq, x_fit, 'logit', stats_seq);
        
        % Plot Shaded CI
        fill([x_fit; flipud(x_fit)], [p_seq - ci_lo; flipud(p_seq + ci_hi)], ...
            col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        % Plot Curve
        plot(x_fit, p_seq, '-', 'Color', col, 'LineWidth', 3, 'DisplayName', leg_names{col_idx});
        
    catch
        warning('Logistic fit failed for Set %d (likely too few points)', ss);
    end
end

% Aesthetics
ylabel('Probability of Response', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Distance from Stimulation Site (µm)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Spatial Probability (Logistic Fit at %.0f µA)', target_amp), 'FontSize', 14);

yline(0.5, '--k', '50% Prob', 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
ax = gca; ax.Box = 'off';
xlim([0 600]); ylim([0 1.05]);
legend('Location', 'northeast', 'Box', 'off');

%% ============================================================
%   5. SAVE RESULTS (Dedicated Spatial File)
% ============================================================
fprintf('\n--- SAVING RESULTS ---\n');

% 1. Setup Directory
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

parts = split(folder_sim, filesep); 
exp_id = parts{end};
% File Name: Result_Spatial_[ExpID].mat
out_filename = fullfile(save_dir, ['Result_Set1_Spatial_5ms_' exp_id '.mat']);

% 2. Structure Data
ResultSpatial = struct();
ResultSpatial.Metadata.Created = datestr(now);
ResultSpatial.Metadata.Source_Sim = folder_sim;
ResultSpatial.Metadata.Source_Seq = folder_seq;
ResultSpatial.Metadata.StimGroups.Sim = stim_chans_sim;
ResultSpatial.Metadata.StimGroups.Seq = uniqueComb_seq(:, uniqueComb_seq(1,:) > 0);
ResultSpatial.Metadata.TargetAmp = target_amp;

% A. Raw Data (For Logistic Regression)
ResultSpatial.Raw.Sim_Dist_Resp = Sim_Results; % [Distance, IsResponsive]
ResultSpatial.Raw.Seq_Dist_Resp = Seq_Results; % Cell array of [Distance, IsResponsive]

% B. Binned Data (For Bar Chart)
ResultSpatial.Binned.Distances = Unique_Dists;
ResultSpatial.Binned.Percentages = Plot_Data; % [Dists x (1 + nSets)]

% 3. Save
save(out_filename, 'ResultSpatial');
fprintf('Success! Spatial Analysis saved to:\n  %s\n', out_filename);

%% ==================== GEOMETRY MAPPING FUNCTION ====================
function Coords = get_electrode_coordinates(type)
    Coords = zeros(64, 2); 
    if type == 1 
        % Single Shank: Ch1=0, Ch2=50...
        for ch = 1:32, Coords(ch, 1) = 0; Coords(ch, 2) = (ch-1)*50; end
    elseif type == 2
        % Four Shank map
        for i = 1:16, ch=i; Coords(ch,1)=0; Coords(ch,2)=(i-1)*50; end
        for i = 1:16, ch=32+i; Coords(ch,1)=200; Coords(ch,2)=(i-1)*50; end
        for i = 1:16, ch=48+i; Coords(ch,1)=400; Coords(ch,2)=(i-1)*50; end
        for i = 1:16, ch=16+i; Coords(ch,1)=600; Coords(ch,2)=(i-1)*50; end
    end
end