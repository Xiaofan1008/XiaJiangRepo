%% ============================================================
%   Spatial Distribution Analysis (with Smooth Curve Fit)
%   - Metric: % of Total Responders located at each Distance
%   - Plot: Grouped Bars + Smooth Trend Lines (Spline Fit)
%   - Output: Saves comprehensive results to 'Result_SpatialDist_... .mat'
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';
Electrode_Type = 1; 
target_amp = 10; 
target_ptd_seq = 5; 
% Distance Bins
dist_edges = 0:50:800; 
dist_centers = dist_edges(1:end-1) + 25; 

%% =================== LOAD DATA ====================
[Rsim, ~, ~, Ssim, QC_Sim] = load_experiment_data(folder_sim);
[Rseq, ~, ~, Sseq, QC_Seq] = load_experiment_data(folder_seq);

amps_all_sim=cell2mat(Ssim.StimParams(2:end,16)); 
trialAmps_sim=amps_all_sim(1:Ssim.simultaneous_stim:end); 
[Amps_sim,~,ampIdx_sim]=unique(trialAmps_sim); Amps_sim(Amps_sim==-1)=0;

simN_seq=Sseq.simultaneous_stim; 
amps_all_seq=cell2mat(Sseq.StimParams(2:end,16)); 
trialAmps_seq=amps_all_seq(1:simN_seq:end); 
[Amps_seq,~,ampIdx_seq]=unique(trialAmps_seq); Amps_seq(Amps_seq==-1)=0;

PTD_all_us=cell2mat(Sseq.StimParams(3:simN_seq:end,6)); PTDs_ms=unique(PTD_all_us/1000); 
stimNames=Sseq.StimParams(2:end,1); [~,idx_all]=ismember(stimNames,Sseq.E_MAP(2:end)); 
if isfield(Sseq,'n_Trials'), nTr_seq=Sseq.n_Trials; else, nTr_seq=(size(Sseq.StimParams,1)-1)/simN_seq; end; 
comb_seq=zeros(nTr_seq,simN_seq); for t=1:nTr_seq, rr=(t-1)*simN_seq+(1:simN_seq); v=idx_all(rr); v=v(v>0); comb_seq(t,1:numel(v))=v(:).'; end; 
[uniqueComb_seq,~,combClass_seq]=unique(comb_seq,'rows','stable'); nSets_seq=size(uniqueComb_seq,1);

%% ================= GEOMETRY & TARGETS =================
Coords = get_electrode_coordinates(Electrode_Type);
if Electrode_Type == 1, valid_channels = 1:32; else, valid_channels = [1:16, 33:48, 49:64, 17:32]; end
idx_sim_amp = find(Amps_sim == target_amp, 1);
idx_seq_amp = find(Amps_seq == target_amp, 1);
idx_ptd = find(abs(PTDs_ms - target_ptd_seq) < 0.1, 1); if isempty(idx_ptd), idx_ptd = 1; end
stim_chans_sim = Rsim.set(1).stimChannels;

%% ================= STEP 1: CALCULATE RAW DISTANCES =================
% A. SIMULTANEOUS
stim_locs_sim  = Coords(stim_chans_sim, :); 
Sim_Results = []; 
chan_structs = Rsim.set(1).amp(idx_sim_amp).ptd(1).channel;
bad_ch_sim = []; if ~isempty(QC_Sim.BadCh), bad_ch_sim = QC_Sim.BadCh{1}; end

for ch = valid_channels
    if ismember(ch, bad_ch_sim), continue; end
    rec_loc = Coords(ch, :);
    dists = sqrt(sum((stim_locs_sim - rec_loc).^2, 2));
    min_dist = round(min(dists));
    if min_dist < 1, continue; end 
    
    is_resp = 0;
    try if chan_structs(ch).is_responsive, is_resp=1; end; catch; end
    Sim_Results = [Sim_Results; min_dist, is_resp];
end

% B. SEQUENTIAL
Seq_Results = cell(nSets_seq, 1);
for ss = 1:nSets_seq
    stim_chans_seq = uniqueComb_seq(ss, :); stim_chans_seq = stim_chans_seq(stim_chans_seq>0);
    stim_locs_seq  = Coords(stim_chans_seq, :);
    chan_structs = Rseq.set(ss).amp(idx_seq_amp).ptd(idx_ptd).channel;
    bad_ch_seq = []; if ~isempty(QC_Seq.BadCh) && ss<=numel(QC_Seq.BadCh), bad_ch_seq = QC_Seq.BadCh{ss}; end
    
    for ch = valid_channels
        if ismember(ch, bad_ch_seq), continue; end
        rec_loc = Coords(ch, :);
        dists = sqrt(sum((stim_locs_seq - rec_loc).^2, 2));
        min_dist = round(min(dists));
        if min_dist < 1, continue; end
        
        is_resp = 0;
        try if chan_structs(ch).is_responsive, is_resp=1; end; catch; end
        Seq_Results{ss} = [Seq_Results{ss}; min_dist, is_resp];
    end
end

%% ================= STEP 2: CALCULATE DISTRIBUTION % =================
% Combine all distances to find unique bins
all_dists = Sim_Results(:, 1);
for ss=1:nSets_seq, all_dists = [all_dists; Seq_Results{ss}(:,1)]; end
Unique_Dists = unique(all_dists);
Unique_Dists = Unique_Dists(Unique_Dists <= 600); 

Plot_Data = zeros(length(Unique_Dists), 1 + nSets_seq);

% Sim
Total_Sim = sum(Sim_Results(:,2));
for i = 1:length(Unique_Dists)
    d = Unique_Dists(i);
    mask = (Sim_Results(:,1) == d);
    if Total_Sim > 0, Plot_Data(i, 1) = sum(Sim_Results(mask, 2)) / Total_Sim * 100; end
end

% Seq
for ss = 1:nSets_seq
    dat = Seq_Results{ss};
    Total_Seq = sum(dat(:,2));
    for i = 1:length(Unique_Dists)
        d = Unique_Dists(i);
        mask = (dat(:,1) == d);
        if Total_Seq > 0, Plot_Data(i, 1+ss) = sum(dat(mask, 2)) / Total_Seq * 100; end
    end
end

prof_colors = [0.25 0.45 0.65; 0.80 0.45, 0.30; 0.55, 0.40, 0.65];

%% ================= STEP 3: PLOT BARS + SMOOTH CURVES =================
figure('Color','w', 'Position',[200 600 800 500]); hold on;

% 1. Plot Bars
b = bar(Unique_Dists, Plot_Data, 'grouped');
for i = 1:length(b)
    if i <= size(prof_colors,1), b(i).FaceColor = prof_colors(i,:); end
    b(i).EdgeColor = 'none'; b(i).FaceAlpha = 0.6; % Lower alpha to see lines
end

% 2. Plot Smooth Curves (Spline Fit)
x_fine = linspace(min(Unique_Dists), max(Unique_Dists), 200);

% Curve for Sim
y_sim = Plot_Data(:,1);
if length(Unique_Dists) > 3
    yy_sim = interp1(Unique_Dists, y_sim, x_fine, 'pchip'); 
    plot(x_fine, yy_sim, 'Color', prof_colors(1,:), 'LineWidth', 3, 'DisplayName', 'Sim Trend');
end

% Curves for Seq
for ss = 1:nSets_seq
    y_seq = Plot_Data(:, 1+ss);
    col_idx = 1+ss; 
    if length(Unique_Dists) > 3
        yy_seq = interp1(Unique_Dists, y_seq, x_fine, 'pchip');
        if col_idx<=size(prof_colors,1), col=prof_colors(col_idx,:); else, col='k'; end
        plot(x_fine, yy_seq, 'Color', col, 'LineWidth', 3, 'DisplayName', sprintf('Seq Set %d', ss));
    end
end

ylabel('% of Total Responding Population', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Distance from Stimulation Site (µm)', 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('Spatial Distribution (at %.0f µA)', target_amp), 'FontSize', 14);
ax = gca; ax.Box = 'off'; xticks(Unique_Dists); xtickangle(0);
legend('Location', 'northeast', 'Box', 'off');

%% ============================================================
%   5. SAVE RESULTS (Comprehensive)
% ============================================================
fprintf('\n--- SAVING RESULTS ---\n');
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

parts = split(folder_sim, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_Set1_SpatialDist_5ms_' exp_id '.mat']);

ResultSpatial = struct();
ResultSpatial.Metadata.Created = datestr(now);
ResultSpatial.Metadata.Source_Sim = folder_sim;
ResultSpatial.Metadata.Source_Seq = folder_seq;

% 1. Parameters
ResultSpatial.Parameters.ElectrodeType = Electrode_Type;
ResultSpatial.Parameters.TargetAmp = target_amp;
ResultSpatial.Parameters.TargetPTD = target_ptd_seq;
ResultSpatial.Parameters.DistEdges = dist_edges;

% 2. Geometry
ResultSpatial.Geometry.StimGroups.Sim = stim_chans_sim;
ResultSpatial.Geometry.StimGroups.Seq = uniqueComb_seq(:, uniqueComb_seq(1,:) > 0);
ResultSpatial.Geometry.Coords = Coords;

% 3. Raw Data (For re-analysis)
% Sim_Results: [Distance, IsResponsive]
% Seq_Results: Cell Array of [Distance, IsResponsive]
ResultSpatial.Raw.Sim_Dist_Resp = Sim_Results; 
ResultSpatial.Raw.Seq_Dist_Resp = Seq_Results; 
ResultSpatial.Raw.Amps_Sim = Amps_sim;
ResultSpatial.Raw.Amps_Seq = Amps_seq;

% 4. Processed Data (The Distribution)
ResultSpatial.Distribution.Distances = Unique_Dists;
ResultSpatial.Distribution.Percentages = Plot_Data; % [Dists x (1+nSets)]
ResultSpatial.Distribution.TotalResponders.Sim = Total_Sim;

% 5. QC Data
ResultSpatial.QC.Sim = QC_Sim; 
ResultSpatial.QC.Seq = QC_Seq;

save(out_filename, 'ResultSpatial');
fprintf('Success! Spatial Distribution saved to:\n  %s\n', out_filename);

%% ==================== HELPERS ====================
function Coords = get_electrode_coordinates(type)
    Coords = zeros(64, 2); 
    if type == 1 
        for ch = 1:32, Coords(ch, 1) = 0; Coords(ch, 2) = (ch-1)*50; end
    elseif type == 2
        for i = 1:16, ch=i; Coords(ch,1)=0; Coords(ch,2)=(i-1)*50; end
        for i = 1:16, ch=32+i; Coords(ch,1)=200; Coords(ch,2)=(i-1)*50; end
        for i = 1:16, ch=48+i; Coords(ch,1)=400; Coords(ch,2)=(i-1)*50; end
        for i = 1:16, ch=16+i; Coords(ch,1)=600; Coords(ch,2)=(i-1)*50; end
    end
end
function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder);
    f = dir('*RespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;
    S = load(dir('*_exp_datafile_*.mat').name);
    sp=[]; trig=[];
    QC.BadCh = []; QC.BadTrials = [];
    f_bc = dir('*.BadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
end