%% ============================================================
%   SPATIAL ANALYSIS: Single Dataset (Separate Sim/Seq Files)
%   - Metric: Spatial Spread (Cumulative Distribution Function)
%   - Logic: 
%       1. Load Sim and Seq datasets separately.
%       2. Determine Electrode Geometry (Single vs 4-Shank).
%       3. For each Ch: Calculate MIN DISTANCE to nearest Stim Ch.
%       4. Analyze: Is the channel Responsive? (1/0).
%          -> Sim data from Sim File (PTD=0)
%          -> Seq data from Seq File (PTD=Target)
%       5. Plot: Histogram (Count) & CDF (Cumulative Fraction).
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
% Define separate folders for Sim and Seq
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Sim'; % Sim Path
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Seq'; % Seq path

Electrode_Type = 2; % 1 = Single Shank (32ch), 2 = 4-Shank (64ch)

% 1. Analysis Parameters
dist_bin_edges = 0 : 50 : 1300; % Distance bins (microns)
Seq_PTD_Target = 5.5;             % The sequential delay to compare (ms)

% 2. Plotting
save_figure = false;
save_dir    = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX005/';

%% ================= 1. LOAD DATA (DUAL) =================
fprintf('Loading Sim Data...\n');
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);

fprintf('Loading Seq Data...\n');
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

% --- Extract Stim Params (Use Sim as Master) ---
Stim_sim = Ssim.StimParams; 
simN_sim = Ssim.simultaneous_stim; 
if isfield(Ssim, 'n_Trials'), nTr_sim = Ssim.n_Trials; else, nTr_sim = (size(Stim_sim, 1) - 1) / simN_sim; end

% Amplitudes (from Sim file)
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps,~,ampIdx_sim] = unique(trialAmps_sim); Amps(Amps==-1) = 0;

% PTDs (Sim File)
if simN_sim > 1, PTD_all_sim = cell2mat(Stim_sim(3:simN_sim:end,6)); else, PTD_all_sim = zeros(nTr_sim,1); end
PTD_ms_sim = PTD_all_sim / 1000; [PTDs_sim,~,ptdIdx_sim] = unique(PTD_ms_sim);

% PTDs (Seq File)
Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim;
if isfield(Sseq, 'n_Trials'), nTr_seq = Sseq.n_Trials; else, nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq; end
if simN_seq > 1, PTD_all_seq = cell2mat(Stim_seq(3:simN_seq:end,6)); else, PTD_all_seq = zeros(nTr_seq,1); end
PTD_ms_seq = PTD_all_seq / 1000; [PTDs_seq,~,ptdIdx_seq] = unique(PTD_ms_seq);

% Sets (Use Sim geometry as reference)
stimNames_sim = Stim_sim(2:end,1); [~, idx_all_sim] = ismember(stimNames_sim, Ssim.E_MAP(2:end));
comb_sim = zeros(nTr_sim, simN_sim);
for t = 1:nTr_sim, rr = (t-1)*simN_sim + (1:simN_sim); v = idx_all_sim(rr); v = v(v>0); comb_sim(t,1:numel(v)) = v(:).'; end
[uniqueComb_sim,~,combClass_sim] = unique(comb_sim,'rows','stable'); 
nSets = size(uniqueComb_sim,1);

% Identify Indices
% We look for 0ms in the SIM file
ptd_sim_idx = find(abs(PTDs_sim - 0) < 0.001);
% We look for Target (e.g. 5ms) in the SEQ file
ptd_seq_idx = find(abs(PTDs_seq - Seq_PTD_Target) < 0.001);

if isempty(ptd_sim_idx), warning('No Sim (0ms) data found in Sim File'); end
if isempty(ptd_seq_idx), warning('No Seq (%.1fms) data found in Seq File', Seq_PTD_Target); end

d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

%% ================= 2. PROCESS SPATIAL DATA =================
% Output Structure: SpatialData.Set(s).Amp(a).Sim/Seq
% Stores raw vectors: [Distance, IsResp, IsBad]
SpatialData = struct();

fprintf('Processing Spatial Data (%d Sets)...\n', nSets);

for ss = 1:nSets
    % 1. Identify Stimulation Electrodes (Using Sim Geometry)
    stim_ch_idx = uniqueComb_sim(ss,:); 
    stim_ch_idx = stim_ch_idx(stim_ch_idx > 0); 
    
    SpatialData.Set(ss).StimCh = stim_ch_idx;
    
    % 2. Calculate Distances for ALL channels
    dist_vector = nan(nCh_Total, 1);
    
    % Get Coords of Stim Electrodes
    stim_coords = [];
    for k = 1:length(stim_ch_idx)
        stim_coords(k,:) = get_electrode_coords(stim_ch_idx(k), Electrode_Type);
    end
    
    % Calc Distance for every recording channel
    for ch = 1:nCh_Total
        ch_coords = get_electrode_coords(ch, Electrode_Type);
        % Euclidean distance to nearest stim electrode
        dists = sqrt(sum((stim_coords - ch_coords).^2, 2));
        dist_vector(ch) = min(dists);
    end
    
    % 3. Extract Response Status per Amplitude
    for ai = 1:length(Amps)
        curr_amp = Amps(ai);
        
        % --- A. SIMULTANEOUS (From Sim File) ---
        resp_vec_sim = zeros(nCh_Total, 1); 
        valid_vec_sim = true(nCh_Total, 1); 
        
        if ~isempty(ptd_sim_idx)
            try
                % Access Rsim
                chn_list = Rsim.set(ss).amp(ai).ptd(ptd_sim_idx).channel;
                for c = 1:min(length(chn_list), nCh_Total)
                    % Check Bad Channel (QC_Sim)
                    if ~isempty(QC_Sim.BadCh) && ss <= length(QC_Sim.BadCh) && ismember(c, QC_Sim.BadCh{ss})
                        valid_vec_sim(c) = false;
                        continue;
                    end
                    % Check Responsiveness
                    if isfield(chn_list(c), 'is_responsive') && chn_list(c).is_responsive
                        resp_vec_sim(c) = 1;
                    end
                end
            catch; end
        end
        
        % --- B. SEQUENTIAL (From Seq File) ---
        resp_vec_seq = zeros(nCh_Total, 1);
        valid_vec_seq = true(nCh_Total, 1);
        
        if ~isempty(ptd_seq_idx)
            try
                % Access Rseq
                chn_list = Rseq.set(ss).amp(ai).ptd(ptd_seq_idx).channel;
                for c = 1:min(length(chn_list), nCh_Total)
                    % Check Bad Channel (QC_Seq)
                    if ~isempty(QC_Seq.BadCh) && ss <= length(QC_Seq.BadCh) && ismember(c, QC_Seq.BadCh{ss})
                        valid_vec_seq(c) = false;
                        continue;
                    end
                    % Check Responsiveness
                    if isfield(chn_list(c), 'is_responsive') && chn_list(c).is_responsive
                        resp_vec_seq(c) = 1;
                    end
                end
            catch; end
        end
        
        % --- STORE RAW DATA ---
        SpatialData.Set(ss).Amp(ai).Val = curr_amp;
        SpatialData.Set(ss).Amp(ai).Dist = dist_vector;
        
        SpatialData.Set(ss).Amp(ai).Sim_Resp = resp_vec_sim;
        SpatialData.Set(ss).Amp(ai).Sim_Valid = valid_vec_sim;
        
        SpatialData.Set(ss).Amp(ai).Seq_Resp = resp_vec_seq;
        SpatialData.Set(ss).Amp(ai).Seq_Valid = valid_vec_seq;
    end
end

%% ================= 3. PLOT (Histogram + CDF) =================
bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;

for ss = 1:nSets
    stimCh = SpatialData.Set(ss).StimCh;
    
    for ai = 1:length(Amps)
        Data = SpatialData.Set(ss).Amp(ai);
        curr_amp = Data.Val;
        if curr_amp == 0, continue; end
        
        % Initialize Counts
        Count_Sim = zeros(size(bin_centers));
        Count_Seq = zeros(size(bin_centers));
        
        for b = 1:length(dist_bin_edges)-1
            % Find channels in this distance bin
            mask = Data.Dist >= dist_bin_edges(b) & Data.Dist < dist_bin_edges(b+1);
            
            % Sim Counts
            valid_mask = mask & Data.Sim_Valid;
            if sum(valid_mask) > 0
                Count_Sim(b) = sum(Data.Sim_Resp(valid_mask));
            end
            
            % Seq Counts
            valid_mask = mask & Data.Seq_Valid;
            if sum(valid_mask) > 0
                Count_Seq(b) = sum(Data.Seq_Resp(valid_mask));
            end
        end
        
        % --- Calculate CDF (Cumulative Distribution) ---
        CDF_Sim = cumsum(Count_Sim);
        if max(CDF_Sim) > 0, CDF_Sim = CDF_Sim / max(CDF_Sim); end % Normalize 0-1
        
        CDF_Seq = cumsum(Count_Seq);
        if max(CDF_Seq) > 0, CDF_Seq = CDF_Seq / max(CDF_Seq); end % Normalize 0-1
        
        % --- Plotting ---
        figTitle = sprintf('Spatial Spread - Set %d (Ch%s) - %.1fuA', ss, num2str(stimCh), curr_amp);
        figure('Color','w', 'Position', [100 100 1000 400], 'Name', figTitle);
        
        % Subplot 1: Raw Count Histogram
        subplot(1,2,1); hold on;
        b = bar(bin_centers, [Count_Sim; Count_Seq]', 'grouped');
        b(1).FaceColor = 'b'; b(1).EdgeColor = 'none'; b(1).DisplayName = 'Simultaneous';
        b(2).FaceColor = 'r'; b(2).EdgeColor = 'none'; b(2).DisplayName = 'Sequential';
        xlabel('Distance (\mum)'); ylabel('# Responding Channels');
        title('Count Histogram'); legend('Location','best');
        
        % Subplot 2: CDF Curve (Spatial Spread)
        subplot(1,2,2); hold on;
        plot(bin_centers, CDF_Sim, '-o', 'Color', 'b', 'LineWidth', 2, 'MarkerFaceColor','b', 'DisplayName','Sim CDF');
        plot(bin_centers, CDF_Seq, '-s', 'Color', 'r', 'LineWidth', 2, 'MarkerFaceColor','r', 'DisplayName','Seq CDF');
        
        % Add 90% Threshold line
        yline(0.9, '--k', '90% Recruitment', 'HandleVisibility','off');
        
        xlabel('Distance (\mum)'); ylabel('Cumulative Fraction');
        ylim([0 1.05]); title('Cumulative Spatial Spread'); legend('Location','best');
        
        sgtitle(figTitle);
        
        if save_figure
            fName = sprintf('SpatialCDF_Set%d_Amp%.1f.png', ss, curr_amp);
            saveas(gcf, fullfile(save_dir, fName));
        end
        % close(gcf); % Uncomment if generating many plots
    end
end
fprintf('>>> All Plots Generated.\n');

%% ================= 4. SAVE RESULTS =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end}; % Use Sim ID
out_filename = fullfile(save_dir, ['Result_SpatialRaw_' exp_id '.mat']);
save(out_filename, 'SpatialData', 'Amps', 'dist_bin_edges', 'Electrode_Type');
fprintf('>>> Spatial Data Saved: %s\n', out_filename);

%% ================= HELPER: GEOMETRY MAP =================
function coords = get_electrode_coords(ch_idx, type)
    % Returns [x, y] in microns
    % ch_idx: 1-based index (Depth Index)
    
    if type == 1 
        % --- Single Shank (1x32) ---
        x = 0;
        y = (ch_idx - 1) * 50; 
        
    elseif type == 2
        % --- 4-Shank (4x16) ---
        % Mapping: Shank 1(0), 2(200), 3(400), 4(600)
        % Note: Ch Mapping provided by user description
        if ch_idx >= 1 && ch_idx <= 16
            shank_x = 0;   local_idx = ch_idx;
        elseif ch_idx >= 33 && ch_idx <= 48
            shank_x = 200; local_idx = ch_idx - 32;
        elseif ch_idx >= 49 && ch_idx <= 64
            shank_x = 400; local_idx = ch_idx - 48;
        elseif ch_idx >= 17 && ch_idx <= 32
            shank_x = 600; local_idx = ch_idx - 16;
        else
            shank_x = 0; local_idx = 1;
        end
        y = (local_idx - 1) * 50;
        x = shank_x;
    end
    coords = [x, y];
end

function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder);
    f = dir('*RespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;
    f = dir('*sp_xia_SSD.mat'); if isempty(f), f=dir('*sp_xia.mat'); end
    if isempty(f), error('No Spike file in %s', folder); end
    S_sp = load(f(1).name);
    if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr; elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD; else, sp = S_sp.sp_in; end
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig = loadTrig(0);
    S = load(dir('*_exp_datafile_*.mat').name);
    QC.BadCh = []; QC.BadTrials = [];
    f_bc = dir('*.BadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
    f_bt = dir('*.BadTrials.mat'); if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end