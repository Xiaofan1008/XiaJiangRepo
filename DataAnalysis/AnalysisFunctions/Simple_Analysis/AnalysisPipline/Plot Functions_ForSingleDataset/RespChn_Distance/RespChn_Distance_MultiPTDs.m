%% ============================================================
%   SPATIAL ANALYSIS: Single Dataset (Response Prob vs Distance)
%   - Logic: 
%       1. Determine Electrode Geometry (Single vs 4-Shank).
%       2. For each Set: Identify Stimulation Electrodes.
%       3. For each Ch: Calculate MIN DISTANCE to nearest Stim Ch.
%       4. Analyze: Is the channel Responsive? (1/0).
%       5. Binning: Calculate Prob(Response) per distance bin.
%       6. Save Raw Data & Plot.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Seq_Sim7';
Electrode_Type = 2; % 1 = Single Shank (32ch), 2 = 4-Shank (64ch)

% 1. Analysis Parameters
dist_bin_edges = 0 : 50 : 1000; % Distance bins (microns) e.g. 0-100, 100-200...
Seq_PTD_Target = 5; % The sequential delay to compare (ms)

% 2. Plotting
save_figure = false;
save_dir    = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/';

%% ================= 1. LOAD DATA =================
fprintf('Loading data...\n');
[R, sp, trig, S, QC] = load_experiment_data(data_folder);

% --- Extract Stim Params ---
Stim = S.StimParams; 
simN = S.simultaneous_stim; 
if isfield(S, 'n_Trials'), nTr = S.n_Trials; else, nTr = (size(Stim, 1) - 1) / simN; end

% Amplitudes
amps_all  = cell2mat(Stim(2:end,16)); trialAmps = amps_all(1:simN:end);
[Amps,~,ampIdx] = unique(trialAmps); Amps(Amps==-1) = 0;

% PTDs
if simN > 1, PTD_all = cell2mat(Stim(3:simN:end,6)); else, PTD_all = zeros(nTr,1); end
PTD_ms = PTD_all / 1000; [PTDs,~,ptdIdx] = unique(PTD_ms);

% Sets
stimNames = Stim(2:end,1); [~, idx_all] = ismember(stimNames, S.E_MAP(2:end));
comb = zeros(nTr, simN);
for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
[uniqueComb,~,combClass] = unique(comb,'rows','stable'); 
nSets = size(uniqueComb,1);

% Identify Sim and Seq Indices
ptd_sim_idx = find(abs(PTDs - 0) < 0.001);
ptd_seq_idx = find(abs(PTDs - Seq_PTD_Target) < 0.001);

if isempty(ptd_sim_idx), warning('No Sim (0ms) data found'); end
if isempty(ptd_seq_idx), warning('No Seq (%.1fms) data found', Seq_PTD_Target); end

d = Depth_s(Electrode_Type); % Mapping from Depth_idx -> Intan_idx
nCh_Total = length(d);

%% ================= 2. PROCESS SPATIAL DATA =================
% Output Structure: SpatialData.Set(s).Amp(a).Sim/Seq
% Stores raw vectors: [Distance, IsResp, IsBad]
SpatialData = struct();

fprintf('Processing Spatial Data (%d Sets)...\n', nSets);

for ss = 1:nSets
    % 1. Identify Stimulation Electrodes for this Set
    stim_ch_idx = uniqueComb(ss,:); 
    stim_ch_idx = stim_ch_idx(stim_ch_idx > 0); % These are Depth Indices (1-64)
    
    SpatialData.Set(ss).StimCh = stim_ch_idx;
    
    % 2. Calculate Distances for ALL channels
    % dist_vector(i) = distance from Channel i to nearest Stim Electrode
    dist_vector = nan(nCh_Total, 1);
    
    % Get Coords of Stim Electrodes
    stim_coords = [];
    for k = 1:length(stim_ch_idx)
        stim_coords(k,:) = get_electrode_coords(stim_ch_idx(k), Electrode_Type);
    end
    
    % Calc Distance for every recording channel
    for ch = 1:nCh_Total
        ch_coords = get_electrode_coords(ch, Electrode_Type);
        
        % Euclidean distance to all stim electrodes
        dists = sqrt(sum((stim_coords - ch_coords).^2, 2));
        
        % Take Minimum Distance (Nearest Source)
        dist_vector(ch) = min(dists);
    end
    
    % 3. Extract Response Status per Amplitude
    for ai = 1:length(Amps)
        curr_amp = Amps(ai);
        
        % --- SIMULTANEOUS ---
        resp_vec_sim = zeros(nCh_Total, 1); 
        valid_vec_sim = true(nCh_Total, 1); % To exclude Bad Channels
        
        if ~isempty(ptd_sim_idx)
            try
                chn_list = R.set(ss).amp(ai).ptd(ptd_sim_idx).channel;
                for c = 1:min(length(chn_list), nCh_Total)
                    % Check Bad Channel
                    if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(c, QC.BadCh{ss})
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
        
        % --- SEQUENTIAL ---
        resp_vec_seq = zeros(nCh_Total, 1);
        valid_vec_seq = true(nCh_Total, 1);
        
        if ~isempty(ptd_seq_idx)
            try
                chn_list = R.set(ss).amp(ai).ptd(ptd_seq_idx).channel;
                for c = 1:min(length(chn_list), nCh_Total)
                    if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(c, QC.BadCh{ss})
                        valid_vec_seq(c) = false;
                        continue;
                    end
                    if isfield(chn_list(c), 'is_responsive') && chn_list(c).is_responsive
                        resp_vec_seq(c) = 1;
                    end
                end
            catch; end
        end
        
        % --- STORE RAW DATA ---
        SpatialData.Set(ss).Amp(ai).Val = curr_amp;
        
        % Store as Table-like struct columns
        % Rows = Channels 1..N
        SpatialData.Set(ss).Amp(ai).Dist = dist_vector;
        
        SpatialData.Set(ss).Amp(ai).Sim_Resp = resp_vec_sim;
        SpatialData.Set(ss).Amp(ai).Sim_Valid = valid_vec_sim;
        
        SpatialData.Set(ss).Amp(ai).Seq_Resp = resp_vec_seq;
        SpatialData.Set(ss).Amp(ai).Seq_Valid = valid_vec_seq;
    end
end


%% ================= 3. PLOT (Per Amp, Per Set) =================
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
        
        % Add 90% Threshold line (Useful for comparing radius)
        yline(0.9, '--k', '90% Recruitment', 'HandleVisibility','off');
        
        xlabel('Distance (\mum)'); ylabel('Cumulative Fraction');
        ylim([0 1.05]); title('Cumulative Spatial Spread'); legend('Location','best');
        
        sgtitle(figTitle);
        
        if save_figure
            fName = sprintf('SpatialCDF_Set%d_Amp%.1f.png', ss, curr_amp);
            saveas(gcf, fullfile(save_dir, fName));
        end
        % close(gcf); 
    end
end
fprintf('>>> All Plots Generated.\n');

%% ================= 4. SAVE RESULTS =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_SpatialRaw_' exp_id '.mat']);
save(out_filename, 'SpatialData', 'Amps', 'dist_bin_edges', 'Electrode_Type');
fprintf('>>> Spatial Data Saved: %s\n', out_filename);

%% ================= HELPER: GEOMETRY MAP =================
function coords = get_electrode_coords(ch_idx, type)
    % Returns [x, y] in microns
    % ch_idx: 1-based index (Depth Index)
    
    if type == 1 
        % --- Single Shank (1x32) ---
        % All in one column
        x = 0;
        y = (ch_idx - 1) * 50; 
        
    elseif type == 2
        % --- 4-Shank (4x16) ---
        % Shank Pitch = 200um, Electrode Pitch = 50um
        % Mapping based on user description:
        % Shank 1 (X=0):   Ch 1-16
        % Shank 2 (X=200): Ch 33-48
        % Shank 3 (X=400): Ch 49-64
        % Shank 4 (X=600): Ch 17-32
        
        if ch_idx >= 1 && ch_idx <= 16
            shank_x = 0;
            local_idx = ch_idx;
        elseif ch_idx >= 33 && ch_idx <= 48
            shank_x = 200;
            local_idx = ch_idx - 32;
        elseif ch_idx >= 49 && ch_idx <= 64
            shank_x = 400;
            local_idx = ch_idx - 48;
        elseif ch_idx >= 17 && ch_idx <= 32
            shank_x = 600;
            local_idx = ch_idx - 16;
        else
            shank_x = 0; local_idx = 1; % Fallback
        end
        
        % Y-coordinate: Assuming linear mapping 1->0, 2->50... 
        % (Adjust if your map is flipped)
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