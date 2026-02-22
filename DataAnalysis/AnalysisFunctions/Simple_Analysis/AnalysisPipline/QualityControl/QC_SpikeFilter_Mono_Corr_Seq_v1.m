%% ============================================================
% SpikeFiltering_Evoked_Cleanup.m
% Refined: Start Slope -> Morphology (Width) -> Z-Cross -> Dynamic Correlation
% ============================================================
clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_1';
FS = 30000;                 

% 1. FILTERING PARAMETERS
Start_Slope_Thresh_uV = 100; % Max allowed voltage swing in 0.2ms window

% [MODIFIED] Replaced SSD with Morphology Width limits (in milliseconds)
min_trough_peak_ms    = 0.15; % Too fast = artifact
max_trough_peak_ms    = 1.0;  % Too slow = low frequency noise
corr_thresh           = 0.6;  % Correlation strictness 

% 2. WINDOWS
baseline_window_ms    = [-100 -5]; % Fallback template window
late_template_window  = [2 20];     % Base Evoked window (Will be made dynamic per ISI)
cleanup_window_ms     = [0 100];   % Where to apply the filter
do_corr_filter        = 1;

%% ================= FOLDER & BASE NAME =================
if ~isfolder(data_folder), error('Folder does not exist.'); end
cd(data_folder);
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) > 4, base_name = last_folder(1 : underscores(end-1)-1);
else, base_name = last_folder; end

%% ================= LOAD SPIKES & STIM PARAMS =================
fname_sp = [base_name '.sp_xia_FirstPulse.mat'];
assert(isfile(fname_sp), 'Cannot find %s.', fname_sp);
S_in = load(fname_sp);
if isfield(S_in,'sp_seq')
    sp_in = S_in.sp_seq;
elseif isfield(S_in,'sp')
    sp_in = S_in.sp;
else
    error('No sp or sp_clipped found in %s', fname_sp);
end
nCh = numel(sp_in);
if isempty(dir('*.trig.dat')), cur=pwd; cleanTrig_sabquick; cd(cur); end
trig = loadTrig(0);
nTrials = numel(trig);

% [MODIFIED] Load StimParams to get the ISI per trial for the dynamic window
S_exp = load(dir('*_exp_datafile_*.mat').name,'StimParams','simultaneous_stim');
simN = S_exp.simultaneous_stim;
if simN > 1
    PTD_all_us = cell2mat(S_exp.StimParams(2:end,6)); 
    trial_PTDs_ms = PTD_all_us(2:simN:end) / 1000; 
else
    trial_PTDs_ms = zeros(nTrials, 1);
end

%% ================= 1) START SLOPE FILTER (Windowed Swing Check) =================
% Detects impossible vertical jumps (Blanking Artifacts) in the first ~0.8ms
sp_slope = sp_in;
check_dur_ms = 0.8;      
slide_win_ms = 0.2;      
n_check = round(check_dur_ms * FS / 1000); 
n_slide = round(slide_win_ms * FS / 1000); 
fprintf('\nStart Slope Filtering (Swing > %d uV in 0.2ms window, scanning first %.1fms)...\n', ...
    Start_Slope_Thresh_uV, check_dur_ms);
for ch = 1:nCh
    if isempty(sp_in{ch}), continue; end
    wfs = sp_in{ch}(:,2:end);
    
    limit_idx = min(n_check, size(wfs, 2));
    if limit_idx > n_slide
        early_wfs = wfs(:, 1:limit_idx);
        max_swing = zeros(size(early_wfs,1), 1);
        
        for k = 1 : (limit_idx - n_slide + 1)
            win_data = early_wfs(:, k : k+n_slide-1);
            swing = max(win_data, [], 2) - min(win_data, [], 2);
            max_swing = max(max_swing, swing);
        end
        
        valid_mask = max_swing < Start_Slope_Thresh_uV;
        sp_slope{ch} = sp_in{ch}(valid_mask, :);
        
        removed = sum(~valid_mask);
        if removed > 0
            fprintf('Ch %2d: Removed %d artifact spikes (rapid swing detected).\n', ch, removed);
        end
    else
        sp_slope{ch} = sp_in{ch};
    end
end

%% ================= 2) MORPHOLOGY FILTER (Trough-to-Peak Width) =================
% [MODIFIED] Replaced SSD filter with a biological width enforcer
sp_morph = sp_slope; 
fprintf('\nMorphology Filtering (Width must be %.2f to %.2f ms)...\n', min_trough_peak_ms, max_trough_peak_ms);
for ch = 1:nCh
    if isempty(sp_slope{ch}), continue; end
    wfs = sp_slope{ch}(:,2:end);
    
    % Find indices of absolute minimum and maximum for each spike
    [~, min_idx] = min(wfs, [], 2);
    [~, max_idx] = max(wfs, [], 2);
    
    % Calculate absolute time difference in milliseconds
    width_ms = abs(max_idx - min_idx) / (FS / 1000);
    
    % Keep only spikes with biologically possible widths
    valid_idx = (width_ms >= min_trough_peak_ms) & (width_ms <= max_trough_peak_ms);
    sp_morph{ch} = sp_slope{ch}(valid_idx, :);
    
    removed = sum(~valid_idx);
    if removed > 0
        fprintf('Ch %2d: Removed %d spikes (Impossible Trough-to-Peak width).\n', ch, removed);
    end
end

%% ================= 3) ZERO-CROSSING CHECK =================
% [MODIFIED] Changed input from sp_SSD to sp_morph
sp_zcross = sp_morph;
fprintf('\nZero-Crossing Check (Removing purely negative waveforms)...\n');
for ch = 1:nCh
    if isempty(sp_zcross{ch}), continue; end
    wfs = sp_zcross{ch}(:,2:end);
    
    has_positive_phase = max(wfs, [], 2) > 0;
    
    removed = sum(~has_positive_phase);
    if removed > 0
        sp_zcross{ch} = sp_zcross{ch}(has_positive_phase, :);
        fprintf('Ch %2d: Removed %d spikes (No positive phase).\n', ch, removed);
    end
end

%% ================= 4) TEMPLATE CORRELATION (Dynamic & Hierarchical) =================
sp_corr = sp_zcross; 
if do_corr_filter
    fprintf('\nCorrelation Filtering (Dynamic Evoked Template with Baseline Fallback)...\n');
    
    for ch = 1:nCh
        if isempty(sp_corr{ch}), continue; end
        spt = sp_corr{ch}(:,1);     
        wfs = sp_corr{ch}(:,2:end); 
        nSpikes = size(wfs,1);
        
        % --- A. BUILD DYNAMIC TEMPLATE ---
        mask_evoked = false(nSpikes,1);
        mask_base   = false(nSpikes,1);
        
        for tr = 1:nTrials
            t0 = trig(tr)/FS*1000;
            isi = trial_PTDs_ms(tr);
            
            % [MODIFIED] Calculate safe dynamic window boundary to avoid Pulse 2
            if isi > 0
                dyn_max = min(late_template_window(2), isi - 0.5);
            else
                dyn_max = late_template_window(2);
            end
            
            % Collect Evoked Spikes (Only if window is valid)
            if dyn_max > late_template_window(1)
                in_evoked = (spt >= t0 + late_template_window(1) & spt <= t0 + dyn_max);
                mask_evoked = mask_evoked | in_evoked;
            end
            
            % Collect Baseline Spikes
            in_base = (spt >= t0 + baseline_window_ms(1) & spt <= t0 + baseline_window_ms(2));
            mask_base = mask_base | in_base;
        end
        
        wfs_evoked = wfs(mask_evoked,:);
        wfs_base   = wfs(mask_base,:);
        
        % [MODIFIED] Strict Hierarchy: Use pure evoked first, fallback to baseline only if empty
        if size(wfs_evoked, 1) >= 5
            template = mean(wfs_evoked, 1);
            source_str = 'Dynamic_Evoked';
        elseif size(wfs_base, 1) >= 5
            template = mean(wfs_base, 1);
            source_str = 'Baseline_Fallback';
        else
            fprintf('Ch %2d: Not enough safe spikes for template. Skipping.\n', ch);
            continue; 
        end
        
        % --- B. SCORE & FILTER ---
        corr_vals = zeros(nSpikes,1);
        for i = 1:nSpikes
            corr_vals(i) = corr(template(:), wfs(i,:)');
        end
        
        bad_mask = false(nSpikes,1);
        for tr = 1:nTrials
            t0 = trig(tr)/FS*1000;
            in_window = (spt >= t0 + cleanup_window_ms(1) & spt <= t0 + cleanup_window_ms(2));
            
            if any(in_window) 
                bad_mask(in_window & (corr_vals < corr_thresh)) = true;
            end
        end
        
        removed_count = sum(bad_mask);
        sp_corr{ch} = sp_corr{ch}(~bad_mask, :);
        
        fprintf('Ch %2d (%s Temp): Removed %d artifact spikes.\n', ...
            ch, source_str, removed_count);
    end
end

%% ================= SAVE OUTPUT =================
QC_params = struct();
QC_params.Start_Slope_Thresh = Start_Slope_Thresh_uV;
% [MODIFIED] Replaced SSD save params with Morphology params
QC_params.min_trough_peak_ms = min_trough_peak_ms;
QC_params.max_trough_peak_ms = max_trough_peak_ms;
QC_params.baseline_window    = baseline_window_ms;
QC_params.late_template_win  = late_template_window;
QC_params.cleanup_window     = cleanup_window_ms;
QC_params.corr_thresh        = corr_thresh;

out_name = [base_name '.sp_xia_SSD.mat'];
if isfile(out_name), delete(out_name); end 

% [MODIFIED] Replaced sp_SSD with sp_morph in the saved output
save(out_name, 'sp_in','sp_slope','sp_morph','sp_zcross','sp_corr','QC_params','-v7.3');
fprintf('\nSaved cleaned spikes to %s\n', out_name);