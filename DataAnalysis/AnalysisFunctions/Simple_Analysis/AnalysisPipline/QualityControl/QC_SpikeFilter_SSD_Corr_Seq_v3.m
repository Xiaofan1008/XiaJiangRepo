%% ============================================================
% SpikeFiltering_Evoked_Cleanup.m
% Refined: Start Slope Filter -> SSD (12) -> Z-Cross -> Dual-Window Correlation
% ============================================================
clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Seq_Sim7';
FS = 30000;                 

% 1. FILTERING PARAMETERS
Start_Slope_Thresh_uV = 100; % Max allowed voltage swing in 0.2ms window
SSD_threshold_factor  = 12;  % SSD Limit(was 16)
corr_thresh           = 0.6; % Correlation strictness 

% 2. WINDOWS
baseline_window_ms    = [-100 -5]; % Window 1 for template (Pre-stim)
late_template_window  = [0 0];  % Window 2 for template (Post-recovery)
cleanup_window_ms     = [0 100];  % Where to apply the filter

do_SSD_filter         = 1;
do_corr_filter        = 1;

%% ================= FOLDER & BASE NAME =================
if ~isfolder(data_folder), error('Folder does not exist.'); end
cd(data_folder);
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4, base_name = last_folder(1 : underscores(end-1)-1);
else, base_name = last_folder; end

%% ================= LOAD SPIKES =================
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

%% ================= 2) SSD-BASED FILTERING (Pre-Clean) =================
sp_SSD = sp_slope; 
if do_SSD_filter
    fprintf('\nSSD Filtering (Threshold = %d)... \n', SSD_threshold_factor);
    for ch = 1:nCh
        if isempty(sp_slope{ch}), continue; end
        wfs = sp_slope{ch}(:,2:end);
        if size(wfs,1) < 5, sp_SSD{ch} = sp_slope{ch}; continue; end
        
        mean_wave = mean(wfs,1);
        SSD       = sum((wfs - mean_wave).^2, 2);
        valid_idx = SSD <= (SSD_threshold_factor * mean(SSD));
        sp_SSD{ch} = sp_slope{ch}(valid_idx,:);
    end
else
    fprintf('Skipping SSD filter.\n');
end

%% ================= 3) ZERO-CROSSING CHECK =================
sp_zcross = sp_SSD;
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

%% ================= 4) TEMPLATE CORRELATION (Dual-Window) =================
sp_corr = sp_zcross; 
if do_corr_filter
    fprintf('\nCorrelation Filtering (Dual Window Template: %d-%dms & %d-%dms)...\n', ...
        baseline_window_ms, late_template_window);
    
    for ch = 1:nCh
        if isempty(sp_corr{ch}), continue; end
        spt = sp_corr{ch}(:,1);     
        wfs = sp_corr{ch}(:,2:end); 
        nSpikes = size(wfs,1);
        
        % --- A. BUILD DUAL TEMPLATE ---
        mask_safe = false(nSpikes,1);
        for tr = 1:nTrials
            t0 = trig(tr)/FS*1000;
            
            % Check Baseline Window
            in_base = (spt >= t0 + baseline_window_ms(1) & spt <= t0 + baseline_window_ms(2));
            
            % Check Late Evoked Window
            in_late = (spt >= t0 + late_template_window(1) & spt <= t0 + late_template_window(2));
            
            mask_safe = mask_safe | in_base | in_late;
        end
        wfs_safe = wfs(mask_safe,:);
        
        if size(wfs_safe, 1) >= 5
            template = mean(wfs_safe, 1);
            source_str = 'Dual_Safe';
        else
            % Fallback: Just Baseline if Late is empty or vice versa
            if sum(mask_safe) > 0
                 template = mean(wfs_safe, 1);
                 source_str = 'Scant_Safe';
            else
                 fprintf('Ch %2d: Not enough safe spikes for template. Skipping.\n', ch);
                 continue; 
            end
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
QC_params.SSD_threshold      = SSD_threshold_factor;
QC_params.baseline_window    = baseline_window_ms;
QC_params.late_template_win  = late_template_window;
QC_params.cleanup_window     = cleanup_window_ms;
QC_params.corr_thresh        = corr_thresh;
out_name = [base_name '.sp_xia_SSD.mat'];
if isfile(out_name), delete(out_name); end 
save(out_name, 'sp_in','sp_slope','sp_SSD','sp_zcross','sp_corr','QC_params','-v7.3');
fprintf('\nSaved cleaned spikes to %s\n', out_name);