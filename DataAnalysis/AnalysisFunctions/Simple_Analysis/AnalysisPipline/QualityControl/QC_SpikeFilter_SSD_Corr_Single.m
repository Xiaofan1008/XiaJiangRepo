%% ============================================================
% SpikeFiltering_Evoked_Cleanup.m
% Refined: Uses Baseline Template + Zero-Crossing Check
% ============================================================
clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Single7';
FS = 30000;                 

% 1. FILTERING PARAMETERS
SSD_threshold_factor = 16;  % First pass: Remove massive outliers
corr_thresh          = 0.6; % Correlation strictness 

% 2. WINDOWS
baseline_window_ms   = [-100 -5]; 
late_evoked_window   = [10 50];   
cleanup_window_ms    = [0 100];   

do_SSD_filter        = 1;
do_corr_filter       = 1;

%% ================= FOLDER & BASE NAME =================
if ~isfolder(data_folder), error('Folder does not exist.'); end
cd(data_folder);
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4, base_name = last_folder(1 : underscores(end-1)-1);
else, base_name = last_folder; end

%% ================= LOAD SPIKES =================
fname_sp = [base_name '.sp_xia.mat'];
assert(isfile(fname_sp), 'Cannot find %s', fname_sp);
S_in = load(fname_sp);
if isfield(S_in,'sp_clipped'), sp_in = S_in.sp_clipped;
elseif isfield(S_in,'sp'), sp_in = S_in.sp;
else, error('No sp or sp_clipped found'); end
nCh = numel(sp_in);

if isempty(dir('*.trig.dat')), cur=pwd; cleanTrig_sabquick; cd(cur); end
trig = loadTrig(0);
nTrials = numel(trig);

%% ================= 1) SSD-BASED FILTERING (Pre-Clean) =================
sp_SSD = sp_in;
if do_SSD_filter
    fprintf('\nSSD Filtering... \n');
    for ch = 1:nCh
        if isempty(sp_in{ch}), continue; end
        wfs = sp_in{ch}(:,2:end);
        if size(wfs,1) < 5, sp_SSD{ch} = sp_in{ch}; continue; end
        
        mean_wave = mean(wfs,1);
        SSD       = sum((wfs - mean_wave).^2, 2);
        valid_idx = SSD <= (SSD_threshold_factor * mean(SSD));
        sp_SSD{ch} = sp_in{ch}(valid_idx,:);
    end
else
    fprintf('Skipping SSD filter.\n');
end

%% ================= 2) ZERO-CROSSING CHECK (New) =================
% A real spike is biphasic (goes down then up). 
% If Max(waveform) <= 0, it is purely negative (likely artifact/noise).
sp_zcross = sp_SSD;
fprintf('\nZero-Crossing Check (Removing purely negative waveforms)...\n');

for ch = 1:nCh
    if isempty(sp_zcross{ch}), continue; end
    wfs = sp_zcross{ch}(:,2:end);
    
    % Check: Does the waveform ever go above 0?
    % logical vector: true if max > 0
    has_positive_phase = max(wfs, [], 2) > 0;
    
    removed = sum(~has_positive_phase);
    if removed > 0
        sp_zcross{ch} = sp_zcross{ch}(has_positive_phase, :);
        fprintf('Ch %2d: Removed %d spikes (No positive phase/Zero Crossing).\n', ch, removed);
    end
end

%% ================= 3) TEMPLATE CORRELATION (Refined) =================
sp_corr = sp_zcross; % Use the Zero-Crossed data as input
if do_corr_filter
    fprintf('\nCorrelation Filtering...\n');
    
    for ch = 1:nCh
        if isempty(sp_corr{ch}), continue; end
        spt = sp_corr{ch}(:,1);     
        wfs = sp_corr{ch}(:,2:end); 
        nSpikes = size(wfs,1);
        
        % --- A. BUILD TEMPLATE ---
        mask_base = false(nSpikes,1);
        for tr = 1:nTrials
            t0 = trig(tr)/FS*1000;
            mask_base = mask_base | (spt >= t0 + baseline_window_ms(1) & spt <= t0 + baseline_window_ms(2));
        end
        wfs_base = wfs(mask_base,:);
        
        if size(wfs_base, 1) >= 5
            template = mean(wfs_base, 1);
            source_str = 'Baseline';
        else
            mask_late = false(nSpikes,1);
            for tr = 1:nTrials
                t0 = trig(tr)/FS*1000;
                mask_late = mask_late | (spt >= t0 + late_evoked_window(1) & spt <= t0 + late_evoked_window(2));
            end
            wfs_late = wfs(mask_late,:);
            
            if size(wfs_late, 1) >= 5
                template = mean(wfs_late, 1);
                source_str = 'Late_Evoked';
            else
                fprintf('Ch %2d: Not enough spikes for template. Skipping.\n', ch);
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
            
            if any(in_window) % Optimization
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
QC_params.SSD_threshold      = SSD_threshold_factor;
QC_params.baseline_window    = baseline_window_ms;
QC_params.cleanup_window     = cleanup_window_ms;
QC_params.corr_thresh        = corr_thresh;
QC_params.zero_crossing_check= 1;

out_name = [base_name '.sp_xia_SSD.mat'];
if isfile(out_name), delete(out_name); end 
save(out_name, 'sp_in','sp_SSD','sp_zcross','sp_corr','QC_params','-v7.3');
fprintf('\nSaved cleaned spikes to %s\n', out_name);