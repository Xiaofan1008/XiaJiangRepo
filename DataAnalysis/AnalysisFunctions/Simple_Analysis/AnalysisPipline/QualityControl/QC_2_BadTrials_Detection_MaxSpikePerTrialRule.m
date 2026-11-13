%% ============================================
%   BAD TRIAL DETECTION (Channel-Specific)
%   - Criteria 1: Sliding Window Density (Artifacts)
%   - Criteria 2: Burst Duration (Sustained Noise)
% ============================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ---------------- USER SETTINGS ----------------
data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq3';
Electrode_Type = 1;       
FS = 30000;               

% 1. BURST DETECTION SETTINGS 
detect_win_ms = [0 50];   
bin_ms        = 1;        
smooth_win_ms = 3;        
SD_factor             = 3;
min_burst_duration_ms = 15;

% 2. SLIDING WINDOW ARTIFACT CHECK 
% Instead of a fixed window, we slide a small window across detect_win_ms
slide_win_size_ms    = 10;     % Look at 10ms chunks
% If > 4 spikes in ANY 10ms chunk, it is an artifact (approx >300Hz)
max_spikes_in_slide  = 4;      

% ---------------- CHECK FOLDER ----------------
if ~isfolder(data_folder), error('Folder not found'); end
cd(data_folder);
fprintf('\nRunning Bad Trial Detection in:\n%s\n\n', data_folder);

% ---------------- EXTRACT BASE NAME ----------------
parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');
if numel(u) >= 4, base_name = last_folder(1 : u(end-1)-1);
else, base_name = last_folder; end

% ---------------- LOAD SPIKE DATA ----------------
ssd_file  = [base_name '.sp_xia_SSD.mat'];
base_file = [base_name '.sp_xia.mat'];
sp_use = [];

if isfile(ssd_file)
    fprintf('Loading SSD file: %s\n', ssd_file);
    S = load(ssd_file);
    if isfield(S,'sp_corr'), sp_use = S.sp_corr;
    elseif isfield(S,'sp_SSD'), sp_use = S.sp_SSD;
    elseif isfield(S,'sp_in'), sp_use = S.sp_in; end
elseif isfile(base_file)
    fprintf('Loading base spike file: %s\n', base_file);
    S = load(base_file);
    if isfield(S,'sp_clipped'), sp_use = S.sp_clipped;
    elseif isfield(S,'sp'), sp_use = S.sp; end
else, error('No spike file found'); end

if ~iscell(sp_use), error('Spike data must be a cell array.'); end
nCh = numel(sp_use);

% ---------------- LOAD TRIGGERS ----------------
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);
nTrials = numel(trig);

% ---------------- ELECTRODE MAP ----------------
d = Depth_s(Electrode_Type); 

% ---------------- INIT OUTPUT ----------------
BadTrials  = cell(nCh,1);
GoodTrials = cell(nCh,1);

edges   = detect_win_ms(1) : bin_ms : detect_win_ms(2);
bin_s   = bin_ms / 1000;
min_bins = ceil(min_burst_duration_ms / bin_ms);
smooth_N = max(1, round(smooth_win_ms / bin_ms));

fprintf('Checking Range: [%d %d] ms\nMethod: Sliding Window (%d ms)\nMax Spikes Allowed: %d (per window)\n', ...
    detect_win_ms(1), detect_win_ms(2), slide_win_size_ms, max_spikes_in_slide);

%% ---------------- MAIN LOOP ----------------
fprintf('\nDetecting bad trials per channel...\n');
total_bad = 0;   

% Define sliding steps (1ms steps)
slide_starts = detect_win_ms(1) : 1 : (detect_win_ms(2) - slide_win_size_ms);

for ich = 1:nCh
    ch = d(ich);
    if isempty(sp_use{ch})
        BadTrials{ich}  = [];
        GoodTrials{ich} = 1:nTrials;
        fprintf('Ch %2d: EMPTY\n', ich);
        continue;
    end
    
    sp_times = sp_use{ch}(:,1);
    bad_list = [];
    
    for tr = 1:nTrials
        t0 = trig(tr)/FS*1000;
        
        % Get all spikes in the broad detection range (e.g. 0-50ms)
        mask_broad = sp_times >= (t0 + detect_win_ms(1)) & ...
                     sp_times <  (t0 + detect_win_ms(2));
        tt_broad = sp_times(mask_broad) - t0;
        
        if isempty(tt_broad), continue; end
        
        % ---- CHECK 1: SLIDING WINDOW DENSITY (New) ----
        is_artifact = false;
        
        % Optimization: Only slide if total spikes exceed threshold
        if numel(tt_broad) > max_spikes_in_slide
            for t_start = slide_starts
                % Count spikes in this specific 10ms window
                % (tt_broad is already relative to t0)
                n_in_win = sum(tt_broad >= t_start & tt_broad <= (t_start + slide_win_size_ms));
                
                if n_in_win > max_spikes_in_slide
                    is_artifact = true;
                    break; % Found artifact, stop checking this trial
                end
            end
        end
        
        if is_artifact
            bad_list(end+1) = tr;
            continue; % Skip to next trial
        end
        
        % ---- CHECK 2: BURST DURATION (Original) ----
        counts = histcounts(tt_broad, edges);
        FR = counts / bin_s;
        
        if smooth_N > 1, FR_s = movmean(FR, smooth_N);
        else, FR_s = FR; end
        
        mu = mean(FR_s);
        sd = std(FR_s);
        
        if sd == 0, continue; end
        
        thr = mu + SD_factor * sd;
        above = FR_s > thr;
        
        if any(above)
            diffA = diff([0 above 0]);
            srt = find(diffA == 1);
            stp = find(diffA == -1) - 1;
            dur = stp - srt + 1;
            if any(dur >= min_bins)
                bad_list(end+1) = tr; 
            end
        end
    end
    
    BadTrials{ich}  = unique(bad_list);
    GoodTrials{ich} = setdiff(1:nTrials, BadTrials{ich});
    
    if ~isempty(BadTrials{ich})
        fprintf('Ch %2d: %d bad trials (Artifact/Burst)\n', ich, numel(BadTrials{ich}));
        total_bad = total_bad + numel(BadTrials{ich});
    end
end

%% ---------------- SUMMARY ----------------
fprintf('\n===============================\n');
fprintf('BAD TRIAL DETECTION SUMMARY\n');
fprintf('Sliding Window: %d ms | Threshold: >%d spikes\n', ...
    slide_win_size_ms, max_spikes_in_slide);
fprintf('TOTAL BAD TRIALS FLAGGED: %d\n', total_bad);

%% ---------------- SAVE RESULTS ----------------
save_file = sprintf('%s.BadTrials.mat', base_name);
save(save_file, 'BadTrials', 'GoodTrials', ...
    'detect_win_ms', 'slide_win_size_ms', ... 
    'bin_ms', 'smooth_win_ms', ...
    'SD_factor', 'min_burst_duration_ms', 'max_spikes_in_slide');
fprintf('Saved bad trial info to: %s\n\n', save_file);