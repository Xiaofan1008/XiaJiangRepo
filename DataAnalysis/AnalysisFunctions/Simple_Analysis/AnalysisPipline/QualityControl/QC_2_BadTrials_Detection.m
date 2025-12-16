%% ============================================
%   BAD TRIAL DETECTION 
% ============================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ---------------- USER SETTINGS ----------------
data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim1';

Electrode_Type = 1;       % 0 rigid, 1 flex single-shank, 2 flex 4-shank
FS = 30000;               % sampling rate (Hz)

% Time window to look for bursts (relative to trigger, ms)
detect_win_ms = [0 50];   % e.g. [0 50] ms

% Firing-rate computation
bin_ms        = 1;        % FR bin size (ms)
smooth_win_ms = 3;        % moving-average window (ms)

% Burst definition: FR > mean + SD_factor*SD for >= min_burst_duration_ms
SD_factor             = 3;
min_burst_duration_ms = 15;

% ---------------- CHECK FOLDER ----------------
if ~isfolder(data_folder)
    error('Folder not found');
end
cd(data_folder);
fprintf('\nRunning Bad Trial Detection in:\n%s\n\n', data_folder);

% ---------------- EXTRACT BASE NAME ----------------
parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');
if numel(u) >= 4
    base_name = last_folder(1 : u(end-1)-1);
else
    base_name = last_folder;
end

% ---------------- LOAD SPIKE DATA ----------------
ssd_file  = [base_name '.sp_xia_SSD.mat'];
base_file = [base_name '.sp_xia.mat'];

sp_use = [];

if isfile(ssd_file)
    fprintf('Loading SSD file: %s\n', ssd_file);
    S = load(ssd_file);

    if isfield(S,'sp_corr')
        sp_use = S.sp_corr;
    elseif isfield(S,'sp_SSD')
        sp_use = S.sp_SSD;
    elseif isfield(S,'sp_in')
        sp_use = S.sp_in;
    end

elseif isfile(base_file)
    fprintf('Loading base spike file: %s\n', base_file);
    S = load(base_file);

    if isfield(S,'sp_clipped')
        sp_use = S.sp_clipped;
    elseif isfield(S,'sp')
        sp_use = S.sp;
    end
else
    error('No spike file found');
end

if ~iscell(sp_use)
    error('Spike data must be a cell array.');
end

nCh = numel(sp_use);

% ---------------- LOAD TRIGGERS ----------------
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
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

fprintf('Detect window: [%d %d] ms, bin = %d ms, smooth = %d bins, burst ≥ %d ms\n', ...
    detect_win_ms(1), detect_win_ms(2), bin_ms, smooth_N, min_burst_duration_ms);

%% ---------------- MAIN LOOP ----------------
fprintf('\nDetecting burst trials per channel...\n\n');

total_bad = 0;   % count all bad trials across channels

for ich = 1:nCh
    ch = d(ich);
    if isempty(sp_use{ch})
        BadTrials{ch}  = [];
        GoodTrials{ch} = 1:nTrials;
        fprintf('Ch %2d: EMPTY → 0 bad / %d total trials\n', ch, nTrials);
        continue;
    end

    sp_times = sp_use{ch}(:,1);
    bad_list = [];

    for tr = 1:nTrials
        t0 = trig(tr)/FS*1000;

        mask = sp_times >= (t0 + detect_win_ms(1)) & ...
               sp_times <  (t0 + detect_win_ms(2));
        tt = sp_times(mask) - t0;

        % ---- empty trial = GOOD ----
        if isempty(tt)
            continue;
        end
        % ---- 1ms-binned firing rate ----
        counts = histcounts(tt, edges);
        FR = counts / bin_s;
        % ---- smoothing ----
        if smooth_N > 1
            FR_s = movmean(FR, smooth_N);
        else
            FR_s = FR;
        end
        % ---- mean + SD*3 ----
        mu = mean(FR_s);
        sd = std(FR_s);

        if sd == 0
            continue;
        end
        thr = mu + SD_factor * sd;
        % ---- find bursts ----
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

    BadTrials{ch}  = unique(bad_list);
    GoodTrials{ch} = setdiff(1:nTrials, BadTrials{ch});

    fprintf('Ch %2d: %d bad / %d total trials\n', ch, numel(BadTrials{ch}), nTrials);
    total_bad = total_bad + numel(BadTrials{ch});
end

%% ---------------- SUMMARY ----------------
fprintf('\n===============================\n');
fprintf('BAD TRIAL DETECTION SUMMARY\n');
fprintf('Total channels: %d\n', nCh);
fprintf('Total trials per channel: %d\n', nTrials);
fprintf('TOTAL BAD TRIALS DETECTED ACROSS DATASET: %d\n', total_bad);

%% ---------------- SAVE RESULTS ----------------
save_file = sprintf('%s.BadTrials.mat', base_name);
save(save_file, 'BadTrials', 'GoodTrials', ...
    'detect_win_ms', 'bin_ms', 'smooth_win_ms', ...
    'SD_factor', 'min_burst_duration_ms');

fprintf('Saved bad trial info to: %s\n\n', save_file);