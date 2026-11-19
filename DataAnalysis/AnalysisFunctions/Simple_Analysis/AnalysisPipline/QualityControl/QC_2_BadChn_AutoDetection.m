%% ============================================================
%        MODE 1 — Automatic Bad Channel Detection
%        Zero-cross based QC + firing rate + amplitude tests
% =============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================= USER SETTINGS =================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Single1_251125_110714';

Electrode_Type = 1;      % 0 rigid, 1 flex 1shank, 2 flex 4shank

% ---- QC thresholds ----
amp_limit             = 300;     % μV (amplitude beyond this = bad spike)
max_pct_amp_fail      = 20;      % (% spikes beyond limit)
max_baseline_FR       = 50;      % sp/s
max_zero_cross_fail   = 40;      % (% waveforms failing zero-cross)

% ---- analysis windows ----
baseline_win_ms = [-60 -5];     % pre-stim baseline window
analysis_win_ms = [0 50];       % QC window for post-stim spikes

FS = 30000;   % sampling rate

% ================= CHECK FOLDER =================
if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end
cd(data_folder);
fprintf('QC Mode 1 Running in folder:\n%s\n\n', data_folder);

% ================= EXTRACT BASE NAME =================
parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');
if numel(u) >= 4
    base_name = last_folder(1 : u(end-1)-1);
else
    base_name = last_folder;
end

% ===== Load spike data =====

fprintf('Searching for spike files...\n');

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
    else
        error('SSD file found but contains no usable spike variable.');
    end

elseif isfile(base_file)
    fprintf('Loading base spike file: %s\n', base_file);
    S = load(base_file);

    if isfield(S,'sp_clipped')
        sp_use = S.sp_clipped;
    elseif isfield(S,'sp')
        sp_use = S.sp;
    else
        error('Base spike file contains no usable spike variable.');
    end
else
    error('No spike files found for this dataset.');
end

% Validate spike data
if ~iscell(sp_use)
    error('Loaded spike data is not a cell array. Got: %s', class(sp_use));
end
nCh = numel(sp_use);
fprintf('Successfully loaded spike data: %d channels.\n\n', nCh);

% ================= LOAD TRIGGERS =================
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);

% ================= ELECTRODE MAP =================
d = Depth_s(Electrode_Type);

%% ======================================================
%              QC COMPUTATION PER CHANNEL
%% ======================================================

QC_metrics = struct();
SuggestedBadChannels = [];
fprintf('Running QC...\n\n');

for ich = 1:length(d)
    ch = d(ich);

    if ch > nCh || isempty(sp_use{ch})
        fprintf('Ch %2d → BAD (empty channel)\n', ich);
        SuggestedBadChannels(end+1) = ich;
        QC_metrics(ich).reason = {'Empty channel'};
        QC_metrics(ich).is_bad = true;
        continue;
    end

    data = sp_use{ch};
    sp_times = data(:,1);          
    sp_wave  = data(:,2:end);       

    % 1) Amplitude QC
    max_vals = max(sp_wave, [], 2);
    min_vals = min(sp_wave, [], 2);
    amp_bad = (max_vals > amp_limit) | (min_vals < -amp_limit);
    pct_over = mean(amp_bad) * 100;

    % 2) Baseline firing QC
    FR_baseline = 0;
    if ~isempty(trig)
        bl_sp = 0;
        for t = 1:length(trig)
            t0 = trig(t) / FS * 1000;
            bl_sp = bl_sp + sum(sp_times > (t0 + baseline_win_ms(1)) & ...
                                sp_times <= (t0 + baseline_win_ms(2)));
        end
        FR_baseline = bl_sp / length(trig) / ((baseline_win_ms(2)-baseline_win_ms(1))/1000);
    end

    % 3) Zero-cross consistency QC
    has_pos = any(sp_wave > 0, 2);
    has_neg = any(sp_wave < 0, 2);
    zero_cross = has_pos & has_neg;
    pct_no_zero_cross = mean(~zero_cross) * 100;

    % 4) Evaluate BAD channel
    is_bad = false;
    reason = {};

    if pct_over > max_pct_amp_fail
        is_bad = true;
        reason{end+1} = sprintf('%.1f%% spikes exceed amplitude limit (%d μV)', pct_over, amp_limit);
    end

    if FR_baseline > max_baseline_FR
        is_bad = true;
        reason{end+1} = sprintf('Baseline firing too high: %.1f sp/s', FR_baseline);
    end

    if pct_no_zero_cross > max_zero_cross_fail
        is_bad = true;
        reason{end+1} = sprintf('%.1f%% waveforms fail zero-crossing test', pct_no_zero_cross);
    end

    %% Store metrics
    QC_metrics(ich).channel           = ich;
    QC_metrics(ich).pct_amp_exceed    = pct_over;
    QC_metrics(ich).baseline_FR       = FR_baseline;
    QC_metrics(ich).pct_no_zero_cross = pct_no_zero_cross;
    QC_metrics(ich).is_bad            = is_bad;
    QC_metrics(ich).reason            = reason;

    %% PRINT BAD CHANNELS
    if is_bad
        SuggestedBadChannels(end+1) = ich;
        fprintf('Ch %2d → BAD\n', ich);

        for r = 1:length(reason)
            fprintf('   - %s\n', reason{r});
        end

        fprintf('\n');
    end
end

%% ======================================================
%                SAVE QC RESULT
%% ======================================================

save_file = sprintf('%s.ChannelQC_Auto.mat', base_name);
save(save_file, 'QC_metrics', 'SuggestedBadChannels');

fprintf('QC complete.\n');
fprintf('Suggested BAD channels:\n');
disp(SuggestedBadChannels);
fprintf('Saved QC results to:\n  %s\n', save_file);
