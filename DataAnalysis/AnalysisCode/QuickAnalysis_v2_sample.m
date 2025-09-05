%% ============================================================
%% Quick Analysis (Clean Refactor) — In-Experiment Check
%% ============================================================

clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Sabrina/Experimental_Design'));

%% -------------------- CONFIG -------------------- %%
% Window around each trigger (ms)
WINDOW_MS        = [-100 100];

% Artifact blank (ms)
PRE_BLANK_MS     = 1.0;      % before stim
POST_BLANK_MS    = 2.0;      % after  stim

% Cross-fade guards (ms) + detection guard around the blank
TAPER_PRE_MS     = 0.8;
TAPER_POST_MS    = 0.8;
GUARD_PRE_MS     = max(TAPER_PRE_MS, 0.5);
GUARD_POST_MS    = max(TAPER_POST_MS, 0.5);

% Spike waveform extraction (ms)  ~1 ms pre, 2 ms post
WF_PRE_MS        = 1.0;
WF_POST_MS       = 2.0;

% Spike shape criteria
PP_MIN_uV        = 30;
PP_MAX_uV        = 500;
TP_MIN_MS        = 0.30;
TP_MAX_MS        = 1.20;

% Filtering (band-pass)
F_HP             = 300;      % Hz
F_LP             = 3000;     % Hz
FIR_ORDER        = 512;      % linear-phase, causal

% Thresholding
THR_DUR_SEC      = 10;       % seconds used to estimate sigma
THR_MULT         = -4.5;     % x MAD/0.6745

% Trials
TRIAL_START      = 1;
NTRIALS          = 1000;

% PSTH
BIN_MS           = 1;
SMOOTH_MS        = 3;        % causal Gaussian FIR
PSTH_RED_THR     = 200;      % sp/s

% Tuning / response window (ms)
RESP_MS          = [0 20];

% Figure layout
GRID_NCOLS       = 4;        % 8x4 grid for 32 ch (change if needed)

%% -------------------- LOAD META -------------------- %%
filepath = pwd;
[amplifier_channels, freq_params] = read_Intan_RHS2000_file;
nChn = numel(amplifier_channels);
FS   = freq_params.amplifier_sample_rate;

dName = 'amplifier';
vFID  = fopen(fullfile(filepath, [dName '.dat']), 'r');

% triggers
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_all = loadTrig(0);
d        = Depth_s(0);       % 0: single-shank rigid (your setting)

% crop to requested trials
trig = trig_all;
if NTRIALS <= numel(trig)
    trig = trig(TRIAL_START : TRIAL_START+NTRIALS-1);
end
NTRIG = numel(trig);

% trial params
TrialParams = loadTrialParams;
trialIDs    = cell2mat(TrialParams(:,2));
trialIDs    = trialIDs(1:NTRIG);

% stim params
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
S = load(fileDIR(1).name, 'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');
StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
E_MAP             = S.E_MAP;
n_Trials_total    = S.n_Trials; %#ok<NASGU>

% amplitudes per trial (μA)
trialAmps_all = cell2mat(StimParams(2:end, 16));
trialAmps_all = trialAmps_all(1:simultaneous_stim:end);
trialAmps     = trialAmps_all(TRIAL_START : TRIAL_START+NTRIG-1);

[uniqAmps, ~, ampIdx] = unique(trialAmps(:));
Amps = uniqAmps;
Amps(Amps == -1) = 0;                  % treat -1 as 0 μA
n_AMP = numel(Amps);

% stim sets (which electrodes per trial)
E_NAME = E_MAP(2:end);
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_NAME);

stimChPerTrial = cell(NTRIG,1);
for t = 1:NTRIG
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v  = unique(idx_all(rr)); v = v(v>0).';
    stimChPerTrial{t} = v;
end
comb = zeros(NTRIG, simultaneous_stim);
for t = 1:NTRIG
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end
[uniqueComb, ~, combClass] = unique(comb, 'rows');
combClass_win = combClass;           % matched NTRIG
uniqueStimCh  = unique([stimChPerTrial{:}]);   % channels that were ever stimulated

%% -------------------- DERIVED CONSTS -------------------- %%
win_samp   = round(WINDOW_MS/1000 * FS);
nSampWin   = diff(win_samp) + 1;
stim_idx   = -win_samp(1) + 1;

preBlank   = round(PRE_BLANK_MS/1000 * FS);
postBlank  = round(POST_BLANK_MS/1000 * FS);
guard_pre  = max(1, round(GUARD_PRE_MS/1000  * FS));
guard_post = max(1, round(GUARD_POST_MS/1000 * FS));

preS       = round(WF_PRE_MS /1000 * FS);
postS      = round(WF_POST_MS/1000 * FS);
refrac_smp = max(1, round(1e-3 * FS));       % 1 ms refractory

tp_min_smp = round(TP_MIN_MS/1000 * FS);
tp_max_smp = round(TP_MAX_MS/1000 * FS);

% causal FIR band-pass
bpFIR  = fir1(FIR_ORDER, [F_HP F_LP]/(FS/2), 'bandpass', hamming(FIR_ORDER+1), 'scale');
GD     = FIR_ORDER/2;   % constant group delay

% PSTH smoothing kernel (causal Gaussian-like)
[gPSTH, bin_s, edges, ctrs] = make_causal_kernel(BIN_MS, SMOOTH_MS, WINDOW_MS);

%% -------------------- THRESHOLDS -------------------- %%
thr = zeros(nChn,1);
num_samp_thr = min( floor(THR_DUR_SEC*FS), floor(dir([dName '.dat']).bytes/(nChn*2)) );
assert(num_samp_thr > 5*FS, 'Recording <5 s; not enough for thresholds.');

fseek(vFID, 0, 'bof');
th_blk_uV = fread(vFID, [nChn, num_samp_thr], 'int16') * 0.195;

for ch = 1:nChn
    x = double(th_blk_uV(d(ch),:));
    xf = filter(bpFIR,1,x);
    if GD>0, xf = [xf(GD+1:end), zeros(1,GD)]; end
    sigma = median(abs(xf))/0.6745;
    thr(ch) = THR_MULT * sigma;
end

%% -------------------- MAIN LOOP -------------------- %%
Spikes  = cell(nChn, NTRIG);                % per channel/trial
MUA_all = zeros(nChn, nSampWin, NTRIG);     % filtered traces (for quick looks)

for tr = 1:NTRIG
    fprintf('Reading trial %d/%d\n', tr, NTRIG);
    start_index = trig(tr) + win_samp(1);
    fseek(vFID, start_index * nChn * 2, 'bof');
    data = fread(vFID, [nChn, nSampWin], 'int16') * 0.195; % μV
    
    % artifact blank
    fill_idx = max(1, stim_idx - preBlank) : min(nSampWin, stim_idx + postBlank);
    keep_idx = setdiff(1:nSampWin, fill_idx);

    for ch = 1:nChn
        raw = data(d(ch),:);
        raw_i = raw;
        raw_i(fill_idx) = interp1(keep_idx, raw(keep_idx), fill_idx, 'linear', 'extrap');

        % store filtered trace (for display if needed)
        xf = filter(bpFIR, 1, raw_i);
        if GD>0, xf = [xf(GD+1:end), zeros(1,GD)]; end
        MUA_all(ch,:,tr) = xf;

        % ---- spike detect ----
        cross = find(xf(2:end) < thr(ch) & xf(1:end-1) >= thr(ch)) + 1;

        % mask away small guard around interpolated region
        gA = max(1, (stim_idx - preBlank) - guard_pre);
        gB = min(nSampWin, (stim_idx + postBlank) + guard_post);
        keep = ~(cross >= gA & cross <= gB);
        cross = cross(keep);

        % refractory de-dupe
        if numel(cross)>1
            cross = cross([true, diff(cross) > refrac_smp]);
        end

        % trough refine in full [preS, postS] window
        trough = zeros(0,1);
        for k = 1:numel(cross)
            a = max(1, cross(k)-preS);
            b = min(nSampWin, cross(k)+postS);
            if a+preS > b-postS, continue; end
            [~, rel] = min(xf(a:b));
            trough(end+1,1) = a + rel - 1;
        end
        if isempty(trough)
            Spikes{d(ch),tr} = struct('idx',[],'t_ms',[],'wf',zeros(0,preS+postS+1));
            continue;
        end

        % enforce refractory via keep-stronger within window
        [trough, ord] = sort(trough(:));
        vals = xf(trough);
        keep = true(size(trough));
        last = 1;
        for m = 2:numel(trough)
            if trough(m) - trough(last) <= refrac_smp
                if vals(m) < vals(last)
                    keep(last) = false; last = m;
                else
                    keep(m) = false;
                end
            else
                last = m;
            end
        end
        trough = trough(keep);

        % extract waveforms + shape filters
        wfs = zeros(0, preS+postS+1);
        idx = zeros(0,1);
        for k = 1:numel(trough)
            ip = trough(k);
            if ip-preS < 1 || ip+postS > nSampWin, continue; end
            wf = xf(ip-preS : ip+postS);
            if max(abs(wf)) > 500, continue; end

            [~, i_tr] = min(wf);
            [~, i_pk_rel] = max(wf(i_tr:end));
            i_pk = i_tr + i_pk_rel - 1;

            is_biphasic = (min(wf) < 0) && (max(wf(i_tr:end)) > 0);
            if ~is_biphasic, continue; end

            pp = max(wf) - min(wf);
            if ~(pp >= PP_MIN_uV && pp <= PP_MAX_uV), continue; end

            tpw = i_pk - i_tr;
            if ~(tpw >= tp_min_smp && tpw <= tp_max_smp), continue; end

            wfs(end+1,:) = wf;
            idx(end+1,1) = ip; 
        end

        t_ms = (idx - stim_idx)/FS*1000;
        Spikes{d(ch),tr} = struct('idx',idx,'t_ms',t_ms,'wf',wfs);
    end
end
fclose(vFID);

%% -------------------- AVERAGE WF PER AMP -------------------- %%
wf_len   = preS + postS + 1;
AvgWFamp = nan(nChn, wf_len, n_AMP);
NspkAmp  = zeros(nChn, n_AMP);

for ch = 1:nChn
    for k = 1:n_AMP
        rows = [];
        for tr = 1:NTRIG
            if ampIdx(tr) ~= k, continue; end
            S = Spikes{ch,tr};
            if isempty(S) || isempty(S.wf), continue; end
            rows = [rows; S.wf];
        end
        if ~isempty(rows)
            AvgWFamp(ch,:,k) = mean(rows,1,'omitnan');
            NspkAmp(ch,k)    = size(rows,1);
        end
    end
end

%% -------------------- QUICK PLOTS -------------------- %%
