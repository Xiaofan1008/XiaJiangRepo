%% MISSING_TRIG_1B_CHECK
%
% Extra check between step 1 (Diagnose) and step 2 (Recover).
% Run AFTER missing_trig_1_diagnose.m on the same folder.
%
% Tests three ideas that came out of reading RunMultiElect_exp.m:
%
%   A) ANCHORS - the run script pauses right after trial 499 (999, ...),
%      and the very last pulse is cut short (~0.5 ms) because line 171
%      sets the line low straight after the final trial. Pauses are found
%      by their signature (very wide pulse followed by a long gap), not by
%      trial number. Between anchors we then know exactly how many trials
%      are missing.
%
%   B) DROPPED SAMPLES - the script waits >= 5 s across each pause, but the
%      recordings show only ~2.7 s. time.dat stores a timestamp for every
%      sample, so any jump there means samples were not recorded.
%
%   C) PARAMETER-UPDATE ARTIFACT - each trial sends its parameters ('SET'),
%      waits r = PARAM..PARAM+125 ms, then fires the TTL. If applying new
%      settings makes a small artifact, every trial should show one
%      ~PARAM..PARAM+125 ms before its trigger. This scans every interval
%      between triggers (and before the first one), finds those events,
%      and counts them inside gaps: a gap with m missing trials should
%      contain m + 1 of them.
%
%   D) WHEN DOES IT APPEAR - for trials whose number is certain (from the
%      anchors), checks whether the event appears only when the electrodes
%      (or any parameter) change from the previous trial.
%
% Output: printed summary, figures in missing_trig_figs/ (1b_check_*.png),
%         results saved as report.check in missing_trig_report.mat.
%
% Read-only -- does NOT modify .trig.dat or any recording file.

clear; clc;

%% ======================= SETTINGS =======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX039/Xia_Linearity_S3_100_400um_SimSeq1_261009_190544';   % <-- same folder as step 1

noise_k            = 8;    % event threshold = median + noise_k x robust SD of pre-trigger noise
noise_ms           = 50;   % pre-trigger noise window (the TTL comes >= PARAM ms after the SET)
guard_fall_ms      = 5;    % start scanning this long after the previous pulse falls...
guard_rise_ms      = 100;  % ...and at least this long after the previous pulse rises (stim artifact)
pre_gap_ms         = 2;    % stop scanning this long before the next trigger
refr_ms            = 150;  % threshold crossings closer than this are one event
img_ms             = 600;  % how far before each trigger the image/profile shows
chunk_s            = 5;    % read long windows in chunks of this many seconds

pause_every        = 500;  % rem(count,500) in RunMultiElect_exp.m
pause_width_factor = 2.5;  % pause pulse: width > this x median width ...
pause_iti_factor   = 4;    % ... AND the next interval > this x median interval
short_last_ms      = 5;    % last pulse shorter than this = final trial (line 171)

%% ======================= LOAD =======================
assert(isfolder(data_folder), 'data_folder not found:\n%s', data_folder);
cd(data_folder);
report_file = fullfile(data_folder, 'missing_trig_report.mat');
assert(isfile(report_file), 'No missing_trig_report.mat here -- run missing_trig_1_diagnose.m first.');
load(report_file, 'report');
assert(isfield(report, 'diagnose'), 'Report has no step 1 results -- run missing_trig_1_diagnose.m first.');
info = report.info;  dg = report.diagnose;

FS    = info.FS;
nSam  = info.nSam;
n_exp = info.n_expected;
rise  = double(dg.rise(:));
fall  = double(dg.fall(:));
n     = numel(rise);
width_ms = (fall - rise + 1) / FS * 1000;
iti_ms   = diff(rise) / FS * 1000;
med_w    = median(width_ms);
med_iti  = median(iti_ms);

exp_path = fullfile(data_folder, info.exp_file);
vlist = whos('-file', exp_path);
want  = intersect({'PARAM', 'PULSE', 'StimParams', 'simultaneous_stim'}, {vlist.name});
S = load(exp_path, want{:});
PARAM = NaN;  PULSE = NaN;
if isfield(S, 'PARAM'), PARAM = double(S.PARAM); end
if isfield(S, 'PULSE'), PULSE = double(S.PULSE); end

fprintf('Folder: %s\n', data_folder);
fprintf('Expected %d trials, detected %d triggers (missing %d)\n', n_exp, n, n_exp - n);
fprintf('From experiment file: PARAM = %g ms (SET -> TTL wait is PARAM..PARAM+125 ms), PULSE = %g ms\n', PARAM, PULSE);
fprintf('Median pulse width %.1f ms, median interval %.1f ms\n\n', med_w, med_iti);

%% ======================= A) ANCHORS =======================
fprintf('--- A) Anchors ---\n');
is_pause = false(n, 1);
is_pause(1:n-1) = width_ms(1:n-1) > pause_width_factor * med_w & iti_ms > pause_iti_factor * med_iti;
pause_idx  = find(is_pause);
last_short = width_ms(end) < short_last_ms;

anchors = zeros(0, 2);   % [trigger index, trial number]
for j = 1:numel(pause_idx)
    anchors(end+1, :) = [pause_idx(j), pause_every*j - 1]; %#ok<SAGROW>
    fprintf('Pause %d: trig(%d) at %.2f s, width %.1f ms, then %.0f ms gap  ->  trial %d\n', ...
        j, pause_idx(j), rise(pause_idx(j))/FS, width_ms(pause_idx(j)), iti_ms(pause_idx(j)), pause_every*j - 1);
end
if last_short
    anchors(end+1, :) = [n, n_exp];
    fprintf('Last pulse trig(%d) is %.2f ms wide  ->  final trial %d\n', n, width_ms(end), n_exp);
else
    fprintf('Last pulse trig(%d) is %.1f ms wide -- NOT the short final pulse; trials may be missing at the end.\n', ...
        n, width_ms(end));
end
expected_pauses = floor(n_exp / pause_every);   % pause happens before trial 500, 1000, ...
if numel(pause_idx) ~= expected_pauses
    warning('Found %d pause(s) but %d trials should give %d. Check pause_width_factor / pause_iti_factor.', ...
        numel(pause_idx), n_exp, expected_pauses);
end

trial_of = nan(n, 1);          % trial number of each trigger, where certain
seg = struct('from_trig', {}, 'to_trig', {}, 'missing', {});
prev = [0 0];
for j = 1:size(anchors, 1) + 1
    if j <= size(anchors, 1)
        a = anchors(j, :);
    else
        if last_short, break; end
        a = [n, n_exp];        % open end: whatever is left
    end
    m = (a(2) - prev(2)) - (a(1) - prev(1));
    seg(end+1) = struct('from_trig', prev(1), 'to_trig', a(1), 'missing', m); %#ok<SAGROW>
    if prev(1) == 0
        fprintf('  before trig(%d) [trial %d]: %d missing\n', a(1), a(2), m);
    else
        fprintf('  trig(%d) -> trig(%d): %d missing\n', prev(1), a(1), m);
    end
    if m < 0
        warning('Negative count between anchors -- an anchor is wrong (or there are extra pulses).');
    elseif m == 0
        trial_of(prev(1)+1:a(1)) = prev(2) + (1:(a(1) - prev(1)));
    end
    prev = a;
end
if ~last_short
    fprintf('  (the last segment includes the end of the recording, so missing trials could also be at the end)\n');
end
fprintf('Triggers with a certain trial number: %d of %d\n\n', sum(~isnan(trial_of)), n);

%% ======================= B) DROPPED SAMPLES (time.dat) =======================
fprintf('--- B) time.dat ---\n');
time_jumps = zeros(0, 2);   % [sample index, samples skipped]
tfile = fullfile(data_folder, 'time.dat');
if isfile(tfile)
    fid = fopen(tfile, 'r');
    ts  = fread(fid, inf, 'int32=>double');
    fclose(fid);
    fprintf('%d timestamps (amplifier has %d samples)\n', numel(ts), nSam);
    d  = diff(ts);
    jx = find(d ~= 1);
    time_jumps = [jx, d(jx) - 1];
    clear ts d
    if isempty(jx)
        fprintf('No jumps: every sample was recorded, none dropped.\n');
    else
        fprintf('%d jump(s) in the timestamps:\n', numel(jx));
        for i = 1:min(20, numel(jx))
            s = jx(i);
            where = '';
            p = find(rise(pause_idx) < s & rise(pause_idx + 1) > s, 1);
            if ~isempty(p), where = sprintf('  <-- inside pause %d', p); end
            fprintf('  at %.3f s: %+d samples (%.1f ms)%s\n', s/FS, time_jumps(i, 2), time_jumps(i, 2)/FS*1000, where);
        end
    end
else
    fprintf('No time.dat in this folder -- skipped.\n');
end
if ~isnan(PARAM) && ~isnan(PULSE)
    exp_min = PULSE + 5000 + PARAM;
    for j = 1:numel(pause_idx)
        k = pause_idx(j);
        skipped = 0;
        if ~isempty(time_jumps)
            in = time_jumps(:, 1) > rise(k) & time_jumps(:, 1) < rise(k+1);
            skipped = sum(max(time_jumps(in, 2), 0)) / FS * 1000;   % forward jumps only; negative = counter reset
        end
        fprintf('Pause %d: recorded %.0f ms (+%.0f ms of forward jumps in time.dat) vs >= %.0f ms expected from the code\n', ...
            j, iti_ms(k), skipped, exp_min);
    end
end
fprintf('\n');

%% ======================= C) PARAMETER-UPDATE ARTIFACT =======================
fprintf('--- C) Scanning between triggers ---\n');
[amp_channels, ~] = read_Intan_RHS2000_file;
nChn = numel(amp_channels);
fid = fopen(fullfile(data_folder, 'amplifier.dat'), 'r');
assert(fid > 0, 'Could not open amplifier.dat');
cleaner = onCleanup(@() fclose_if_open(fid));

% Threshold from the quiet stretch just before each trigger
ns = round(noise_ms / 1000 * FS);  pg = max(1, round(pre_gap_ms / 1000 * FS));
B = nan(n, 1);
for k = 1:n
    e = artifact_trace(read_block(fid, rise(k) - pg - ns, ns, nChn, nSam));
    B(k) = max(e);
end
noise_med = median(B, 'omitnan');
noise_rsd = 1.4826 * median(abs(B - noise_med), 'omitnan');
thr = noise_med + noise_k * noise_rsd;
fprintf('Threshold %.2f uV/sample (noise median %.2f, robust SD %.2f, k = %g)\n', thr, noise_med, noise_rsd, noise_k);

% Window k ends just before trigger k; window 1 is the start of the recording
gf = round(guard_fall_ms / 1000 * FS);  gr = round(guard_rise_ms / 1000 * FS);
refr = round(refr_ms / 1000 * FS);      chunk = round(chunk_s * FS);
win_a = zeros(n, 1);  win_b = rise - pg;
win_a(1) = 1;
win_a(2:end) = max(fall(1:end-1) + gf, rise(1:end-1) + gr);

IMG    = nan(n, img_ms, 'single');   % max per 1 ms bin, img_ms before each trigger
events = cell(n, 1);
n_ev   = zeros(n, 1);
lead   = nan(n, 1);                  % ms from the event closest to the trigger, to the trigger
peak   = nan(n, 1);
for k = 1:n
    a = win_a(k);  b = win_b(k);
    if b - a < 3, continue; end
    cross = [];  cval = [];
    for c0 = a:chunk:b
        L = min(chunk + 1, b - c0 + 1);
        if L < 3, continue; end
        e  = artifact_trace(read_block(fid, c0, L, nChn, nSam));
        sx = c0 + (1:numel(e))' - 1;                     % sample index of each value
        ix = find(e(:) > thr);
        cross = [cross; sx(ix)];  cval = [cval; e(ix)']; %#ok<AGROW>
        rel = (sx - rise(k)) / FS * 1000;                % ms relative to trigger k
        keep = rel >= -img_ms & rel < 0;
        if any(keep)
            bins = floor(rel(keep)) + img_ms + 1;
            row  = accumarray(bins, double(e(keep)'), [img_ms 1], @max, NaN);
            cur  = double(IMG(k, :)');
            IMG(k, :) = single(max(cur, row));           % max ignores NaN
        end
    end
    if ~isempty(cross)
        st = [1; find(diff(cross) > refr) + 1];
        en = [st(2:end) - 1; numel(cross)];
        on = cross(st);
        pk = arrayfun(@(i) max(cval(st(i):en(i))), (1:numel(st))');
        events{k} = [on, pk];
        n_ev(k)   = numel(on);
        lead(k)   = (rise(k) - on(end)) / FS * 1000;
        peak(k)   = pk(end);
    end
    if mod(k, 200) == 0, fprintf('  scanned %d / %d\n', k, n); end
end
clear cleaner;

% Window types
wtype = repmat({'normal'}, n, 1);
wtype{1} = 'start';
for k = 2:n
    if is_pause(k-1)
        wtype{k} = 'pause';
    elseif iti_ms(k-1) > dg.settings.gap_factor * med_iti
        wtype{k} = 'gap';
    end
end
isN = strcmp(wtype, 'normal');

fprintf('\nNormal intervals (%d):\n', sum(isN));
fprintf('  with >= 1 event: %d (%.1f%%), with >= 2 events: %d\n', ...
    sum(n_ev(isN) > 0), 100*mean(n_ev(isN) > 0), sum(n_ev(isN) > 1));
ld = lead(isN & n_ev > 0);
if ~isempty(ld)
    fprintf('  event -> trigger: median %.0f ms, 5-95%% range %.0f - %.0f ms\n', ...
        median(ld), prctile(ld, 5), prctile(ld, 95));
    if ~isnan(PARAM)
        fprintf('  within PARAM..PARAM+125 (+25 ms slack for TCP checks): %.1f%%\n', ...
            100*mean(ld >= PARAM & ld <= PARAM + 150));
    end
end

fprintf('\nStart, gaps and pauses (a window with m missing trials should hold m + 1 events):\n');
for k = find(~isN)'
    if k == 1
        lbl = sprintf('start (0 - %.2f s)', rise(1)/FS);
    else
        lbl = sprintf('%-5s trig(%d)->trig(%d) %.0f ms', wtype{k}, k-1, k, iti_ms(k-1));
    end
    fprintf('  %-38s %d event(s)', lbl, n_ev(k));
    if n_ev(k) > 0
        fprintf('  at %s ms before trig(%d)', mat2str(round((rise(k) - events{k}(:,1)')/FS*1000)), k);
    end
    fprintf('  -> suggests %d missing\n', max(n_ev(k) - 1, 0));
end
fprintf('\n');

%% ======================= D) WHEN DOES THE EVENT APPEAR? =======================
fprintf('--- D) Event vs change from the previous trial ---\n');
crosstab_elec = nan(2);  crosstab_any = nan(2);
if isfield(S, 'StimParams') && isfield(S, 'simultaneous_stim') && any(~isnan(trial_of))
    SP  = S.StimParams;  sim = S.simultaneous_stim;
    nT  = floor((size(SP, 1) - 1) / sim);
    elec_key = cell(nT, 1);  full_key = cell(nT, 1);
    for t = 1:nT
        blk = SP(1 + (t-1)*sim + (1:sim), :);
        nm  = cellfun(@tostr, blk(:, 1), 'UniformOutput', false);
        elec_key{t} = strjoin(sort(nm(:))', '|');
        full_key{t} = strjoin(cellfun(@tostr, blk(:)', 'UniformOutput', false), '|');
    end
    use = find(isN & ~isnan(trial_of) & trial_of > 1 & trial_of <= nT);
    t   = trial_of(use);
    ch_e = ~strcmp(elec_key(t), elec_key(t - 1));
    ch_a = ~strcmp(full_key(t), full_key(t - 1));
    has  = n_ev(use) > 0;
    crosstab_elec = [sum(ch_e & has), sum(ch_e & ~has); sum(~ch_e & has), sum(~ch_e & ~has)];
    crosstab_any  = [sum(ch_a & has), sum(ch_a & ~has); sum(~ch_a & has), sum(~ch_a & ~has)];
    fprintf('%d trials with a certain number, in normal intervals:\n', numel(use));
    fprintf('                          event    no event\n');
    fprintf('  electrodes changed      %5d    %5d\n', crosstab_elec(1, :));
    fprintf('  electrodes same         %5d    %5d\n', crosstab_elec(2, :));
    fprintf('  any parameter changed   %5d    %5d\n', crosstab_any(1, :));
    fprintf('  all parameters same     %5d    %5d\n', crosstab_any(2, :));
else
    fprintf('Skipped (no StimParams, or no trials with a certain number).\n');
end

%% ======================= FIGURES =======================
fig_dir = fullfile(data_folder, 'missing_trig_figs');
if ~isfolder(fig_dir), mkdir(fig_dir); end
tx = (-img_ms:-1) + 0.5;

f1 = figure('Color', 'w', 'Position', [80 60 1100 850], 'Name', 'Check: before each trigger');
subplot(4, 1, 1:3);
imagesc(tx, 1:n, log10(max(double(IMG), eps)), 'AlphaData', ~isnan(IMG));
set(gca, 'Color', [0.85 0.85 0.85]);
clim = log10([noise_med, max(thr * 4, prctile(double(IMG(:)), 99.5))]);
if all(isfinite(clim)) && diff(clim) > 0, caxis(clim); end
colormap(gca, parula); cb = colorbar; cb.Label.String = 'log_{10} median |jump| (\muV)';
hold on
if ~isnan(PARAM), xline(-PARAM, 'r--', 'PARAM'); xline(-(PARAM + 125), 'r--', 'PARAM+125'); end
for k = find(~isN)', yline(k, 'm-'); end
ylabel('Trigger index (pink = start/gap/pause)'); title('Signal before each trigger (grey = previous trial / outside window)');
subplot(4, 1, 4); hold on
plot(tx, mean(double(IMG(isN, :)), 1, 'omitnan'), 'k');
yline(thr, 'b--', 'threshold');
if ~isnan(PARAM), xline(-PARAM, 'r--'); xline(-(PARAM + 125), 'r--'); end
xlabel('Time before trigger (ms)'); ylabel('Mean (\muV)'); title('Average over normal intervals'); grid on
saveas(f1, fullfile(fig_dir, '1b_check_before_trigger.png'));

f2 = figure('Color', 'w', 'Position', [120 120 1100 450], 'Name', 'Check: event timing');
subplot(1, 2, 1); hold on
if ~isempty(ld), histogram(ld, 0:10:max(600, max(ld) + 10)); end
if ~isnan(PARAM), xline(PARAM, 'r--', 'PARAM'); xline(PARAM + 125, 'r--', 'PARAM+125'); end
xlabel('Event -> trigger (ms)'); ylabel('Count'); title('Normal intervals'); grid on
subplot(1, 2, 2); hold on
plot(find(isN), n_ev(isN), '.', 'Color', [0.5 0.5 0.5]);
k = find(~isN);
plot(k, n_ev(k), 'mo', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Window before trigger index'); ylabel('Events found');
title('Events per window (pink = start/gap/pause)'); grid on
saveas(f2, fullfile(fig_dir, '1b_check_event_timing.png'));

%% ======================= SAVE =======================
report.check = struct( ...
    'run_time', datetime('now'), 'diagnose_run_time', dg.run_time, ...
    'settings', struct('noise_k', noise_k, 'guard_fall_ms', guard_fall_ms, 'guard_rise_ms', guard_rise_ms, ...
        'refr_ms', refr_ms, 'pause_every', pause_every, 'pause_width_factor', pause_width_factor, ...
        'pause_iti_factor', pause_iti_factor, 'short_last_ms', short_last_ms), ...
    'PARAM', PARAM, 'PULSE', PULSE, ...
    'pause_idx', pause_idx, 'last_short', last_short, 'anchors', anchors, 'segments', seg, ...
    'trial_of', trial_of, 'time_jumps', time_jumps, ...
    'threshold', thr, 'window_type', {wtype}, 'n_events', n_ev, 'events', {events}, ...
    'lead_ms', lead, 'peak', peak, ...
    'crosstab_elec', crosstab_elec, 'crosstab_any', crosstab_any);
save(report_file, 'report');
fprintf('\nSaved to report.check in %s\nFigures in %s\n', report_file, fig_dir);

%% ======================= LOCAL FUNCTIONS =======================
function X = read_block(fid, s0, L, nChn, nSam)
% L samples of all channels from sample s0 (1-based), in uV; NaN outside the file.
    X = nan(nChn, L);
    a = max(s0, 1);  b = min(s0 + L - 1, nSam);
    if b < a, return; end
    fseek(fid, (a - 1) * nChn * 2, 'bof');
    D = fread(fid, [nChn, b - a + 1], 'int16=>double') * 0.195;
    X(:, (a - s0 + 1):(a - s0 + size(D, 2))) = D;
end

function e = artifact_trace(X)
% Median across channels of the absolute sample-to-sample jump.
    e = median(abs(diff(X, 1, 2)), 1, 'omitnan');
end

function s = tostr(x)
    if ischar(x),        s = x;
    elseif isstring(x),  s = char(x);
    elseif isnumeric(x) || islogical(x), s = mat2str(x);
    else,                s = class(x);
    end
end

function fclose_if_open(fid)
    if ~isempty(fopen(fid)), fclose(fid); end
end