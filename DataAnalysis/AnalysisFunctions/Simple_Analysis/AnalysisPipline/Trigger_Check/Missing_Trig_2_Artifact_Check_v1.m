%% MISSING_TRIG_2_RECOVER  (version 2 -- reset counting)
%
% Step 2 of 3 in the missing-trigger repair workflow:
%   1) missing_trig_1_diagnose.m   (recount pulses + locate gaps)
%   2) missing_trig_2_recover.m    <-- this script (work out which trials are missing)
%   3) missing_trig_3_apply.m      (back up and rewrite .trig.dat)
%
% How it works
%   RunMultiElect_exp.m sends every trial's parameters to the Intan ('SET'),
%   which makes the Intan stop and restart acquisition. Each restart resets
%   the sample counter in time.dat. So time.dat holds exactly one reset per
%   trial -- including trials whose TTL never fired. Counting the resets
%   before each detected trigger gives its trial number. A reset with no
%   trigger after it is a missing trial: the parameters were sent but the
%   TTL (which also drives the stimulator) never fired, so that trial was
%   not stimulated and gets -500.
%
%   The recording is split into segments by fixed anchors from the run
%   script: the last pulse before each pause is trial 499, 999, ... (the
%   pause pulse is very wide and followed by a long gap), and a very short
%   final pulse is the last trial. Each segment is solved on its own, so a
%   problem in one place (e.g. an extra reset from a TCP resend) cannot
%   shift trials anywhere else.
%
% Checks
%   - resets in each segment = trials in that segment
%   - at most one trigger per reset, last trigger of a segment = its anchor
%   - reset -> trigger time is roughly PARAM..PARAM+125 ms
%   - nothing left over after the final trial
%
% Also saved: for every trial, how much continuous recording follows the
%   stimulus before the next restart (relevant for analysis windows).
%
% Output: report.recover (incl. proposed trig) in missing_trig_report.mat,
%         figure missing_trig_figs/2_recover_summary.png.
% Read-only -- does NOT modify .trig.dat or any recording file.

clear; clc;

%% ======================= SETTINGS =======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX039/Xia_Linearity_S3_100_400um_SimSeq1_261009_190544';   % <-- same folder as step 1

pause_every        = 500;  % rem(count,500) in RunMultiElect_exp.m
pause_width_factor = 2.5;  % pause pulse: width > this x median width ...
pause_iti_factor   = 4;    % ... AND the next interval > this x median interval
short_last_ms      = 5;    % last pulse shorter than this = final trial
lead_tol_ms        = 60;   % tolerance on PARAM..PARAM+125 (time.dat resets come in ~34 ms blocks)
lead_bad_frac      = 0.02; % more than this fraction outside the range -> needs review

% Manual decisions, applied last (leave empty normally):
% rows of [trial, raw_rise_sample] or [trial, -500]
manual_trials = [];

%% ======================= LOAD =======================
assert(isfolder(data_folder), 'data_folder not found:\n%s', data_folder);
cd(data_folder);
report_file = fullfile(data_folder, 'missing_trig_report.mat');
assert(isfile(report_file), 'No missing_trig_report.mat here -- run missing_trig_1_diagnose.m first.');
load(report_file, 'report');
assert(isfield(report, 'diagnose'), 'Report has no step 1 results -- run missing_trig_1_diagnose.m first.');
if isfield(report, 'apply')
    warning('Step 3 has already rewritten .trig.dat for this folder. Results below are for review only.');
end
info = report.info;  dg = report.diagnose;

FS    = info.FS;
nSam  = info.nSam;
n_exp = info.n_expected;
rise  = double(dg.rise(:));
fall  = double(dg.fall(:));
edge  = double(dg.new_trig(:));        % edge step 1 matched to the old trig file
n     = numel(rise);
width_ms = (fall - rise + 1) / FS * 1000;
iti_ms   = diff(rise) / FS * 1000;
med_w    = median(width_ms);
med_iti  = median(iti_ms);
amps     = info.trialAmps(:);
if numel(amps) < n_exp, amps(end+1:n_exp) = NaN; end

exp_path = fullfile(data_folder, info.exp_file);
vlist = whos('-file', exp_path);
want  = intersect({'PARAM', 'PULSE', 'StimParams', 'simultaneous_stim'}, {vlist.name});
S = load(exp_path, want{:});
PARAM = NaN;  PULSE = NaN;
if isfield(S, 'PARAM'), PARAM = double(S.PARAM); end
if isfield(S, 'PULSE'), PULSE = double(S.PULSE); end
elec = repmat({''}, n_exp, 1);
if isfield(S, 'StimParams') && isfield(S, 'simultaneous_stim')
    sim = S.simultaneous_stim;
    for t = 1:min(n_exp, floor((size(S.StimParams, 1) - 1) / sim))
        nm = cellfun(@tostr, S.StimParams(1 + (t-1)*sim + (1:sim), 1), 'UniformOutput', false);
        elec{t} = strjoin(nm(:)', ' > ');
    end
end

fprintf('Folder: %s\n', data_folder);
fprintf('Expected %d trials, detected %d triggers (missing %d). PARAM = %g ms, PULSE = %g ms\n\n', ...
    n_exp, n, n_exp - n, PARAM, PULSE);

%% ======================= TRIG OFFSET (match the old file) =======================
old_trig = double(info.old_trig(:));
old_real = old_trig(old_trig > 0);
if numel(old_real) >= 2
    offs        = round(nearest_offset(old_real, edge));
    trig_offset = mode(offs);
    offset_frac = mean(offs == trig_offset);
    fprintf('Old .trig.dat = %s edge + %d samples for %.1f%% of its triggers -> same offset used here\n', ...
        dg.edge_used, trig_offset, offset_frac*100);
    if offset_frac < 0.95
        warning('Old trig offset is not consistent (%.1f%%). Check before applying.', offset_frac*100);
    end
else
    trig_offset = 0;  offset_frac = NaN;
    warning('No usable old .trig.dat -- using the raw edge (offset 0).');
end

%% ======================= RESETS FROM time.dat =======================
tfile = fullfile(data_folder, 'time.dat');
have_resets = isfile(tfile);
resets = zeros(0, 1);  n_drop_jumps = 0;
if have_resets
    fid = fopen(tfile, 'r');
    ts  = fread(fid, inf, 'int32=>double');
    fclose(fid);
    if numel(ts) ~= nSam
        warning('time.dat has %d samples but amplifier has %d.', numel(ts), nSam);
    end
    d = diff(ts);
    resets = find(d < 0) + 1;          % first sample after each restart
    n_drop_jumps = sum(d > 1);
    clear ts d
    fprintf('time.dat: %d counter resets (expected one per trial = %d)', numel(resets), n_exp);
    if n_drop_jumps > 0, fprintf(', plus %d forward jump(s) (dropped samples)', n_drop_jumps); end
    fprintf('\n');
else
    warning('No time.dat -- only the anchors can be used.');
end

%% ======================= ANCHORS =======================
is_pause = false(n, 1);
is_pause(1:n-1) = width_ms(1:n-1) > pause_width_factor * med_w & iti_ms > pause_iti_factor * med_iti;
pause_idx  = find(is_pause);
last_short = width_ms(end) < short_last_ms;
anchors = [pause_idx(:), pause_every * (1:numel(pause_idx))' - 1];
if last_short, anchors(end+1, :) = [n, n_exp]; end
fprintf('Anchors: %d pause(s)%s\n', numel(pause_idx), ...
    ternary(last_short, sprintf(' + short final pulse (trig %d = trial %d)', n, n_exp), ' (last pulse is NOT short)'));
for j = 1:numel(pause_idx)
    fprintf('  trig(%d) at %.2f s = trial %d\n', pause_idx(j), rise(pause_idx(j))/FS, anchors(j, 2));
end
if numel(pause_idx) ~= floor(n_exp / pause_every)
    warning('Found %d pause(s) but %d trials should give %d.', numel(pause_idx), n_exp, floor(n_exp / pause_every));
end
fprintf('\n');

%% ======================= SOLVE EACH SEGMENT =======================
% Boundaries: [trigger index, trial number, sample]. Segment j covers
% triggers b(j)+1 .. b(j+1) and trials t(j)+1 .. t(j+1).
B = [0, 0, 0];
for j = 1:size(anchors, 1)
    B(end+1, :) = [anchors(j, 1), anchors(j, 2), rise(anchors(j, 1))]; %#ok<SAGROW>
end
B(end+1, :) = [n, n_exp, nSam];       % tail: rest of the recording (empty if last pulse is an anchor)

trial_of     = nan(n, 1);
reset_sample = nan(n_exp, 1);          % restart that belongs to each trial
seg_info = struct('first_trig', {}, 'last_trig', {}, 'first_trial', {}, 'last_trial', {}, ...
                  'missing', {}, 'resets', {}, 'method', {}, 'ok', {});
problems = {};
fprintf('--- Segments ---\n');
for j = 1:size(B, 1) - 1
    a0 = B(j, 1);  t0 = B(j, 2);  s0 = B(j, 3);
    a1 = B(j+1, 1);  t1 = B(j+1, 2);  s1 = B(j+1, 3);
    if j + 1 == size(B, 1) && last_short
        % after the short final pulse: no triggers, no trials -- just check for leftover resets
        s0 = rise(n);  a0 = n;  t0 = n_exp;
    end
    trig_ix  = (a0 + 1):a1;
    nTrig    = numel(trig_ix);
    nTrial   = t1 - t0;
    m        = nTrial - nTrig;
    R        = resets(resets > s0 & resets <= s1);
    method   = '';  ok = true;

    if nTrig == 0 && nTrial == 0
        if ~isempty(R)
            problems{end+1} = sprintf('%d reset(s) after the final trial', numel(R)); %#ok<SAGROW>
            ok = false;
        end
        seg_info(end+1) = struct('first_trig', NaN, 'last_trig', NaN, 'first_trial', NaN, 'last_trial', NaN, ...
            'missing', 0, 'resets', numel(R), 'method', 'after final trial', 'ok', ok); %#ok<SAGROW>
        fprintf('  after final trial: %d reset(s) %s\n', numel(R), ternary(ok, '(ok)', '<-- unexpected'));
        continue
    end

    if m < 0
        problems{end+1} = sprintf('trials %d-%d: %d more triggers than trials -- an anchor is wrong', t0+1, t1, -m); %#ok<SAGROW>
        ok = false;  method = 'impossible';
    elseif have_resets && numel(R) == nTrial
        % map by counting resets
        local = arrayfun(@(k) sum(R <= rise(k)), trig_ix(:));
        mapped = t0 + local;
        if any(local < 1) || numel(unique(local)) < nTrig || (a1 > a0 && j + 1 < size(B, 1) && mapped(end) ~= t1)
            ok = false;  method = 'resets (inconsistent)';
            problems{end+1} = sprintf('trials %d-%d: reset counting gives an inconsistent order', t0+1, t1); %#ok<SAGROW>
        else
            trial_of(trig_ix) = mapped;
            reset_sample(t0 + (1:nTrial)) = R;
            method = 'resets';
        end
    elseif m == 0
        % nothing missing here; resets only used as a check
        trial_of(trig_ix) = t0 + (1:nTrig)';
        method = 'anchors (no missing)';
        if have_resets
            problems{end+1} = sprintf('trials %d-%d: %d resets for %d trials (nothing missing, so mapping is still certain)', ...
                t0+1, t1, numel(R), nTrial); %#ok<SAGROW>
        end
    else
        ok = false;  method = 'unresolved';
        if have_resets
            problems{end+1} = sprintf('trials %d-%d: %d missing but %d resets for %d trials', t0+1, t1, m, numel(R), nTrial); %#ok<SAGROW>
        else
            problems{end+1} = sprintf('trials %d-%d: %d missing and no time.dat', t0+1, t1, m); %#ok<SAGROW>
        end
    end
    seg_info(end+1) = struct('first_trig', a0 + 1, 'last_trig', a1, 'first_trial', t0 + 1, 'last_trial', t1, ...
        'missing', m, 'resets', numel(R), 'method', method, 'ok', ok); %#ok<SAGROW>
    fprintf('  trig %4d-%4d = trials %4d-%4d: %d missing, %d resets -> %s\n', ...
        a0 + 1, a1, t0 + 1, t1, m, numel(R), method);
end

% Problems flagged on segments with nothing missing are notes, not failures
hard = ~[seg_info.ok];

%% ======================= TIMING CHECKS =======================
known = find(~isnan(trial_of));
lead_ms = nan(n, 1);
has_rs  = known(~isnan(reset_sample(trial_of(known))));
lead_ms(has_rs) = (rise(has_rs) - reset_sample(trial_of(has_rs))) / FS * 1000;
bad_lead = [];
if ~isnan(PARAM) && ~isempty(has_rs)
    lo = PARAM - lead_tol_ms;  hi = PARAM + 125 + lead_tol_ms;
    bad_lead = has_rs(lead_ms(has_rs) < lo | lead_ms(has_rs) > hi);
    fprintf('\nReset -> trigger: median %.0f ms (5-95%%: %.0f - %.0f ms); expected %g-%g ms (+/-%g)\n', ...
        median(lead_ms(has_rs)), prctile(lead_ms(has_rs), 5), prctile(lead_ms(has_rs), 95), PARAM, PARAM + 125, lead_tol_ms);
    fprintf('Outside the expected range: %d of %d\n', numel(bad_lead), numel(has_rs));
    for k = bad_lead(1:min(end, 10))'
        fprintf('  trig(%d) = trial %d: %.0f ms\n', k, trial_of(k), lead_ms(k));
    end
    if numel(bad_lead) > lead_bad_frac * numel(has_rs)
        problems{end+1} = sprintf('%d triggers have an unusual reset -> trigger time', numel(bad_lead));
        hard(end+1) = true;
    end
end

% Continuous recording after each stimulus (until the next restart)
cont_ms = nan(n_exp, 1);
if have_resets
    for k = known(:)'
        nx = resets(find(resets > rise(k), 1));
        if isempty(nx), nx = nSam + 1; end
        cont_ms(trial_of(k)) = (nx - rise(k)) / FS * 1000;
    end
    c = cont_ms(~isnan(cont_ms));
    fprintf('\nContinuous recording after each stimulus: median %.0f ms, 5th percentile %.0f ms, minimum %.0f ms\n', ...
        median(c), prctile(c, 5), min(c));
    fprintf('(analysis windows longer than this cross a restart -- see report.recover.cont_ms per trial)\n');
end

%% ======================= RESULT =======================
proposed = nan(1, n_exp);
proposed(trial_of(known)) = edge(known) + trig_offset;
missing_trials = find(isnan(proposed));
resolved = all(~hard);
if resolved
    proposed(missing_trials) = -500;
end

% Manual decisions
for i = 1:size(manual_trials, 1)
    t = manual_trials(i, 1);  v = manual_trials(i, 2);
    assert(t >= 1 && t <= n_exp, 'manual_trials: trial %d out of range.', t);
    if v == -500, proposed(t) = -500; else, proposed(t) = v + trig_offset; end
end
if ~resolved && ~isempty(manual_trials) && ~any(isnan(proposed))
    fprintf('\nManual decisions fill every unresolved trial.\n');
    resolved = true;
end

% Missing-trial table
mt = missing_trials(:);
prev_trig = nan(numel(mt), 1);  next_trig = nan(numel(mt), 1);
for i = 1:numel(mt)
    p = find(trial_of < mt(i), 1, 'last');  q = find(trial_of > mt(i), 1, 'first');
    if ~isempty(p), prev_trig(i) = p; end
    if ~isempty(q), next_trig(i) = q; end
end
missing = table(mt, prev_trig, next_trig, reset_sample(mt) / FS, amps(mt), elec(mt), ...
    'VariableNames', {'trial', 'after_trig', 'before_trig', 'reset_time_s', 'amp', 'electrodes'});

pos = proposed(proposed > 0);
if resolved && ~all(diff(pos) > 0)
    problems{end+1} = 'proposed trigger times are not in increasing order';
    resolved = false;
end
if resolved, status = 'pass'; else, status = 'needs_review'; end

fprintf('\n--- %s (%d) ---\n', ternary(resolved, 'Missing trials', 'Trials without a trigger (NOT all resolved)'), numel(mt));
for i = 1:height(missing)
    fprintf('  trial %4d  between trig(%d) and trig(%d)', missing.trial(i), missing.after_trig(i), missing.before_trig(i));
    if ~isnan(missing.reset_time_s(i)), fprintf(', restart at %.3f s', missing.reset_time_s(i)); end
    fprintf(', amp %g, %s\n', missing.amp(i), missing.electrodes{i});
end
fprintf('\nSTATUS: %s\n', upper(status));
for i = 1:numel(problems), fprintf('  - %s\n', problems{i}); end
if strcmp(status, 'pass')
    fprintf('Missing trials will be written as -500. Check the figure, then run missing_trig_3_apply.m.\n');
else
    fprintf('Not ready to apply. Send the output; manual_trials can fill unresolved trials if needed.\n');
end

%% ======================= FIGURE =======================
fig_dir = fullfile(data_folder, 'missing_trig_figs');
if ~isfolder(fig_dir), mkdir(fig_dir); end
f = figure('Color', 'w', 'Position', [100 100 1300 420], 'Name', 'Recover v2');
subplot(1, 3, 1); hold on
plot(1:n, trial_of - (1:n)', 'k.');
for j = 1:size(anchors, 1), xline(anchors(j, 1), 'b:'); end
xlabel('Detected trigger index'); ylabel('Trial number - trigger index');
title('Each step up = a missing trial (blue = anchors)'); grid on
subplot(1, 3, 2); hold on
if ~isempty(has_rs), histogram(lead_ms(has_rs), 0:5:max(400, max(lead_ms(has_rs)) + 5)); end
if ~isnan(PARAM), xline(PARAM, 'r--', 'PARAM'); xline(PARAM + 125, 'r--', 'PARAM+125'); end
xlabel('Restart -> trigger (ms)'); ylabel('Count'); title('Timing check'); grid on
subplot(1, 3, 3);
if any(~isnan(cont_ms)), histogram(cont_ms(~isnan(cont_ms)), 50); end
xlabel('Continuous data after stimulus (ms)'); ylabel('Trials'); title('Until the next restart'); grid on
saveas(f, fullfile(fig_dir, '2_recover_summary.png'));

%% ======================= SAVE =======================
report.recover = struct( ...
    'version', 2, 'run_time', datetime('now'), 'diagnose_run_time', dg.run_time, ...
    'settings', struct('pause_every', pause_every, 'pause_width_factor', pause_width_factor, ...
        'pause_iti_factor', pause_iti_factor, 'short_last_ms', short_last_ms, ...
        'lead_tol_ms', lead_tol_ms, 'manual_trials', manual_trials), ...
    'PARAM', PARAM, 'PULSE', PULSE, 'trig_offset', trig_offset, 'offset_frac', offset_frac, ...
    'resets', resets, 'anchors', anchors, 'segments', seg_info, ...
    'trial_of', trial_of, 'reset_sample', reset_sample, 'lead_ms', lead_ms, 'bad_lead', bad_lead, ...
    'cont_ms', cont_ms, 'missing', missing, ...
    'proposed_trig', proposed, 'status', status);
report.recover.problems = problems;
save(report_file, 'report');
fprintf('\nSaved report: %s\nFigure:       %s\n', report_file, fullfile(fig_dir, '2_recover_summary.png'));

%% ======================= LOCAL FUNCTIONS =======================
function off = nearest_offset(x, edges)
    edges = double(edges(:));
    nn  = interp1(edges, edges, double(x(:)), 'nearest', 'extrap');
    off = double(x(:)) - nn;
end

function s = tostr(x)
    if ischar(x),        s = x;
    elseif isstring(x),  s = char(x);
    elseif isnumeric(x) || islogical(x), s = mat2str(x);
    else,                s = class(x);
    end
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end