%% MISSING_TRIG_1_DIAGNOSE
%
% Step 1 of 3 in the missing-trigger repair workflow:
%   1) missing_trig_1_diagnose.m   <-- this script (recount + locate gaps)
%   2) missing_trig_2_recover.m    (find real times in amplifier.dat + check)
%   3) missing_trig_3_apply.m      (back up and rewrite .trig.dat)
%
% Only run this on recordings where the trigger count is wrong.
%
% What it does:
%   A) RECOUNT - re-reads digitalin.dat and detects pulses using ONLY the
%      stim bit. cleanTrig_sabquick uses (digital_in == 1), which misses
%      any pulse where another digital pin happens to be high at the same
%      time. If the recount finds all expected pulses, we're done.
%   B) LOCATE - if pulses are still missing, finds unusually long
%      intervals between triggers, estimates how many pulses each one
%      hides, and checks whether that accounts for the whole shortfall
%      (if not, some are probably missing at the start or end).
%
% Output: missing_trig_report.mat in the data folder (read by step 2),
%         plus figures in data_folder/missing_trig_figs/.
%
% Diagnostic only -- does NOT modify .trig.dat or any recording file.
% Rerunning this script starts a fresh report (any step 2 results are
% cleared, so step 3 can never apply a fix based on outdated gaps).

clear; clc;

%% ======================= SETTINGS =======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX039/Xia_Linearity_S3_100_400um_SimSeq1_261009_190544';   % <-- set this

stim_bit          = 0;    % digital pin carrying the stim TTL (bit 0 in cleanTrig_sabquick)
min_pulse_samples = 2;    % high stretches shorter than this are treated as noise glitches
gap_factor        = 1.5;  % interval > gap_factor x median ITI  -> gap (missing pulse/s)
short_factor      = 0.5;  % interval < short_factor x median ITI -> suspicious extra pulse
pause_every       = 500;  % run script pauses every N trials (pause(5)); 0 = ignore
pause_window      = 5;    % flag a gap as a possible pause if within +/- this many trials

%% ======================= LOAD EXPERIMENT INFO =======================
assert(isfolder(data_folder), 'data_folder not found:\n%s', data_folder);
cd(data_folder);
fprintf('Data folder:\n  %s\n\n', data_folder);

report_file = fullfile(data_folder, 'missing_trig_report.mat');
if isfile(report_file)
    prev = load(report_file);
    if isfield(prev, 'report') && isfield(prev.report, 'apply')
        warning(['A previous report shows step 3 (Apply) already rewrote .trig.dat ' ...
                 'for this folder. The "original" trig recorded below is the PATCHED ' ...
                 'version. Restore the backup first if you want to start from scratch.']);
    end
    fprintf('Existing report found -- it will be replaced (step 2/3 results cleared).\n\n');
end

[~, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found in data_folder.');
S = load(fullfile(data_folder, fileDIR(1).name), 'StimParams', 'simultaneous_stim', 'n_Trials');
n_expected        = S.n_Trials;
simultaneous_stim = S.simultaneous_stim;

trialAmps_all = cell2mat(S.StimParams(2:end, 16));
trialAmps     = trialAmps_all(1:simultaneous_stim:end);
trialAmps     = trialAmps(:);
is_catch      = (trialAmps == -1);
if numel(trialAmps) ~= n_expected
    warning('StimParams gives %d trials but n_Trials = %d. Using n_Trials as expected count.', ...
        numel(trialAmps), n_expected);
end
fprintf('Expected triggers (n_Trials): %d   (%d catch trials)\n', n_expected, sum(is_catch));

% Existing trig file (what the rest of the pipeline currently uses)
trigDIR  = dir(fullfile(data_folder, '*.trig.dat'));
old_trig = [];
if ~isempty(trigDIR)
    old_trig = double(loadTrig(0));
    old_trig = old_trig(:);
    fprintf('Existing %s: %d entries (%d are -500 placeholders, %d other non-positive)\n', ...
        trigDIR(1).name, numel(old_trig), sum(old_trig == -500), ...
        sum(old_trig <= 0 & old_trig ~= -500));
else
    fprintf('No .trig.dat found -- comparing against expected count only.\n');
end
old_real = old_trig(old_trig > 0);

%% ======================= A) RECOUNT FROM digitalin.dat =======================
fid = fopen(fullfile(data_folder, 'digitalin.dat'), 'r');
assert(fid > 0, 'Could not open digitalin.dat');
digital_in = fread(fid, inf, 'uint16=>uint16');
fclose(fid);
nSam = numel(digital_in);
fprintf('Recording length: %.1f s (%d samples)\n\n', nSam/FS, nSam);

stim_line = bitand(digital_in, uint16(2^stim_bit)) > 0;
d    = diff([false; stim_line; false]);
rise = find(d == 1);          % first high sample of each pulse
fall = find(d == -1) - 1;     % last high sample of each pulse
width = fall - rise + 1;

keep     = width >= min_pulse_samples;
n_glitch = sum(~keep);
rise = rise(keep);  fall = fall(keep);  width = width(keep);

% Compare with the old exact-match method
exact_line = (digital_in == 1);
n_exact    = sum(diff([false; exact_line]) == 1);
mixed      = stim_line & ~exact_line;          % stim high but another pin also high
cs0        = [0; cumsum(double(mixed))];
mixed_cnt  = cs0(fall + 1) - cs0(rise);        % mixed samples within each pulse
n_mixed_pulses  = sum(mixed_cnt > 0);
n_fully_masked  = sum(mixed_cnt == width);     % invisible to (digital_in == 1)

fprintf('--- A) Recount ---\n');
fprintf('Old method (digital_in == 1):   %d pulses\n', n_exact);
fprintf('Stim bit only (bit %d):          %d pulses', stim_bit, numel(rise));
if n_glitch > 0, fprintf('  (%d glitches < %d samples ignored)', n_glitch, min_pulse_samples); end
fprintf('\n');
fprintf('Pulses overlapping another active pin: %d (of which %d fully hidden from old method)\n', ...
    n_mixed_pulses, n_fully_masked);
fprintf('Pulse width: median %.2f ms, range %.2f - %.2f ms\n', ...
    median(width)/FS*1000, min(width)/FS*1000, max(width)/FS*1000);

% Which edge does the existing trig file use? Match it so values stay compatible.
edge_used = 'rise';
edge_offset_samples = NaN;
if numel(old_real) >= 2 && numel(rise) >= 2
    dr = nearest_offset(old_real, rise);
    df = nearest_offset(old_real, fall);
    if median(abs(df)) < median(abs(dr))
        edge_used = 'fall';  edge_offset_samples = median(df);
    else
        edge_offset_samples = median(dr);
    end
    fprintf('Existing trig matches the %s edge (median offset %g samples)\n', ...
        upper(edge_used), edge_offset_samples);
    if abs(edge_offset_samples) > 2
        warning(['Existing trig values do not sit on either edge of the stim pulses ' ...
                 '(offset %g samples). Check before applying any fix.'], edge_offset_samples);
    end
else
    fprintf('No usable existing trig to match -- using RISING edge.\n');
end
if strcmp(edge_used, 'fall'), new_trig = fall; else, new_trig = rise; end
new_trig  = double(new_trig(:));
n_found   = numel(new_trig);
shortfall = n_expected - n_found;
fprintf('\nExpected %d, found %d  ->  shortfall = %d\n\n', n_expected, n_found, shortfall);

%% ======================= B) LOCATE GAPS =======================
iti_ms  = diff(new_trig) / FS * 1000;
med_iti = median(iti_ms);

gap_idx     = find(iti_ms > gap_factor * med_iti);     % gap is between trigger k and k+1
est_missing = max(round(iti_ms(gap_idx) / med_iti) - 1, 1);
short_idx   = find(iti_ms < short_factor * med_iti);

% Provisional trial number of each detected trigger, ASSUMING nothing is
% missing before the first detected trigger (step 2 will test that).
add = zeros(n_found, 1);
add(gap_idx + 1) = est_missing;
trial_of = (1:n_found)' + cumsum(add);

nG = numel(gap_idx);
possible_pause  = false(nG, 1);
missing_trials  = cell(nG, 1);
missing_catch   = cell(nG, 1);
for g = 1:nG
    k  = gap_idx(g);
    tk = trial_of(k);
    if pause_every > 0 && tk >= pause_every - pause_window
        possible_pause(g) = abs(tk - pause_every*round(tk/pause_every)) <= pause_window;
    end
    mt = tk + (1:est_missing(g));
    missing_trials{g} = mt;
    mt_valid = mt(mt <= numel(is_catch));
    missing_catch{g}  = mt_valid(is_catch(mt_valid));
end

gaps = table((1:nG)', gap_idx, new_trig(gap_idx)/FS, iti_ms(gap_idx), est_missing, ...
    possible_pause, missing_trials, missing_catch, ...
    'VariableNames', {'gap', 'after_trig', 'time_s', 'interval_ms', 'est_missing', ...
                      'possible_pause', 'provisional_trials', 'catch_trials'});

lead_ms = new_trig(1) / FS * 1000;
tail_ms = (nSam - new_trig(end)) / FS * 1000;
residual = shortfall - sum(est_missing);

fprintf('--- B) Locate ---\n');
fprintf('Median ITI: %.1f ms   (normal range %.1f - %.1f ms, excluding gaps)\n', med_iti, ...
    min(iti_ms(iti_ms <= gap_factor*med_iti)), max(iti_ms(iti_ms <= gap_factor*med_iti)));
fprintf('Time before first trigger: %.2f s   after last trigger: %.2f s\n', lead_ms/1000, tail_ms/1000);

if nG == 0
    fprintf('No gaps found.\n');
else
    fprintf('\n%d gap(s) found:\n', nG);
    for g = 1:nG
        fprintf('  Gap %d: after trig(%d) at %.2f s, interval %.0f ms (~%.1fx median) -> ~%d missing', ...
            g, gap_idx(g), gaps.time_s(g), gaps.interval_ms(g), gaps.interval_ms(g)/med_iti, est_missing(g));
        if possible_pause(g), fprintf('  [near the %d-trial pause]', pause_every); end
        fprintf('\n      provisional trial(s): %s', mat2str(missing_trials{g}));
        if ~isempty(missing_catch{g}), fprintf('   catch: %s', mat2str(missing_catch{g})); end
        fprintf('\n');
    end
end
if ~isempty(short_idx)
    fprintf('\nWARNING: %d suspiciously SHORT interval(s) (< %.0f ms) -- possible extra/double pulses:\n', ...
        numel(short_idx), short_factor*med_iti);
    for s = short_idx(:)'
        fprintf('  between trig(%d) and trig(%d): %.0f ms at %.2f s\n', s, s+1, iti_ms(s), new_trig(s)/FS);
    end
end

% Status
if shortfall == 0 && isempty(short_idx)
    status = 'resolved';
    msg = 'Stim-bit recount finds every expected pulse. Go straight to step 3 (Apply).';
elseif shortfall < 0
    status = 'extra_pulses';
    msg = sprintf('Found %d MORE pulses than expected -- check the short intervals / glitch setting.', -shortfall);
elseif residual == 0
    status = 'gaps_located';
    msg = 'Gaps account for the whole shortfall. Continue with step 2 (Recover).';
elseif residual > 0
    status = 'count_mismatch';
    msg = sprintf(['Gaps explain %d of %d missing; %d unaccounted for -- likely at the START ' ...
                   'or END of the recording. Step 2 will search there.'], sum(est_missing), shortfall, residual);
else
    status = 'count_mismatch';
    msg = sprintf(['Gaps suggest %d missing but only %d are -- some gaps are probably deliberate ' ...
                   'pauses rather than lost pulses (check the [pause] flags).'], sum(est_missing), shortfall);
end
fprintf('\nSTATUS: %s\n%s\n', upper(status), msg);

%% ======================= FIGURES =======================
fig_dir = fullfile(data_folder, 'missing_trig_figs');
if ~isfolder(fig_dir), mkdir(fig_dir); end

fh = figure('Color', 'w', 'Position', [100 100 1200 650], 'Name', 'Diagnose: intervals');
subplot(2, 1, 1); hold on
plot(1:numel(iti_ms), iti_ms, '.', 'Color', [0.4 0.4 0.4]);
yline(med_iti, 'b-', 'median');
yline(gap_factor*med_iti, 'r--', 'gap threshold');
if ~isempty(gap_idx),   plot(gap_idx, iti_ms(gap_idx), 'ro', 'MarkerSize', 9, 'LineWidth', 1.5); end
if ~isempty(short_idx), plot(short_idx, iti_ms(short_idx), 'ms', 'MarkerSize', 9, 'LineWidth', 1.5); end
xlabel('Interval index k (between trig(k) and trig(k+1))'); ylabel('Interval (ms)');
title(sprintf('Inter-trigger intervals  |  expected %d, found %d, shortfall %d  |  %s', ...
    n_expected, n_found, shortfall, status), 'Interpreter', 'none');
set(gca, 'YScale', 'log'); grid on

subplot(2, 1, 2);
histogram(iti_ms(iti_ms <= gap_factor*med_iti), 50);
xlabel('Interval (ms)'); ylabel('Count'); title('Normal interval distribution (gaps excluded)'); grid on
saveas(fh, fullfile(fig_dir, '1_diagnose_intervals.png'));

%% ======================= SAVE REPORT =======================
report = struct();
report.info = struct( ...
    'data_folder', data_folder, 'exp_file', fileDIR(1).name, 'FS', FS, 'nSam', nSam, ...
    'n_expected', n_expected, 'simultaneous_stim', simultaneous_stim, ...
    'trialAmps', trialAmps, 'is_catch', is_catch, ...
    'trig_file', '', 'old_trig', old_trig, 'n_old_entries', numel(old_trig), ...
    'created', datetime('now'));
if ~isempty(trigDIR), report.info.trig_file = trigDIR(1).name; end

report.diagnose = struct( ...
    'run_time', datetime('now'), ...
    'settings', struct('stim_bit', stim_bit, 'min_pulse_samples', min_pulse_samples, ...
        'gap_factor', gap_factor, 'short_factor', short_factor, ...
        'pause_every', pause_every, 'pause_window', pause_window), ...
    'n_exact_method', n_exact, 'n_found', n_found, 'n_glitch', n_glitch, ...
    'n_mixed_pulses', n_mixed_pulses, 'n_fully_masked', n_fully_masked, ...
    'edge_used', edge_used, 'edge_offset_samples', edge_offset_samples, ...
    'rise', rise, 'fall', fall, 'new_trig', new_trig, ...
    'iti_ms', iti_ms, 'median_iti_ms', med_iti, ...
    'gaps', gaps, 'short_idx', short_idx, ...
    'lead_ms', lead_ms, 'tail_ms', tail_ms, ...
    'shortfall', shortfall, 'residual', residual, ...
    'status', status, 'message', msg);

save(report_file, 'report');
fprintf('\nSaved report:  %s\nSaved figure:  %s\n', report_file, fullfile(fig_dir, '1_diagnose_intervals.png'));

%% ======================= LOCAL FUNCTIONS =======================
function off = nearest_offset(x, edges)
% Signed offset from each value in x to the nearest value in edges.
    edges = double(edges(:));
    nn  = interp1(edges, edges, double(x(:)), 'nearest', 'extrap');
    off = double(x(:)) - nn;
end