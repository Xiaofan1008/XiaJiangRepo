%% ============================================================
% SpikeFiltering_Evoked_Cleanup.m
% Amplitude ceiling -> Zero-Crossing -> PCA + Mahalanobis-distance
% outlier rejection, applied to the waveform file and carried over to
% the spike-times file.
% ============================================================
% Loads BOTH the spike-times file (*.sp_xia.mat / sp_clipped) and the
% spike-waveforms file (*.sp_xia_waveforms.mat / sp_waveforms) for this
% dataset. Filtering decisions are made from the waveform shapes (since
% that's the only one of the two with the shape information), then the
% SAME per-spike keep/reject decision is applied to both, so the two
% outputs stay in sync with each other -- exactly like the original
% pair they're derived from.
%
% Neither original file is modified. Two NEW files are written:
%   <base_name>.sp_xia_filtered.mat           -- sp_clipped only
%   <base_name>.sp_xia_waveforms_filtered.mat -- sp_waveforms only
% Both also carry 'fs' and 'QC_params' (the settings used to produce
% them) for traceability. Downstream pipeline scripts are untouched for
% now -- point them at the _filtered files later, once you're happy
% with the result.
% ============================================================
clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
% Enter the two input files directly -- no auto-detection/guessing from
% the folder name, since that broke on folder-name patterns that don't
% match the assumed convention (e.g. base names ending up truncated to
% "Xia_Linearity_S4_100" instead of the real file prefix).
times_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1/Xia_Linearity_S4_100_400um_Single1.sp_xia.mat';
wf_file    = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1/Xia_Linearity_S4_100_400um_Single1.sp_xia_waveforms.mat';

% Must match the Electrode_Type used by the detection script for THIS
% dataset -- needed here only to translate probe_elecs (map positions)
% into physical-channel indices, since sp_clipped/sp_waveforms are both
% indexed by physical recording channel.
Electrode_Type = 3; % 0:single shank rigid; 1:single shank flex; 2:four shank flex; 3:64chn+32chn hybrid

% -------- PROBES (map positions) -- same convention as the detection
% script. Update to match whatever probe layout this dataset used.
probe_names = {'4-shank (Ports A+B)', 'Single shank (Port C)'};
probe_elecs = {1:64, 65:96};

do_amp_filter  = 1;
do_zc_filter   = 1;
do_pca_filter  = 1;

% 1. AMPLITUDE CEILING CHECK
% Flat, fixed cutoff -- not adaptive per channel. Anything whose
% waveform reaches above this many uV (absolute value, anywhere in the
% snippet) is treated as artifact, not a real action potential, full
% stop. Same idea as abs_max_probe already used at detection time, just
% re-checked here.
amp_ceiling_uV = 500;

% 2. ZERO-CROSSING CHECK -- no threshold to set; a spike must show some
% positive-going phase, not just a trough, to count as a plausible
% extracellular action potential.

% 3. PCA + MAHALANOBIS DISTANCE (per probe)
% pca_nComp:  number of principal components of the waveform SHAPE to
%             keep.
% pca_alpha:  significance level for the chi-squared cutoff, e.g. 0.001
%             means "flag anything outside the 99.9% confidence
%             boundary of this channel's own spike population in PCA
%             space."
pca_nComp_probe = [3, 3];
pca_alpha_probe = [0.001, 0.001];
pca_min_spikes  = 30;   % need at least this many spikes on a channel to
                        % fit PCA/covariance meaningfully; channels with
                        % fewer are left unfiltered at this stage (with
                        % a warning) rather than guessed at.

%% ================= LOAD BOTH INPUT FILES =================
assert(isfile(times_file), 'Cannot find %s.', times_file);
S_t = load(times_file);
assert(isfield(S_t, 'sp_clipped'), 'Variable "sp_clipped" not found in %s.', times_file);
sp_clipped_in = S_t.sp_clipped;

assert(isfile(wf_file), 'Cannot find %s.', wf_file);
S_w = load(wf_file);
assert(isfield(S_w, 'sp_waveforms'), 'Variable "sp_waveforms" not found in %s.', wf_file);
sp_waveforms_in = S_w.sp_waveforms;

if isfield(S_w, 'fs')
    fs = S_w.fs;
elseif isfield(S_t, 'fs')
    fs = S_t.fs;
else
    fs = 30000;
    warning('fs not found in either input file -- defaulting to 30000 Hz.');
end

nCh = numel(sp_waveforms_in);
assert(numel(sp_clipped_in) == nCh, ...
    'sp_clipped (%d channels) and sp_waveforms (%d channels) disagree -- check these came from the same run.', ...
    numel(sp_clipped_in), nCh);

% Sanity check: the time column inside sp_waveforms should exactly match
% sp_clipped, since both came from the same detection run on the same
% spikes, in the same order.
for ch = 1:nCh
    if isempty(sp_waveforms_in{ch}) && isempty(sp_clipped_in{ch}), continue; end
    if ~isequal(sp_waveforms_in{ch}(:,1), sp_clipped_in{ch}(:))
        error(['Ch %d: sp_waveforms'' time column does not match sp_clipped -- ' ...
            '%s and %s are out of sync (not from the same detection run?).'], ch, times_file, wf_file);
    end
end

% -------- Map probe_elecs (map positions) onto physical channels --------
% Uses ChnMap instead of Depth_s: identical electrode-to-channel mapping
% (same ProbeMAP/E_MAP lookup), but takes the channel count directly
% (nCh, already known from the loaded spike files) instead of reading it
% from the raw Intan recording files -- so this works even when
% amplifier.dat etc. have since been deleted. Remember to use ChnMap
% (not Depth_s) anywhere else in the pipeline that gets rewritten to
% work from the saved spike files alone.
map_nums_plus = ChnMap(Electrode_Type, nCh);   % physical channel per map position
nMapPos = numel(map_nums_plus);
if nMapPos ~= nCh
    error(['ChnMap(%d, %d) returned %d map positions but the spike files have %d channels -- ' ...
        'Electrode_Type does not match this dataset.'], Electrode_Type, nCh, nMapPos, nCh);
end

all_probe_elecs = sort(horzcat(probe_elecs{:}));
if numel(all_probe_elecs) ~= nMapPos || ~isequal(all_probe_elecs, 1:nMapPos)
    error(['probe_elecs (spanning %d positions) does not match Electrode_Type %d ' ...
        '(%d map positions) -- update probe_names/probe_elecs to match.'], ...
        numel(all_probe_elecs), Electrode_Type, nMapPos);
end

probe_of_mappos = zeros(1, nMapPos);
for p = 1:numel(probe_elecs)
    probe_of_mappos(probe_elecs{p}) = p;
end
probe_of = zeros(1, nCh);                 % indexed by PHYSICAL channel
probe_of(map_nums_plus) = probe_of_mappos;
if any(probe_of == 0)
    error('Every physical channel must map to a probe (check probe_elecs/Electrode_Type).');
end

%% ================= FILTERING (one keep-mask per channel) =================
sp_waveforms = cell(1, nCh);
sp_clipped   = cell(1, nCh);

% reject_log{ch} records, for EVERY spike in the raw input (same row
% order as sp_waveforms_in{ch}), which stage removed it ('kept' if it
% survived all stages). This is what lets the QC companion script
% (QC_FilterCheck_Waveforms.m) color rejected waveforms by the reason
% they were removed, instead of lumping them all together.
reject_log = cell(1, nCh);

fprintf('\nFiltering %d channels (amplitude ceiling -> zero-crossing -> PCA/Mahalanobis)...\n', nCh);
for ch = 1:nCh
    wf = sp_waveforms_in{ch};
    if isempty(wf)
        sp_waveforms{ch} = wf;
        sp_clipped{ch}   = sp_clipped_in{ch};
        reject_log{ch}   = struct('time', [], 'reason', {{}});
        continue;
    end

    spt    = wf(:,1);
    wfs    = wf(:,2:end);
    keep   = true(size(spt));
    reason = repmat({'kept'}, size(spt));

    % ----- 1) Amplitude ceiling -----
    if do_amp_filter
        valid = all(abs(wfs) <= amp_ceiling_uV, 2);
        removed_idx = keep & ~valid;
        removed = sum(removed_idx);
        reason(removed_idx) = {'amplitude'};
        keep = keep & valid;
        if removed > 0
            fprintf('Ch %3d: removed %d spikes (exceeded %g uV).\n', ch, removed, amp_ceiling_uV);
        end
    end

    % ----- 2) Zero-crossing -----
    if do_zc_filter
        valid = max(wfs, [], 2) > 0;
        removed_idx = keep & ~valid;
        removed = sum(removed_idx);
        reason(removed_idx) = {'zero_crossing'};
        keep = keep & valid;
        if removed > 0
            fprintf('Ch %3d: removed %d spikes (no positive phase).\n', ch, removed);
        end
    end

    % ----- 3) PCA + Mahalanobis distance (on whatever survived so far) -----
    if do_pca_filter
        idx = find(keep);
        n   = numel(idx);
        if n < pca_min_spikes
            fprintf('Ch %3d: only %d spikes left (< %d) -- PCA/Mahalanobis skipped for this channel.\n', ...
                ch, n, pca_min_spikes);
        else
            p     = probe_of(ch);
            nComp = min([pca_nComp_probe(p), n-1, size(wfs,2)-1]);
            alpha = pca_alpha_probe(p);

            [~, score] = pca(wfs(idx,:), 'NumComponents', nComp);
            d2     = mahal(score, score);
            thresh = chi2inv(1 - alpha, nComp);
            valid_sub = d2 <= thresh;

            removed = sum(~valid_sub);
            reason(idx(~valid_sub)) = {'pca'};
            keep(idx(~valid_sub)) = false;
            if removed > 0
                fprintf('Ch %3d: removed %d spikes (Mahalanobis d2 > %.1f, %d PCs, alpha=%.4f).\n', ...
                    ch, removed, thresh, nComp, alpha);
            end
        end
    end

    sp_waveforms{ch} = wf(keep, :);
    sp_clipped{ch}   = spt(keep);
    reject_log{ch}   = struct('time', spt, 'reason', {reason});
end

%% ================= SAVE OUTPUT (two separate files) =================
QC_params = struct();
QC_params.probe_names      = probe_names;
QC_params.probe_elecs      = probe_elecs;
QC_params.Electrode_Type   = Electrode_Type;
QC_params.amp_ceiling_uV   = amp_ceiling_uV;
QC_params.pca_nComp_probe  = pca_nComp_probe;
QC_params.pca_alpha_probe  = pca_alpha_probe;
QC_params.pca_min_spikes   = pca_min_spikes;
QC_params.do_amp_filter    = do_amp_filter;
QC_params.do_zc_filter     = do_zc_filter;
QC_params.do_pca_filter    = do_pca_filter;

% Output names derived directly from the two input paths given above --
% no folder-name guessing. If a file somehow doesn't end in the
% expected suffix, fall back to inserting "_filtered" before ".mat".
times_out = strrep(times_file, '.sp_xia.mat', '.sp_xia_filtered.mat');
if strcmp(times_out, times_file)
    [p,n,e] = fileparts(times_file); times_out = fullfile(p, [n '_filtered' e]);
end
save(times_out, 'sp_clipped', 'fs', 'QC_params', '-v7.3');
fprintf('\nSaved filtered spike times to %s\n', times_out);

wf_out = strrep(wf_file, '.sp_xia_waveforms.mat', '.sp_xia_waveforms_filtered.mat');
if strcmp(wf_out, wf_file)
    [p,n,e] = fileparts(wf_file); wf_out = fullfile(p, [n '_filtered' e]);
end
save(wf_out, 'sp_waveforms', 'fs', 'QC_params', '-v7.3');
fprintf('Saved filtered spike waveforms to %s\n', wf_out);

% Third file: per-spike reject reason (amplitude/zero_crossing/pca/kept),
% aligned row-for-row with sp_waveforms_in{ch} (the RAW waveforms file).
% Used by QC_FilterCheck_Waveforms.m to color-split rejected waveforms by
% which stage removed them; not needed by any other downstream script.
rejectlog_out = strrep(wf_file, '.sp_xia_waveforms.mat', '.sp_xia_rejectlog.mat');
if strcmp(rejectlog_out, wf_file)
    [p,n,e] = fileparts(wf_file); rejectlog_out = fullfile(p, [n '_rejectlog' e]);
end
save(rejectlog_out, 'reject_log', 'fs', 'QC_params', '-v7.3');
fprintf('Saved per-spike reject log to %s\n', rejectlog_out);