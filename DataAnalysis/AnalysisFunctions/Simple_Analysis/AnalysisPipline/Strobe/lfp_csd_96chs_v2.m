function R = lfp_csd_analysis()
% LFP depth profile + CSD for Intan RHS recordings (one-file-per-signal-type)
%
% Probes (picked automatically from the number of channels in amplifier.dat):
%   32 ch : single-shank flexible ('flex1') or rigid ('rigid1')  -> cfg.probe32
%   64 ch : 4-shank flexible ('flex4')
%   96 ch : hybrid = 4-shank flexible (ports A+B) + single-shank flexible (port C)
%
% Self-contained: probe maps and the info.rhs header reader are built in, so
% ProbeMAP.m, Depth_s.m and read_Intan_RHS2000_file.m are NOT needed.
% (loadTrig.m is only needed if cfg.trigSource = 'loadTrig'.)
%
% Per channel:
%   raw int16 -> µV -> remove DC -> interpolate stimulus artifacts on RAW data
%   -> resample to fs_lfp (anti-aliased) -> 50 Hz notches -> band-pass (zero-phase)
%   -> epoch -> baseline (artifact-free pre-stim samples only)
% Then:
%   noisy-trial rejection -> trial average -> interpolate bad channels
%   -> light spatial smoothing -> CSD (Vaknin boundary) -> plots
%
% Conventions
%   Depth axis: top contact at the top.
%   Red = negative, blue = positive in every plot.
%   CSD: -sigma * d2(phi)/dz2, units µA/mm^3.  Sink = negative = RED,
%        source = positive = BLUE.
%
% Requires: Signal Processing Toolbox, MATLAB R2020b or newer.

%% ============================ CONFIG ==================================
cfg.folder      = '/Volumes/MACData/Data/Data_Xia/DX037/strobe_261002_103144';
                               % must contain amplifier.dat, time.dat, info.rhs (+ digitalin.dat)
cfg.gain        = 0.195;       % µV per bit (Intan)

% Probe. 'auto' picks by channel count: 96 -> 'hybrid96', 64 -> 'flex4',
%        32 -> cfg.probe32.  Or set one of these explicitly:
%   'hybrid96' : 4-shank flex + single-shank flex (ProbeMAP column 7)
%   'flex4'    : 4-shank flexible, 64 ch (ProbeMAP column 6)
%   'flex1'    : single-shank flexible, 32 ch (ProbeMAP column 5)
%   'rigid1'   : single-shank rigid, 32 ch (ProbeMAP column 3)
cfg.probeType   = 'auto';
cfg.probe32     = 'flex1';     % which 32-ch probe 'auto' should assume
cfg.badElec     = [];          % electrode numbers (E#) to interpolate, e.g. [7 40 70]

% Geometry, separately for the 4-shank part and the single shank
% (the hybrid probe has both; other probes use only one of these).
%   pitch_um  : contact spacing along the shank
%   angle_deg : insertion angle
%   angleRef  : 'vertical' = angle from the cortical normal -> dz = pitch*cos(angle)
%               'surface'  = angle from the cortical surface -> dz = pitch*sin(angle)
%   offset_um : depth of the top contact (0 = depths relative to the top contact)
%   flipDepth : true if electrode 1 of each shank is at the TIP
cfg.geom.shank4 = struct('pitch_um', 50, 'angle_deg', 30, 'angleRef', 'vertical', ...
                         'offset_um', 0, 'flipDepth', false);
cfg.geom.single = struct('pitch_um', 50, 'angle_deg', 30, 'angleRef', 'vertical', ...
                         'offset_um', 0, 'flipDepth', false);

% Triggers
cfg.trigSource  = 'digitalin'; % 'digitalin' (reads digitalin.dat) or 'loadTrig'
cfg.trigBit     = [];          % digital input bit (0 = DIGITAL-IN-01); [] = auto-detect
cfg.trigEdge    = 'rising';    % 'rising' or 'falling'
cfg.trigMinGap_ms = 5;         % ignore edges closer than this (debounce)

% LFP
cfg.fs_lfp      = 2000;        % Hz after downsampling
cfg.band        = [1 300];     % Hz
cfg.filtOrder   = 3;           % Butterworth order (band-pass -> 2x poles)
cfg.notchHz     = 50:50:300;   % mains line noise + harmonics; [] = off
cfg.notchBW     = 2;           % Hz, total width of each notch

% Epoch (ms relative to trigger)
cfg.pre_ms      = 50;
cfg.post_ms     = 300;

% Stimulus artifacts (ms relative to trigger); assumed time-locked to the trigger
cfg.blankWin    = [0 10];      % strobe window; [] to disable
cfg.artPeriod   = 34;          % periodic artifact every N ms; [] to disable
cfg.artHalf     = 1.5;         % half-width of each periodic artifact
cfg.artPad_ms   = 100;         % also clean this far beyond the epoch (filter edges)

% Trial rejection: drop trials whose peak |LFP| > median + k * robust SD
cfg.rejectK     = 6;           % Inf = off

% CSD
cfg.smoothW     = [0.23 0.54 0.23];   % spatial (Hamming) smoothing of LFP; [] = none
cfg.sigma       = 0.3;         % extracellular conductivity, S/m
cfg.plotUpsample = 4;          % depth interpolation for CSD DISPLAY only (1 = off)

% Display
cfg.lfpGain     = 1.0;         % trace size: 1 = biggest deflection fills one channel gap
cfg.troughWin_ms = [40 150];   % window to find each channel's main negative peak
cfg.blankDisplay = 'marks';    % 'marks': draw through blanked windows (interpolated),
                               %          grey ticks at the top show where they are
                               % 'gaps' : leave blanked windows empty
cfg.scaleMode   = 'shared';    % 'shared': same µV / CSD scale on every shank (comparable)
                               % 'perShank': each shank scaled to its own range

cfg.saveFile    = 'lfp_csd_results.mat';
cfg.useSaved    = false;       % true = skip processing, just re-plot cfg.saveFile

%% ===================== RE-PLOT FROM SAVED RESULTS =====================
if cfg.useSaved && isfile(cfg.saveFile)
    R = load(cfg.saveFile);
    if ~iscell(R.depths)                        % file from the older 4-shank script
        R.depths = repmat({R.depths}, 1, numel(R.lfp));
        R.dz     = repmat(R.dz, 1, numel(R.lfp));
    end
    R.cfg.lfpGain      = cfg.lfpGain;           % display settings from this run
    R.cfg.troughWin_ms = cfg.troughWin_ms;
    R.cfg.plotUpsample = cfg.plotUpsample;
    R.cfg.blankDisplay = cfg.blankDisplay;
    R.cfg.scaleMode    = cfg.scaleMode;
    iv = artifact_intervals(R.cfg);
    plot_lfp(R, iv);
    plot_csd(R, iv);
    return;
end

%% ======================== FILES AND HEADER ============================
ampFile  = fullfile(cfg.folder, 'amplifier.dat');
ia = dir(ampFile);
it = dir(fullfile(cfg.folder, 'time.dat'));
assert(~isempty(ia) && ~isempty(it), 'amplifier.dat or time.dat not found in cfg.folder.');
nSamples      = it.bytes / 4;                        % time.dat: one int32 per sample
cfg.nChanFile = ia.bytes / (2 * nSamples);
assert(mod(cfg.nChanFile, 1) == 0, 'amplifier.dat size does not match time.dat.');

% Probe
if strcmpi(cfg.probeType, 'auto')
    switch cfg.nChanFile
        case 96, cfg.probeType = 'hybrid96';
        case 64, cfg.probeType = 'flex4';
        case 32, cfg.probeType = cfg.probe32;
        otherwise, error('No automatic probe for %d channels: set cfg.probeType.', cfg.nChanFile);
    end
end
P = probe_definition(cfg.probeType);
for s = 1:numel(P.shanks)
    if cfg.geom.(P.part{s}).flipDepth
        P.shanks{s} = fliplr(P.shanks{s});
    end
end
cfg.shanks = P.shanks;  cfg.shankNames = P.shankNames;  cfg.shankPart = P.part;

% Recorded channel names and settings from info.rhs
H = [];
try
    H = read_rhs_header(fullfile(cfg.folder, 'info.rhs'));
catch err
    warning('Could not read info.rhs (%s). Assuming full 32-ch ports in order.', err.message);
end
if ~isempty(H)
    assert(numel(H.ampNames) == cfg.nChanFile, ...
        'info.rhs lists %d amplifier channels but amplifier.dat holds %d.', numel(H.ampNames), cfg.nChanFile);
    cfg.fs   = H.fs;
    ampNames = H.ampNames;
    fprintf('info.rhs: fs %g Hz, hardware band %.2f-%.0f Hz, DSP %s (%.2f Hz)\n', ...
        H.fs, H.lower_bw, H.upper_bw, onoff(H.dsp_enabled), H.dsp_cutoff);
    lowest = H.lower_bw;
    if H.dsp_enabled, lowest = max(lowest, H.dsp_cutoff); end
    if lowest > cfg.band(1)
        warning('Recording was high-passed at %.2f Hz, above cfg.band(1) = %g Hz.', lowest, cfg.band(1));
    end
else
    cfg.fs = 30000;
    ports  = unique(cellfun(@(s) s(1), P.labels));
    assert(cfg.nChanFile == 32*numel(ports), 'Cannot guess channel names; fix info.rhs reading.');
    ampNames = {};
    for q = ports, ampNames = [ampNames, arrayfun(@(c) sprintf('%c-%03d', q, c), 0:31, 'UniformOutput', false)]; end %#ok<AGROW>
end
fprintf('amplifier.dat: %d channels, %.1f min, probe ''%s''\n', ...
    cfg.nChanFile, nSamples/cfg.fs/60, cfg.probeType);

elec2row = labels_to_rows(P.labels, ampNames);   % E# -> row in amplifier.dat

%% ============================ SETUP ===================================
nS   = numel(cfg.shanks);
nChS = cellfun(@numel, cfg.shanks);
depths = cell(1, nS);  dz = zeros(1, nS);
for s = 1:nS
    G = cfg.geom.(cfg.shankPart{s});
    switch lower(G.angleRef)
        case 'vertical', dz(s) = G.pitch_um * cosd(G.angle_deg);
        case 'surface',  dz(s) = G.pitch_um * sind(G.angle_deg);
        otherwise, error('angleRef must be ''vertical'' or ''surface''.');
    end
    depths{s} = G.offset_um + (0:nChS(s)-1) * dz(s);
    fprintf('%-14s %2d contacts, %.1f µm/channel, span %.0f µm\n', ...
        cfg.shankNames{s}, nChS(s), dz(s), depths{s}(end) - depths{s}(1));
end
shankRows = mat2cell(1:sum(nChS), 1, nChS);      % rows of each shank in ep/mu

trig    = get_triggers(cfg, nSamples);
nTrials = numel(trig);

mm = memmapfile(ampFile, 'Format', {'int16', [cfg.nChanFile nSamples], 'x'});

% Artifact intervals and epoch time axis
iv     = artifact_intervals(cfg);
pre_s  = round(cfg.pre_ms  * cfg.fs_lfp / 1000);
post_s = round(cfg.post_ms * cfg.fs_lfp / 1000);
rel    = -pre_s:post_s;                    % samples relative to trigger
t_ms   = rel / cfg.fs_lfp * 1000;          % t = 0 exactly at the trigger
blank  = false(size(t_ms));
for k = 1:size(iv,1)
    blank = blank | (t_ms >= iv(k,1) & t_ms <= iv(k,2));
end
base = t_ms < 0 & ~blank;                  % artifact-free baseline
assert(any(base), 'No artifact-free baseline samples: widen pre_ms or check artifact settings.');

% Triggers in the downsampled time base
nLfp   = ceil(nSamples * cfg.fs_lfp / cfg.fs);
trig_l = round((trig - 1) * cfg.fs_lfp / cfg.fs) + 1;
valid  = (trig_l + rel(1) >= 1) & (trig_l + rel(end) <= nLfp);
fprintf('%d triggers, %d fully inside the recording\n', nTrials, sum(valid));

% Filters (second-order sections for numerical stability)
[z, p, k] = butter(cfg.filtOrder, cfg.band / (cfg.fs_lfp/2), 'bandpass');
[sos, g]  = zp2sos(z, p, k);
notch = cell(0, 2);
for f0 = cfg.notchHz(cfg.notchHz + cfg.notchBW/2 < cfg.fs_lfp/2)
    [zn, pn, kn] = butter(2, [f0 - cfg.notchBW/2, f0 + cfg.notchBW/2] / (cfg.fs_lfp/2), 'stop');
    [sn, gn] = zp2sos(zn, pn, kn);
    notch(end+1, :) = {sn, gn}; %#ok<AGROW>
end

%% ======================= PER-CHANNEL PROCESSING =======================
allE = [cfg.shanks{:}];
nE   = numel(allE);
ep   = nan(nE, numel(rel), nTrials);       % electrodes x time x trials

for i = 1:nE
    row = elec2row(allE(i));
    fprintf('E%-3d (%s, row %2d) ... ', allE(i), ampNames{row}, row);

    x = double(mm.Data.x(row, :)) * cfg.gain;   % µV
    x = x - mean(x);
    x = interp_artifacts(x, trig, iv, cfg.fs);  % before any filtering
    x = resample(x, cfg.fs_lfp, cfg.fs);
    for q = 1:size(notch, 1)
        x = filtfilt(notch{q,1}, notch{q,2}, x);
    end
    x = filtfilt(sos, g, x);

    idx = trig_l(valid)' + rel;                 % [nValid x nTime]
    e   = x(idx);
    e   = e - mean(e(:, base), 2);              % per-trial baseline
    ep(i, :, valid) = permute(e, [3 2 1]);
    fprintf('done\n');
end

%% ============================ QC ======================================
isBad = ismember(allE, cfg.badElec);

% Noisy / dead channel hint: median single-trial baseline SD per electrode,
% compared within each shank (the two probe parts may differ in noise)
bsd   = squeeze(std(ep(:, base, valid), 0, 2));
noise = median(bsd, 2, 'omitnan');
sus   = [];
for s = 1:nS
    r  = shankRows{s};
    mn = median(noise(r(~isBad(r))));
    sus = [sus, allE(r(noise(r) > 3*mn | noise(r) < mn/3))]; %#ok<AGROW>
end
if ~isempty(sus)
    fprintf('Check these electrodes (baseline noise far from their shank''s median): %s\n', mat2str(sus));
end

% Trial rejection on peak |LFP| across good electrodes, artifact-free samples
pk   = squeeze(max(max(abs(ep(~isBad, ~blank, :)), [], 1), [], 2))';
keep = valid;
if isfinite(cfg.rejectK)
    pv   = pk(valid);
    thr  = median(pv) + cfg.rejectK * 1.4826 * mad(pv, 1);
    keep = valid & pk <= thr;
end
fprintf('Trials kept: %d / %d valid\n', sum(keep), sum(valid));

%% ======================= AVERAGE, FIX, CSD ============================
mu = mean(ep(:, :, keep), 3);   % blanked samples hold interpolated values; see R.blank

R = struct();
R.cfg = cfg;  R.t_ms = t_ms;  R.depths = depths;  R.dz = dz;
R.keep = keep;  R.nKept = sum(keep);  R.blank = blank;  R.trig = trig;
R.electrodeRow = elec2row;  R.ampNames = ampNames;
R.lfp = cell(nS,1);  R.csd = cell(nS,1);  R.bad = cell(nS,1);

for s = 1:nS
    rows = shankRows{s};
    bad  = isBad(rows);
    lfp  = fix_bad_channels(mu(rows, :), bad, depths{s});
    R.lfp{s} = lfp;                                               % µV
    R.csd{s} = compute_csd(lfp, dz(s), cfg.smoothW, cfg.sigma);   % µA/mm^3
    R.bad{s} = bad;
end

if ~isempty(cfg.saveFile)
    save(cfg.saveFile, '-struct', 'R');
    fprintf('Saved %s\n', cfg.saveFile);
end

%% ============================ PLOTS ===================================
plot_lfp(R, iv);
plot_csd(R, iv);
end


%% ======================================================================
%  Probe maps (from ProbeMAP.m)
%% ======================================================================

function P = probe_definition(type)
% P.labels{E} = Intan channel name for electrode E
% P.shanks    = electrode numbers per shank (top -> tip assumed)
% P.part      = which cfg.geom entry each shank uses ('shank4' or 'single')
lab = @(port, ch) arrayfun(@(c) sprintf('%c-%03d', port, c), ch(:)', 'UniformOutput', false);

% 4-shank flexible (ProbeMAP NN_flexible_4SHANKA / B)
A4 = [31,7,0,24,30,6,1,25,29,5,2,26,28,4,3,27,  8,16,23,15,9,17,22,14,10,18,21,13,11,19,20,12]; % E1-16 S1 | E17-32 S4
B4 = [8,11,9,15,10,19,12,23,13,22,14,21,16,20,17,18,  4,7,0,6,28,5,24,3,25,2,26,1,27,31,29,30]; % E33-48 S2 | E49-64 S3
% Single-shank flexible (ProbeMAP NN_flexible)
F1 = [8,7,9,6,10,5,12,3,13,2,14,1,23,24,22,25,21,26,19,28,18,29,17,30,16,31,20,27,15,0,11,4];

shanks4 = {1:16, 33:48, 49:64, 17:32};
names4  = {'Shank 1', 'Shank 2', 'Shank 3', 'Shank 4'};

switch lower(type)
    case 'hybrid96'     % ProbeMAP column 7: 4-shank on A+B, single shank on C (E65-96)
        P.labels     = [lab('A', A4), lab('B', B4), lab('C', F1)];
        P.shanks     = [shanks4, {65:96}];
        P.shankNames = [names4, {'Single shank'}];
        P.part       = {'shank4', 'shank4', 'shank4', 'shank4', 'single'};
    case 'flex4'        % ProbeMAP column 6, E1-64
        P.labels     = [lab('A', A4), lab('B', B4)];
        P.shanks     = shanks4;
        P.shankNames = names4;
        P.part       = {'shank4', 'shank4', 'shank4', 'shank4'};
    case 'flex1'        % ProbeMAP column 5
        P.labels     = lab('A', F1);
        P.shanks     = {1:32};
        P.shankNames = {'Single shank (flex)'};
        P.part       = {'single'};
    case 'rigid1'       % ProbeMAP column 3, via connectors
        MOLC_MALE     = [32,30,31,28,29,27,25,22,23,21,17,18,19,20,24,26,1,4,13,14,15,16,12,10,8,6,2,3,5,7,9,11];
        MOLC_FEMALE   = 32:-1:1;
        OMNETICS_MALE = [23,25,27,29,31,19,17,21,11,15,13,1,3,5,7,9,10,8,6,4,2,14,16,12,22,18,20,32,30,28,26,24];
        NN_I          = [1,32,2,31,3,30,4,29,5,28,6,27,7,26,8,25,9,24,10,23,11,22,12,21,13,20,14,19,15,18,16,17];
        T = zeros(1, 32);
        for n = 1:32
            T(n) = find(OMNETICS_MALE == find(MOLC_FEMALE == find(MOLC_MALE == NN_I(n)))) - 1;
        end
        P.labels     = lab('A', T);
        P.shanks     = {1:32};
        P.shankNames = {'Single shank (rigid)'};
        P.part       = {'single'};
    otherwise
        error('Unknown probe type ''%s''.', type);
end
assert(numel(unique(P.labels)) == numel(P.labels), 'Probe map has duplicate channels.');
end

function rows = labels_to_rows(labels, names)
% Row in amplifier.dat for each probe label, by matching info.rhs channel names.
% If the probe was plugged into different port(s) than the map assumes
% (e.g. port B instead of A), port letters are remapped in order.
[tf, rows] = ismember(labels, names);
if ~all(tf)
    pl = unique(cellfun(@(s) s(1), labels));
    pr = unique(cellfun(@(s) s(1), names));
    if ~any(tf) && numel(pl) == numel(pr)
        lab2 = labels;
        for q = 1:numel(pl)
            hit = cellfun(@(s) s(1) == pl(q), labels);
            lab2(hit) = cellfun(@(s) [pr(q) s(2:end)], labels(hit), 'UniformOutput', false);
        end
        warning('Probe map uses port(s) %s but recording has %s: remapping.', pl, pr);
        [tf, rows] = ismember(lab2, names);
    end
end
assert(all(tf), 'These probe channels are not in the recording: %s', strjoin(labels(~tf), ', '));
rows = rows(:)';
end

%% ======================================================================
%  info.rhs header reader (minimal; follows Intan's RHS file format)
%% ======================================================================

function H = read_rhs_header(file)
fid = fopen(file, 'r', 'ieee-le');
assert(fid > 0, 'Cannot open %s', file);
cln = onCleanup(@() fclose(fid));

magic = fread(fid, 1, 'uint32');
assert(magic == hex2dec('d69127ac'), 'Not an Intan RHS header file.');
fread(fid, 2, 'int16');                    % version major, minor
H.fs          = fread(fid, 1, 'single');
H.dsp_enabled = fread(fid, 1, 'int16');
f = fread(fid, 8, 'single');               % actual dsp/lower/lower-settle/upper, desired x4
H.dsp_cutoff  = f(1);
H.lower_bw    = f(2);
H.upper_bw    = f(4);
fread(fid, 1, 'int16');                    % notch filter mode
fread(fid, 2, 'single');                   % impedance test frequency (desired, actual)
fread(fid, 2, 'int16');                    % amp settle mode, charge recovery mode
fread(fid, 3, 'single');                   % stim step size, recovery current, target voltage
for q = 1:3, read_qstring(fid); end        % notes
fread(fid, 2, 'int16');                    % dc amplifier data saved, board mode
read_qstring(fid);                         % reference channel
nGroups = fread(fid, 1, 'int16');

names = {};
for gI = 1:nGroups
    read_qstring(fid);  read_qstring(fid); % group name, prefix
    en = fread(fid, 1, 'int16');
    nc = fread(fid, 1, 'int16');
    fread(fid, 1, 'int16');                % number of amplifier channels
    if nc > 0 && en > 0
        for c = 1:nc
            nat = read_qstring(fid);
            read_qstring(fid);             % custom name
            fread(fid, 2, 'int16');        % native order, custom order
            sigType = fread(fid, 1, 'int16');
            chEn    = fread(fid, 1, 'int16');
            fread(fid, 3, 'int16');        % chip channel, command stream, board stream
            fread(fid, 4, 'int16');        % spike-scope trigger settings
            fread(fid, 2, 'single');       % impedance magnitude, phase
            if chEn && sigType == 0        % enabled amplifier channel
                names{end+1} = nat; %#ok<AGROW>
            end
        end
    end
end
assert(H.fs > 1000 && H.fs < 100000 && ~isempty(names), 'info.rhs parsed to implausible values.');
H.ampNames = names;
end

function s = read_qstring(fid)
len = fread(fid, 1, 'uint32');
if isempty(len) || len == hex2dec('ffffffff') || len == 0
    s = '';  return;
end
s = char(fread(fid, len/2, 'uint16')');
end

function s = onoff(v)
if v, s = 'on'; else, s = 'off'; end
end

%% ======================================================================
%  Triggers
%% ======================================================================

function trig = get_triggers(cfg, nSamples)
switch lower(cfg.trigSource)
    case 'loadtrig'
        trig = loadTrig(0);
        trig = trig(1,:);
        fprintf('Triggers from loadTrig: %d\n', numel(trig));
    case 'digitalin'
        fid = fopen(fullfile(cfg.folder, 'digitalin.dat'), 'r');
        assert(fid > 0, 'digitalin.dat not found.');
        d = fread(fid, Inf, 'uint16=>uint16');
        fclose(fid);
        assert(numel(d) == nSamples, 'digitalin.dat length does not match time.dat.');
        bit = cfg.trigBit;
        if isempty(bit)
            act = [];
            for b = 0:15
                v = bitget(d, b+1);
                if any(v) && ~all(v), act(end+1) = b; end %#ok<AGROW>
            end
            assert(numel(act) == 1, ...
                'Active digital inputs (bits): %s. Set cfg.trigBit to the stimulus line.', mat2str(act));
            bit = act;
        end
        v = double(bitget(d, bit+1));
        if strcmpi(cfg.trigEdge, 'rising')
            on = find(diff(v) ==  1) + 1;
        else
            on = find(diff(v) == -1) + 1;
        end
        gap  = cfg.trigMinGap_ms * cfg.fs / 1000;
        keep = true(size(on));  last = -Inf;
        for i = 1:numel(on)
            if on(i) - last < gap, keep(i) = false; else, last = on(i); end
        end
        trig = on(keep)';
        fprintf('Triggers from DIGITAL-IN-%02d: %d, median interval %.1f ms\n', ...
            bit+1, numel(trig), median(diff(trig)) / cfg.fs * 1000);
    otherwise
        error('cfg.trigSource must be ''digitalin'' or ''loadTrig''.');
end
assert(~isempty(trig), 'No triggers found.');
end

%% ======================================================================
%  Signal processing helpers
%% ======================================================================

function iv = artifact_intervals(cfg)
% Merged [start end] ms intervals to clean, covering the epoch plus padding.
lo = -cfg.pre_ms - cfg.artPad_ms;
hi =  cfg.post_ms + cfg.artPad_ms;
iv = zeros(0, 2);
if ~isempty(cfg.blankWin)
    iv = [iv; cfg.blankWin(:)'];
end
if ~isempty(cfg.artPeriod)
    tk = [fliplr(-cfg.artPeriod:-cfg.artPeriod:lo), 0:cfg.artPeriod:hi];
    iv = [iv; tk(:) - cfg.artHalf, tk(:) + cfg.artHalf];
end
iv = merge_intervals(iv);
end

function m = merge_intervals(iv)
if isempty(iv), m = zeros(0,2); return; end
iv = sortrows(iv, 1);
m  = iv(1,:);
for k = 2:size(iv,1)
    if iv(k,1) <= m(end,2)
        m(end,2) = max(m(end,2), iv(k,2));
    else
        m(end+1,:) = iv(k,:); %#ok<AGROW>
    end
end
end

function x = interp_artifacts(x, trig, iv, fs)
% Linear interpolation across each artifact interval around every trigger,
% on the raw (unfiltered) signal so filters cannot smear the artifact.
n = numel(x);
for k = 1:size(iv,1)
    a = round(iv(k,1) * fs / 1000);
    b = round(iv(k,2) * fs / 1000);
    for t = trig(:)'
        i1 = t + a;  i2 = t + b;
        if i1 < 2 || i2 > n-1, continue; end
        x(i1:i2) = linspace(x(i1-1), x(i2+1), i2 - i1 + 1);
    end
end
end

function lfp = fix_bad_channels(lfp, bad, depths)
% Replace bad channels by interpolation across depth (nearest at the edges).
if ~any(bad), return; end
g = ~bad;
assert(sum(g) >= 2, 'Need at least two good channels per shank.');
v  = interp1(depths(g), lfp(g,:), depths(bad), 'linear');
vn = interp1(depths(g), lfp(g,:), depths(bad), 'nearest', 'extrap');
miss = isnan(v) & ~isnan(vn);
v(miss) = vn(miss);
lfp(bad,:) = v;
end

function csd = compute_csd(lfp, dz, w, sigma)
% lfp [nCh x nT] in µV, dz in µm -> CSD in µA/mm^3 (sink negative).
phi = lfp;
if ~isempty(w)
    w   = w(:) / sum(w);
    h   = (numel(w) - 1) / 2;
    pad = [repmat(phi(1,:), h, 1); phi; repmat(phi(end,:), h, 1)];
    phi = conv2(pad, w, 'valid');                  % same number of channels
end
phi = [phi(1,:); phi; phi(end,:)];                 % Vaknin boundary condition
csd = -sigma * diff(phi, 2, 1) / dz^2 * 1e3;       % µV/µm^2 * S/m -> µA/mm^3
end

%% ======================================================================
%  Plotting
%% ======================================================================

function shade_blanks(ax, iv, t_ms, yl, mode)
% 'gaps' : full-height light band;  'marks': short grey tick at the top edge
for k = 1:size(iv,1)
    a = max(iv(k,1), t_ms(1));  b = min(iv(k,2), t_ms(end));
    if a >= b, continue; end
    if strcmpi(mode, 'gaps')
        y = yl;  col = [0.95 0.95 0.95];
    else
        y = [yl(1), yl(1) + 0.015 * diff(yl)];  col = [0.55 0.55 0.55];
    end
    patch(ax, [a b b a], [y(1) y(1) y(2) y(2)], col, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end

function v = nice_value(x)
p = 10^floor(log10(x));  f = x / p;
if f >= 5, v = 5*p; elseif f >= 2, v = 2*p; else, v = p; end
end

function [showLbl, showBar] = panel_flags(R, s)
% Depth labels on the first panel of each probe part; scale bar / colour bar
% on the last panel of each part (or on every panel when scaled per shank).
nS = numel(R.lfp);
showLbl = s == 1 || ~isequal(R.depths{s}, R.depths{s-1});
showBar = s == nS || ~isequal(R.depths{s}, R.depths{s+1}) || strcmpi(R.cfg.scaleMode, 'perShank');
end

function pk = panel_peaks(R, field, pct)
% Robust peak |value| per shank (blanked samples excluded); shared = max over shanks.
pk = cellfun(@(X) prctile(reshape(abs(X(:, ~R.blank)), [], 1), pct), R.(field));
if ~strcmpi(R.cfg.scaleMode, 'perShank'), pk(:) = max(pk); end
end

function set_depth_axis(ax, showLbl, d)
set(ax, 'YDir', 'reverse', 'YTick', d, 'TickDir', 'out', 'Box', 'off');
lbl = arrayfun(@(x) sprintf('%.0f', x), d, 'UniformOutput', false);
if numel(d) > 16, lbl(2:2:end) = {''}; end         % label every other contact
if showLbl
    ylabel(ax, 'Depth (µm)');
    ax.YTickLabel = lbl;
else
    ax.YTickLabel = [];
end
end

function plot_lfp(R, iv)
% Row 1: stacked traces, negative deflections filled red, positive blue,
%        black dot = each channel's main negative peak (trough).
% Row 2: the same LFP as a colour image (red = negative, blue = positive).
nS    = numel(R.lfp);  t = R.t_ms;
gaps  = strcmpi(R.cfg.blankDisplay, 'gaps');
red   = [0.85 0.15 0.15];  blue = [0.15 0.35 0.85];
tw    = t >= R.cfg.troughWin_ms(1) & t <= R.cfg.troughWin_ms(2) & ~R.blank;
tt    = t(tw);
peaks = panel_peaks(R, 'lfp', 99.5);

figure('Color', 'w', 'Name', 'LFP depth profile', 'Position', [40 40 260 + 330*nS 1000]);
tl = tiledlayout(2, nS, 'TileSpacing', 'compact', 'Padding', 'compact');

for s = 1:nS
    d = R.depths{s};  dz = R.dz(s);  peak = peaks(s);
    scale = R.cfg.lfpGain * dz / peak;              % µm on the plot per µV
    yl    = [d(1) - dz, d(end) + dz];
    [showLbl, showBar] = panel_flags(R, s);
    L = R.lfp{s};
    if gaps, L(:, R.blank) = NaN; end

    % ---- traces ----
    ax = nexttile(tl, s);  hold(ax, 'on');
    shade_blanks(ax, iv, t, yl, R.cfg.blankDisplay);
    fprintf('\n%s: main negative peak per channel (%d-%d ms)\n', R.cfg.shankNames{s}, R.cfg.troughWin_ms);
    fprintf('  depth(µm)  latency(ms)  amplitude(µV)\n');
    for c = 1:numel(d)
        fill_trace(ax, t, L(c,:), d(c), scale, red, blue, 0.35);
        lc = [0.15 0.15 0.15];  if R.bad{s}(c), lc = [0.6 0 0.6]; end   % purple = interpolated
        plot(ax, t, d(c) - L(c,:) * scale, 'Color', lc, 'LineWidth', 0.8);
        [mn, ix] = min(L(c, tw));
        plot(ax, tt(ix), d(c) - mn * scale, 'k.', 'MarkerSize', 10);
        fprintf('  %8.0f  %11.1f  %13.1f\n', d(c), tt(ix), mn);
    end
    xline(ax, 0, 'k--', 'LineWidth', 1);
    set_depth_axis(ax, showLbl, d);
    xlim(ax, [t(1) t(end)]);  ylim(ax, yl);
    title(ax, R.cfg.shankNames{s}, 'FontSize', 12);
    if showBar                                      % amplitude scale bar
        barUV = nice_value(0.5 * peak);
        x0 = t(end) - 0.04 * (t(end) - t(1));
        y0 = d(end) - dz/2;
        plot(ax, [x0 x0], [y0, y0 - barUV * scale], 'k', 'LineWidth', 2.5);
        text(ax, x0, y0 - barUV * scale / 2, sprintf('%g µV ', barUV), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
    end

    % ---- colour image ----
    ax2 = nexttile(tl, nS + s);  hold(ax2, 'on');
    imagesc(ax2, t, d, L, 'AlphaData', ~isnan(L));
    if ~gaps, shade_blanks(ax2, iv, t, [d(1) - dz/2, d(end) + dz/2], 'marks'); end
    xline(ax2, 0, 'k--', 'LineWidth', 1);
    set_depth_axis(ax2, showLbl, d);
    set(ax2, 'CLim', [-peak peak], 'Color', [0.92 0.92 0.92], 'Layer', 'top', 'Box', 'on');
    colormap(ax2, sinksource(256));
    xlim(ax2, [t(1) t(end)]);  ylim(ax2, [d(1) - dz/2, d(end) + dz/2]);
    xlabel(ax2, 'Time (ms)');
    if showBar
        cb = colorbar(ax2);  cb.Label.String = 'LFP (µV)';
    end
end
title(tl, sprintf('LFP  (n = %d trials; red = negative, blue = positive; grey ticks = blanked)', ...
    R.nKept), 'FontSize', 13);
end

function fill_trace(ax, t, v, y0, scale, colNeg, colPos, alpha)
% Fill between a trace and its baseline, skipping NaN gaps.
ok = ~isnan(v);
e  = diff([0 ok 0]);
st = find(e == 1);  en = find(e == -1) - 1;
for k = 1:numel(st)
    ii = st(k):en(k);
    x  = t(ii);  vv = v(ii);  z = y0 * ones(1, numel(ii));
    patch(ax, [x fliplr(x)], [y0 - min(vv,0)*scale, z], colNeg, ...
        'EdgeColor', 'none', 'FaceAlpha', alpha, 'HandleVisibility', 'off');
    patch(ax, [x fliplr(x)], [y0 - max(vv,0)*scale, z], colPos, ...
        'EdgeColor', 'none', 'FaceAlpha', alpha, 'HandleVisibility', 'off');
end
end

function plot_csd(R, iv)
nS   = numel(R.csd);  t = R.t_ms;
gaps = strcmpi(R.cfg.blankDisplay, 'gaps');
cls  = panel_peaks(R, 'csd', 99);                   % colour limits
lpk  = panel_peaks(R, 'lfp', 100);                  % for overlaid traces
up   = max(1, round(R.cfg.plotUpsample));

figure('Color', 'w', 'Name', 'CSD', 'Position', [50 100 260 + 330*nS 650]);
tl = tiledlayout(1, nS, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:nS
    d  = R.depths{s};  dz = R.dz(s);  cl = cls(s);
    scale = 0.8 * dz / lpk(s);
    [showLbl, showBar] = panel_flags(R, s);
    dq = linspace(d(1), d(end), (numel(d)-1) * up + 1);
    C  = R.csd{s};  L = R.lfp{s};
    if gaps, C(:, R.blank) = NaN;  L(:, R.blank) = NaN; end
    if up > 1, C = interp1(d, C, dq, 'linear'); end % display only

    ax = nexttile(tl);  hold(ax, 'on');
    imagesc(ax, t, dq, C, 'AlphaData', ~isnan(C));
    for c = 1:numel(d)                              % LFP traces on top for reference
        plot(ax, t, d(c) - L(c,:) * scale, 'Color', [0 0 0 0.5], 'LineWidth', 0.6);
    end
    if ~gaps, shade_blanks(ax, iv, t, [d(1) - dz/2, d(end) + dz/2], 'marks'); end
    xline(ax, 0, 'k--', 'LineWidth', 1);
    set_depth_axis(ax, showLbl, d);
    set(ax, 'CLim', [-cl cl], 'Color', [0.85 0.85 0.85], 'Layer', 'top', 'Box', 'on');
    colormap(ax, sinksource(256));
    xlim(ax, [t(1) t(end)]);  ylim(ax, [d(1) - dz/2, d(end) + dz/2]);
    xlabel(ax, 'Time (ms)');  title(ax, R.cfg.shankNames{s}, 'FontSize', 12);
    if showBar
        cb = colorbar(ax);  cb.Label.String = 'CSD (µA/mm^3)';
    end
end
title(tl, sprintf('CSD  (n = %d trials, sigma = %.2f S/m; sink = red, source = blue)', ...
    R.nKept, R.cfg.sigma), 'FontSize', 13);
end

function cmap = sinksource(n)
% Negative (sink) -> red, zero -> white, positive (source) -> blue
h = floor(n/2);
r = [ones(h,1), linspace(0,1,h)', linspace(0,1,h)'];
b = [linspace(1,0,h)', linspace(1,0,h)', ones(h,1)];
cmap = [r; b];
end