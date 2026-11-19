%% ============================================================
% NN_Train_FR_CountVsPulses_save_v2.m
% NN probe | pulse TRAIN data | FIRING-RATE analysis A | SAVE the per-trial / per-channel results for cross-dataset averaging | v2
% (v2 = ALL channels stored + channel-set masks, per-channel NORMALISATION result, richer metadata; v1 is not changed. NN data only: NP has its own script)
%
% PURPOSE
%   Runs the count / firing-rate analysis of NN_Train_FR_CountVsPulses_v1 on ONE recovery file and SAVES everything that a
%   later cross-dataset analysis needs. The RAW numbers are the source of truth; the per-channel normalisation of
%   NN_Train_FR_CountVsPulses_Norm_v3 is stored as well (Norm and Pat.rateN), so the averaging script can use it directly.
%   All levels of every chain are saved (half-steps of seq included, flagged), so the later script can choose pair levels only or all.
%
% PER LEVEL (each recovered level = its own condition; chain positions 1 ... L of a pattern)
%   window  = first pulse (+trainStart_ms) ... last pulse of the level + tail_ms
%   CNT     = spikes in the window, BAS = spikes in baseWin_ms          (per good trial and channel)
%   CORR    = CNT - BAS * (window length / baseline length)             (baseline-subtracted count per trial and channel)
%   cm      = mean of CORR over the good trials                         (per channel)
%   rate    = cm / window (spikes/s)          spp = cm / pulses per electrode
%   add     = cm(step) - cm(previous step of the chain); addRate = add / (window(step) - window(previous step))
%   addPair / addRatePair = the same but between consecutive PAIR levels only (sim n; seq [n n]) -> what the plotting script shows
%   ALL channels are computed and stored; the channel set of a level (chSet / chMask) says which ones are responsive.
%   Channel sets come from <folder>.FR_RespChn.mat (else chRange),
%   bad trials from <folder>.FR_BadTrials.mat (else all trials kept): exactly as in NN_Train_FR_CountVsPulses_v1.
%
% OUTPUT FILE  <out_folder>/<out_name>[_<rec_tag>].mat   (out_folder = '' -> the data folder; out_name = '' -> <dataset>.FR_CountVsPulses)   with the struct  FRRes:
%   FRRes.meta    dataset, animal, area, paradigm, note, paths of the files used, date, version, FS, settings used
%   FRRes.nChn, FRRes.patName(1..nPat)
%   FRRes.Pat(p)  .name .kind .finalCond .chain .cond(L) .isPair(L) .pairIdx .pe(L x nElectrodes) .x_n(L) .t_last_ms(L) .win_ms(L)
%                 .nTrialsUsed(L) .nBad(L) .trialRel{L} .trialAbs{L} .badRel{L}   (GOOD trials only are stored; the bad trials as indices)
%                 .chSet{L}  channels of the level's set (from the responsive-channel file, else chRange);  .chMask (L x nChn) the same as a mask
%                 .CNT{L} .BAS{L} .CORR{L}   (nGoodTrials x nChn each; ALL channels are stored, use chMask / Norm.chResp to select)
%                 .cm .rate .spp .add .addRate (L x nChn; ALL channels) .addPair .addRatePair (nPair x nChn)
%                 .rateN (L x nChn) = rate / own reference (NaN for channels not in Norm.keep)
%                 .rateN_mean .rateN_sem .rateN_n (L x 1)   group value of rateN over Norm.keep channels (SEM over channels)
%                 .stim   per electrode: electrode label (StimParams column 1), amp (column 16), period_us, freq_Hz, start_ms, nPulses,
%                         polarity (column 12), phase1_us, phase2_us, interphase_us (columns 13, 14, 15) of the final condition
%   FRRes.Norm    method, refLevels, minPeak_Hz, chResp (responsive channels, union over all levels), keep (mask of channels normalised),
%                 refC (1 x nChn, spikes/s; NaN = none), refPat / refCond / refX (where each channel peaks)
%   FRRes.recInfo   the settings stored in Recovery_params of the recovery file (thinning: rand_seed, percent_k, event_ms_k, sim/seq overrides ...)
%   FRRes.recLog    the recovery log of the file (per condition: nWindow, nKept, percent ... for a thinned file)
%   NORMALISATION: reference of a channel = the highest point of its own raw, baseline-corrected rate over all levels (refLevels: 'pair' =
%   sim levels and seq [n n] levels, 'all' = every level) of all patterns; normalised rate = rate / reference. Only the channels of the
%   responsive sets (union over the levels) are normalised. The rate uses the GOOD trials only.
%   Use:  S = load(file, 'FRRes');  S.FRRes.Pat(1).rateN(S.FRRes.Pat(1).pairIdx, :)  etc.
%   An existing file is NOT overwritten unless overwrite = true (the old file is then copied to a dated backup first).
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_2_100Hz_30uA';
rec_tag     = '';          % '' = *.sp_xia_TrainRecovery.mat (unthinned); 'RandomEvent' etc. = the tagged file
rec_file    = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_2_100Hz_30uA/Xia_Train_15Pulses_2_100Hz_30uA.sp_xia_TrainRecovery_RandomEvent.mat';          % full path (overrides rec_tag)
FS          = 30000;       % used only if the recovery file has no FS

patternConds = [];         % FINAL condition number(s) of the patterns; [] = all sim / AB / BA patterns found
chRange      = 1:32;       % manual channel list: used when no responsive-channel file / pattern entry exists
useRespChn   = true;       % true = use <folder>.FR_RespChn.mat when it exists
chSetMode    = 'union';    % 'pattern' | 'intersect' | 'union'  (only decides chSet / chMask; the values of all channels are stored; Norm always uses the union)
useBadTrials = true;       % true = use <folder>.FR_BadTrials.mat when it exists
resp_file    = '';         % '' = <data_folder>/<folder name>.FR_RespChn.mat
bad_file     = '';         % '' = <data_folder>/<folder name>.FR_BadTrials.mat

% ---- counting window and baseline ----
tail_ms       = 15;        % window = first pulse + trainStart_ms ... last pulse + tail_ms
trainStart_ms = 0;
baseWin_ms    = [-50 -5];  % baseline window relative to the trigger (ms)

% ---- normalisation (per channel: own peak of the raw rate; same as NN_Train_FR_CountVsPulses_Norm_v3) ----
minPeak_Hz = 0;            % channels whose own peak rate is <= this are left out of the normalisation
refLevels  = 'pair';       % 'pair' = reference over sim levels and seq [n n] levels | 'all' = over every level

% ---- free-text information stored with the result (for grouping datasets later) ----
meta_animal = 'CJ287';          % e.g. 'CJ287'
meta_area   = 'V1';          % e.g. 'V1' | 'V2'
meta_paradigm = 'Train_15Pulses';   % e.g. 'Train_Freq' | 'Train_Amp' | 'Train_ISI' | 'Train_15Pulses'
meta_note   = '';          % anything else

% ---- output ----
out_folder = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Marmo_Train/FR_RandomEvent_Analysis';           % folder for the result file (created if missing); '' = the data folder. The file name is made automatically.
out_name  = 'CJ287_15Pulses_2';            % file name WITHOUT .mat, e.g. 'CJ287_Train_Freq1'; '' = <dataset>.FR_CountVsPulses. The tag _<rec_tag> is added when rec_tag is set.
out_file  = '';            % full path WITH file name; overrides out_folder and out_name (nothing is added).
overwrite = false;         % true = replace an existing file (a dated backup of the old one is made first)
do_save   = true;          % false = only analyse and print
Plot_Check = true;         % quick look at the result (rate and added rate of the PAIR levels)

%% ================= CHECK SETTINGS =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
assert(ismember(chSetMode, {'pattern','intersect','union'}), 'chSetMode: pattern, intersect or union.');
assert(numel(baseWin_ms) == 2 && baseWin_ms(2) > baseWin_ms(1), 'baseWin_ms must be [start end].');
assert(tail_ms >= 0, 'tail_ms must be >= 0.');
assert(ismember(refLevels, {'pair','all'}), 'refLevels: pair or all.');
assert(minPeak_Hz >= 0, 'minPeak_Hz must be >= 0.');

%% ================= LOAD =================
if isempty(rec_file)
    if isempty(rec_tag), pat_ = '*.sp_xia_TrainRecovery.mat'; else, pat_ = ['*.sp_xia_TrainRecovery_' rec_tag '.mat']; end
    d = dir(fullfile(data_folder, pat_));  d = d(~startsWith({d.name}, '.'));
    assert(~isempty(d), 'No %s in %s.', pat_, data_folder);
    assert(numel(d) == 1, 'More than one %s in %s - set rec_file.', pat_, data_folder);
    rec_file = fullfile(data_folder, d(1).name);
end
assert(isfile(rec_file), 'Cannot find %s', rec_file);
tok_ = regexp(rec_file, '\.sp_xia_TrainRecovery_([^\\/]+)\.mat$', 'tokens', 'once');
if isempty(rec_tag) && ~isempty(tok_), rec_tag = tok_{1}; end     % the tag of a tagged recovery file goes into the output name
vn_ = who('-file', rec_file);
R = load(rec_file, vn_{ismember(vn_, {'sp_clipped','Recovery_params','FS','isRecovered','recLog'})});
assert(isfield(R, 'sp_clipped') && isfield(R, 'Recovery_params'), '%s is not a TrainRecovery file.', rec_file);
sp = R.sp_clipped;
if isfield(R, 'FS') && ~isempty(R.FS), FS = double(R.FS); end
trialCond = R.Recovery_params.trialCond(:);
condTable = R.Recovery_params.condTable;
nChn = numel(sp);  nTrialsAll = numel(trialCond);  nCond = numel(condTable);
if isfield(R, 'isRecovered') && numel(R.isRecovered) == nTrialsAll
    isRec = logical(R.isRecovered(:));
else
    isRec = true(nTrialsAll, 1);
end
recInfo = struct();  recLogS = [];
RP_ = R.Recovery_params;  fn_ = fieldnames(RP_);
for i_ = 1:numel(fn_)
    if ~ismember(fn_{i_}, {'trialCond','condTable'}), recInfo.(fn_{i_}) = RP_.(fn_{i_}); end
end
if isfield(R, 'recLog'), recLogS = R.recLog; end
clear R RP_ fn_

ft = dir(fullfile(data_folder, '*.trig.dat'));  ft = ft(~startsWith({ft.name}, '.'));
assert(~isempty(ft), 'No *.trig.dat in %s.', data_folder);
fid = fopen(fullfile(data_folder, ft(1).name), 'r');  trig = fread(fid, [1 inf], 'double');  fclose(fid);
t0 = nan(nTrialsAll, 1);  nUse = min(numel(trig), nTrialsAll);  t0(1:nUse) = trig(1:nUse)/FS*1000;

fe = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));  fe = fe(~startsWith({fe.name}, '.'));
assert(~isempty(fe), 'No *_exp_datafile_*.mat in %s.', data_folder);
E = load(fullfile(data_folder, fe(1).name), 'StimParams', 'simultaneous_stim');
nEl = double(E.simultaneous_stim);
peK = zeros(nCond, nEl);  stimK = cell(nCond, 1);
for k = 1:nCond
    t1 = find(trialCond == k, 1);
    st = struct('amp', zeros(1, nEl), 'period_us', nan(1, nEl), 'freq_Hz', nan(1, nEl), 'start_ms', nan(1, nEl), 'nPulses', zeros(1, nEl), ...
        'electrode', {cell(1, nEl)}, 'polarity', nan(1, nEl), 'phase1_us', nan(1, nEl), 'phase2_us', nan(1, nEl), 'interphase_us', nan(1, nEl));
    if ~isempty(t1)
        for e = 1:nEl
            row = 1 + (t1-1)*nEl + e;
            amp = double(E.StimParams{row,16});
            st.amp(e) = amp;
            st.electrode{e} = char(E.StimParams{row,1});  st.polarity(e) = double(E.StimParams{row,12});
            st.phase1_us(e) = double(E.StimParams{row,13});  st.phase2_us(e) = double(E.StimParams{row,14});  st.interphase_us(e) = double(E.StimParams{row,15});
            if amp > 0
                st.nPulses(e) = double(E.StimParams{row,8});  st.period_us(e) = double(E.StimParams{row,9});
                st.start_ms(e) = double(E.StimParams{row,6})/1000;
                if st.period_us(e) > 0, st.freq_Hz(e) = 1e6/st.period_us(e); end
                peK(k,e) = st.nPulses(e);
            end
        end
    end
    stimK{k} = st;
end
clear E

[~, folderName] = fileparts(regexprep(data_folder, '[\\/]+$', ''));
haveResp = false;  RespChn = [];
if useRespChn
    if isempty(resp_file), resp_file = fullfile(data_folder, [folderName '.FR_RespChn.mat']); end
    if isfile(resp_file)
        S_ = load(resp_file, 'RespChn');  RespChn = S_.RespChn;  haveResp = true;  clear S_
        assert(numel(RespChn.respChn) == nCond, 'The responsive-channel file has %d conditions, the recovery file %d.', numel(RespChn.respChn), nCond);
    end
end
haveBad = false;  BadTrials = [];
if useBadTrials
    if isempty(bad_file), bad_file = fullfile(data_folder, [folderName '.FR_BadTrials.mat']); end
    if isfile(bad_file)
        S_ = load(bad_file, 'BadTrials');  BadTrials = S_.BadTrials;  haveBad = true;  clear S_
        assert(isfield(BadTrials, 'badRel'), 'The bad-trial file must hold BadTrials.badRel{cond}.');
    end
end

%% ================= PATTERNS AND CHAINS =================
kinds = arrayfun(@(c) char(c.kind), condTable, 'UniformOutput', false);
isFinal = arrayfun(@(c) isfield(c,'isFinal') && ~isempty(c.isFinal) && logical(c.isFinal), condTable);
if isempty(patternConds)
    patternConds = find(isFinal(:)' & ismember(kinds(:)', {'sim','AB','BA'}));
end
assert(~isempty(patternConds), 'No sim / AB / BA pattern found in the recovery file.');
nPat = numel(patternConds);
chain = cell(nPat,1);  patName = cell(nPat,1);
for p = 1:nPat
    c = patternConds(p);  ch_ = c;
    while isfinite(condTable(ch_(1)).sourceCond), ch_ = [condTable(ch_(1)).sourceCond, ch_]; end %#ok<AGROW>
    chain{p} = ch_;
    switch kinds{c}
        case 'sim', patName{p} = 'sim';
        case 'AB',  patName{p} = 'seq A->B';
        case 'BA',  patName{p} = 'seq B->A';
        otherwise,  patName{p} = kinds{c};
    end
end

%% ================= CHANNEL SETS (every chain position) =================
chSet = cell(nPat, 1);
usedFile = false(nPat, 1);
for p = 1:nPat
    chSet{p} = cell(1, numel(chain{p}));
    for n = 1:numel(chain{p})
        k = chain{p}(n);
        if haveResp && RespChn.isSet(k)
            chSet{p}{n} = intersect(double(RespChn.respChn{k}(:)'), 1:nChn);  usedFile(p) = true;
        else
            chSet{p}{n} = intersect(chRange(:)', 1:nChn);
        end
    end
end
if ismember(chSetMode, {'intersect','union'})
    if strcmp(chSetMode, 'intersect'), common = 1:nChn; else, common = zeros(1, 0); end
    for p = 1:nPat
        for n = 1:numel(chSet{p})
            if strcmp(chSetMode, 'intersect'), common = intersect(common, chSet{p}{n}); else, common = union(common, chSet{p}{n}); end
        end
    end
    common = common(:)';
    for p = 1:nPat, for n = 1:numel(chSet{p}), chSet{p}{n} = common; end, end
    if haveResp && ~all(usedFile)
        warning('SV:notset', 'Pattern(s) %s have no entry in the responsive-channel file: chRange was used for them, which widens the %s set.', ...
            strjoin(patName(~usedFile), ', '), chSetMode);
    end
end

%% ================= ANALYSIS (every level of every chain) =================
baseLen = diff(baseWin_ms);
Pat = struct([]);
for p = 1:nPat
    cc = chain{p};  L = numel(cc);
    t = arrayfun(@(k) double(condTable(k).lastPulse_ms), cc);
    P = struct();
    P.name = patName{p};  P.kind = kinds{patternConds(p)};  P.finalCond = patternConds(p);  P.chain = cc(:)';
    P.cond = cc(:)';
    P.isPair = arrayfun(@(k) all(peK(k,:) == peK(k,1)) && peK(k,1) > 0, cc);
    P.pairIdx = find(P.isPair);
    P.pe = peK(cc, :);
    P.x_n = mean(P.pe, 2)';
    P.t_last_ms = t(:)';
    P.win_ms = nan(1, L);  P.nTrialsUsed = zeros(1, L);  P.nBad = zeros(1, L);
    P.trialRel = cell(1, L);  P.trialAbs = cell(1, L);  P.badRel = cell(1, L);  P.chSet = chSet{p};
    P.chMask = false(L, nChn);
    for n_ = 1:L, P.chMask(n_, chSet{p}{n_}) = true; end
    P.CNT = cell(1, L);  P.BAS = cell(1, L);  P.CORR = cell(1, L);
    P.cm = nan(L, nChn);
    for n = 1:L
        cond = cc(n);
        all_ = find(trialCond == cond);
        isL1 = ~isfinite(condTable(cond).sourceCond);
        ok = isfinite(t0(all_)) & (isRec(all_) | isL1);
        badR = zeros(1, 0);
        if haveBad && numel(BadTrials.badRel) >= cond && ~isempty(BadTrials.badRel{cond})
            br = BadTrials.badRel{cond}(:)';  br = br(br >= 1 & br <= numel(all_));
            badR = br;  P.nBad(n) = sum(ok(br));  ok(br) = false;
        end
        tr = all_(ok);  tt = t0(tr);
        a = t(1) + trainStart_ms;  b = t(n) + tail_ms;  scaleF = (b - a)/baseLen;
        P.win_ms(n) = b - a;  P.nTrialsUsed(n) = numel(tr);
        P.trialAbs{n} = tr(:)';  P.trialRel{n} = find(ok(:))';  P.badRel{n} = badR;
        CNT = nan(numel(tr), nChn);  BAS = nan(numel(tr), nChn);
        for ch = 1:nChn                                % ALL channels; chMask / Norm say which ones are responsive
            s = sp{ch}(:,1);
            for i = 1:numel(tt)
                CNT(i,ch) = lb_(s, tt(i)+b) - lb_(s, tt(i)+a);
                BAS(i,ch) = lb_(s, tt(i)+baseWin_ms(2)) - lb_(s, tt(i)+baseWin_ms(1));
            end
        end
        CORR = CNT - BAS*scaleF;
        P.CNT{n} = CNT;  P.BAS{n} = BAS;  P.CORR{n} = CORR;
        if numel(tr) >= 1, P.cm(n,:) = mean(CORR, 1); end
    end
    P.rate = P.cm ./ (P.win_ms(:)/1000);
    P.spp  = P.cm ./ P.x_n(:);
    dW = [P.win_ms(1), diff(P.win_ms)];
    P.add     = [P.cm(1,:); P.cm(2:end,:) - P.cm(1:end-1,:)];
    P.addRate = P.add ./ (dW(:)/1000);
    pi_ = P.pairIdx;
    cmPair = P.cm(pi_, :);  winPair = P.win_ms(pi_);
    P.addPair     = [cmPair(1,:); cmPair(2:end,:) - cmPair(1:end-1,:)];
    P.addRatePair = P.addPair ./ ([winPair(1), diff(winPair)]'/1000);
    fin = patternConds(p);
    P.stim = stimK{fin};
    Pat = [Pat, P]; %#ok<AGROW>
end

%% ================= NORMALISATION (per channel, own peak; as in NN_Train_FR_CountVsPulses_Norm_v3) =================
chResp = zeros(1, 0);
for p = 1:nPat, for n = 1:numel(chSet{p}), chResp = union(chResp, chSet{p}{n}); end, end
chResp = chResp(:)';
refC = nan(1, nChn);  refPatIdx = nan(1, nChn);  refCond = nan(1, nChn);  refX = nan(1, nChn);
for p = 1:nPat
    P = Pat(p);
    if strcmp(refLevels, 'pair'), idx = P.pairIdx; else, idx = 1:numel(P.cond); end
    [mx, j] = max(P.rate(idx, :), [], 1, 'omitnan');
    upd = isfinite(mx) & (~isfinite(refC) | mx > refC);
    refC(upd) = mx(upd);  refPatIdx(upd) = p;  refCond(upd) = P.cond(idx(j(upd)));  refX(upd) = P.x_n(idx(j(upd)));
end
keepN = ismember(1:nChn, chResp) & isfinite(refC) & refC > minPeak_Hz;
assert(any(keepN), 'No responsive channel has a positive peak rate: cannot normalise.');
for p = 1:nPat
    L = numel(Pat(p).cond);
    Pat(p).rateN = Pat(p).rate ./ refC;
    Pat(p).rateN(:, ~keepN) = NaN;
    Pat(p).rateN_mean = nan(L, 1);  Pat(p).rateN_sem = nan(L, 1);  Pat(p).rateN_n = zeros(L, 1);
    for n = 1:L
        [Pat(p).rateN_mean(n), Pat(p).rateN_sem(n), Pat(p).rateN_n(n)] = ms_(Pat(p).rateN(n, :));
    end
end
Norm = struct('method', 'rate / own peak of the channel (raw baseline-corrected rate, good trials), reference over refLevels of all patterns', ...
    'refLevels', refLevels, 'minPeak_Hz', minPeak_Hz, 'chResp', chResp, 'keep', keepN, 'refC', refC, ...
    'refPat', refPatIdx, 'refCond', refCond, 'refX', refX);

%% ================= RESULT STRUCT =================
meta = struct('dataset', folderName, 'data_folder', data_folder, 'rec_file', rec_file, 'rec_tag', rec_tag, ...
    'date', datestr(now, 'yyyy-mm-dd HH:MM:SS'), 'version', 'NN_Train_FR_CountVsPulses_save_v2', 'FS', FS, ...
    'tail_ms', tail_ms, 'trainStart_ms', trainStart_ms, 'baseWin_ms', baseWin_ms, 'chSetMode', chSetMode, 'chRange', chRange, ...
    'resp_file_used', ternary_(haveResp, resp_file, ''), 'bad_file_used', ternary_(haveBad, bad_file, ''), ...
    'animal', meta_animal, 'area', meta_area, 'paradigm', meta_paradigm, 'note', meta_note, 'refLevels', refLevels, 'minPeak_Hz', minPeak_Hz);
FRRes = struct('meta', meta, 'nChn', nChn, 'patName', {patName}, 'Pat', Pat, 'Norm', Norm, 'recInfo', recInfo, 'recLog', {recLogS});

%% ================= PRINT =================
fprintf('\n============================================================\n');
fprintf('NN TRAIN FR | COUNT / FIRING RATE VS PULSES | SAVE (v2)\n');
fprintf('============================================================\n');
fprintf('Folder   : %s\nRecovery : %s\n', data_folder, rec_file);
fprintf('Window   : first pulse%+g ms ... last pulse +%g ms | baseline [%g, %g] ms\n', trainStart_ms, tail_ms, baseWin_ms);
if haveResp, fprintf('Channels : responsive-channel file %s (chSetMode ''%s'')\n', resp_file, chSetMode);
else,        fprintf('Channels : no responsive-channel file -> chRange %s (chSetMode ''%s'')\n', fmtList_(chRange), chSetMode); end
if haveBad,  fprintf('Bad trials: file %s\n', bad_file); else, fprintf('Bad trials: none (all trials kept)\n'); end
for p = 1:nPat
    P = Pat(p);
    fprintf('\n--- %s (final cond %d): %d levels, %d pair levels | stim: freq %s Hz, amp(col16) %s, start %s ms ---\n', upper(P.name), P.finalCond, ...
        numel(P.cond), numel(P.pairIdx), mat2str(P.stim.freq_Hz, 5), mat2str(P.stim.amp), mat2str(P.stim.start_ms, 4));
    fprintf('  cond  pair   x     trials(bad)  chn    rate(Hz)   addedRate(Hz)   count\n');
    for n = 1:numel(P.cond)
        fprintf('  %3d   %s   %5.1f   %4d (%d)   %3d   %8.1f   %10.1f   %7.2f\n', P.cond(n), ternary_(P.isPair(n), 'Y', '-'), P.x_n(n), ...
            P.nTrialsUsed(n), P.nBad(n), numel(P.chSet{n}), mean(P.rate(n, chResp), 'omitnan'), mean(P.addRate(n, chResp), 'omitnan'), ...
            mean(P.cm(n, chResp), 'omitnan'));
    end
end

fprintf('\n--- NORMALISATION (per channel, own peak over %s levels of sim and seq) ---\n', refLevels);
fprintf('  channels normalised: %d  %s | reference (spikes/s): median %.1f, min %.1f, max %.1f\n', nnz(keepN), fmtList_(find(keepN)), median(refC(keepN)), min(refC(keepN)), max(refC(keepN)));
for p = 1:nPat
    P = Pat(p);  [vN, jN] = max(P.rateN_mean(P.pairIdx));
    fprintf('  %-9s group normalised rate: highest %.3f at x = %g | %d channels\n', P.name, vN, P.x_n(P.pairIdx(jN)), P.rateN_n(P.pairIdx(jN)));
end
fprintf('  stored: all %d channels, chMask per level, Norm, recInfo%s\n', nChn, ternary_(isempty(recLogS), '', ', recLog'));

%% ================= SAVE =================
suffix = '';  if ~isempty(rec_tag), suffix = ['_' rec_tag]; end
if isempty(out_file)
    if isempty(out_folder), out_folder = data_folder; end
    if do_save && ~isfolder(out_folder)
        [okm, msg] = mkdir(out_folder);
        assert(okm, 'Cannot create the output folder %s: %s', out_folder, msg);
        fprintf('Created folder %s\n', out_folder);
    end
    nm = out_name;
    if isempty(nm), nm = [folderName '.FR_CountVsPulses']; end
    nm = regexprep(nm, '\.mat$', '');
    out_file = fullfile(out_folder, [nm suffix '.mat']);
end
fprintf('\nOutput file : %s\n', out_file);
if ~do_save
    fprintf('do_save = false: nothing saved.\n');
elseif isfile(out_file) && ~overwrite
    error('The file already exists and overwrite = false. Nothing was changed. Set overwrite = true to replace it (a dated backup is made first).');
else
    if isfile(out_file)
        bak = strrep(out_file, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
        copyfile(out_file, bak);
        fprintf('Old file copied to %s\n', bak);
    end
    save(out_file, 'FRRes', '-v7.3');
    d_ = dir(out_file);
    fprintf('Saved (%.1f MB).\n', d_(1).bytes/1e6);
end

%% ================= QUICK LOOK =================
if Plot_Check
    cols = [0.00 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556];
    mks = {'o','s','d','^'};
    figure('Color','w','Position',[80 80 1500 430], 'Name', 'Saved result: quick look (pair levels)');
    yl = {'firing rate (spikes/s)', 'added rate of the new stretch (spikes/s)', 'firing rate / own peak of the channel'};
    for k_ = 1:3
        ax = subplot(1, 3, k_);  hold(ax, 'on');  hl = gobjects(1, nPat);
        for p = 1:nPat
            P = Pat(p);  c = cols(mod(p-1, size(cols,1))+1, :);  mk = mks{mod(p-1, numel(mks))+1};
            if k_ == 1, Y = P.rate(P.pairIdx, chResp); elseif k_ == 2, Y = P.addRatePair(:, chResp); else, Y = P.rateN(P.pairIdx, :); end
            m = nan(1, size(Y,1));  se = m;
            for q = 1:size(Y,1), [m(q), se(q)] = ms_(Y(q,:)); end
            hl(p) = errorbar(ax, P.x_n(P.pairIdx), m, se, ['-' mk], 'Color', c, 'MarkerFaceColor', c, 'LineWidth', 1.5, 'CapSize', 4);
        end
        yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
        xlabel(ax, 'pulses per electrode');  ylabel(ax, yl{k_});  legend(ax, hl, patName, 'Box', 'off');  box(ax, 'off');
    end
    sgtitle(sprintf('%s | mean +/- SEM over channels (quick look)', strrep(folderName, '_', '\_')), 'FontSize', 10);
end

%% ================= HELPERS =================
function [m, se, n] = ms_(v)
v = v(isfinite(v));  n = numel(v);
if n == 0, m = NaN;  se = NaN;  return; end
m = mean(v);
if n > 1, se = std(v)/sqrt(n); else, se = NaN; end
end

function i = lb_(s, x)
lo = 1;  hi = numel(s) + 1;
while lo < hi
    m = floor((lo + hi)/2);
    if s(m) < x, lo = m + 1; else, hi = m; end
end
i = lo;
end

function s = fmtList_(v)
v = sort(v(:)');
if isempty(v), s = '[]'; return; end
parts = {};  i = 1;
while i <= numel(v)
    j = i;
    while j < numel(v) && v(j+1) == v(j) + 1, j = j + 1; end
    if j == i,        parts{end+1} = sprintf('%d', v(i));          %#ok<AGROW>
    elseif j == i+1,  parts{end+1} = sprintf('%d %d', v(i), v(j)); %#ok<AGROW>
    else,             parts{end+1} = sprintf('%d:%d', v(i), v(j)); %#ok<AGROW>
    end
    i = j + 1;
end
s = ['[' strjoin(parts, ' ') ']'];
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
