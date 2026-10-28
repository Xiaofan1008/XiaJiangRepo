function c = expTrainISI(subject,varargin)
% EXPTRAINISI  ISI sweep with pulse trains: two electrodes stimulated
% SIMULTANEOUSLY or SEQUENTIALLY with different ISIs (the delay between the
% two trains), interleaved with VISUAL-ONLY trials, in one session. Electrical
% and visual stimulation never occur in the same trial. Default ISIs 5 ms
% (pulses interleaved) and 50 ms (one train, then the other).
%
% STATUS: written but NOT yet run in MATLAB or on the rig. Test with
% 'noIntan',true first (see BENCH SWITCHES), then with Intan.
%
% Built on exp1A.m (rig setup, plugin order, trigger timing, cleanup), the
% multi-electrode plugin estimIntanMulti.m and the event-level recovery logic
% of design_exp_v7_VariableReps.m.
%
%   c = expTrainISI('CJ286','chn',{'A-007','A-012'});
%   c = expTrainISI('CJ286','chn',{'A-007','A-012'},'amp',[60 90], ...
%                    'nPulses',5,'contrast',[0.25 1],'visType','both');
%
% The first input is the SUBJECT (animal ID, a label only, as in exp1A).
%
% STANDALONE: this one file holds the experiment, the condition builder and the
% charge-limit check. It needs only estimIntanMulti.m plus what the
% stimulation PC already has (neurostim, marmolab, the Intan helper functions,
% Psychtoolbox, the DAQ functions).
%
% =========================================================================
% INPUT LABELS  (name-value pairs, e.g. 'amp',[60 90]). Every label has a
% default. Only the SUBJECT and 'chn' must be given: a default channel could
% stimulate the wrong site.
%
%   label          meaning                                   default
%   -------------  ----------------------------------------  -------------
%  ELECTRICAL
%   'chn'          Intan names of electrode A and B          REQUIRED, e.g. {'A-007','A-012'}
%   'amp'          amplitude LEVELS in uA (same on A and B)  75
%                  every electrical condition runs at each level
%   'ampB'         levels for B if different from A          [] (= same as amp)
%   'freq'         pulse frequency within a train, Hz        100
%   'nPulses'      pulses per electrode in the FULL train    5
%   'isi'          ISI LIST in ms: delay of the second       [5 50]
%                  electrode's train after the first one's
%                  (every ISI is run). AVOID whole multiples of the pulse
%                  period (10, 20, ... ms at 100 Hz): the second train's
%                  pulses then coincide with the first's. A warning prints.
%                  At 100 Hz, 5 ms interleaves the pulses; with 5 pulses the
%                  first train lasts 40 ms, so 50 ms means one train, then
%                  the other.
%   'order'        which electrode first: 'both'|'AB'|'BA'   'both'
%   'fpd','spd'    first / second phase width, us            200, 200
%   'ipi'          interphase delay, us                      100
%  VISUAL (the same labels will be used in the other experiment files)
%   'withVisual'   include visual-only trials                true
%   'visType'      'spot' | 'gabor' | 'both'                 'spot'
%                  spot = uniform disc; gabor = grating in the same circular
%                  window (same position, size, onset and duration). The
%                  gabor is used ONLY if you ask for it.
%   'polarity'     'bright' | 'dark' | 'both'                'both'
%                  (applies to the spot; for the gabor it sets the phase)
%   'contrast'     visual contrast LEVELS, 0..1              [0.06 0.12 0.25 0.5 1]
%                  (each contrast level costs 2 conditions with both polarities, 1 with one; an amplitude level costs 5N)
%   'visDur'       visual duration, ms (rounded to frames)   250
%   'diameter'     deg                                       6
%   'rfX','rfY'    deg, centre of the patch. PLACEHOLDER 0, 0:   0, 0
%                  set from RF mapping (a warning prints if both are 0)
%   'ori'          gabor orientation, deg                    0
%   'sf'           gabor spatial frequency, cycles/deg       1
%  TIMING
%   'pulseTime'    ms from trial start to the trigger        50
%   'trialDur'     trial length, ms (keep under ~550)        500
%   'trigDuration' ms the trigger line stays high            100
%   'iti'          inter-trial interval, s (inert, exp1A)    1.5
%   'nReps'        repeats of every condition                20
%   'headerPause'  s of recording before first/after last    25
%  SAFETY (checkChargeLimit)
%   'siteArea'     electrode site area, um^2                 703
%   'capacity'     charge capacity, uC/cm^2 (UNCONFIRMED)    3600
%   'kMax'         Shannon k limit                           1.7
%   'allowOver'    true = warn instead of refuse             false
%  BENCH SWITCHES (none change the science)
%   'noIntan'      no TCP; trigger line and visual only      false
%   'debug'        rig debug mode                            false
%   'skipSync'     Psychtoolbox SkipSyncTests                false
%   'confirm'      print a summary and ask before running    true
%
% =========================================================================
% CONDITIONS (built by the local function buildTrainConds below)
% =========================================================================
% For every amplitude level, with N = nPulses:
%   simultaneous     A and B pulses together. N levels: the first k pulses
%                    of each train, k = 1..N (event-level recovery)
% and for every ISI in 'isi':
%   A -> B           B's train is delayed by the ISI. Levels add one pulse
%                    at a time in TIME ORDER (all pulses of both trains are
%                    sorted). For a short ISI (5 ms, N = 5) that is 2N levels:
%                    [A B] pulses = [1 0] [1 1] [2 1] [2 2] ... [N N]
%                    For a long ISI (50 ms) it is [1 0] [2 0] ... [N 0] [N 1]
%                    ... [N N]. Pulses at exactly the same time share a level.
%   B -> A           the same with B first
% Only the LAST level of each pattern (all N pulses on both) is the
% fully delivered condition you report; the earlier levels exist so the
% response hidden under each pulse's stimulation artifact can be recovered
% (design_exp_v7 trainMode 2). All levels use the same number of repeats.
% Shared by all amplitudes:
%   blank            nothing delivered
%   visual           one condition per (visType x polarity x contrast), only
%                    if 'withVisual' is true. No electrical stimulation.
% Defaults: 1 amplitude x (5 + 2 ISIs x (10 + 10)) = 45 electrical levels, but
% identical levels are listed ONCE (the first pulse of A alone, [1 0], is the
% same condition for every ISI, and likewise [0 1] for B), which leaves 43.
% With 1 blank and 10 visual (bright and dark spot) that is 54 conditions, 20 repeats, about 36 min.
% Each extra ISI adds up to 20 conditions (about 13 min); each extra contrast
% adds 2 (about 80 s, both polarities). A shared level keeps the ISI of its first occurrence in
% 'info'; rebuild the levels from nsp and ptd in the analysis. THE SESSION LENGTH IS PRINTED BEFORE
% YOU CONFIRM: cut nPulses, the number of amplitudes, nReps or order if it is
% too long.
%
% =========================================================================
% TRIAL STRUCTURE
% =========================================================================
%   0 ms        trial starts (first frame)
%   pulseTime   trigger fires; electrode A (or B, for B -> A) pulses at once
%               and the visual patch appears on the same frame   (default 50 ms)
%   + ISI       the second electrode starts, via its post-trigger delay
%   trialDur    trial ends                                       (default 500 ms)
%               (the visual patch runs pulseTime .. pulseTime+visDur)
% The parameters are sent and confirmed (SET + verifyParameters) in
% beforeTrial, BEFORE the trial starts. Align rasters and PSTHs to the TRIGGER
% recorded in Open Ephys (electrical) or to the photodiode (visual). The
% electrical response starts earlier than the visual one, so handle the
% latency difference in the analysis, not by shifting the stimuli.
% Keep trialDur under about 550 ms so the trial finishes before the Intan hole.
%
% SESSION START AND END. estimIntanMulti waits headerPause seconds after
% RECORD and again before STOP so Intan can write its file header and closing
% buffer; both Intan and NP record during it.
%
% =========================================================================
% PROVENANCE TAGS (as in exp1A)
%   [EXP1A] copied from exp1A.m, behaviour unchanged
%   [NEW]   written for this experiment
% =========================================================================

aoutChannel = 0;   % [EXP1A] Intan trigger: MCC ANALOG out, via DaqAOut

% [EXP1A] argument handling
if ~exist('subject','var')
    error('No subject name provided.');
end
validateattributes(subject,{'char'},{'nonempty'},'','subject',1);
p = inputParser();
p.KeepUnmatched = true;

% [NEW] electrical
p.addParameter('chn',{},@(x) iscellstr(x) && numel(x)==2); %#ok<ISCLSTR>
p.addParameter('amp',75,@(x) isnumeric(x) && ~isempty(x) && all(x>0));
p.addParameter('ampB',[],@(x) isnumeric(x) && all(x>0));
p.addParameter('freq',100,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('nPulses',5,@(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('isi',[5 50],@(x) isnumeric(x) && ~isempty(x) && all(x>0));
p.addParameter('order','both',@(x) ismember(lower(x),{'both','ab','ba'}));
p.addParameter('fpd',200,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('spd',200,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('ipi',100,@(x) isnumeric(x) && isscalar(x) && x>=0);
% [NEW] visual
p.addParameter('withVisual',true,@islogical);
p.addParameter('visType','spot',@(x) ismember(lower(x),{'spot','gabor','both'}));
p.addParameter('polarity','both',@(x) ismember(lower(x),{'bright','dark','both'}));
p.addParameter('contrast',[0.06 0.12 0.25 0.5 1],@(x) isnumeric(x) && ~isempty(x) && all(x>0) && all(x<=1));
p.addParameter('visDur',250,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('diameter',6,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('rfX',0,@isnumeric);
p.addParameter('rfY',0,@isnumeric);
p.addParameter('ori',0,@isnumeric);
p.addParameter('sf',1,@(x) isnumeric(x) && isscalar(x) && x>0);
% [NEW] timing
p.addParameter('pulseTime',50,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('trialDur',500,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('trigDuration',100,@isnumeric);
p.addParameter('iti',1.5,@isnumeric);
p.addParameter('nReps',20,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('headerPause',25,@(x) isnumeric(x) && isscalar(x) && x>=0);
% [NEW] safety
p.addParameter('siteArea',703,@isnumeric);
p.addParameter('capacity',3600,@isnumeric);
p.addParameter('kMax',1.7,@isnumeric);
p.addParameter('allowOver',false,@islogical);
% [NEW] bench switches
p.addParameter('noIntan',false,@islogical);
p.addParameter('debug',false,@islogical);
p.addParameter('skipSync',false,@islogical);
p.addParameter('confirm',true,@islogical);

p.parse(varargin{:});
args = p.Results;

if isempty(args.chn)
    error('expTrainISI:chn','Give the two Intan channel names, e.g. ''chn'',{''A-007'',''A-012''}.');
end
if args.trialDur > 550
    warning('expTrainISI:trialDur', ...
        'trialDur %g ms is over ~550 ms: the end of the trial may fall in the Intan hole.',args.trialDur);
end

% [NEW] train timing, checked before anything is opened
period_ms = 1000/args.freq;
period_us = round(1e6/args.freq);
isi_ms = unique(args.isi(:)');
for q = isi_ms
    if args.nPulses > 1 && mod(round(q*1000),period_us) == 0   % irrelevant for single pulses
        warning('expTrainISI:isi','ISI %g ms is a whole number of pulse periods (%g ms): the second train''s pulses coincide with the first''s.',q,period_ms);
    end
end
trainEnd = args.pulseTime + max(isi_ms) + (args.nPulses-1)*period_ms;   % last pulse starts here
if args.withVisual && args.pulseTime + args.visDur > args.trialDur
    error('expTrainISI:timing','The visual patch (%g ms from %g ms) does not fit in trialDur %g ms.', ...
        args.visDur,args.pulseTime,args.trialDur);
end
if args.withVisual && args.rfX == 0 && args.rfY == 0
    warning('expTrainISI:rf','rfX and rfY are both 0: this is a PLACEHOLDER. Set them from RF mapping.');
end
if args.pulseTime + args.trigDuration > args.trialDur || trainEnd + 2 > args.trialDur
    error('expTrainISI:timing', ...
        'The train does not fit: last pulse at %.1f ms, trigger high until %.1f ms, trialDur %g ms.', ...
        trainEnd,args.pulseTime+args.trigDuration,args.trialDur);
end

% [NEW] amplitudes and conditions
ampA = args.amp(:)';
if isempty(args.ampB), ampB = ampA; else, ampB = args.ampB(:)'; end
if numel(ampB) ~= numel(ampA)
    error('expTrainISI:amp','ampB needs the same number of levels as amp (%d).',numel(ampA));
end
vis = struct('on',args.withVisual,'type',lower(args.visType), ...
             'polarity',lower(args.polarity),'contrast',args.contrast(:)','sf',args.sf);
[cond,info] = buildTrainConds(ampA,ampB,args.nPulses,period_us,isi_ms,args.order,vis);
nCond = numel(cond);
nTrials = nCond * args.nReps;

% [NEW] charge limit, checked BEFORE anything is opened
if args.fpd ~= args.spd
    warning('expTrainISI:unbalanced','fpd (%g us) ~= spd (%g us) at equal amplitude: the pulse is NOT charge balanced.',args.fpd,args.spd);
end
checkChargeLimit([ampA ampB],max(args.fpd,args.spd), ...
    'siteArea',args.siteArea,'capacity',args.capacity,'kMax',args.kMax,'allowOver',args.allowOver);

nElec = sum(~strcmp({info.kind},'visual') & ~strcmp({info.kind},'blank'));
fprintf('\nexpTrainISI: %d conditions x %d repeats = %d trials (about %.0f min at ~2 s per trial)\n', ...
    nCond,args.nReps,nTrials,nTrials*2/60);
fprintf('  %d electrical + %d visual + %d blank conditions\n', ...
    nElec,sum(strcmp({info.kind},'visual')),sum(strcmp({info.kind},'blank')));
fprintf('  channels A = %s, B = %s; amplitude A [%s] uA, B [%s] uA\n', ...
    args.chn{1},args.chn{2},num2str(ampA),num2str(ampB));
fprintf('  train: %d pulses at %g Hz (period %.3g ms), ISI [%s] ms, last pulse at %.1f ms; phases %g / %g us, interphase %g us\n', ...
    args.nPulses,args.freq,period_ms,num2str(isi_ms),trainEnd,args.fpd,args.spd,args.ipi);
if args.withVisual
    fprintf('  visual: %s, %s, contrasts [%s], %g ms, %g deg at (%g,%g)\n', ...
        vis.type,vis.polarity,num2str(vis.contrast),args.visDur,args.diameter,args.rfX,args.rfY);
end
for k = 1:nCond
    fprintf('  %3d  %-34s amp %-14s nsp %-10s ptd %s us\n',k,cond(k).label, ...
        mat2str(cond(k).amp),mat2str(cond(k).nsp),mat2str(cond(k).ptd));
end

if args.confirm
    ans_ = input('Run this? (y/n) ','s');
    if ~strcmpi(strtrim(ans_),'y')
        c = [];
        disp('Cancelled.');
        return
    end
end

import neurostim.*
commandwindow;

% [EXP1A] PTB sync workaround
if args.skipSync
    Screen('Preference','SkipSyncTests',1);
    warning('expTrainISI:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
end

% [EXP1A] rigcfg builds cic, opens the screen, and creates the oephys plugin
% and the diode stimulus.
c = marmolab.rigcfg('debug',args.debug,'subject',subject,p.Unmatched);

% [EXP1A] oephys.fakeConnection MUST stay false (its two blocking messages are
% part of the pulse recovery, see exp1A).
c.oephys.fakeConnection = false;

% -------------------------------------------------------------------------
% [NEW] ELECTRICAL STIMULUS
% -------------------------------------------------------------------------
e = estimIntanMulti(c,'estim');
e.chnList  = args.chn(:)';
e.fpd      = args.fpd;
e.spd      = args.spd;
e.ipi      = args.ipi;
e.trainPeriod_us = round(1e6/args.freq);   % shared by every channel's train
e.testMode = args.noIntan;
e.headerPause = args.headerPause;
e.setConditions(cond);

% frame-grid pulse time. estim is "on" from e.on for trigDuration ms; its
% onsetFunction fires the trigger.
ivl = 1000/c.screen.frameRate;
nf  = round(args.pulseTime/ivl);
e.on       = nf*ivl;
e.duration = args.trigDuration;
fprintf('expTrainISI timing: pulse requested at %g ms, on frame %d, delivered at %.2f ms after trial start\n', ...
    args.pulseTime,nf+1,e.on);

% [NEW] the visual patch, one gabor object for spot and gabor alike. Spot =
% spatial frequency 0 (a disc of uniform luminance). The design writes
% contrast, phase and spatial frequency every trial. Contrast 0 hides the
% patch, and the design sets 'on' to Inf on every non-visual trial, so electrical
% and blank trials show nothing AND the photodiode does not flash.
% PHASE (checked in the shader source, GLSLShaders/gabor.frag.txt): luminance =
% background x (1 + contrast x sin(2 pi f x + phase)). A spot has f = 0, so its
% pattern is sin(phase): bright = 90, dark = 270. Phase 0 or 180 would be
% INVISIBLE (sin = 0). The gabor keeps 0 / 180 (a grating and its inverse).
if args.withVisual
    g = stimuli.gabor(c,'flash');
    g.X = args.rfX;
    g.Y = args.rfY;
    g.mask   = 'CIRCLE';
    g.sigma  = args.diameter/2;      % CIRCLE: sigma is the radius
    g.width  = ceil(args.diameter) + 1;
    g.height = ceil(args.diameter) + 1;
    g.orientation = args.ori;
    g.phaseSpeed  = 0;
    g.frequency   = 0;               % the design overrides this per trial
    g.phase       = 0;               % the design overrides this per trial
    g.contrast    = NaN;             % the design writes this every trial
    g.on       = e.on;               % same frame as the pulse
    g.duration = args.visDur;
end

% -------------------------------------------------------------------------
% [EXP1A] Intan TCP, opened before c.run covers the console
% -------------------------------------------------------------------------
if ~args.noIntan
    fprintf('\n---- connect Intan to the TCP/IP server now ----\n');
    e.tcp = local_TCPIP_server();
    fprintf(e.tcp,'NOMCS');
    fprintf('Intan connected.\n');
end

% [EXP1A] Intan trigger. NOT gated on debug.
DAQ = DaqDeviceIndex;
if isempty(DAQ)
    error('expTrainISI: DaqDeviceIndex found no DAQ, so no trigger can be sent.');
end
fprintf('DAQ index %s, trigger on analog out %d\n',mat2str(DAQ),aoutChannel);
DaqAOut(DAQ,aoutChannel,0);          % line starts low

e.daq         = DAQ;
e.aoutChannel = aoutChannel;

% Every trial fires the trigger, including blank and visual-only trials (all
% channels are off then), so timing is the same in every condition.
e.onsetFunction  = @(o,t) firePulse(o);
e.offsetFunction = @(o,t) DaqAOut(DAQ, aoutChannel, 0);

% [EXP1A] session hooks
if ~args.debug
    c.addScript('beforeExperiment', @(c) doBeforeExperiment(c));
    c.addScript('afterExperiment', @(c) doAfterExperiment(c));
end

    function doBeforeExperiment(cc)
        io.forceBlackDiode(cc);
        % [NEW] move the condition table next to the Neurostim file, so the
        % two always stay together: <subject>.<paradigm>.<HHMMSS>.conds.mat
        try
            movefile(condFile,[cc.fullFile '.conds.mat']);
            fprintf('Condition table moved to %s.conds.mat\n',cc.fullFile);
        catch me
            warning('expcond:move','Could not move the condition table (%s). It stays in %s.', ...
                    me.message,pwd);
        end
    end
    function doAfterExperiment(cc) %#ok<INUSD>
        io.beep3;
    end
    function firePulse(o)
        DaqAOut(DAQ,aoutChannel,1);
        o.markPulse();               % stamp after the edge, never before
    end

% [EXP1A] diode flasher: c.diode is disabled; with the visual patch the
% flasher is attached to it so it marks the visual onset.
addDiodeFlasher(c.diode,'enabled',false);
if args.withVisual
    addDiodeFlasher(g,'size',0.075,'location','nw');
end

c.trialDuration = args.trialDur;      % ms
c.iti = args.iti*1e3;                 % s -> ms

% -------------------------------------------------------------------------
% [NEW] DESIGN: one factor, the condition index. The visual patch properties
% are listed in the same factor: contrast 0 (hidden) on every non-visual
% condition.
% -------------------------------------------------------------------------
d = neurostim.design('myDesign');
d.fac1.estim.stimCond = 1:nCond;
if args.withVisual
    d.fac1.flash.contrast  = [info.contrast];
    d.fac1.flash.phase     = [info.phase];
    d.fac1.flash.frequency = [info.sf];
    % The patch, and so its photodiode flash, exists ONLY on visual trials:
    % on = Inf means 'never on' (electrical and blank trials show nothing).
    visOn = inf(1,nCond);
    visOn(strcmp({info.kind},'visual')) = e.on;
    d.fac1.flash.on = visOn;
end

blk = block('myBlock',d);
blk.nrRepeats = args.nReps;

c.subject  = subject;
c.paradigm = 'expTrainISI';
c.addPropsToInform('estim.stimCond');

% [NEW] keep the condition table with the data
stamp = datestr(now,'yyyymmdd_HHMMSS');
condFile = sprintf('%s_expTrainISI_conds_%s.mat',subject,stamp);
chnList = e.chnList;                  %#ok<NASGU>
save(condFile,'cond','info','chnList','args','isi_ms','period_us');
fprintf('Condition table saved to %s\n',fullfile(pwd,condFile));

% [EXP1A] plugin order: oephys first, so its two blocking messages sit between
% the previous pulse and estim's parameter exchange. DO NOT reorder.
if args.withVisual
    c.setPluginOrder('oephys','flash','estim');
else
    c.setPluginOrder('oephys','estim');
end

try
    c.run(blk);
catch me
    cleanup(e);
    rethrow(me);
end
cleanup(e);

end

function cleanup(e)
% [EXP1A] Release everything, in an order that works even when the run died
% early. Each step is wrapped individually so one failure cannot prevent the
% others.
try; sca; end                                      %#ok<TRYNC>
try; ListenChar(0); end                            %#ok<TRYNC>
try
    if ~isempty(e.daq); DaqAOut(e.daq,e.aoutChannel,0); end
end                                                %#ok<TRYNC>
try
    if ~isempty(e.tcp); fclose(e.tcp); end
end                                                %#ok<TRYNC>
end

% =========================================================================
% LOCAL FUNCTIONS (so this file is standalone)
% =========================================================================

function [cond,info] = buildTrainConds(ampA,ampB,N,period_us,isi_ms,order,vis)
% BUILDTRAINCONDS  Train conditions for one electrode pair, with event-level
% recovery levels (design_exp_v7 trainMode 2), plus blank and visual-only
% conditions.
%
%   [cond,info] = buildTrainConds(ampA,ampB,N,period_us,isi_ms,order,vis)
%
% ampA, ampB  vectors of equal length: one entry per amplitude LEVEL (uA)
% N           pulses per electrode in the full train
% period_us   time between pulses of one train (us)
% isi_ms      list of ISIs (ms): B's train starts this long after A's
%             (A -> B), or A's after B's (B -> A)
% order       'both' | 'AB' | 'BA'
% vis         struct: on, type ('spot'|'gabor'|'both'), polarity
%             ('bright'|'dark'|'both'), contrast (vector), sf (cycles/deg)
%
% RECOVERY LEVELS. All pulse times of both trains are sorted. Level 1 holds
% the earliest pulse, and each next level adds the next pulse in time, so a
% level is the pair [pulses of A, pulses of B] delivered so far. A pulse that
% falls at exactly the same time as another one is added at the SAME level
% (this happens when the ISI is a whole number of periods). Normally this
% gives 2N levels per ISI and order; the simultaneous pattern gives N.
%
% cond  struct array for estimIntanMulti (amp, ptd, nsp, label; two channels
%       [A B]). A channel with nsp = 0 is given amp 0 and ptd 0 (switched off).
% info  one entry per condition: kind ('sim'|'AB'|'BA'|'blank'|'visual'),
%       isi (ms, 0 if none), level, nLevels, isFinal, ampA, ampB, and the
%       visual settings contrast, phase (spot 90 bright / 270 dark; gabor 0 / 180), sf.
%       Electrical and blank conditions have contrast 0 (patch hidden).
%
% PTD is in MICROSECONDS (what estimIntanMulti writes into StimParams).

order = lower(order);
assert(ismember(order,{'both','ab','ba'}), ...
    'buildTrainConds: order must be ''both'', ''AB'' or ''BA''.');
isi_ms = unique(isi_ms(:)');

cond = struct('amp',{},'ptd',{},'nsp',{},'label',{});
info = struct('kind',{},'isi',{},'level',{},'nLevels',{},'isFinal',{},'ampA',{},'ampB',{}, ...
              'contrast',{},'phase',{},'sf',{});
keys = {};

    function add(aA,aB,nA,nB,ptdA,ptdB,label,kind,isi,lev,nLev,vc,vp,vsf)
        if nA < 1, aA = 0; ptdA = 0; end        % channel off
        if nB < 1, aB = 0; ptdB = 0; end
        amp = [aA aB];  ptd = [ptdA ptdB];  nsp = [max(nA,0) max(nB,0)];
        key = [mat2str(amp) mat2str(ptd) mat2str(nsp) sprintf('|%g|%g|%g',vc,vp,vsf)];
        if any(strcmp(keys,key)), return; end    % identical condition, keep one
        keys{end+1} = key; %#ok<AGROW>
        k = numel(cond) + 1;
        cond(k).amp = amp;  cond(k).ptd = ptd;  cond(k).nsp = nsp;  cond(k).label = label;
        info(k).kind = kind;  info(k).isi = isi;  info(k).level = lev;  info(k).nLevels = nLev;
        info(k).isFinal = (lev == nLev);
        info(k).ampA = aA;  info(k).ampB = aB;
        info(k).contrast = vc;  info(k).phase = vp;  info(k).sf = vsf;
    end

multi = numel(ampA) > 1;
for L = 1:numel(ampA)
    aA = ampA(L);  aB = ampB(L);
    if multi, tag = sprintf(' @%g/%g uA',aA,aB); else, tag = ''; end

    % simultaneous: both trains together, level k = first k pulses
    lv = seqLevels(N,period_us,0,0);
    for k = 1:size(lv,1)
        add(aA,aB,lv(k,1),lv(k,2),0,0,sprintf('sim %d/%d%s',k,size(lv,1),tag), ...
            'sim',0,k,size(lv,1),0,0,0);
    end
    for isi = isi_ms
        us = round(isi*1000);
        if any(strcmp(order,{'both','ab'}))          % A first, B delayed by isi
            lv = seqLevels(N,period_us,0,us);
            for k = 1:size(lv,1)
                add(aA,aB,lv(k,1),lv(k,2),0,us, ...
                    sprintf('A->B %gms %d/%d [%d %d]%s',isi,k,size(lv,1),lv(k,1),lv(k,2),tag), ...
                    'AB',isi,k,size(lv,1),0,0,0);
            end
        end
        if any(strcmp(order,{'both','ba'}))          % B first, A delayed by isi
            lv = seqLevels(N,period_us,us,0);
            for k = 1:size(lv,1)
                add(aA,aB,lv(k,1),lv(k,2),us,0, ...
                    sprintf('B->A %gms %d/%d [%d %d]%s',isi,k,size(lv,1),lv(k,1),lv(k,2),tag), ...
                    'BA',isi,k,size(lv,1),0,0,0);
            end
        end
    end
end

add(0,0,1,1,0,0,'blank','blank',0,1,1,0,0,0);
if vis.on
    types = {vis.type};  if strcmp(vis.type,'both'), types = {'spot','gabor'}; end
    pols  = {vis.polarity};  if strcmp(vis.polarity,'both'), pols = {'bright','dark'}; end
    for it = 1:numel(types)
        for ip = 1:numel(pols)
            for ic = 1:numel(vis.contrast)
                if strcmp(types{it},'spot')
                    if strcmp(pols{ip},'bright'), ph = 90; else, ph = 270; end   % sin(phase)
                else
                    if strcmp(pols{ip},'bright'), ph = 0; else, ph = 180; end
                end
                if strcmp(types{it},'spot'), sf = 0; else, sf = vis.sf; end
                lab = sprintf('visual %s %s c%g',types{it},pols{ip},vis.contrast(ic));
                add(0,0,1,1,0,0,lab,'visual',0,1,1,vis.contrast(ic),ph,sf);
            end
        end
    end
end
end

function lv = seqLevels(N,period_us,dA,dB)
% SEQLEVELS  Cumulative [nA nB] after each distinct pulse time. A pulses are
% at dA + (0..N-1)*period, B pulses at dB + (0..N-1)*period (all in us).
tA = dA + (0:N-1)*period_us;
tB = dB + (0:N-1)*period_us;
t  = unique([tA tB]);
lv = zeros(numel(t),2);
for i = 1:numel(t)
    lv(i,:) = [sum(tA <= t(i)) sum(tB <= t(i))];
end
end

function lim = checkChargeLimit(amps_uA,pw_us,varargin)
% CHECKCHARGELIMIT  Refuse stimulation amplitudes above the charge limit.
%
%   lim = checkChargeLimit(amps_uA, pw_us)
%   lim = checkChargeLimit(amps_uA, pw_us, 'siteArea',703, 'capacity',3600, 'kMax',1.7)
%
% amps_uA  every amplitude the run will deliver (any shape; values <= 0 are
%          ignored because they mean "channel off")
% pw_us    phase duration in us (first = second phase, charge balanced)
%
% Two limits, and the lower one applies (NeuroNexus microstimulation guide):
%   electrode limit : capacity [uC/cm^2] x site area [cm^2]       (uC per phase)
%   tissue limit    : Shannon criterion, Q <= sqrt(A x 10^k)      (uC per phase)
%
% Defaults match the marmoset probe plan: 703 um^2 site, capacity 3600 uC/cm^2
% (3 x iridium oxide from the Z-Coat table, UNCONFIRMED, so use 1200 to be
% cautious), k = 1.7.
%
% Errors if any amplitude is over the limit, unless 'allowOver' is true, in
% which case it warns. Returns a struct with the limits.

p = inputParser();
p.addParameter('siteArea',703);        % um^2
p.addParameter('capacity',3600);       % uC/cm^2
p.addParameter('kMax',1.7);
p.addParameter('allowOver',false);
p.addParameter('verbose',true);
p.parse(varargin{:});
a = p.Results;

A  = a.siteArea * 1e-8;                % cm^2
Qe = a.capacity * A;                   % uC per phase
Qs = sqrt(A * 10^a.kMax);              % uC per phase
Q1 = min(Qe,Qs);
if Qe < Qs, binds = 'electrode'; else, binds = 'tissue (Shannon)'; end
maxI = Q1 / pw_us * 1e6;               % uC / us = A  ->  uA

lim = struct('electrode_nC',Qe*1e3,'tissue_nC',Qs*1e3,'binds',binds, ...
             'maxAmp_uA',maxI,'pw_us',pw_us);

amps = amps_uA(:);
amps = amps(amps > 0);

if a.verbose
    fprintf(['charge limit: site %.0f um^2, capacity %.0f uC/cm^2, k %.2f\n' ...
             '  electrode limit %.1f nC/phase, tissue limit %.1f nC/phase -> %s binds\n' ...
             '  at %g us per phase the maximum amplitude is %.1f uA\n'], ...
        a.siteArea,a.capacity,a.kMax,Qe*1e3,Qs*1e3,binds,pw_us,maxI);
end

over = amps(amps * pw_us * 1e-6 > Q1*(1+1e-9));    % uA x us = 1e-6 uC
if ~isempty(over)
    msg = sprintf('Amplitude %.1f uA at %g us per phase = %.1f nC, over the %.1f nC limit (max %.1f uA).', ...
        max(over),pw_us,max(over)*pw_us*1e-3,Q1*1e3,maxI);
    if a.allowOver
        warning('checkChargeLimit:overLimit','%s (allowOver is ON)',msg);
    else
        error('checkChargeLimit:overLimit', ...
            '%s Lower the amplitude or the pulse width.',msg);
    end
elseif a.verbose
    fprintf('  all requested amplitudes are within the limit.\n');
end
end
