function c = expSinglePulseISI(subject,varargin)
% EXPSINGLEPULSEISI  Single-pulse mechanism check: simultaneous vs sequential,
% ISI sweep, both orders, with the single-electrode recovery conditions.
%
% STATUS: written but NOT yet run in MATLAB or on the rig. Test with
% 'noIntan',true first (see BENCH SWITCHES), then with Intan.
%
% Built on exp1A.m (rig setup, plugin order, trigger timing, cleanup) and the
% multi-electrode stimulation plugin estimIntanMulti.m. Every value below is
% an input with a default; change anything when you call it:
%
%   c = expSinglePulseISI('CJ286','chn',{'A-007','A-012'});
%   c = expSinglePulseISI('CJ286','chn',{'A-007','A-012'},'amp',60, ...
%                         'isi',[2 4 6 8 10],'nReps',20);
%
% The first input is the SUBJECT: the animal's ID, as a string. It is only a
% label (it names the neurostim data file and the Open Ephys recording
% folder), exactly as in exp1A (e.g. 'CJ286'). It changes nothing else.
%
% STANDALONE: this one file holds the experiment, the condition builder and the
% charge-limit check. It needs only estimIntanMulti.m (a MATLAB class has to be
% its own file) plus what the stimulation PC already has (neurostim, marmolab,
% the Intan helper functions, Psychtoolbox, the DAQ functions).
%
% =========================================================================
% INPUT LABELS  (name-value pairs, e.g. 'amp',60). Every label has a default,
% so leave out any you do not want to change. Only the SUBJECT (first input)
% and 'chn' must be given: a default channel could stimulate the wrong site.
%
%   label          meaning                                   default
%   -------------  ----------------------------------------  -------------
%  ELECTRICAL
%   'chn'          Intan names of electrode A and B          REQUIRED, e.g. {'A-007','A-012'}
%   'amp'          amplitude LEVELS in uA, for A (and B)     75
%                  e.g. [25 50 75 90] runs every condition at each level
%   'ampB'         levels for B if different from A          [] (= same as amp)
%   'isi'         ISIs A->B in ms (simultaneous always in)  [3 6 9 12 15]
%   'order'        which electrode first: 'both'|'AB'|'BA'   'both'
%   'fpd'          first-phase width, us                     200
%   'spd'          second-phase width, us                    200
%                  (amplitude is the same in both phases, so keep
%                  fpd = spd for a charge-balanced pulse; a warning
%                  is printed if they differ)
%   'ipi'          interphase delay, us                      100
%  TIMING
%   'pulseTime'    ms from trial start to the trigger        50
%   'trialDur'     trial length, ms (keep under ~550)        400
%   'trigDuration' ms the trigger line stays high            100
%   'iti'          inter-trial interval, s (inert, exp1A)    1.5
%   'nReps'        repeats of every condition                30
%   'headerPause'  s of recording before first/after last    25
%  SAFETY (checkChargeLimit)
%   'siteArea'     electrode site area, um^2                 703
%   'capacity'     charge capacity, uC/cm^2 (UNCONFIRMED)    3600
%   'kMax'         Shannon k limit                           1.7
%   'allowOver'    true = warn instead of refuse             false
%  VISUAL ARM (off by default)
%   'withVisual'   add a visual-only condition               false
%   'rfX','rfY'    spot centre, deg                          0, 0
%   'diameter'     spot diameter, deg                        6
%   'flashDur'     ms                                        100
%   'flashPolarity' 'bright'|'dark'     'bright'
%   'flashContrast' 0..1                                     1
%  BENCH SWITCHES (none change the science)
%   'noIntan'      no TCP; trigger line and visual only      false
%   'debug'        rig debug mode                            false
%   'skipSync'     Psychtoolbox SkipSyncTests                false
%   'confirm'      print a summary and ask before running    true
%
% Examples:  order 'AB' only:   expSinglePulseISI('CJ286','chn',{'A-007','A-012'},'order','AB')
%            different amplitude on B: expSinglePulseISI('CJ286','chn',{'A-007','A-012'},'amp',40,'ampB',60)
% =========================================================================

% =========================================================================
% CONDITIONS (built by the local function buildPairConds below). A and B are
% chn{1} and chn{2}.
% =========================================================================
%   simultaneous     A and B together
%   A -> B, ISI g    A first, B g ms later          (one per ISI)
%   B -> A, ISI g    B first, A g ms later          (one per ISI)
%   A only, B only   single-electrode recovery conditions, no delay
%   blank            nothing delivered (trigger still fires, all channels off)
%   visual only      only if 'withVisual' is true
% With 5 ISIs, both orders and ONE amplitude: 1 + 10 + 2 + 1 = 14 conditions.
% With several amplitude levels, every level repeats the first three lines
% (blank and visual only are added once).
%
% =========================================================================
% TRIAL STRUCTURE
% =========================================================================
%   0 ms        trial starts (first frame)
%   pulseTime   trigger fires; electrode A pulses at once      (default 50 ms)
%   + ISI       second electrode of a sequential pair, via its post-trigger
%               delay (ptd), e.g. 5 ms later
%   trialDur    trial ends                                   (default 400 ms)
%
% WHY pulseTime = 50. The Intan parameters are sent and confirmed (SET plus
% verifyParameters) in beforeTrial, before the trial starts. The 50 ms is
% extra margin so the stimulator is certainly ready, and it gives the Intan
% trial segment a short pre-pulse baseline. pulseTime 0 (trigger on the first
% frame) has NOT been tested; try it in the lab if you want it, because exp1A
% always fires at 242 ms or later.
%
% ALIGNMENT. Align rasters and PSTHs to the TRIGGER recorded in Open Ephys,
% not to trial start. Time zero is then the first pulse, and the pre-trigger
% period is the baseline. Neuropixels records continuously, so its baseline
% can also come from the ITI before the trigger. Keep trialDur under about
% 550 ms so the trial finishes before the Intan hole (see estimIntan.m). The
% pulse time is rounded to the frame grid; the delivered value is printed.
%
% SESSION START AND END. estimIntanMulti waits headerPause seconds (default 25)
% after RECORD and again before STOP. That wait is only so Intan can write its
% file header and closing buffer (RunMultiElect_exp.m used 20 s). I found no
% code that uses it for alignment or as a baseline, but both Intan and NP are
% recording during it, so it is spontaneous activity if you want it.
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
p.addParameter('isi',[3 6 9 12 15],@(x) isnumeric(x) && all(x>=0));
p.addParameter('order','both',@(x) ismember(lower(x),{'both','ab','ba'}));
p.addParameter('fpd',200,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('spd',200,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('ipi',100,@(x) isnumeric(x) && isscalar(x) && x>=0);
% [NEW] timing
p.addParameter('pulseTime',50,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('trialDur',400,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('trigDuration',100,@isnumeric);
p.addParameter('iti',1.5,@isnumeric);
p.addParameter('nReps',30,@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('headerPause',25,@(x) isnumeric(x) && isscalar(x) && x>=0);
% [NEW] safety
p.addParameter('siteArea',703,@isnumeric);
p.addParameter('capacity',3600,@isnumeric);
p.addParameter('kMax',1.7,@isnumeric);
p.addParameter('allowOver',false,@islogical);
% [NEW] visual arm
p.addParameter('withVisual',false,@islogical);
p.addParameter('rfX',0,@isnumeric);
p.addParameter('rfY',0,@isnumeric);
p.addParameter('diameter',6,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('flashDur',100,@isnumeric);
p.addParameter('flashPolarity','bright',@(x) ismember(lower(x),{'bright','dark'}));
p.addParameter('flashContrast',1,@(x) isnumeric(x) && x>=0 && x<=1);
% [NEW] bench switches
p.addParameter('noIntan',false,@islogical);
p.addParameter('debug',false,@islogical);
p.addParameter('skipSync',false,@islogical);
p.addParameter('confirm',true,@islogical);

p.parse(varargin{:});
args = p.Results;

if isempty(args.chn)
    error('expSinglePulseISI:chn','Give the two Intan channel names, e.g. ''chn'',{''A-007'',''A-012''}.');
end
if args.trialDur > 550
    warning('expSinglePulseISI:trialDur', ...
        'trialDur %g ms is over ~550 ms: the end of the trial may fall in the Intan hole.',args.trialDur);
end
if args.pulseTime + args.trigDuration > args.trialDur
    error('expSinglePulseISI:timing','pulseTime + trigDuration (%g ms) must fit inside trialDur (%g ms).', ...
        args.pulseTime+args.trigDuration,args.trialDur);
end

% [NEW] amplitudes and conditions
% 'amp' is a LIST of levels (uA). Electrode A uses amp; electrode B uses ampB,
% or the same levels when ampB is not given.
ampA = args.amp(:)';
if isempty(args.ampB), ampB = ampA; else, ampB = args.ampB(:)'; end
if numel(ampB) ~= numel(ampA)
    error('expSinglePulseISI:amp','ampB needs the same number of levels as amp (%d).',numel(ampA));
end
[cond,isVisual] = buildPairConds(ampA,ampB,args.isi,args.order,args.withVisual);
nCond = numel(cond);
nTrials = nCond * args.nReps;

% [NEW] charge limit, checked BEFORE anything is opened
if args.fpd ~= args.spd
    warning('expSinglePulseISI:unbalanced','fpd (%g us) ~= spd (%g us) at equal amplitude: the pulse is NOT charge balanced.',args.fpd,args.spd);
end
% same amplitude in both phases, so the wider phase carries the most charge
checkChargeLimit([ampA ampB],max(args.fpd,args.spd), ...
    'siteArea',args.siteArea,'capacity',args.capacity,'kMax',args.kMax,'allowOver',args.allowOver);

fprintf('\nexpSinglePulseISI: %d conditions x %d repeats = %d trials (about %.0f min at ~2 s per trial)\n', ...
    nCond,args.nReps,nTrials,nTrials*2/60);
fprintf('  channels A = %s, B = %s; amplitude A [%s] uA, B [%s] uA; phase widths %g / %g us, interphase %g us\n', ...
    args.chn{1},args.chn{2},num2str(ampA),num2str(ampB),args.fpd,args.spd,args.ipi);
for k = 1:nCond
    fprintf('  %2d  %-16s amp [%g %g]  ptd [%g %g] us\n',k,cond(k).label, ...
        cond(k).amp(1),cond(k).amp(2),cond(k).ptd(1),cond(k).ptd(2));
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
    warning('expSinglePulseISI:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
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
e.testMode = args.noIntan;
e.headerPause = args.headerPause;
e.setConditions(cond);

% frame-grid pulse time. estim is "on" from e.on for trigDuration ms; its
% onsetFunction fires the trigger.
ivl = 1000/c.screen.frameRate;
nf  = round(args.pulseTime/ivl);
e.on       = nf*ivl;
e.duration = args.trigDuration;
fprintf('expSinglePulseISI timing: pulse requested at %g ms, on frame %d, delivered at %.2f ms after trial start\n', ...
    args.pulseTime,nf+1,e.on);

% [NEW] the visual arm (optional). A uniform spot made with the same gabor
% class as exp1A: spatial frequency 0, so it is a disc of uniform luminance.
% PHASE (checked in the shader source, GLSLShaders/gabor.frag.txt): luminance =
% background x (1 + contrast x sin(2 pi f x + phase)). The spot has f = 0, so its
% pattern is sin(phase): bright = 90, dark = 270. Phase 0 or 180 would be
% INVISIBLE (sin = 0). Contrast 0 hides the spot, so
% the electrical trials show nothing (the photodiode patch flashes
% only on the visual-only condition).
if args.withVisual
    g = stimuli.gabor(c,'flash');
    g.X = args.rfX;
    g.Y = args.rfY;
    g.mask   = 'CIRCLE';
    g.sigma  = args.diameter/2;      % CIRCLE: sigma is the radius
    g.width  = ceil(args.diameter) + 1;
    g.height = ceil(args.diameter) + 1;
    g.frequency   = 0;
    g.orientation = 0;
    g.phaseSpeed  = 0;
    if strcmpi(args.flashPolarity,'bright'), g.phase = 90; else, g.phase = 270; end   % sin(phase)
    g.contrast = NaN;                % the design writes this every trial
    g.on       = e.on;               % same frame as the pulse
    g.duration = args.flashDur;
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
    error('expSinglePulseISI: DaqDeviceIndex found no DAQ, so no trigger can be sent.');
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

% [EXP1A] diode flasher: c.diode is disabled; with the visual arm the patch is
% attached to the spot so it marks the visual onset.
addDiodeFlasher(c.diode,'enabled',false);
if args.withVisual
    addDiodeFlasher(g,'size',0.075,'location','nw');
end

c.trialDuration = args.trialDur;      % ms
c.iti = args.iti*1e3;                 % s -> ms

% -------------------------------------------------------------------------
% [NEW] DESIGN: one factor, the condition index
% -------------------------------------------------------------------------
d = neurostim.design('myDesign');
d.fac1.estim.stimCond = 1:nCond;
if args.withVisual
    % same factor, listed together: contrast is flashContrast on the
    % visual-only condition and 0 (hidden) on all the others
    d.fac1.flash.contrast = args.flashContrast * double(isVisual(:)');
    % The spot, and so its photodiode flash, exists ONLY on the visual-only
    % condition: on = Inf means 'never on'.
    visOn = inf(1,nCond);
    visOn(isVisual(:)') = e.on;
    d.fac1.flash.on = visOn;
end

blk = block('myBlock',d);
blk.nrRepeats = args.nReps;

c.subject  = subject;
c.paradigm = 'expSinglePulseISI';
c.addPropsToInform('estim.stimCond');

% [NEW] keep the condition table with the data
stamp = datestr(now,'yyyymmdd_HHMMSS');
condFile = sprintf('%s_expSinglePulseISI_conds_%s.mat',subject,stamp);
chnList = e.chnList;                  %#ok<NASGU>
save(condFile,'cond','isVisual','chnList','args');
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

function [cond,isVisual] = buildPairConds(ampA,ampB,isi_ms,order,withVisual)
% BUILDPAIRCONDS  Single-pulse conditions for one electrode pair (A and B).
%
%   [cond,isVisual] = buildPairConds(ampA,ampB,isi_ms,order,withVisual)
%
% ampA, ampB  uA for channel A and B: vectors of equal length, one entry per
%             amplitude LEVEL. Every level gets its own full set of conditions
%             (simultaneous, sequential, A only, B only); blank and visual
%             only are added once.
% isi_ms     positive ISIs in ms between the two pulses. Simultaneous
%             (ISI 0) is ALWAYS added, so do not list 0 here.
% order       'both' | 'AB' | 'BA'   (which electrode fires first)
% withVisual  true adds a visual-only condition
%
% Returns the condition table for estimIntanMulti (struct array with fields
% amp, ptd, nsp, label; one entry per channel in chnList = {A,B}) and a
% logical vector marking the visual-only condition.
%
% Conditions, for ISI g:
%   simultaneous      amp [A B], ptd [0 0]
%   A -> B, ISI g     amp [A B], ptd [0 g]    (A at 0, B delayed by g)
%   B -> A, ISI g     amp [A B], ptd [g 0]
%   A only            amp [A 0] (B off), ptd [0 0]   recovery condition for A -> B
%   B only            amp [0 B] (A off), ptd [0 0]   recovery condition for B -> A
%   blank             nothing delivered (all channels off)
%   visual only       nothing delivered electrically, visual on (if withVisual)
% PTD is in MICROSECONDS here (what estimIntanMulti writes into StimParams).

isi_ms = unique(isi_ms(:)');
isi_ms = isi_ms(isi_ms > 0);
order = lower(order);
assert(ismember(order,{'both','ab','ba'}), ...
    'buildPairConds: order must be ''both'', ''AB'' or ''BA''.');

cond = struct('amp',{},'ptd',{},'nsp',{},'label',{});
isVisual = false(0,1);

    function add(amp,ptd,label,vis)
        k = numel(cond) + 1;
        cond(k).amp = amp;
        cond(k).ptd = ptd;
        cond(k).nsp = [1 1];
        cond(k).label = label;
        isVisual(k,1) = vis;
    end

multi = numel(ampA) > 1;          % add the amplitude to the label only if >1 level
for L = 1:numel(ampA)
    aA = ampA(L);  aB = ampB(L);
    if multi, tag = sprintf(' @%g/%g uA',aA,aB); else, tag = ''; end

    add([aA aB],[0 0],['simultaneous' tag],false);

    for g = isi_ms
        us = round(g*1000);
        if any(strcmp(order,{'both','ab'}))
            add([aA aB],[0 us],sprintf('A->B %g ms%s',g,tag),false);
        end
        if any(strcmp(order,{'both','ba'}))
            add([aA aB],[us 0],sprintf('B->A %g ms%s',g,tag),false);
        end
    end

    % Recovery (single-electrode) conditions at THIS level. A->B needs A only;
    % B->A needs B only.
    if any(strcmp(order,{'both','ab'}))
        add([aA 0],[0 0],['A only' tag],false);
    end
    if any(strcmp(order,{'both','ba'}))
        add([0 aB],[0 0],['B only' tag],false);
    end
end

add([0 0],[0 0],'blank',false);

if withVisual
    add([0 0],[0 0],'visual only',true);
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
