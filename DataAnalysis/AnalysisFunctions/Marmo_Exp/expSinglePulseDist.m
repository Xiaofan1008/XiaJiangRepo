function c = expSinglePulseDist(subject,varargin)
% EXPSINGLEPULSEDIST  Single-pulse mechanism check: simultaneous vs sequential,
% DISTANCE sweep at one ISI (default 5 ms, 75 uA). Electrode A is fixed; each
% partner electrode B1, B2, ... sits at a different distance from A. Both
% orders, with the single-electrode recovery conditions for every pair.
%
% STATUS: written but NOT yet run in MATLAB or on the rig. Test with
% 'noIntan',true first (see BENCH SWITCHES), then with Intan.
%
% Built on exp1A.m (rig setup, plugin order, trigger timing, cleanup) and the
% multi-electrode stimulation plugin estimIntanMulti.m. Every value below is
% an input with a default; change anything when you call it:
%
%   c = expSinglePulseDist('CJ286','chnA','A-007','chnB',{'A-012','A-020','A-030'}, ...
%                          'dist',[100 300 600]);
%   c = expSinglePulseDist('CJ286','chnA','A-007','chnB',{'A-012','A-020'}, ...
%                          'dist',[100 300],'amp',60,'isi',5,'nReps',20);
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
% and 'chnA' and 'chnB' must be given: a default channel could stimulate the
% wrong site.
%
%   label          meaning                                   default
%   -------------  ----------------------------------------  -------------
%  ELECTRICAL
%   'chnA'         Intan name of the fixed electrode A       REQUIRED, e.g. 'A-007'
%   'chnB'         Intan names of the partner electrodes,    REQUIRED, e.g. {'A-012','A-020'}
%                  one per distance (cell array of strings)
%   'dist'         distance A to each B in um (labels only;  [] (not recorded)
%                  not used to pick channels). Same length as chnB
%   'amp'          amplitude LEVELS in uA, for A (and B)     75
%                  e.g. [50 75] runs every condition at each level
%   'ampB'         levels for B if different from A          [] (= same as amp)
%   'isi'         ISIs A->B in ms (simultaneous always in)  5
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
% Examples:  order 'AB' only:   expSinglePulseDist('CJ286','chnA','A-007','chnB',{'A-012','A-020'},'order','AB')
%            different amplitude on B: expSinglePulseDist('CJ286','chnA','A-007','chnB',{'A-012'},'amp',40,'ampB',60)
% =========================================================================

% =========================================================================
% CONDITIONS (built by the local function buildDistConds below). A is chnA and
% B1..Bn are the chnB entries. The stimulator holds all n+1 channels; in every
% trial only the ones in the condition are on, the rest are switched off.
% =========================================================================
% For each partner Bj (one pair A+Bj per distance):
%   simultaneous     A and Bj together
%   A -> Bj, ISI g   A first, Bj g ms later         (one per ISI)
%   Bj -> A, ISI g   Bj first, A g ms later         (one per ISI)
%   Bj only          single-electrode recovery condition, no delay
% Shared by all pairs (listed once):
%   A only           recovery condition for A -> Bj (identical for every j)
%   blank            nothing delivered (trigger still fires, all channels off)
%   visual only      only if 'withVisual' is true
% With 3 partners, 1 ISI, both orders, one amplitude:
%   3 x (1 + 2 + 1) + 1 (A only) + 1 (blank) = 14 conditions.
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
p.addParameter('chnA','',@(x) ischar(x));
p.addParameter('chnB',{},@(x) iscellstr(x) && ~isempty(x)); %#ok<ISCLSTR>
p.addParameter('dist',[],@(x) isnumeric(x) && all(x>=0));
p.addParameter('amp',75,@(x) isnumeric(x) && ~isempty(x) && all(x>0));
p.addParameter('ampB',[],@(x) isnumeric(x) && all(x>0));
p.addParameter('isi',5,@(x) isnumeric(x) && all(x>=0));
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

if isempty(args.chnA) || isempty(args.chnB)
    error('expSinglePulseDist:chn','Give ''chnA'' (one name) and ''chnB'' (a cell array of names), e.g. ''chnA'',''A-007'',''chnB'',{''A-012'',''A-020''}.');
end
args.chnB = args.chnB(:)';
if any(strcmp(args.chnA,args.chnB)) || numel(unique(args.chnB)) ~= numel(args.chnB)
    error('expSinglePulseDist:chn','chnA and the chnB entries must all be different channels.');
end
if ~isempty(args.dist) && numel(args.dist) ~= numel(args.chnB)
    error('expSinglePulseDist:dist','dist needs one value per chnB entry (%d).',numel(args.chnB));
end
chnAll = [{args.chnA} args.chnB];      % channel order used by every condition
if args.trialDur > 550
    warning('expSinglePulseDist:trialDur', ...
        'trialDur %g ms is over ~550 ms: the end of the trial may fall in the Intan hole.',args.trialDur);
end
if args.pulseTime + args.trigDuration > args.trialDur
    error('expSinglePulseDist:timing','pulseTime + trigDuration (%g ms) must fit inside trialDur (%g ms).', ...
        args.pulseTime+args.trigDuration,args.trialDur);
end

% [NEW] amplitudes and conditions
% 'amp' is a LIST of levels (uA). Electrode A uses amp; electrode B uses ampB,
% or the same levels when ampB is not given.
ampA = args.amp(:)';
if isempty(args.ampB), ampB = ampA; else, ampB = args.ampB(:)'; end
if numel(ampB) ~= numel(ampA)
    error('expSinglePulseDist:amp','ampB needs the same number of levels as amp (%d).',numel(ampA));
end
[cond,isVisual] = buildDistConds(numel(args.chnB),ampA,ampB,args.isi,args.order,args.withVisual,args.chnB,args.dist);
nCond = numel(cond);
nTrials = nCond * args.nReps;

% [NEW] charge limit, checked BEFORE anything is opened
if args.fpd ~= args.spd
    warning('expSinglePulseDist:unbalanced','fpd (%g us) ~= spd (%g us) at equal amplitude: the pulse is NOT charge balanced.',args.fpd,args.spd);
end
% same amplitude in both phases, so the wider phase carries the most charge
checkChargeLimit([ampA ampB],max(args.fpd,args.spd), ...
    'siteArea',args.siteArea,'capacity',args.capacity,'kMax',args.kMax,'allowOver',args.allowOver);

fprintf('\nexpSinglePulseDist: %d conditions x %d repeats = %d trials (about %.0f min at ~2 s per trial)\n', ...
    nCond,args.nReps,nTrials,nTrials*2/60);
fprintf('  channels (order in the amp/ptd vectors): %s\n',strjoin(chnAll,', '));
fprintf('  amplitude A [%s] uA, B [%s] uA; phase widths %g / %g us, interphase %g us\n', ...
    num2str(ampA),num2str(ampB),args.fpd,args.spd,args.ipi);
for k = 1:nCond
    fprintf('  %2d  %-28s amp %s  ptd %s us\n',k,cond(k).label, ...
        mat2str(cond(k).amp),mat2str(cond(k).ptd));
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
    warning('expSinglePulseDist:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
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
e.chnList  = chnAll;
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
fprintf('expSinglePulseDist timing: pulse requested at %g ms, on frame %d, delivered at %.2f ms after trial start\n', ...
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
    error('expSinglePulseDist: DaqDeviceIndex found no DAQ, so no trigger can be sent.');
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
c.paradigm = 'expSinglePulseDist';
c.addPropsToInform('estim.stimCond');

% [NEW] keep the condition table with the data
stamp = datestr(now,'yyyymmdd_HHMMSS');
condFile = sprintf('%s_expSinglePulseDist_conds_%s.mat',subject,stamp);
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

function [cond,isVisual] = buildDistConds(nB,ampA,ampB,isi_ms,order,withVisual,nameB,dist) %#ok<INUSL>
% BUILDDISTCONDS  Single-pulse conditions for one fixed electrode A and nB
% partner electrodes (the distance sweep).
%
%   [cond,isVisual] = buildDistConds(nB,ampA,ampB,isi_ms,order,withVisual,nameB,dist)
%
% Channel order in every amp/ptd vector: [A B1 B2 ... BnB], the same order as
% chnList passed to estimIntanMulti. A channel that is not in the condition has
% amp 0 (the plugin switches it off).
%
% ampA, ampB  vectors of equal length, one entry per amplitude LEVEL
% isi_ms     positive ISIs in ms. Simultaneous (ISI 0) is ALWAYS added
% order       'both' | 'AB' | 'BA'  (A first / Bj first)
% nameB, dist partner names and distances in um (labels only; dist may be [])
%
% Duplicate conditions (the A-only recovery condition is the same for every
% partner) are listed once. PTD is in MICROSECONDS.

isi_ms = unique(isi_ms(:)');
isi_ms = isi_ms(isi_ms > 0);
order = lower(order);
assert(ismember(order,{'both','ab','ba'}), ...
    'buildDistConds: order must be ''both'', ''AB'' or ''BA''.');
nCh = nB + 1;

cond = struct('amp',{},'ptd',{},'nsp',{},'label',{});
isVisual = false(0,1);
keys = {};

    function add(amp,ptd,label,vis)
        key = [mat2str(amp) mat2str(ptd)];
        if any(strcmp(keys,key)), return; end    % identical condition, keep one
        keys{end+1} = key; %#ok<AGROW>
        k = numel(cond) + 1;
        cond(k).amp = amp;
        cond(k).ptd = ptd;
        cond(k).nsp = ones(1,nCh);
        cond(k).label = label;
        isVisual(k,1) = vis;
    end

multi = numel(ampA) > 1;
for L = 1:numel(ampA)
    aA = ampA(L);  aB = ampB(L);
    if multi, tag = sprintf(' @%g/%g uA',aA,aB); else, tag = ''; end

    for j = 1:nB
        col = j + 1;                              % column of Bj
        if isempty(dist), nm = sprintf('B%d',j); else, nm = sprintf('B%d(%gum)',j,dist(j)); end
        amp0 = zeros(1,nCh);  amp0([1 col]) = [aA aB];

        add(amp0,zeros(1,nCh),[nm ' simult' tag],false);
        for g = isi_ms
            us = round(g*1000);
            if any(strcmp(order,{'both','ab'}))
                p = zeros(1,nCh);  p(col) = us;   % A at 0, Bj delayed
                add(amp0,p,sprintf('A->%s %g ms%s',nm,g,tag),false);
            end
            if any(strcmp(order,{'both','ba'}))
                p = zeros(1,nCh);  p(1) = us;     % Bj at 0, A delayed
                add(amp0,p,sprintf('%s->A %g ms%s',nm,g,tag),false);
            end
        end
        % Recovery: Bj only recovers Bj -> A
        if any(strcmp(order,{'both','ba'}))
            a1 = zeros(1,nCh);  a1(col) = aB;
            add(a1,zeros(1,nCh),[nm ' only' tag],false);
        end
    end
    % Recovery: A only recovers A -> Bj (same for all j)
    if any(strcmp(order,{'both','ab'}))
        a1 = zeros(1,nCh);  a1(1) = aA;
        add(a1,zeros(1,nCh),['A only' tag],false);
    end
end

add(zeros(1,nCh),zeros(1,nCh),'blank',false);
if withVisual
    add(zeros(1,nCh),zeros(1,nCh),'visual only',true);
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
