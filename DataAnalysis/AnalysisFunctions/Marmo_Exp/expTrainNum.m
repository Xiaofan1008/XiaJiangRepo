function c = expTrainNum(subject,varargin)
% EXPTRAINNUM  ELECTRODE-NUMBER experiment with pulse trains: the same TOTAL
% current is delivered through 1, 2 or 3 electrodes, simultaneously or
% sequentially, interleaved with VISUAL-ONLY trials, in one session.
% Electrical and visual stimulation never occur in the same trial.
%
% STATUS: written but NOT yet run in MATLAB or on the rig. Test with
% 'noIntan',true first (see BENCH SWITCHES), then with Intan.
%
% =========================================================================
% >>> DESIGN REMINDER: KEEP THE OUTER SPAN FIXED <<<
% =========================================================================
% Adding electrodes along a shank also widens the area they cover, so number
% and spread are confounded. This file avoids that:
%     E1 ---------- E2 ---------- E3        (E2 is the CENTRE site)
%   1 electrode   : E2
%   2 electrodes  : E1 + E3     (the OUTER pair)
%   3 electrodes  : E1 + E2 + E3
% The 2- and 3-electrode conditions share the SAME outer span (E1 to E3), so
% only the NUMBER changes. Going from 1 to 2 electrodes does change the extent;
% say so when you report it.
% WHEN YOU CHOOSE THE CHANNELS:
%   * Choose E1, E2, E3 from the sites that SCREENED well (a poor site wastes
%     its share of the current).
%   * Keep the spacing the same in every session of this experiment and
%     write it down ('span', um). Spacing itself is tested on its own in
%     expTrainDist.m. This is the ONE-FACTOR-AT-A-TIME plan (as in Kunigk et
%     al. 2022): distance in expTrainDist, number here with the span fixed.
%   * 'chn' = {E1, E2, E3} in that order: outer, CENTRE, outer.
%
% Built on exp1A.m (rig setup, plugin order, trigger timing, cleanup), the
% multi-electrode plugin estimIntanMulti.m and the event-level recovery logic
% of design_exp_v7_VariableReps.m.
%
%   c = expTrainNum('CJ286','chn',{'A-007','A-012','A-017'},'span',200);
%   c = expTrainNum('CJ286','chn',{'A-007','A-012','A-017'},'amp',[60 90], ...
%                   'order','both','match','both');
%
% The first input is the SUBJECT (animal ID, a label only, as in exp1A).
%
% STANDALONE: this one file holds the experiment, the condition builder and the
% charge-limit check. It needs only estimIntanMulti.m plus what the
% stimulation PC already has (neurostim, marmolab, the Intan helper functions,
% Psychtoolbox, the DAQ functions).
%
% =========================================================================
% INPUT LABELS  (name-value pairs, e.g. 'amp',60). Every label has a default.
% Only the SUBJECT and 'chn' must be given: a default channel could stimulate
% the wrong site.
%
%   label          meaning                                   default
%   -------------  ----------------------------------------  -------------
%  ELECTRICAL
%   'chn'          Intan names {E1,E2,E3}: outer, CENTRE,    REQUIRED, e.g. {'A-007','A-012','A-017'}
%                  outer (see the design reminder)
%   'span'         distance E1 to E3 in um (stored only)     [] (not recorded)
%   'amp'          TOTAL current levels in uA                75
%                  (every electrical condition runs at each level)
%   'match'        'total': each electrode gets total/k      'total'
%                  'each' : every electrode gets the full amount (unmatched
%                           control); 'both': run both sets
%   'freq'         pulse frequency within a train, Hz        100
%   'nPulses'      pulses per electrode in the FULL train    5
%   'isi'          stagger in ms between neighbouring        [] = period/k
%                  electrodes in a sequential pattern.       (3 el.: 3.33 ms,
%                  [] = automatic, so pulses alternate        2 el.: 5 ms)
%                  evenly. A number sets it for all counts.
%   'order'        sequential order: 'fwd' (E1 first) |      'fwd'
%                  'rev' (last electrode first) | 'both'
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
% CONDITIONS (built by the local function buildTrainNumConds below)
% =========================================================================
% For every TOTAL-current level, with N = nPulses, and for k = 1, 2, 3 electrodes
% (sets E2 | E1+E3 | E1+E2+E3, each electrode at total/k when match = 'total'):
%   k = 1            one train on E2: N levels (the first j pulses, j = 1..N)
%   k >= 2, sim      all k electrodes together: N levels
%   k >= 2, fwd      staggered, E1 first: levels add one pulse at a time in
%                    TIME ORDER, about k*N levels
%   k >= 2, rev      the same with the last electrode first (if order allows)
% Only the LAST level of each pattern (all N pulses on every electrode) is
% the fully delivered condition you report; the earlier levels exist so the
% response hidden under each pulse's stimulation artifact can be recovered
% (design_exp_v7 trainMode 2). All levels use the same number of repeats.
% Shared by all levels:
%   blank            nothing delivered
%   visual           one condition per (visType x polarity x contrast), only
%                    if 'withVisual' is true. No electrical stimulation.
% Counts (N = 5, one total-current level, match 'total'):
%   order 'fwd':  5 + (5+10) + (5+15) = 40 electrical + 1 blank + 10 visual (bright and dark spot)
%                 = 51 conditions, 20 repeats, about 34 min
%   order 'both': 40 + 10 + 15 = 65 electrical, 76 conditions, about 51 min
%   match 'both' adds the 'each' set (+34 for 'fwd'); identical conditions
%                 are listed once.
% Each extra contrast adds 2 conditions with both polarities (about 80 s).
% THE SESSION LENGTH IS PRINTED BEFORE YOU CONFIRM: cut nPulses, the number of
% levels, nReps or order if it is too long.
%
% =========================================================================
% TRIAL STRUCTURE
% =========================================================================
%   0 ms        trial starts (first frame)
%   pulseTime   trigger fires; the first electrode pulses at once
%               and the visual patch appears on the same frame   (default 50 ms)
%   + stagger   the next electrodes start via their post-trigger delays
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
p.addParameter('chn',{},@(x) iscellstr(x) && numel(x)==3); %#ok<ISCLSTR>
p.addParameter('span',[],@(x) isnumeric(x) && (isempty(x) || (isscalar(x) && x>=0)));
p.addParameter('amp',75,@(x) isnumeric(x) && ~isempty(x) && all(x>0));
p.addParameter('match','total',@(x) ismember(lower(x),{'total','each','both'}));
p.addParameter('freq',100,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('nPulses',5,@(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('isi',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('order','fwd',@(x) ismember(lower(x),{'fwd','rev','both'}));
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
    error('expTrainNum:chn','Give ''chn'' = {E1,E2,E3}: outer, CENTRE, outer, e.g. ''chn'',{''A-007'',''A-012'',''A-017''}.');
end
if numel(unique(args.chn)) ~= 3
    error('expTrainNum:chn','The three channels in chn must all be different.');
end
chnAll = args.chn(:)';                 % channel order used by every condition
if args.trialDur > 550
    warning('expTrainNum:trialDur', ...
        'trialDur %g ms is over ~550 ms: the end of the trial may fall in the Intan hole.',args.trialDur);
end

% [NEW] train timing, checked before anything is opened
period_ms = 1000/args.freq;
period_us = round(1e6/args.freq);
if isempty(args.isi)
    stagger_ms = period_ms./[1 2 3];                  % per electrode count k
else
    stagger_ms = args.isi*[1 1 1];
    for kq = 2:3
        if args.nPulses > 1 && mod(round((kq-1)*args.isi*1000),period_us) == 0
            warning('expTrainNum:isi','isi %g ms puts some pulses of different electrodes at the same time (a whole number of periods, %g ms).',args.isi,period_ms);
        end
    end
end
maxDelay_ms = 2*stagger_ms(3);                        % the last electrode of the 3-electrode set
trainEnd = args.pulseTime + maxDelay_ms + (args.nPulses-1)*period_ms;   % last pulse starts here
if args.withVisual && args.pulseTime + args.visDur > args.trialDur
    error('expTrainNum:timing','The visual patch (%g ms from %g ms) does not fit in trialDur %g ms.', ...
        args.visDur,args.pulseTime,args.trialDur);
end
if args.withVisual && args.rfX == 0 && args.rfY == 0
    warning('expTrainNum:rf','rfX and rfY are both 0: this is a PLACEHOLDER. Set them from RF mapping.');
end
if args.pulseTime + args.trigDuration > args.trialDur || trainEnd + 2 > args.trialDur
    error('expTrainNum:timing', ...
        'The train does not fit: last pulse at %.1f ms, trigger high until %.1f ms, trialDur %g ms.', ...
        trainEnd,args.pulseTime+args.trigDuration,args.trialDur);
end

% [NEW] amplitudes and conditions
ampTot = args.amp(:)';
vis = struct('on',args.withVisual,'type',lower(args.visType), ...
             'polarity',lower(args.polarity),'contrast',args.contrast(:)','sf',args.sf);
[cond,info] = buildTrainNumConds(ampTot,args.nPulses,period_us,args.isi,lower(args.order),lower(args.match),vis);
nCond = numel(cond);
nTrials = nCond * args.nReps;

% [NEW] charge limit, checked BEFORE anything is opened
if args.fpd ~= args.spd
    warning('expTrainNum:unbalanced','fpd (%g us) ~= spd (%g us) at equal amplitude: the pulse is NOT charge balanced.',args.fpd,args.spd);
end
% the largest single-electrode amplitude equals the total current
checkChargeLimit(ampTot,max(args.fpd,args.spd), ...
    'siteArea',args.siteArea,'capacity',args.capacity,'kMax',args.kMax,'allowOver',args.allowOver);

nElec = sum(~strcmp({info.kind},'visual') & ~strcmp({info.kind},'blank'));
fprintf('\nexpTrainNum: %d conditions x %d repeats = %d trials (about %.0f min at ~2 s per trial)\n', ...
    nCond,args.nReps,nTrials,nTrials*2/60);
fprintf('  %d electrical + %d visual + %d blank conditions\n', ...
    nElec,sum(strcmp({info.kind},'visual')),sum(strcmp({info.kind},'blank')));
fprintf('  channels E1,E2,E3 (outer, CENTRE, outer; order in the amp/ptd/nsp vectors): %s\n',strjoin(chnAll,', '));
fprintf('  total current [%s] uA, match %s; per electrode (match total): 1 el. = total, 2 el. = total/2, 3 el. = total/3\n',num2str(ampTot),lower(args.match));
fprintf('  train: %d pulses at %g Hz (period %.3g ms), stagger 2 el. %.3g ms, 3 el. %.3g ms, last pulse at %.1f ms; phases %g / %g us, interphase %g us\n', ...
    args.nPulses,args.freq,period_ms,stagger_ms(2),stagger_ms(3),trainEnd,args.fpd,args.spd,args.ipi);
fprintf('  REMINDER: E1 and E3 are the OUTER pair (fixed span%s), E2 is the centre.\n', ...
    ternary(isempty(args.span),'',sprintf(' %g um',args.span)));
% [NEW] total charge per phase of the FULL train, per electrode count, so the
% matching can be read off before you confirm (nC = uA * us / 1000)
pwMax = max(args.fpd,args.spd);
for t = ampTot
    fprintf('  total %g uA, full train charge per phase (all electrodes, %d pulses):', t, args.nPulses);
    for kq = 1:3
        if strcmp(lower(args.match),'each'), q = kq*t; else, q = t; end
        fprintf('  %d el. %.1f nC', kq, q*pwMax*args.nPulses/1000);
    end
    if strcmp(lower(args.match),'both')
        fprintf('  (set ''each'' injects %d / %d / %d times as much)', 1, 2, 3);
    end
    fprintf('\n');
end
% the stimulator may not take every fractional amplitude: warn if a per-electrode
% amplitude is not a whole number of uA (NOT verified on the rig)
perEl = [];
for t = ampTot
    perEl = [perEl, t, t/2, t/3]; %#ok<AGROW>
end
if any(abs(perEl - round(perEl)) > 1e-9)
    warning('expTrainNum:amp','Some per-electrode amplitudes are not whole uA (%s). Check on the rig what the stimulator delivers, or pick a total divisible by 6 (60, 72, 90).',num2str(unique(perEl(abs(perEl-round(perEl))>1e-9))));
end
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
    warning('expTrainNum:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
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
fprintf('expTrainNum timing: pulse requested at %g ms, on frame %d, delivered at %.2f ms after trial start\n', ...
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
    error('expTrainNum: DaqDeviceIndex found no DAQ, so no trigger can be sent.');
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
c.paradigm = 'expTrainNum';
c.addPropsToInform('estim.stimCond');

% [NEW] keep the condition table with the data
stamp = datestr(now,'yyyymmdd_HHMMSS');
condFile = sprintf('%s_expTrainNum_conds_%s.mat',subject,stamp);
chnList = e.chnList;                  %#ok<NASGU>
save(condFile,'cond','info','chnList','args','period_us');
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

function [cond,info] = buildTrainNumConds(ampTot,N,period_us,isiArg_ms,order,match,vis)
% BUILDTRAINNUMCONDS  Electrode-number conditions with event-level recovery
% levels (design_exp_v7 trainMode 2), plus blank and visual-only conditions.
%
%   [cond,info] = buildTrainNumConds(ampTot,N,period_us,isiArg_ms,order,match,vis)
%
% Three channels [E1 E2 E3]; E1 and E3 are the OUTER pair, E2 the centre.
% The electrode sets are fixed so that the OUTER SPAN stays the same:
%   1 electrode   E2
%   2 electrodes  E1 + E3       (outer pair)
%   3 electrodes  E1 + E2 + E3
%
% ampTot      vector of TOTAL current levels (uA)
% N           pulses per electrode in the full train
% period_us   time between pulses of one train (us)
% isiArg_ms   [] = automatic stagger of period/k for k electrodes (pulses
%             alternate evenly), or a fixed stagger in ms between neighbours
% order       'fwd' | 'rev' | 'both'   (E1 first / last electrode first)
% match       'total' (each electrode gets total/k) | 'each' (each gets the
%             full amount) | 'both'
% vis         struct: on, type, polarity, contrast, sf
%
% RECOVERY LEVELS. For each pattern, all pulse times of all electrodes in the
% set are sorted; level i holds the pulses delivered up to the i-th distinct
% time. Pulses at exactly the same time share a level. Identical conditions
% are listed once. A channel with nsp 0 is switched off (amp 0, ptd 0).
%
% info  per condition: kind ('single'|'sim'|'fwd'|'rev'|'blank'|'visual'),
%       nElec, match, ampEach, level, nLevels, isFinal, and the visual
%       settings contrast, phase (spot 90 bright / 270 dark; gabor 0 / 180), sf. PTD in us.

order = lower(order);
assert(ismember(order,{'fwd','rev','both'}), ...
    'buildTrainNumConds: order must be ''fwd'', ''rev'' or ''both''.');
match = lower(match);
assert(ismember(match,{'total','each','both'}), ...
    'buildTrainNumConds: match must be ''total'', ''each'' or ''both''.');
if strcmp(match,'both'), matches = {'total','each'}; else, matches = {match}; end
sets = {2,[1 3],[1 2 3]};
nCh = 3;

cond = struct('amp',{},'ptd',{},'nsp',{},'label',{});
info = struct('kind',{},'nElec',{},'match',{},'ampEach',{},'level',{},'nLevels',{}, ...
              'isFinal',{},'contrast',{},'phase',{},'sf',{});
keys = {};

    function add(amp,ptd,nsp,label,kind,k,mt,aE,lev,nLev,vc,vp,vsf)
        off = nsp < 1;
        amp(off) = 0;  ptd(off) = 0;  nsp(off) = 0;       % channel off
        key = [mat2str(amp) mat2str(ptd) mat2str(nsp) sprintf('|%g|%g|%g',vc,vp,vsf)];
        if any(strcmp(keys,key)), return; end              % identical condition, keep one
        keys{end+1} = key; %#ok<AGROW>
        n = numel(cond) + 1;
        cond(n).amp = amp;  cond(n).ptd = ptd;  cond(n).nsp = nsp;  cond(n).label = label;
        info(n).kind = kind;  info(n).nElec = k;  info(n).match = mt;  info(n).ampEach = aE;
        info(n).level = lev;  info(n).nLevels = nLev;  info(n).isFinal = (lev == nLev);
        info(n).contrast = vc;  info(n).phase = vp;  info(n).sf = vsf;
    end

multi = numel(ampTot) > 1;
for L = 1:numel(ampTot)
    tot = ampTot(L);
    if multi, tag = sprintf(' @%guA',tot); else, tag = ''; end
    for im = 1:numel(matches)
        mt = matches{im};
        if strcmp(mt,'each'), mtag = ' each'; else, mtag = ''; end
        for ki = 1:3
            S = sets{ki};  k = numel(S);
            if strcmp(mt,'total'), aE = tot/k; else, aE = tot; end
            if isempty(isiArg_ms), stag_us = round(period_us/k); else, stag_us = round(isiArg_ms*1000); end

            pats = {};                                    % {kind, delays_us}
            if k == 1
                pats(end+1,:) = {'single',0}; %#ok<AGROW>
            else
                pats(end+1,:) = {'sim',zeros(1,k)}; %#ok<AGROW>
                if any(strcmp(order,{'fwd','both'})), pats(end+1,:) = {'fwd',(0:k-1)*stag_us}; end %#ok<AGROW>
                if any(strcmp(order,{'rev','both'})), pats(end+1,:) = {'rev',(k-1:-1:0)*stag_us}; end %#ok<AGROW>
            end
            for ip = 1:size(pats,1)
                kind = pats{ip,1};  dl = pats{ip,2};
                lv = seqLevelsN(N,period_us,dl);
                for i = 1:size(lv,1)
                    amp = zeros(1,nCh);  amp(S) = aE;
                    nsp = zeros(1,nCh);  nsp(S) = lv(i,:);
                    ptd = zeros(1,nCh);  ptd(S) = dl;
                    lab = sprintf('%de %s %d/%d [%s]%s%s',k,kind,i,size(lv,1),num2str(lv(i,:)),mtag,tag);
                    add(amp,ptd,nsp,lab,kind,k,mt,aE,i,size(lv,1),0,0,0);
                end
            end
        end
    end
end

z = zeros(1,nCh);
add(z,z,ones(1,nCh),'blank','blank',0,'',0,1,1,0,0,0);
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
                add(z,z,ones(1,nCh),lab,'visual',0,'',0,1,1,vis.contrast(ic),ph,sf);
            end
        end
    end
end
end

function lv = seqLevelsN(N,period_us,delays)
% SEQLEVELSN  Cumulative pulse counts [n1 n2 ...] after each distinct pulse
% time. Electrode j fires at delays(j) + (0..N-1)*period (all in us).
nE = numel(delays);
tt = cell(1,nE);
for j = 1:nE
    tt{j} = delays(j) + (0:N-1)*period_us;
end
t  = unique([tt{:}]);
lv = zeros(numel(t),nE);
for i = 1:numel(t)
    for j = 1:nE
        lv(i,j) = sum(tt{j} <= t(i));
    end
end
end

function out = ternary(c,a,b)
if c, out = a; else, out = b; end
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
