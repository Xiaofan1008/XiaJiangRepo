function c = expTrainScreen(subject,varargin)
% EXPTRAINSCREEN  Electrode SCREENING: stimulate each candidate channel ALONE
% with one pulse train, plus blank trials, so you can see online which sites
% give a response before choosing the electrodes for the main experiments.
% Optional visual-only trials (off by default). Short: with 8 candidates and
% 10 repeats it is about 3 minutes.
%
% STATUS: written but NOT yet run in MATLAB or on the rig. Test with
% 'noIntan',true first (see BENCH SWITCHES), then with Intan.
%
% WHY. A poor electrode wastes its share of the current (in expTrainNum,
% expTrainDist and the others). Run this first on the sites you are thinking of
% using, look at the responses, and keep the sites that work.
%
% HOW TO READ IT. Nothing is analysed here. In Open Ephys / the Intan viewer
% compare the activity in the window after the trigger on stimulation trials
% with the blank trials, for the recording sites you care about (the NP sites
% and the Intan channels near the candidates). The condition list printed at
% the start and saved with the data tells you which trials were which channel:
% write down the channels that gave a clear response. Also check IMPEDANCE of
% the candidates in Intan first: a very high or very low value suggests a
% broken or shorted site.
% One caution: a stimulation artifact also appears on the recording; judge the
% response AFTER the artifact, not the artifact itself.
%
%   c = expTrainScreen('CJ286','chn',{'A-003','A-007','A-012','A-017'});
%   c = expTrainScreen('CJ286','chn',{'A-003','A-007','A-012','A-017'}, ...
%                      'amp',[50 75],'nReps',15);
%
% The first input is the SUBJECT (animal ID, a label only, as in exp1A).
%
% STANDALONE: this one file holds the experiment, the condition builder and the
% charge-limit check. It needs only estimIntanMulti.m plus what the
% stimulation PC already has (neurostim, marmolab, the Intan helper functions,
% Psychtoolbox, the DAQ functions).
%
% =========================================================================
% INPUT LABELS  (name-value pairs, e.g. 'amp',50). Every label has a default.
% Only the SUBJECT and 'chn' must be given: a default channel could stimulate
% the wrong site.
%
%   label          meaning                                   default
%   -------------  ----------------------------------------  -------------
%  ELECTRICAL
%   'chn'          Intan names of ALL candidate channels     REQUIRED, e.g. {'A-003','A-007','A-012'}
%                  (one cell array, 1 to about 12 names)
%   'amp'          amplitude LEVELS in uA, one channel       75
%                  stimulated at a time
%   'freq'         pulse frequency within the train, Hz      100
%   'nPulses'      pulses in the train                       5
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
%                  (a contrast level costs 1 condition; an amplitude level costs 5N)
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
% CONDITIONS (built by the local function buildScreenConds below)
% =========================================================================
% For every amplitude level and every candidate channel: ONE train on that
% channel alone (all other channels off). There are NO recovery levels here.
%   blank            nothing delivered
%   visual           one condition per (visType x polarity x contrast), only
%                    if 'withVisual' is true (default FALSE here)
% Count: (channels x amplitude levels) + 1 blank (+ visual conditions).
% Example: 8 channels, 1 amplitude: 9 conditions x 10 repeats x ~2 s = 3 min.
% THE SESSION LENGTH IS PRINTED BEFORE YOU CONFIRM.
%
% =========================================================================
% TRIAL STRUCTURE
% =========================================================================
%   0 ms        trial starts (first frame)
%   pulseTime   trigger fires; the train on the one stimulated electrode starts
%               and the visual patch appears on the same frame   (default 50 ms)
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
p.addParameter('chn',{},@(x) iscellstr(x) && ~isempty(x)); %#ok<ISCLSTR>
p.addParameter('amp',75,@(x) isnumeric(x) && ~isempty(x) && all(x>0));
p.addParameter('freq',100,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('nPulses',5,@(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('fpd',200,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('spd',200,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('ipi',100,@(x) isnumeric(x) && isscalar(x) && x>=0);
% [NEW] visual
p.addParameter('withVisual',false,@islogical);
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
p.addParameter('nReps',10,@(x) isnumeric(x) && isscalar(x) && x>=1);
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
    error('expTrainScreen:chn','Give ''chn'' = a cell array of candidate channel names, e.g. ''chn'',{''A-003'',''A-007'',''A-012''}.');
end
args.chn = args.chn(:)';
if numel(unique(args.chn)) ~= numel(args.chn)
    error('expTrainScreen:chn','The channels in chn must all be different.');
end
if numel(args.chn) > 12
    warning('expTrainScreen:chn','%d candidate channels: every channel is sent every trial, so the parameter exchange gets slower.',numel(args.chn));
end
chnAll = args.chn;                     % channel order used by every condition
if args.trialDur > 550
    warning('expTrainScreen:trialDur', ...
        'trialDur %g ms is over ~550 ms: the end of the trial may fall in the Intan hole.',args.trialDur);
end

% [NEW] train timing, checked before anything is opened
period_ms = 1000/args.freq;
period_us = round(1e6/args.freq);
trainEnd = args.pulseTime + (args.nPulses-1)*period_ms;   % last pulse starts here
if args.withVisual && args.pulseTime + args.visDur > args.trialDur
    error('expTrainScreen:timing','The visual patch (%g ms from %g ms) does not fit in trialDur %g ms.', ...
        args.visDur,args.pulseTime,args.trialDur);
end
if args.withVisual && args.rfX == 0 && args.rfY == 0
    warning('expTrainScreen:rf','rfX and rfY are both 0: this is a PLACEHOLDER. Set them from RF mapping.');
end
if args.pulseTime + args.trigDuration > args.trialDur || trainEnd + 2 > args.trialDur
    error('expTrainScreen:timing', ...
        'The train does not fit: last pulse at %.1f ms, trigger high until %.1f ms, trialDur %g ms.', ...
        trainEnd,args.pulseTime+args.trigDuration,args.trialDur);
end

% [NEW] amplitudes and conditions
ampTot = args.amp(:)';
vis = struct('on',args.withVisual,'type',lower(args.visType), ...
             'polarity',lower(args.polarity),'contrast',args.contrast(:)','sf',args.sf);
[cond,info] = buildScreenConds(ampTot,args.nPulses,numel(chnAll),vis);
nCond = numel(cond);
nTrials = nCond * args.nReps;

% [NEW] charge limit, checked BEFORE anything is opened
if args.fpd ~= args.spd
    warning('expTrainScreen:unbalanced','fpd (%g us) ~= spd (%g us) at equal amplitude: the pulse is NOT charge balanced.',args.fpd,args.spd);
end
% every stimulated channel runs at one of these amplitudes
checkChargeLimit(ampTot,max(args.fpd,args.spd), ...
    'siteArea',args.siteArea,'capacity',args.capacity,'kMax',args.kMax,'allowOver',args.allowOver);

nElec = sum(~strcmp({info.kind},'visual') & ~strcmp({info.kind},'blank'));
fprintf('\nexpTrainScreen: %d conditions x %d repeats = %d trials (about %.0f min at ~2 s per trial)\n', ...
    nCond,args.nReps,nTrials,nTrials*2/60);
fprintf('  %d electrical + %d visual + %d blank conditions\n', ...
    nElec,sum(strcmp({info.kind},'visual')),sum(strcmp({info.kind},'blank')));
fprintf('  candidate channels (order in the amp/ptd/nsp vectors): %s\n',strjoin(chnAll,', '));
fprintf('  amplitude levels [%s] uA, one channel at a time\n',num2str(ampTot));
fprintf('  train: %d pulses at %g Hz (period %.3g ms), last pulse at %.1f ms; phases %g / %g us, interphase %g us\n', ...
    args.nPulses,args.freq,period_ms,trainEnd,args.fpd,args.spd,args.ipi);
fprintf('  charge per phase per train: %s nC\n',num2str(ampTot*max(args.fpd,args.spd)*args.nPulses/1000));
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
    warning('expTrainScreen:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
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
fprintf('expTrainScreen timing: pulse requested at %g ms, on frame %d, delivered at %.2f ms after trial start\n', ...
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
    error('expTrainScreen: DaqDeviceIndex found no DAQ, so no trigger can be sent.');
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
c.paradigm = 'expTrainScreen';
c.addPropsToInform('estim.stimCond');

% [NEW] keep the condition table with the data
stamp = datestr(now,'yyyymmdd_HHMMSS');
condFile = sprintf('%s_expTrainScreen_conds_%s.mat',subject,stamp);
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

function [cond,info] = buildScreenConds(ampLevels,N,nCh,vis)
% BUILDSCREENCONDS  Screening conditions: ONE candidate channel at a time, one
% full train (no recovery levels), plus blank and visual-only conditions.
%
%   [cond,info] = buildScreenConds(ampLevels,N,nCh,vis)
%
% ampLevels   vector of amplitude levels (uA)
% N           pulses in the train
% nCh         number of candidate channels (all sent every trial; only the
%             stimulated one is on, the rest are off)
% vis         struct: on, type, polarity, contrast, sf
%
% info  per condition: kind ('screen'|'blank'|'visual'), chan (index of the
%       stimulated channel, 0 if none), ampEach, level, nLevels, isFinal,
%       and the visual settings contrast, phase (spot 90 bright / 270 dark; gabor 0 / 180), sf.

cond = struct('amp',{},'ptd',{},'nsp',{},'label',{});
info = struct('kind',{},'chan',{},'ampEach',{},'level',{},'nLevels',{},'isFinal',{}, ...
              'contrast',{},'phase',{},'sf',{});

    function add(amp,nsp,label,kind,ch,aE,vc,vp,vsf)
        n = numel(cond) + 1;
        cond(n).amp = amp;  cond(n).ptd = zeros(1,nCh);  cond(n).nsp = nsp;  cond(n).label = label;
        info(n).kind = kind;  info(n).chan = ch;  info(n).ampEach = aE;
        info(n).level = 1;  info(n).nLevels = 1;  info(n).isFinal = true;
        info(n).contrast = vc;  info(n).phase = vp;  info(n).sf = vsf;
    end

multi = numel(ampLevels) > 1;
for a = ampLevels
    if multi, tag = sprintf(' @%guA',a); else, tag = ''; end
    for i = 1:nCh
        amp = zeros(1,nCh);  amp(i) = a;
        nsp = ones(1,nCh);   nsp(i) = N;
        add(amp,nsp,sprintf('chn%d train%s',i,tag),'screen',i,a,0,0,0);
    end
end

z = zeros(1,nCh);
add(z,ones(1,nCh),'blank','blank',0,0,0,0,0);
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
                add(z,ones(1,nCh),lab,'visual',0,0,vis.contrast(ic),ph,sf);
            end
        end
    end
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
