function c = exp1B(subject,varargin)
% EXP1B  Delay. WHEN in the visual response does the artificial volley get
% through? Separates the onset transient from the sustained state.
%
% Same machinery as exp1A -- read that file first. Only the design and the
% grating onset differ.
%
% 14 conditions x 50 repeats = 700 trials, ~23 min.
%
%   c = exp1B('CJ286','rfX',-4,'rfY',-3,'diameter',6,'orientation',45, ...
%             'sf',1,'intanChn','A-007');
%
% =========================================================================
% TRIAL STRUCTURE
% =========================================================================
%   0      trial start, screen blank
%   0      pulse if delay is -100
%   100    grating on
%   150    pulse if delay is  +50
%   200    pulse if delay is +100
%   350    pulse if delay is +250
%   450    pulse if delay is +350
%   600    grating off
%   650    trial ends
%
% Delays in the parameter list are RELATIVE TO GRATING ONSET, which is why
% the grating starts at 100 -- it is what lets -100 put the pulse before any
% visual drive inside an otherwise identical trial.
%
% Trial length is the same in every condition. If it tracked the delay, trial
% length would be confounded with delay and so would adaptation state.
%
% =========================================================================
% WHAT CHANGED FROM THE PROTOCOL, AND WHY
% =========================================================================
% Protocol: contrast {0,25,100} x delay {-100,50,100,250,450} @50 uA = 15,
% plus 3 visual-alone = 18, 600 ms stimulus, 900 trials.
%
% 1. Stimulus 600 -> 500 ms, longest delay 450 -> 350. The Intan hole opens
%    560-650 ms after the trial's first frame (see estimIntan.m), so a
%    720 ms trial loses its tail. At 350 ms into a 500 ms grating the pulse
%    sits at the same fraction of the stimulus as 450 into 600.
%
% 2. The contrast-0 row drops from 6 cells to 2. With no visual stimulus a
%    pulse delay of 50 vs 450 ms measures the same thing, so five of those
%    six cells were duplicates. What remains is one ICMS-in-darkness cell
%    and one empty cell. 18 -> 14 conditions, 900 -> 700 trials.
%
% The 14 conditions are built as parallel vectors in a single factor rather
% than as a cross, because the design is deliberately not rectangular.

aoutChannel = 0;

if ~exist('subject','var')
    error('No subject name provided.');
end
validateattributes(subject,{'char'},{'nonempty'},'','subject',1);
p = inputParser();
p.KeepUnmatched = true;

p.addParameter('rfX',0,@(x) validateattributes(x,{'numeric'},{'nonempty','scalar'}));
p.addParameter('rfY',0,@(x) validateattributes(x,{'numeric'},{'nonempty','scalar'}));
p.addParameter('diameter',6,@(x) validateattributes(x,{'numeric'},{'nonempty','scalar','positive'}));

% [NEW] Aperture edge.
%   'hard'   CIRCLE mask. Contrast is full out to diameter/2 and then drops to
%            zero in one pixel. This is what every run up to 231026 used.
%   'gauss'  GAUSS3 mask. Gaussian envelope truncated (and made transparent)
%            at 3 sigma, with sigma = diameter/6 so the OUTER EXTENT is
%            unchanged: the patch still vanishes at diameter/2.
%
% The gabor shader (GLSLShaders/gabor.frag.txt) implements exactly five masks
% -- GAUSS, CIRCLE, ANNULUS, GAUSS3, GAUSS2D. There is NO raised-cosine, so a
% cosine taper on the outer edge is not reachable without adding a branch to
% the shader. GAUSS3 is the only smooth-edged option with a defined extent.
%
% READ THIS BEFORE USING 'gauss': a Gaussian is not a disc with a softened
% rim. It changes the whole contrast profile, not just the edge. With
% diameter 6 (sigma 1):
%        r =  0 deg   contrast x 1.00
%        r =  1 deg   contrast x 0.61
%        r =  2 deg   contrast x 0.14
%        r =  3 deg   contrast x 0.01  (cut to 0 here)
% The old disc delivered x1.00 everywhere out to r = 3. So 'gauss' drives the
% RF considerably less, which matters for a contrast-gating measurement. The
% startup printout lists the profile for whatever diameter is passed.
p.addParameter('edge','gauss',@(x) ismember(lower(x),{'hard','gauss'}));

p.addParameter('orientation',0,@(x) validateattributes(x,{'numeric'},{'nonempty','scalar'}));
p.addParameter('intanChn','A-007',@ischar);

% [25 Oct] BIPOLAR. Leave returnChn empty for the monopolar behaviour every
% run so far has used. Set it and the SAME pulse goes to that second site
% with the polarity inverted, so current flows between the pair.
%
% returnSep is the centre-to-centre distance in um and YOU must supply it:
% the NeuroNexus channel map is not monotonic, so a distance derived from the
% two channel names would be a guess. estimIntan refuses to start below
% 400 um, and refuses a bipolar run with returnSep unset.
%
% Cost at the rig: one extra parameter_update per trial. Still ONE SET and so
% still one board restart and one Intan hole -- this is RunMultiElect_exp.m's
% multi-electrode pattern (lines 109-127), not a second exchange.
p.addParameter('returnChn','',@ischar);         % '' = monopolar
p.addParameter('returnSep',NaN,@isnumeric);     % [um] >= 400, required if bipolar
p.addParameter('returnPolarity',[], ...         % [] = opposite of stimPolarity
    @(x) isempty(x) || isnumeric(x));

p.addParameter('sf',1,@isnumeric);
p.addParameter('tf',4,@isnumeric);
p.addParameter('gratingOn',100);            % [ms] so the -108.14 control fits
p.addParameter('stimDuration',500);         % [ms] grating 100-600
% [25 Oct] The pre-grating control was -100. At 60 Hz that needs an 8.135 ms
% busy-wait, and the wait has to be timed from a flip the grating has not made
% yet -- the one anchor in this file that no run has ever exercised. -108.14 is
% the nearest whole flip (nf = -5), so its residual is zero and the pulse fires
% straight on a measured flip with no wait and no anchor. The control only has
% to land well before the visual response (~50-65 ms AFTER the grating); 8 ms
% earlier changes nothing physiological. Set it back to -100 if you ever need
% the round number -- the anchor below is fixed and will handle it.
p.addParameter('delays',[-108.14 50 100 250 350]);  % [ms] RELATIVE to grating on
p.addParameter('contrasts',[0.25 1.00]);    % the two driven contrasts
p.addParameter('current',50);               % [uA]
p.addParameter('blankDelay',250);           % [ms] the one delay used for the
                                            %      two contrast-0 cells
p.addParameter('onsetCorrection',0);        % [ms] SUPERSEDED -- keep at 0. See
                                            %      exp1A. Replaced by K, applied
                                            %      in frames.
p.addParameter('K',-24.8);                  % [ms] display-path constant.
                                            %      measured diode->trigger
                                            %        = (Fe-Fg)*flip + K
                                            %      -24.8 +/- 1.1 ms. Full
                                            %      derivation in exp1A. This is
                                            %      the protocol that proved it:
                                            %      one constant fits all five
                                            %      delays from -4 to +23 frames
                                            %      to within 2 ms.
p.addParameter('nominalFrameRate',60);      % chooses the frame counts only
p.addParameter('delayFrames',[]);           % pin the frames, bypassing ms.
                                            %      Must be numel(delays) long.
p.addParameter('subFrame',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
p.addParameter('rMax',10);   % [ms] see the note above the frame block                  % [ms] longest busy-wait allowed in
                                            %      afterFlip. ~16.7 ms until the
                                            %      next draw deadline, ~3 ms
                                            %      already spent, so 12 is the
                                            %      safe ceiling. A longer
                                            %      residual goes on the next
                                            %      flip instead: at most
                                            %      (ivl-rMax) late, but drops no
                                            %      frame. A dropped frame would
                                            %      be a 24 deg gabor phase stall
                                            %      locked to every pulse.
p.addParameter('postPulsePause',0.35);      % [s] estimIntan's PULSE recovery.
                                            %      0.35, not the 0.2 default:
                                            %      exp1B is the only protocol
                                            %      where the Intan hole reaches
                                            %      into the post-pulse window.
                                            %      At 0.2 the hole opened as
                                            %      close as 179 ms after the
                                            %      trigger and 4 of 140 trials
                                            %      lost part of a 200 ms
                                            %      window; 0.35 moves the
                                            %      minimum to ~330. The hole
                                            %      closes ~770 ms before the
                                            %      next trial's first frame, so
                                            %      there is room for +600 ms
                                            %      before it clips the NEXT
                                            %      trial -- do not exceed 0.5.
p.addParameter('trigDuration',100);
p.addParameter('iti',1.5);
p.addParameter('nrRepeats',50);

p.addParameter('noIntan',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('skipSync',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));

p.parse(varargin{:});
args = p.Results;

% [25 Oct] This used to read min(delays) + gratingOn < 0, which ignored K and
% so rejected any delay more negative than -gratingOn. The pulse does not fire
% at the requested delay; it fires K earlier in flip time, i.e. at
% gratingOn + delay - K after the trial starts. With K = -24.8 that is 24.8 ms
% of headroom the old check threw away, and it would have refused -108.14 even
% though that lands at +16.7 ms. onAbs >= 0 below is the real guard; this one
% only exists to fail early with a readable message.
earliest = min(args.delays) + args.onsetCorrection - args.K + args.gratingOn;
if earliest < 0
    error(['exp1B: delay %.1f with gratingOn %d fires %.1f ms before the trial ' ...
           'starts. Raise gratingOn by at least %.0f.'], ...
           min(args.delays),args.gratingOn,-earliest,ceil(-earliest));
end
if max(args.delays) + 150 > args.stimDuration
    error(['exp1B: delay %.1f leaves only %.1f ms of grating after the pulse. Keep ' ...
           'at least 150, or that condition measures the response to grating ' ...
           'OFFSET rather than to a later delay.'], ...
           max(args.delays),args.stimDuration-max(args.delays));
end

% -------------------------------------------------------------------------
% [NEW] THE 14 CONDITIONS, as parallel vectors
% -------------------------------------------------------------------------
% For each driven contrast: the five delays with current, plus one
% visual-alone. Then two contrast-0 cells: ICMS in darkness, and nothing.
nD   = numel(args.delays);
cnd  = struct('contrast',[],'on',[],'fpa',[]);
for k = 1:numel(args.contrasts)
    cnd.contrast = [cnd.contrast, repmat(args.contrasts(k),1,nD+1)];
    cnd.on       = [cnd.on,       args.delays, args.blankDelay];
    cnd.fpa      = [cnd.fpa,      repmat(args.current,1,nD), -1];
end
cnd.contrast = [cnd.contrast, 0,                0];
cnd.on       = [cnd.on,       args.blankDelay,  args.blankDelay];
cnd.fpa      = [cnd.fpa,      args.current,    -1];

% onAbs is NOT computed here any more. It needs c.screen.frameRate, so it is
% built after rigcfg -- see the frame block below.

nCond = numel(cnd.fpa);
fprintf('exp1B: %d conditions x %d repeats = %d trials\n', ...
        nCond,args.nrRepeats,nCond*args.nrRepeats);
fprintf('       grating %d-%d ms; delays (rel) %s\n', ...
        args.gratingOn,args.gratingOn+args.stimDuration,mat2str(args.delays));

import neurostim.*
commandwindow;

if args.skipSync
    Screen('Preference','SkipSyncTests',1);
    warning('exp1B:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
end

c = marmolab.rigcfg('debug',args.debug,'subject',subject,p.Unmatched);
c.oephys.fakeConnection = false;   % MUST stay false -- see exp1A

% -------------------------------------------------------------------------
% [NEW] SUB-FRAME TIMING
% -------------------------------------------------------------------------
% The grating's diode patch is drawn in the same pass as the gabor
% (stimulus.m:461-466), so the diode marks the grating's flip exactly. The
% DAQ line, though, is a separate wire -- nothing ties the pulse to a flip
% except that onsetFunction is called from afterFlip (stimulus.m:566). So
% rather than quantising, anchor to the grating's MEASURED flip timestamp:
%
%   T = gratingFlipPtb + (wanted - K)
%
% and the delivered delay is the requested one. estim is still scheduled on
% the last flip before T; its onsetFunction busy-waits the residual rW and
% fires. rW is capped at rMax so the frame loop is never blocked long enough
% to miss the next draw deadline.
%
% The grating's own onset is SNAPPED to a flip first, so the offset is
% measured from a real flip rather than a requested millisecond.
%
% nf and rW are VECTORS here, one entry per condition, because every delay
% has its own residual. firePulse looks up the current trial's entry by the
% estim.on the design wrote.
ivl = 1000/c.screen.frameRate;

% The grating's own onset is snapped to a flip first, so the pulse offset is
% measured from a real flip rather than from a requested millisecond.
mG  = round(args.gratingOn/ivl);
gOn = mG*ivl;

% off(k) = how long after the grating's FLIP condition k should fire, so that
% it lands cnd.on(k) after the LIGHT. K carries the whole display path.
off = cnd.on + args.onsetCorrection - args.K;

if isempty(args.delayFrames)
    if args.subFrame
        nf = floor(off/ivl);
        rW = off - nf*ivl;                   % residual per condition, [0,ivl)
        bump = rW > args.rMax;               % too long to block in afterFlip
        nf(bump) = nf(bump) + 1;  rW(bump) = 0;
        if any(rW > args.rMax - 2)
            warning('exp1B:rMaxMargin', ...
                ['%d condition(s) wait within 2 ms of rMax %.1f. 8.1 ms drops ' ...
                 'no frame, 11.47 ms drops one (231023).'], ...
                 sum(rW > args.rMax - 2), args.rMax);
        end
    else
        nf = round(off/(1000/args.nominalFrameRate));  rW = zeros(size(off));
    end
else
    assert(numel(args.delayFrames)==numel(args.delays), ...
        'exp1B: delayFrames must have %d entries, one per delay.', numel(args.delays));
    nf = nan(size(cnd.on));
    for k = 1:numel(args.delays)
        nf(abs(cnd.on - args.delays(k)) < 1e-9) = args.delayFrames(k);
    end
    nf(abs(cnd.on - args.blankDelay) < 1e-9 & isnan(nf)) = ...
        args.delayFrames(find(abs(args.delays-args.blankDelay)<1e-9,1));
    assert(~any(isnan(nf)),'exp1B: delayFrames did not cover every condition.');
    rW = zeros(size(nf));
end

onAbs   = gOn + nf*ivl;        % commanded estim.on, per condition
gOffVec = nf*ivl + rW;         % fire this long after the grating's flip
gPtb = NaN; gTrial = -1;       % shared with the nested functions below

% -100 with gratingOn 100 puts the pulse on frame 1 or earlier; guard it.
assert(all(onAbs >= 0), ...
    'exp1B: condition(s) %s command a negative estim.on. Raise gratingOn.', ...
    mat2str(cnd.on(onAbs < 0)));

fprintf('exp1B timing: flip %.4f ms (%.4f Hz measured), K %.1f ms\n', ...
        ivl, c.screen.frameRate, args.K);
fprintf('  grating on frame %d, i.e. %.3f ms (asked %g)\n', mG+1, gOn, args.gratingOn);
fprintf('  %6s %7s %7s %10s %11s\n','asked','frames','wait','commanded','delivered');
[uq,iq] = unique(cnd.on,'stable');
for q = 1:numel(uq)
    fprintf('  %6g %7d %7.2f %10.3f %11.2f\n', ...
            uq(q), nf(iq(q)), rW(iq(q)), onAbs(iq(q)), gOffVec(iq(q)) + args.K);
end
fprintf('  (delivered = frames*flip + wait + K)\n');
if args.subFrame && any(rW == 0 & off ~= 0)
    fprintf('  NOTE %d condition(s) exceeded rMax %.1f ms and went on the next flip.\n', ...
            sum(rW == 0 & off ~= 0), args.rMax);
end

%% grating
g = stimuli.gabor(c,'grating');
g.X           = args.rfX;
g.Y           = args.rfY;
switch lower(args.edge)
    case 'hard'
        g.mask  = 'CIRCLE';               % hard-edged disc
        g.sigma = args.diameter/2;        % CIRCLE: sigma is the RADIUS
    case 'gauss'
        g.mask  = 'GAUSS3';               % Gaussian, transparent beyond 3 sigma
        g.sigma = args.diameter/6;        % 3*sigma = diameter/2, same extent
end
g.width       = ceil(args.diameter) + 1;   % integer -- see exp1A
g.height      = ceil(args.diameter) + 1;
fprintf('exp1B aperture: %s mask, sigma %.3f deg, extent +/-%.2f deg\n', ...
        g.mask, g.sigma, args.diameter/2);
if strcmpi(args.edge,'gauss')
    rr = 0:args.diameter/6:args.diameter/2;
    fprintf('  contrast profile r =');  fprintf(' %5.2f',rr);  fprintf(' deg\n');
    fprintf('                  x    ');
    fprintf(' %5.2f',exp(-0.5*(rr/g.sigma).^2));  fprintf('\n');
end
g.frequency   = args.sf;
g.orientation = args.orientation;
g.phaseSpeed  = 360*args.tf/c.screen.frameRate;
g.contrast    = NaN;                  % design writes it
g.on          = gOn;                  % ~100, NOT 0 -- see header. Snapped to a
                                      % flip by the frame block above.
g.duration    = args.stimDuration;
g.phase = plugins.jitter(c,{0,360},'distribution','uniform');

%% electrical stimulus
e = estimIntan(c,'estim');
e.on          = NaN;                  % design writes it
e.duration    = args.trigDuration;
e.intanChn    = args.intanChn;
e.returnChn      = args.returnChn;       % '' = monopolar
e.returnSep      = args.returnSep;       % [um], checked against the 400 floor
e.returnPolarity = args.returnPolarity;  % [] = inverted automatically
e.testMode    = args.noIntan;
e.postPulsePause = args.postPulsePause;   % 0.35 here, not the 0.2 default --
                                          % see the parameter comment above.

e.fpa  = NaN;      e.spa  = NaN;
e.fpd  = 100;      e.spd  = 100;      e.ipi  = 60;
e.pot  = 0;        e.nsp  = 1;        e.pulseTrainPeriod = 10000;
e.ptr  = 1000;     e.stSH = 1;        e.stimPolarity = 0;
e.enAS = 1;        e.prAS = 200;      e.poAS = 160000;
e.maAS = 1;        e.enCR = 1;        e.prCR = 100000;   e.poCR = 150000;
e.triggerSource = 0; e.triggerEdgeOrLevel = 0; e.triggerHighOrLow = 0;
e.postTriggerDelay = 0;

%% Intan TCP
if ~args.noIntan
    fprintf('\n---- connect Intan to the TCP/IP server now ----\n');
    e.tcp = local_TCPIP_server();
    fprintf(e.tcp,'NOMCS');
    fprintf('Intan connected.\n');
end

%% trigger
DAQ = DaqDeviceIndex;
if isempty(DAQ)
    error('exp1B: DaqDeviceIndex found no DAQ, so no trigger can be sent.');
end
fprintf('DAQ index %s, trigger on analog out %d\n',mat2str(DAQ),aoutChannel);
DaqAOut(DAQ,aoutChannel,0);

e.daq            = DAQ;
e.aoutChannel    = aoutChannel;

% [NEW] The grating stamps its own flip time; estim fires gOffVec(k) after it.
% Both run from afterFlip (stimulus.m:566), which passes ptbTime in ms on the
% GetSecs clock. The grating turns on no later than estim and the plugin order
% is ('oephys','grating','estim'), so even when both land on the same flip the
% grating's callback is registered first and runs first -- gPtb is always set
% before firePulse reads it. DO NOT reorder the plugins.
g.onsetFunction  = @(o,t) captureGrating(o,t);
e.onsetFunction  = @(o,t) firePulse(o,t);
e.offsetFunction = @(o,t) DaqAOut(DAQ, aoutChannel, 0);

if ~args.debug
    c.addScript('beforeExperiment', @(c) doBeforeExperiment(c));
    c.addScript('afterExperiment', @(c) doAfterExperiment(c));
end

    function doBeforeExperiment(cc)
        io.forceBlackDiode(cc);
    end
    function doAfterExperiment(cc) %#ok<INUSD>
        io.beep3;
    end

    % [NEW] Nested, so they share gPtb/gTrial/gOffVec/rW with the body above.
    % Local functions get their own workspace and could not see the timestamp.
    function captureGrating(o,t)
        gPtb   = t;                 % ms, GetSecs clock, this frame's flip
        gTrial = o.cic.trial;
    end

    function firePulse(o,tE)
        % nf/rW are per condition, so find this trial's entry by the estim.on
        % the design wrote. Conditions that share an on-value share a residual,
        % so the first match is always the right one.
        k = find(abs(onAbs - o.on) < 1e-6, 1);
        if isempty(k)
            DaqAOut(DAQ,aoutChannel,1); o.markPulse(); return;
        end
        % Anchor. Prefer the grating's MEASURED flip, but for the -100
        % condition estim is on an earlier frame than the grating, so the
        % grating has not flipped yet. Fall back to estim's OWN measured flip,
        % which is the same instant: onAbs = gOn + nf*ivl (line 263), so
        % estim flips exactly nf*ivl after the grating, and
        %     tE + rW  ==  gPtb + nf*ivl + rW  ==  gPtb + gOffVec.
        % [25 Oct] This replaced an anchor on cic.firstFrame, which is NOT an
        % absolute GetSecs time -- it logs 0.3-0.7 in test7 -- so the
        % WaitSecs('UntilTime') target was already long past and the residual
        % was silently skipped on exactly the condition this branch exists for.
        T = NaN;
        if gTrial == o.cic.trial && ~isnan(gPtb)
            T = gPtb + gOffVec(k);
        elseif ~isempty(tE) && isfinite(tE)
            T = tE + rW(k);
        end
        if rW(k) > 0 && ~isnan(T)
            WaitSecs('UntilTime',T/1000);           % ms -> s
        end
        DaqAOut(DAQ,aoutChannel,1);
        o.markPulse();              % stamp after the edge, never before
    end

    % Patch on the grating, so the burst marks the stimulus. It starts
    % 100 ms into the trial, so diode->trigger reads NEGATIVE for the -100
    % delay. That is correct, not a fault.
    addDiodeFlasher(c.diode,'enabled',false);
    addDiodeFlasher(g,'size',0.075,'location','nw');

c.trialDuration = '@grating.stopTime+50';
c.iti = args.iti*1e3;

%% design -- 14 conditions in one factor, fully interleaved
d = neurostim.design('myDesign');
d.fac1.grating.contrast = cnd.contrast;
d.fac1.estim.on         = onAbs;
d.fac1.estim.fpa        = cnd.fpa;
d.fac1.estim.spa        = cnd.fpa;

blk = block('myBlock',d);
blk.nrRepeats = args.nrRepeats;

c.subject = subject;
c.paradigm = 'exp1B';
c.addPropsToInform('grating.contrast','estim.on','estim.fpa');

c.setPluginOrder('oephys','grating','estim');

try
    c.run(blk);
catch me
    cleanup(e);
    rethrow(me);
end
cleanup(e);

end

function cleanup(e)
try; sca; end                                      %#ok<TRYNC>
try; ListenChar(0); end                            %#ok<TRYNC>
try
    if ~isempty(e.daq); DaqAOut(e.daq,e.aoutChannel,0); end
end                                                %#ok<TRYNC>
try
    if ~isempty(e.tcp); fclose(e.tcp); end
end                                                %#ok<TRYNC>
end
