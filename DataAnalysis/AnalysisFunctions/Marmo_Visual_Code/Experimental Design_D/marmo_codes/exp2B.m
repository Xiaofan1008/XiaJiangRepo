function c = exp2B(subject,varargin)
% EXP2B  Pairing, ICMS LEADS. The other sign of exp2A: the artificial V1
% volley reaches V2 BEFORE the visual drive does.
%
% Same machinery as exp1A -- read that file first.
%
% 2 conditions x 750 repeats = 1500 presentations, 750 pulses.
%
%   c = exp2B('CJ286','rfX',-4,'rfY',-3,'diameter',6,'orientation',45, ...
%             'sf',1,'intanChn','A-007','visLatency',58);
%
% Use the SAME oriOffsets and pairCurrent as exp2A, and the matching
% exp2A sites, or the two arms are not comparable.
%
% =========================================================================
% WHAT THIS FILE IS
% =========================================================================
% The INDUCTION only. Two orientations interleaved 50/50:
%
%   theta_A = theta_stim + 45    WITH ICMS      750 trials
%   theta_B = theta_A  + 90      WITHOUT ICMS   750 trials  (no trigger)
%
% theta_A and theta_B are 90 deg apart and sit symmetrically about
% theta_stim, so anything driven by the stimulated column's preference
% lands equally on both. The result is therefore a WITHIN-UNIT ASYMMETRY of
% tuning shift, not a pre/post change -- adaptation and non-specific
% excitability are symmetric and cancel.
%
% Run the tuning blocks yourself, either side, without the stimulator:
%
%   1. PRE  tuning (24 directions, your protocol)
%   2. exp2B
%   3. POST tuning at t = 0, 30, 90 min
%
% ONE SIGN PER PENETRATION. After pairing the site is no longer naive, so
% exp2B and exp2B cannot both be run at one penetration. Sites 1 and 2 in
% each animal, one protocol each, order counterbalanced across animals.
%
% =========================================================================
% WHY theta_B HAS NO TRIGGER
% =========================================================================
% The Intan parameters never change in this run: one current, set once.
% theta_B trials simply do not fire the trigger (estim.on = Inf), so no
% current is delivered. That means parameter_update/SET runs ONCE, on trial
% 1, and never again -- so the Intan board never restarts and the recording
% stays continuous. estimIntan.skipIfUnchanged is what allows that, and it
% is on only here.
%
% The cost: no digital edge on theta_B trials, so those 750 trials are
% located by photodiode, not by trigger. chk_intan will report 750 triggers
% against 1500 trials; that is correct, not a failure.
%
% =========================================================================
% THE PULSE TIME -- the number that decides the whole experiment
% =========================================================================
%   estim.on is derived from pulseTime, which is HARDCODED at 33 ms after
%   the grating becomes VISIBLE. This is the EARLY condition (ICMS leads:
%   the pulse arrives BEFORE the visual response begins).
%
% WHY IT IS HARDCODED, AND NOT visLatency + deltaT.
% deltaT is defined against the visual RESPONSE onset, so computing the pulse
% time from it makes the stimulus depend on a latency you have estimated. A
% 10 ms error in that estimate then mistimes every pulse of the induction.
% Hardcoding the absolute time removes the estimate from the stimulus
% entirely: the pulse lands at 33 ms whatever the latency turns out to be,
% and deltaT becomes something you REPORT afterwards rather than something
% you command. visLatency is still a parameter, but it only labels the
% result -- a wrong value mislabels, it cannot mistime.
%
% WHY 33 AND 76 SPECIFICALLY (chosen 241024).
% The pulse can only sit a few ms after a flip before the busy-wait starts
% costing a frame (8.1 ms is safe, 11.5 ms dropped one on 231023). That
% leaves usable windows at 25-33, 42-50, 59-67 and 75-83 ms. The gap from 67
% to 75 needs an 8.5-14.5 ms wait and is unusable.
%
% Within those windows, 33/76 stays at least 4 SD clear of the visual
% response for ANY latency between 44 and 65 ms, and is symmetric (+/-21.5)
% at a latency of 54.5. Schuett, Bonhoeffer & Hubener 2001 used +15/-15
% around a ~50 ms cat V1 latency, i.e. 35/65 -- but 35/65 only survives a
% latency of 46-54, and above 54 the late pulse falls inside the response and
% silently turns pairing trials into anti-pairing ones.
%
% IF YOU MEASURE THE LATENCY and it is 55 or less, move exp2A's late
% condition from 76 to 66: that gives +/-18, closer to the STDP optimum.
% Above 55, keep 76. This file's 33 is unaffected either way.
%
% The other condition of this pair is exp2A at 76 ms; the two protocols must
% use the same pair.
%
% =========================================================================
% TIMING
% =========================================================================
% =========================================================================
%   0      grating on, 500 ms, 100% contrast, randomised phase
%   ~70    pulse (theta_A trials only)
%   500    grating off
%   550    trial ends
%
% Protocol asks for 1 Hz. With the Open Ephys messages left on, beforeTrial
% is ~1000 ms and the cycle is ~1520 ms, so 1500 presentations take ~38 min
% rather than 25. Cut nrPairing to 1000 if the 25 min matters more than the
% pairing count.

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
p.addParameter('stimDuration',500);         % [ms] protocol pairing stimulus
p.addParameter('contrast',1.00);            % 100%, protocol

% ---- the three numbers that set the pulse time -------------------------
p.addParameter('pulseTime',33);           % [ms] after the grating becomes
                                            %      VISIBLE. HARDCODED on purpose.
                                            %      This is the EARLY condition.
                                            %
                                            % 33 and 76 are the pair chosen on
                                            % 241024. They are NOT free numbers:
                                            % the pulse can only sit a few ms
                                            % after a flip without a long
                                            % busy-wait, so the usable windows
                                            % are 25-33, 42-50, 59-67, 75-83 ms.
                                            % Everything from 67 to 75 needs an
                                            % 8.5-14.5 ms wait, which is the
                                            % range that dropped a frame on
                                            % 231023.
                                            %
                                            % 33/76 stays >=4 SD clear of the
                                            % visual response for ANY latency
                                            % from 44 to 65 ms, and is symmetric
                                            % (+/-21.5) at a latency of 54.5.
                                            % Schuett's 35/65 only survives a
                                            % latency of 46-54.
                                            %
                                            % If you MEASURE the latency and it
                                            % is <= 55, move the late condition
                                            % to 66: that gives +/-18 and sits
                                            % closer to the STDP optimum.
p.addParameter('visLatency',60);            % [ms] REPORTING ONLY. It does not
                                            %      touch the stimulus any more:
                                            %      the pulse time is hardcoded
                                            %      above. This is used solely to
                                            %      print the delivered deltaT,
                                            %      so a wrong value here mislabels
                                            %      the result but cannot mistime
                                            %      the pulse. That is the whole
                                            %      point of hardcoding.
p.addParameter('onsetCorrection',0);        % [ms] SUPERSEDED -- keep at 0. See
                                            %      exp1A. Replaced by K, applied
                                            %      in frames.
p.addParameter('K',-24.8);                  % [ms] display-path constant.
                                            %      measured diode->trigger
                                            %        = (Fe-Fg)*flip + K
                                            %      Light arrives ~25 ms AFTER
                                            %      the flip; the trigger follows
                                            %      its flip by 2-4 ms, because
                                            %      onsetFunction runs inside
                                            %      afterFlip (stimulus.m:566).
                                            %      -24.8 +/- 1.1 ms. Derivation
                                            %      in exp1A.
                                            %      exp2B IS THE DANGEROUS ONE:
                                            %      on 21 Oct, commanding 70 ms
                                            %      with onsetCorrection 25 was
                                            %      numerically identical to
                                            %      exp2A with correction 0. One
                                            %      wrong constant and 2B fires
                                            %      where 2A should.
p.addParameter('nominalFrameRate',60);      % chooses the frame count only
p.addParameter('pulseFrame',[]);            % pin the frame, bypassing ms
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

p.addParameter('oriOffsetA',45);            % theta_A = orientation + 45
p.addParameter('oriOffsetB',135);           % theta_B = orientation + 135
p.addParameter('pairCurrent',50);           % [uA]. -1 = no-current control.
p.addParameter('nrPairing',1500);           % total presentations, half paired

p.addParameter('trigDuration',100);
p.addParameter('iti',0.5);                  % [s] as short as the rig allows;
                                            %     beforeTrial dominates anyway
p.addParameter('noIntan',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('skipSync',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));

p.parse(varargin{:});
args = p.Results;

if rem(args.nrPairing,2) ~= 0
    error('exp2B: nrPairing must be even so the two orientations are 50/50.');
end

oriA    = mod(args.orientation + args.oriOffsetA, 360);
oriB    = mod(args.orientation + args.oriOffsetB, 360);
% Nominal target, for the sanity checks only. The frame-exact command is built
% after rigcfg, where c.screen.frameRate exists.
wantOn = args.pulseTime + args.onsetCorrection;   % absolute, hardcoded

if wantOn < 0
    error('exp2B: pulse would be at %g ms, before the trial starts.',wantOn);
end
if wantOn > args.stimDuration - 100
    error('exp2B: pulse at %g ms leaves under 100 ms of grating after it.',wantOn);
end

fprintf('exp2B: %d presentations, %d pulses, ~%.0f min\n', ...
        args.nrPairing,args.nrPairing/2,args.nrPairing*1.52/60);
fprintf('       theta_A %g deg WITH %g uA, theta_B %g deg WITHOUT\n', ...
        oriA,args.pairCurrent,oriB);
if args.pairCurrent < 0
    fprintf('       *** NO-CURRENT CONTROL RUN ***\n');
end

import neurostim.*
commandwindow;

if args.skipSync
    Screen('Preference','SkipSyncTests',1);
    warning('exp2B:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
end

c = marmolab.rigcfg('debug',args.debug,'subject',subject,p.Unmatched);
c.oephys.fakeConnection = false;   % MUST stay false -- see exp1A

% -------------------------------------------------------------------------
% [NEW] FRAME-BASED PULSE TIME -- the part of this protocol that matters most
% -------------------------------------------------------------------------
% The pulse time is hardcoded (pulseTime), so nothing here depends on the
% latency estimate. The DAQ is a separate wire from the display, so the pulse
% is anchored to the grating's MEASURED flip timestamp and fired at
%   T = gratingFlipPtb + (pulseTime - K)
% which delivers the requested time exactly, with no frame quantisation. estim is scheduled on the last
% flip before T and busy-waits the residual rW inside afterFlip. The wait is
% capped at rMax so the frame loop is never blocked long enough to drop a
% frame -- a dropped frame here would be a 24 deg gabor phase stall at a fixed
% time relative to every pulse, i.e. an artifact locked to the manipulation.
% If the residual exceeds rMax the pulse goes on the NEXT flip instead, at
% most (flip - rMax) late; the printout says when that happens.
% The sign of the DELIVERED value is still checked below: getting it wrong
% turns exp2B into exp2A.
ivl = 1000/c.screen.frameRate;
off = wantOn - args.K;                 % ms after the grating's FLIP
if args.subFrame
    nf = floor(off/ivl);
    rW = off - nf*ivl;                 % residual to busy-wait, [0,ivl)
    if rW > args.rMax
        nf = nf + 1;  rW = 0;          % too long to block; next flip
    elseif rW > args.rMax - 2
        warning('exp2B:rMaxMargin', ...
            ['busy-wait is %.2f ms against an rMax of %.1f. Measured ' ...
             '231023: 8.1 ms drops no content frame, 11.47 ms drops one, ' ...
             'and the dropped-frames check CANNOT see it -- it reads the ' ...
             'panel 120 Hz ripple, not the content rate. Verify the ' ...
             'trigger width against its prediction in chk_exp2.'], ...
             rW, args.rMax);
    end
else
    % pure frame scheduling: use the NOMINAL rate so the frame count is
    % deterministic across sessions (see K / x.488 note in exp1A).
    nf = round(off/(1000/args.nominalFrameRate));  rW = 0;
end
if ~isempty(args.pulseFrame), nf = args.pulseFrame; rW = 0; end
assert(nf >= 0, ...
    'exp2B: pulse would be %d frames BEFORE the grating.', nf);
pulseOn = nf*ivl;
gOffMs  = nf*ivl + rW;                 % fire this long after the flip
gOn  = 0;                     % grating is frame 1 in this protocol
gPtb = NaN; gTrial = -1;               % shared with nested functions

fprintf('\nexp2B timing: flip %.4f ms (%.4f Hz measured), K %.1f ms\n', ...
        ivl, c.screen.frameRate, args.K);
fprintf('  pulse hardcoded at %g ms after VISIBLE onset\n', args.pulseTime);
fprintf('  visLatency %g ms (reporting only) -> nominal deltaT %+.1f\n', ...
        args.visLatency, args.pulseTime - args.visLatency);
fprintf('  pulse on frame %d + %.2f ms wait, commanded estim.on = %.3f ms\n', ...
        nf+1, rW, pulseOn);
fprintf('  DELIVERED deltaT = %+.2f ms   <-- the number for the methods\n', ...
        gOffMs + args.K - args.visLatency);
if args.subFrame && rW == 0
    fprintf(['  NOTE residual exceeded rMax %.1f ms, so the pulse went on the\n' ...
             '       next flip, %.2f ms from the request. Raise rMax towards\n' ...
             '       %.1f for exactness, at the cost of dropped frames.\n'], ...
            args.rMax, gOffMs + args.K - args.pulseTime, ivl);
end

if sign(gOffMs + args.K - args.visLatency) ~= -1
    warning('exp2B:sign', ...
        ['DELIVERED deltaT is %+.1f ms -- NOT negative. exp2B is the ' ...
         'ICMS-leads protocol. Pick an earlier frame.'], ...
         gOffMs + args.K - args.visLatency);
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
fprintf('exp2B aperture: %s mask, sigma %.3f deg, extent +/-%.2f deg\n', ...
        g.mask, g.sigma, args.diameter/2);
if strcmpi(args.edge,'gauss')
    rr = 0:args.diameter/6:args.diameter/2;
    fprintf('  contrast profile r =');  fprintf(' %5.2f',rr);  fprintf(' deg\n');
    fprintf('                  x    ');
    fprintf(' %5.2f',exp(-0.5*(rr/g.sigma).^2));  fprintf('\n');
end
g.frequency   = args.sf;
g.orientation = NaN;                  % design writes it
g.phaseSpeed  = 360*args.tf/c.screen.frameRate;
g.contrast    = args.contrast;
g.on          = 0;
g.duration    = args.stimDuration;

% Randomised starting phase. Without it every pairing lands at the same
% point in the 4 Hz response cycle, so the induction happens at one phase
% rather than averaged over the cycle.
g.phase = plugins.jitter(c,{0,360},'distribution','uniform');

%% electrical stimulus
e = estimIntan(c,'estim');
e.on          = NaN;                  % design writes it: pulseOn, or
                                      % NOPULSE on theta_B trials
e.duration    = args.trigDuration;
e.intanChn    = args.intanChn;
e.returnChn      = args.returnChn;       % '' = monopolar
e.returnSep      = args.returnSep;       % [um], checked against the 400 floor
e.returnPolarity = args.returnPolarity;  % [] = inverted automatically
e.testMode    = args.noIntan;

% The current is FIXED for the whole run and never re-sent. theta_B trials
% withhold the trigger instead of changing the current -- see the header.
e.fpa  = args.pairCurrent;
e.spa  = args.pairCurrent;
e.skipIfUnchanged = true;             % ON here only. One SET, no restarts.

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
    error('exp2B: DaqDeviceIndex found no DAQ, so no trigger can be sent.');
end
fprintf('DAQ index %s, trigger on analog out %d\n',mat2str(DAQ),aoutChannel);
DaqAOut(DAQ,aoutChannel,0);

e.daq            = DAQ;
e.aoutChannel    = aoutChannel;

% [NEW] The grating stamps its own flip time; estim fires gOffMs after it.
% Both run from afterFlip (stimulus.m:566), which passes ptbTime in ms on the
% GetSecs clock. The grating turns on no later than estim and the plugin order
% is ('oephys','grating','estim'), so even when both land on the same flip the
% grating's callback is registered first and runs first -- gPtb is always set
% before firePulse reads it. DO NOT reorder the plugins.
g.onsetFunction  = @(o,t) captureGrating(o,t);
e.onsetFunction  = @(o,t) firePulse(o);
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

    % [NEW] Nested, so they share gPtb/gTrial/gOffMs/rW with the body above.
    % Local functions get their own workspace and could not see the timestamp.
    function captureGrating(o,t)
        gPtb   = t;                 % ms, GetSecs clock, this frame's flip
        gTrial = o.cic.trial;
    end

    function firePulse(o)
        % Anchor. Prefer the grating's MEASURED flip, but it does not exist
        % yet when the pulse precedes the grating (exp1B's -100 condition
        % fires on frame 2 while the grating comes on at frame 7). In that
        % case anchor to the trial's first flip: the grating's flip is
        % firstFrame + gOn by construction, and with zero dropped frames the
        % two agree to the microsecond.
        T = NaN;
        if gTrial == o.cic.trial && ~isnan(gPtb)
            T = gPtb + gOffMs;
        else
            ff = o.cic.firstFrame;
            if ~isempty(ff) && isfinite(ff)
                T = ff + gOn + gOffMs;
            end
        end
        if rW > 0 && ~isnan(T)
            WaitSecs('UntilTime',T/1000);           % ms -> s
        end
        DaqAOut(DAQ,aoutChannel,1);
        o.markPulse();              % stamp after the edge, never before
    end

    % Patch on the grating. It marks EVERY trial, paired or not, which is
    % what locates the theta_B trials in both recordings.
    addDiodeFlasher(c.diode,'enabled',false);
    addDiodeFlasher(g,'size',0.075,'location','nw');

c.trialDuration = '@grating.stopTime+50';
c.iti = args.iti*1e3;

%% design -- 2 conditions, 50/50
% NOPULSE is a trial time the trial never reaches, so estim never turns on,
% onsetFunction never runs and no trigger is sent. A large finite number
% rather than Inf: stimulus.baseBeforeFrame computes
% round(on*frameRate/1000)+1, and Inf there is asking for trouble.
% The current stays loaded in the Intan and is simply not triggered.
NOPULSE = 1e6;   % [ms]
d = neurostim.design('pairing');
d.fac1.grating.orientation = [oriA  oriB];
d.fac1.estim.on            = [pulseOn NOPULSE];

blk = block('pairing',d);
blk.nrRepeats = args.nrPairing/2;

c.subject = subject;
c.paradigm = 'exp2B';
c.addPropsToInform('grating.orientation','estim.on','estim.fpa');

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
