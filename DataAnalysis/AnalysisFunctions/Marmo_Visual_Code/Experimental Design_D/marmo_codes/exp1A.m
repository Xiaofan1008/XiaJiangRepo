function c = exp1A(subject,varargin)
% EXP1A  Contrast x current. Does visual drive gate V1->V2 transmission?
%
% Your noisegridPassive.m, with the noise grid swapped for a drifting grating
% and an estimIntan stimulus added. rigcfg, the inputParser block, the
% addScript hooks, the diode flasher and the design/block/run sequence are
% unchanged.
%
% 24 conditions (6 contrasts x 4 currents) x 50 repeats = 1200 trials.
% Cycle ~2000 ms measured, so ~40 min.
%
% Set the Intan recording directory and file format by hand in the GUI first.
% Nothing about the save path is sent over TCP.
%
%   c = exp1A('CJ286','rfX',-4,'rfY',-3,'diameter',6,'orientation',45, ...
%             'sf',1,'intanChn','A-007');
%
% =========================================================================
% TRIAL STRUCTURE -- shared by exp1A, exp1C and (shifted) exp1B, exp2A, exp2B
% =========================================================================
%   0      grating on            (frame 1)
%   ~242   stim pulse            (frame 17; see the frame block below)
%   500    grating off
%   550    trial ends
%
% The pulse time is quantised to the flip grid and then shifted by K, the
% display-path constant, so "250 ms" is delivered as ~242 ms. The file prints
% the delivered value at startup -- that is the number for the methods.
%
% 550 ms is not arbitrary. The Intan board restarts once per SET and the hole
% opens ~560-650 ms after the trial's first frame (see estimIntan.m). A 550 ms
% trial finishes 10-100 ms before that. Every protocol in this set keeps the
% same 550 ms so the margin is the same everywhere and so the protocols are
% comparable. If "bursts complete" starts failing, drop stimDuration to 400
% rather than changing anything else.
%
% PROVENANCE TAGS
%   [NOISEGRID] your noisegridPassive.m, behaviour unchanged
%   [RUNMULTI]  your RunMultiElect_exp.m, behaviour unchanged
%   [NEW]       written for this experiment
%   [GLUE]      new, exists only so the two systems can meet

aoutChannel = 0;   % [RUNMULTI] Intan trigger: MCC ANALOG out, via DaqAOut

% [NOISEGRID] ARGUMENT HANDLING -- lifted wholesale
if ~exist('subject','var')
    error('No subject name provided.');
end
validateattributes(subject,{'char'},{'nonempty'},'','subject',1);
p = inputParser();
p.KeepUnmatched = true;

% [NEW] Site-specific. All five MUST be set per penetration.
%   rfX/rfY      centre of the aggregate NeuroNexus RF, from the noise map
%   diameter     2 x the aggregate RF, floor 4 deg (protocol)
%   orientation  PREFERRED orientation of the stimulated column, from the
%                tuning run. The default of 0 is a placeholder -- at the wrong
%                orientation the high-contrast conditions do not drive the
%                column and the contrast axis collapses.
%   intanChn     Intan's own name for the stimulating site
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

% [NEW] Grating and timing. Defaults are the protocol values.
p.addParameter('sf',1,@isnumeric);          % cycles/deg. From the tuning run.
p.addParameter('tf',4,@isnumeric);          % Hz. Drift speed is tf/sf deg/s.
p.addParameter('stimDuration',500);         % [ms] grating on
p.addParameter('pulseTime',250);            % [ms] after grating onset. In the
                                            %      sustained response, past the
                                            %      onset transient.
p.addParameter('onsetCorrection',0);        % [ms] SUPERSEDED -- keep at 0.
                                            %      It used to carry a hand
                                            %      correction for the display
                                            %      path. That is now K, applied
                                            %      in frames. Non-zero values
                                            %      still shift the REQUEST, so
                                            %      it is harmless but pointless.

% [NEW] THE ONE TIMING CONSTANT.
%   measured diode->trigger = (F_estim - F_grating)*flip + K
%   F = round(on*frameRate/1000) + 1                        stimulus.m:408
%
% Both events are flip-locked. The trigger follows its flip by the DaqAOut
% write, 2-4 ms, because onsetFunction runs inside afterFlip (stimulus.m:566),
% AFTER Screen('Flip') has returned. The diode crosses threshold ~28 ms after
% the flip -- panel input lag plus pixel rise. K is the difference, so light
% arrives ~25 ms after the flip and the pulse has to be pushed that much later
% to meet it.
%
% K = -24.8 +/- 1.1 ms, from 8 conditions over 3 sessions and 4 protocols.
% Independent confirmation: trigger high measures 99.97 ms against 6 flips =
% 99.995, so both trigger edges carry the same latency.
%
% chk_core reports K on every run. If it moves, something in the display path
% changed -- do not absorb it into onsetCorrection again.
p.addParameter('K',-24.8);

% [NEW] Nominal rate used ONLY to choose the frame count, so that count is
% identical every session. The command that places the pulse uses the rate
% actually measured at runtime. See the frame block below for why.
p.addParameter('nominalFrameRate',60);

% [NEW] Pin the pulse frame explicitly, bypassing the ms -> frame conversion.
% Empty means compute it from pulseTime.
p.addParameter('pulseFrame',[]);

% [NEW] Sub-frame pulse timing. true anchors the pulse to the grating's
% measured flip timestamp and busy-waits the residual, so the delivered time
% is the requested one rather than the nearest flip. false reverts to pure
% frame scheduling.
p.addParameter('subFrame',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

% [NEW] Longest busy-wait allowed inside afterFlip, in ms. After a flip there
% is ~16.7 ms before the next draw deadline and ~3 ms is already spent on
% DaqAOut and the hooks, so ~12 ms is the safe ceiling. A residual longer than
% this puts the pulse on the next flip instead -- at most (ivl - rMax) late,
% but with no risk of dropping a frame. Raise towards 16.7 for exactness at
% the cost of dropped frames; set 0 to never wait.
p.addParameter('rMax',10);   % [ms] see the note above the frame block
p.addParameter('trigDuration',100);         % [ms] trigger line high. Only has
                                            %      to outlast Intan's edge
                                            %      detection; pulse is 260 us.
p.addParameter('iti',1.5);                  % [s] Inert in practice: the
                                            %     postPulsePause plus the two
                                            %     oephys HTTP messages plus the
                                            %     Intan exchange make
                                            %     beforeTrial longer than this,
                                            %     and neurostim subtracts
                                            %     elapsed time from the ITI.
                                            %     Raising it does NOT move the
                                            %     Intan hole -- tested 261016.

% [NEW] The two design axes, exposed so they can be changed without editing.
%
% CURRENTS: 50 and 100 uA were INDISTINGUISHABLE on the Neuropixels in every
% test run where the current actually varied (261015, 261016, 261019), and
% non-monotonic in two of them. The Intan saturates at both so it cannot
% arbitrate. If that is still true at this penetration, use [-1 10 25 50] --
% same 24 conditions, same trial count, but a current axis that spans a range
% instead of two levels that are the same.
p.addParameter('contrasts',[0 0.06 0.12 0.25 0.50 1.00]);   % Michelson
p.addParameter('currents',[-1 10 50 100]);                  % uA, -1 = catch
p.addParameter('nrRepeats',50);             % 50 x 24 conditions = 1200 trials

% [NEW] Bench switches. None change the science.
p.addParameter('noIntan',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('skipSync',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));

p.parse(varargin{:});
args = p.Results;

nCond = numel(args.contrasts)*numel(args.currents);
fprintf('exp1A: %d conditions x %d repeats = %d trials\n', ...
        nCond,args.nrRepeats,nCond*args.nrRepeats);

import neurostim.*
commandwindow;

% [NEW] PTB SYNC WORKAROUND. This machine's beamposition queries are broken,
% so PTB falls back to a timing estimate that sometimes fails its own
% threshold and refuses to open the window. skipSync forces it through.
% Timing guarantees are then weaker -- tolerable only because the photodiode
% is the ground truth. For real data fix the display output mapping instead
% (help DisplayOutputMappings).
if args.skipSync
    Screen('Preference','SkipSyncTests',1);
    warning('exp1A:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
end

% [NOISEGRID] rigcfg builds cic, opens the screen, and creates the oephys
% plugin, the eScript plugin and the diode stimulus.
c = marmolab.rigcfg('debug',args.debug,'subject',subject,p.Unmatched);

% oephys.fakeConnection MUST stay false. Its two blocking HTTP messages sit
% between the stim pulse and estimIntan's parameter_update and form part of
% the PULSE recovery. With them disabled the exchange stalled from 513 ms to
% 1849 ms and the run dropped triggers (261016, 360-trial run).
c.oephys.fakeConnection = false;

% -------------------------------------------------------------------------
% [NEW] THE GRATING
% -------------------------------------------------------------------------
% Screen units on this rig are DEGREES (c.screen.width is 66.15), so every
% spatial field below is already in degrees.
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
fprintf('exp1A aperture: %s mask, sigma %.3f deg, extent +/-%.2f deg\n', ...
        g.mask, g.sigma, args.diameter/2);
if strcmpi(args.edge,'gauss')
    rr = 0:args.diameter/6:args.diameter/2;
    fprintf('  contrast profile r =');  fprintf(' %5.2f',rr);  fprintf(' deg\n');
    fprintf('                  x    ');
    fprintf(' %5.2f',exp(-0.5*(rr/g.sigma).^2));  fprintf('\n');
end
g.frequency   = args.sf;
g.orientation = args.orientation;

% Degrees of PHASE per FRAME. Nothing to do with visual angle.
% 360 deg = one cycle; tf cycles/s / frameRate frames/s = cycles/frame.
% At 4 Hz this is 24, because c.screen.frameRate is the CONTENT rate, 60 Hz.
% The panel refreshes at 120 and repaints each frame twice, but that does not
% enter here: new content only appears once per 16.664 ms. Using 120 would
% halve the drift to 2 Hz.
% Do NOT use gabor.tfToPhaseSpeed; it computes something else.
g.phaseSpeed  = 360*args.tf/c.screen.frameRate;

g.contrast    = NaN;                  % the design writes this every trial.
                                      % NaN so an unassigned value gives a
                                      % blank screen you notice, not a
                                      % plausible wrong contrast.
g.on          = 0;                    % [ms] grating on at trial start
g.duration    = args.stimDuration;    % [ms]

% [NEW] Randomise starting spatial phase every trial. Not cosmetic: at 4 Hz
% one cycle is 250 ms, so a fixed pulse delay would always land at the same
% point in the response cycle. Modulation depth grows with contrast, so that
% offset would grow with contrast too -- manufacturing a contrast-dependent
% interaction that looks exactly like the gating effect.
g.phase = plugins.jitter(c,{0,360},'distribution','uniform');

% -------------------------------------------------------------------------
% [GLUE] THE ELECTRICAL STIMULUS
% -------------------------------------------------------------------------
e = estimIntan(c,'estim');

% -------------------------------------------------------------------------
% [NEW] SUB-FRAME PULSE TIME
% -------------------------------------------------------------------------
% The DAQ line is a separate wire. Nothing forces the pulse onto a flip --
% it sat on the flip grid only because onsetFunction is called from
% afterFlip. So instead of quantising, anchor to the grating's MEASURED flip
% timestamp and fire at an absolute time:
%
%   T = gratingFlipPtb + (wanted - K)
%
% K is daqLatency - panelLag, so using it cancels both and the measured
% diode->trigger comes out at `wanted` exactly, whatever is asked for.
%
% estim is still scheduled on the last flip BEFORE T; its onsetFunction then
% busy-waits the residual rW and fires. The wait blocks the frame loop, and
% after a flip there is only ~16.7 ms before the next draw deadline with ~3 ms
% already spent, so rW is capped at rMax. Past that the pulse goes on the NEXT
% flip instead, which is at most (ivl - rMax) late but drops no frame. A
% dropped frame would be a 24 deg phase stall in the gabor at a fixed time
% relative to every pulse -- a stimulus artifact locked to the manipulation,
% which is worse than a few ms of error.
%
%   subFrame=false  -> old pure-frame behaviour, nf = round(off/ivl)
%   rMax=0          -> always use a flip, never wait
ivl = 1000/c.screen.frameRate;
off = args.pulseTime + args.onsetCorrection - args.K;   % ms after grating FLIP
if args.subFrame
    nf = floor(off/ivl);
    rW = off - nf*ivl;                  % residual in [0,ivl)
    if rW > args.rMax
        nf = nf + 1;  rW = 0;           % too long to block; take the next flip
    elseif rW > args.rMax - 2
        warning('exp1A:rMaxMargin', ...
            ['busy-wait is %.2f ms against an rMax of %.1f. Measured on ' ...
             '231023: 8.1 ms drops no frame, 11.47 ms drops one. A dropped ' ...
             'content frame is a gabor phase stall locked to the pulse, and ' ...
             'the dropped-frames check CANNOT see it (it reads the panel''s ' ...
             '120 Hz ripple, not the content rate). Check the trigger width ' ...
             'against its prediction in the check script.'], rW, args.rMax);
    end
else
    nf = round(off/(1000/args.nominalFrameRate));  rW = 0;
end
if ~isempty(args.pulseFrame), nf = args.pulseFrame; rW = 0; end
assert(nf >= 0, ...
    'exp1A: pulse would be %d frames BEFORE the grating. Raise pulseTime.', nf);

e.on   = nf*ivl;              % g.on = 0, so the grating is frame 1
gOffMs = nf*ivl + rW;         % fire this long after the grating's flip

% Shared with the nested functions below. gTrial guards against firing
% against a stale timestamp if the grating never came on this trial.
gOn  = 0;                     % grating is frame 1 in this protocol
gPtb = NaN; gTrial = -1;

fprintf(['exp1A timing: flip %.4f ms (%.4f Hz measured), K %.1f ms\n' ...
         '  pulse on frame %d + %.2f ms wait, commanded estim.on = %.3f ms\n' ...
         '  lands %.2f ms after VISIBLE onset (asked for %g)\n'], ...
        ivl, c.screen.frameRate, args.K, nf, rW, e.on, gOffMs + args.K, args.pulseTime);
if rW == 0 && args.subFrame
    fprintf('  (residual exceeded rMax %.1f ms, so it went on the next flip)\n', args.rMax);
end

e.duration    = args.trigDuration;
e.intanChn    = args.intanChn;
e.returnChn      = args.returnChn;       % '' = monopolar
e.returnSep      = args.returnSep;       % [um], checked against the 400 floor
e.returnPolarity = args.returnPolarity;  % [] = inverted automatically
e.testMode    = args.noIntan;

% [RUNMULTI] StimParams values, by name instead of column index. These are
% your defaults out of 2_pulse_exp_datafile_.mat.
e.fpa  = NaN;      % varied in the design; -1 = catch. NaN so an unassigned
e.spa  = NaN;      % value fails loudly rather than delivering a wrong current.
e.fpd  = 100;      % us  \ 100 us per phase, 60 us interphase.
e.spd  = 100;      % us  | fpd = spd and fpa = spa is what makes it
e.ipi  = 60;       % us  / charge balanced.
e.pot  = 0;        % single pulse, not a train
e.nsp  = 1;
e.pulseTrainPeriod = 10000;   % us, ignored while pot = 0
e.ptr  = 1000;     % REFRACTORY PERIOD column. parameter_update ignores it and
                   % derives the real one as poAS + 100.
e.stSH = 1;        % biphasic with interphase delay
e.stimPolarity = 0;% cathodic first
e.enAS = 1;        % amp settle on
e.prAS = 200;      % us before the pulse. Also feeds postTriggerDelay.
e.poAS = 160000;   % us after = 160 ms. Kills Intan LFP for that window; MUA
                   % is unaffected.
e.maAS = 1;
e.enCR = 1;        % charge recovery on
e.prCR = 100000;   % us = 100 ms
e.poCR = 150000;   % us = 150 ms
e.triggerSource      = 0;   % \
e.triggerEdgeOrLevel = 0;   %  | column alignment only; parameter_update
e.triggerHighOrLow   = 0;   % /  never sends these three
e.postTriggerDelay   = 0;   % parameter_update adds prAS + 10, so Intan
                            % receives 210 us.

%% [GLUE] Intan TCP -- opened HERE, before c.run covers the console.
if ~args.noIntan
    fprintf('\n---- connect Intan to the TCP/IP server now ----\n');
    e.tcp = local_TCPIP_server();
    fprintf(e.tcp,'NOMCS');
    fprintf('Intan connected.\n');
end

%% [RUNMULTI] Intan trigger. NOT gated on debug -- debug only shrinks the
% window; it must not disable the trigger.
DAQ = DaqDeviceIndex;
if isempty(DAQ)
    error('exp1A: DaqDeviceIndex found no DAQ, so no trigger can be sent.');
end
fprintf('DAQ index %s, trigger on analog out %d\n',mat2str(DAQ),aoutChannel);
DaqAOut(DAQ,aoutChannel,0);          % line starts low

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

% [NOISEGRID] SESSION HOOKS -- unchanged
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

    % [NEW] Nested, so they share gPtb/gTrial/gOffMs with the body above.
    % Nested rather than local on purpose: local functions get their own
    % workspace and could not see the timestamp.
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

    %% Diode Flasher
    % [NOISEGRID] The patch is attached to the GRATING, so the burst marks
    % the stimulus, not the trial. c.diode (from MU00011699) draws to the
    % same corner and has duration Inf, so it must be disabled or the two
    % overlap and the burst length stops meaning anything.
    addDiodeFlasher(c.diode,'enabled',false);
    addDiodeFlasher(g,'size',0.075,'location','nw');

% [NOISEGRID] trialDuration is an @-expression evaluated at runtime.
% stopTime is when the grating goes off (500 ms); the +50 tail gives
% offsetFunction a frame to fire on before the trial is torn down.
c.trialDuration = '@grating.stopTime+50';
c.iti = args.iti*1e3;                 % s -> ms

% -------------------------------------------------------------------------
% [NEW] THE DESIGN -- contrasts x currents, fully crossed
% -------------------------------------------------------------------------
d = neurostim.design('myDesign');
d.fac1.grating.contrast = args.contrasts;

% fpa and spa get the SAME vector, listed twice. fpa = spa keeps the waveform
% charge balanced, so the two must move together. Both in fac2 makes them one
% factor with four levels, not a 4 x 4 cross.
%
% -1 is the catch level: parameter_update turns it into enabled = 0, so the
% trigger still goes out and no current is delivered. Better than 0 uA, which
% would still fire the stimulator and leave a switching artifact.
d.fac2.estim.fpa        = args.currents;
d.fac2.estim.spa        = args.currents;

blk = block('myBlock',d);
blk.nrRepeats = args.nrRepeats;

c.subject = subject;
c.paradigm = 'exp1A';
c.addPropsToInform('grating.contrast','estim.fpa');

% [NEW] setPluginOrder, NOT c.order -- c.order exists in the marmolab fork,
% not in this neurostim (cic.m line 585).
%
% oephys FIRST. Its two blocking HTTP messages then sit between the previous
% trial's pulse and estim's parameter_update. That separation is part of the
% PULSE recovery; putting estim first removes it and the exchange stalls.
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
% [GLUE] Release everything, in an order that works even when the run died
% early. Each step is wrapped individually so one failure cannot prevent the
% others -- if sca throws, the socket must still close.
try; sca; end                                      %#ok<TRYNC>
try; ListenChar(0); end                            %#ok<TRYNC>
try
    if ~isempty(e.daq); DaqAOut(e.daq,e.aoutChannel,0); end
end                                                %#ok<TRYNC>
try
    if ~isempty(e.tcp); fclose(e.tcp); end
end                                                %#ok<TRYNC>
end
