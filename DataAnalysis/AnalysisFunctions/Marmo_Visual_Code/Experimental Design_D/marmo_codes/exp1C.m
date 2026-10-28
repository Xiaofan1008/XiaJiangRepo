function c = exp1C(subject,varargin)
% EXP1C  Spatial extent. Is the gate local to the stimulated column, or is
% it the extended normalisation pool? This is the block that decides H2.
%
% Same machinery and the same trial as exp1A -- read that file first.
%
% 8 conditions (4 apertures x 2 currents) x 50 repeats = 400 trials, ~13 min.
%
%   c = exp1C('CJ286','rfX',-4,'rfY',-3,'diameter',6,'orientation',45, ...
%             'sf',1,'intanChn','A-007');
%
% =========================================================================
% THE FOUR APERTURES
% =========================================================================
%   disc      centred on the RF, summation diameter
%   annulus   inner 1.5 x diameter, outer full screen -- NO light in the
%             stimulated RF at all
%   full      whole screen
%   blank     no visual stimulus
%
% Divisive normalisation pools over a region far larger than the classical
% RF. If visual drive suppresses transmission only when it falls INSIDE the
% stimulated RF, the gate is local -- shunting at the column, not
% normalisation. If an annulus that puts no light in that RF still
% suppresses, the gate is the extended pool. That is the distinguishing
% test, and it is the reason this block exists.
%
%   disc suppresses, annulus does not     -> local gate, shunting
%   annulus suppresses nearly as much     -> extended pool, H2 proper
%   full suppresses more than either      -> pool sums centre and surround
%
% blank x 50 uA is the pure electrical baseline: ICMS in darkness.
%
% =========================================================================
% HOW THE APERTURES ARE DRAWN
% =========================================================================
% neurostim's gabor has a native ANNULUS mask (gabor.m line 83: sigma is
% [Inner Outer]; GLSLShaders/gabor.frag.txt mask==3 keeps the texture where
% inner <= r <= outer). So all four apertures are one mask type with sigma
% varied, and no occluder or second stimulus is needed:
%
%   disc      sigma = [0                 diameter/2]
%   annulus   sigma = [0.75 x diameter   screenRadius]
%   full      sigma = [0                 screenRadius]
%   blank     sigma = [0                 screenRadius], contrast 0
%
% Inner 0.75 x diameter is the protocol's "inner 1.5 x diameter" expressed
% as a RADIUS, since sigma is a radius and diameter is a diameter.
%
% width and height are the texture window, not the stimulus, and are fixed
% at the SCREEN SIZE for every trial so the texture is never rebuilt
% mid-run. They must be WHOLE NUMBERS and must not exceed the screen:
% gabor passes them to Screen('SetOpenGLTexture') as glWidth/glHeight,
% which rejects a fractional or oversized value.

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

% [NEW] Aperture edge, for the DISC only.
%   'hard'   every condition uses ANNULUS. disc = [0 dR], a one-pixel edge.
%            This is what every run up to 231026 used.
%   'gauss'  the disc uses GAUSS3 with sigma = dR/3, so it fades to nothing at
%            dR instead of stopping dead there. The other three conditions are
%            unchanged.
%
% Only the DISC has a visible outer edge in this protocol. The annulus outer
% radius and the full-field radius are both the screen-corner distance, which
% is off-screen, so there is nothing there to smooth. The annulus INNER edge
% stays hard -- it is an inner edge, and because innerMult is 1.5 there is
% already a deliberate gap between disc and annulus, so the two were never
% complementary and nothing has to tile.
%
% mask is set in gabor.beforeTrial (gabor.m:140), not only at texture
% creation, so it can vary per trial as a design factor.
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
p.addParameter('stimDuration',500);
p.addParameter('pulseTime',250);
p.addParameter('onsetCorrection',0);        % [ms] SUPERSEDED -- keep at 0. See
                                            %      exp1A. Replaced by K, which
                                            %      is applied in frames.
p.addParameter('K',-24.8);                  % [ms] display-path constant.
                                            %      measured diode->trigger
                                            %        = (Fe-Fg)*flip + K
                                            %      Light arrives ~25 ms after
                                            %      the flip (panel lag + pixel
                                            %      rise); the trigger follows
                                            %      its flip by 2-4 ms because
                                            %      onsetFunction runs in
                                            %      afterFlip (stimulus.m:566).
                                            %      -24.8 +/- 1.1 ms over 8
                                            %      conditions. Full derivation
                                            %      in exp1A.
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
p.addParameter('trigDuration',100);
p.addParameter('iti',1.5);
p.addParameter('innerMult',1.5);            % annulus inner edge, x diameter
p.addParameter('screenRadius',[]);          % [deg] outer edge of "full screen".
                                            %      [] = computed from c.screen.
                                            %      Check the printout: if it is
                                            %      not plausible in DEGREES for
                                            %      your display, set it by hand.
p.addParameter('contrast',1.00);            % 100% for every non-blank cell
p.addParameter('current',50);               % [uA]
p.addParameter('nrRepeats',50);             % 50 x 8 = 400 trials

p.addParameter('noIntan',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));
p.addParameter('skipSync',false,@(x) validateattributes(x,{'logical'},{'scalar','nonempty'}));

p.parse(varargin{:});
args = p.Results;

import neurostim.*
commandwindow;

if args.skipSync
    Screen('Preference','SkipSyncTests',1);
    warning('exp1C:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
end

c = marmolab.rigcfg('debug',args.debug,'subject',subject,p.Unmatched);
c.oephys.fakeConnection = false;   % MUST stay false -- see exp1A

% THE TEXTURE. gabor.createProcGabor passes width and height into
% Screen('SetOpenGLTexture') as glWidth/glHeight, so they must be WHOLE
% NUMBERS and must not exceed the window. exp1A gets away with it because
% diameter+1 = 7; any fractional size fails with "Create Procedural Gabor
% Failed ... fractional or outside of allowed range".
%
% So: the texture is the screen rounded down, and the ANNULUS outer radius is
% the farthest screen corner FROM THE RF -- the grating is centred on the
% RF, not on the screen, so that distance is what "out to the edge" means
% here. Anything past the rect is clipped by the texture, which is the
% same thing the screen edge would do.
% INTEGERS. gabor.createProcGabor passes width/height straight into
% Screen('SetOpenGLTexture',...) as glWidth/glHeight, which must be whole
% numbers. 66.1527 x 37.2 is what killed exp1C twice on 261021. floor, not
% round, so the rect never exceeds the screen either.
texW = floor(c.screen.width);
texH = floor(c.screen.height);
if isempty(args.screenRadius)
    cx = [-1 -1 1 1]*texW/2;
    cy = [-1 1 -1 1]*texH/2;
    scrR = max(hypot(cx - args.rfX, cy - args.rfY));
else
    scrR = args.screenRadius;
end
if ~isfinite(scrR) || scrR <= 0 || scrR > 200
    error(['exp1C: computed screen radius %.1f deg is not plausible. ' ...
           'c.screen.width = %g, height = %g. Pass screenRadius by hand.'], ...
           scrR, c.screen.width, c.screen.height);
end
inR  = 0.5*args.innerMult*args.diameter;    % annulus inner RADIUS
dR   = 0.5*args.diameter;                   % disc RADIUS

if inR >= scrR
    error('exp1C: annulus inner radius %.1f deg is outside the screen (%.1f).',inR,scrR);
end

% ---- the 8 conditions, as parallel vectors ------------------------------
% sigma is [inner outer] per trial, so it goes in as a cell array.
apName = {'disc','annulus','full','blank'};
apCon  = [args.contrast, args.contrast, args.contrast, 0];
switch lower(args.edge)
    case 'hard'
        apSig  = {[0 dR], [inR scrR], [0 scrR], [0 scrR]};
        apMask = {'ANNULUS','ANNULUS','ANNULUS','ANNULUS'};
    case 'gauss'
        % GAUSS3 reads sigma1 only and ignores sigma2, but sigma is logged as
        % a pair for every condition so the log stays a clean n x 2 and the
        % check script can read it without ragged rows.
        % NOTE this softens the DISC only. The annulus keeps a hard inner
        % rim: no mask in this shader has a soft inner radius, and a
        % contrast-0 occluder cannot help because alpha is binary in every
        % mask. The annulus inner edge is hard, full stop.
        apSig  = {[dR/3 0], [inR scrR], [0 scrR], [0 scrR]};
        apMask = {'GAUSS3','ANNULUS','ANNULUS','ANNULUS'};
end

sig = {}; con = []; fpa = []; nm = {}; msk = {};
for a = 1:4
    for f = [-1 args.current]
        sig{end+1} = apSig{a};          %#ok<AGROW>
        msk{end+1} = apMask{a};         %#ok<AGROW>
        con(end+1) = apCon(a);          %#ok<AGROW>
        fpa(end+1) = f;                 %#ok<AGROW>
        nm{end+1}  = apName{a};         %#ok<AGROW>
    end
end
fprintf('exp1C aperture edge: %s\n', lower(args.edge));
for a = 1:4
    fprintf('  %-8s %-8s sigma %s\n', apName{a}, apMask{a}, mat2str(apSig{a},4));
end
nCond = numel(fpa);

fprintf('exp1C: %d conditions x %d repeats = %d trials\n', ...
        nCond,args.nrRepeats,nCond*args.nrRepeats);
fprintf('       texture %.1f x %.1f deg (screen); disc r %.1f; annulus r %.1f-%.1f\n', ...
        texW,texH,dR,inR,scrR);

%% grating
g = stimuli.gabor(c,'grating');
g.X           = args.rfX;
g.Y           = args.rfY;
g.mask        = apMask{1};            % the design writes mask AND sigma per
g.sigma       = apSig{1};             % trial. These two must still be VALID
                                      % here: beforeExperiment builds the
                                      % procedural texture before any design
                                      % value has been written, and that is
                                      % where exp1C already failed once.
                                      % Taking them from the disc condition
                                      % keeps the pair self-consistent.
g.width       = texW;                 % texture window = the screen, fixed.
g.height      = texH;                 % NOT 2*scrR -- see above.
g.frequency   = args.sf;
g.orientation = args.orientation;
g.phaseSpeed  = 360*args.tf/c.screen.frameRate;
g.contrast    = NaN;                  % design writes it
g.on          = 0;
g.duration    = args.stimDuration;
g.phase = plugins.jitter(c,{0,360},'distribution','uniform');

%% electrical stimulus
e = estimIntan(c,'estim');

% [NEW] SUB-FRAME PULSE TIME. Same scheme as exp1A -- see that file for the
% full reasoning. The pulse is anchored to the grating's MEASURED flip
% timestamp and fired at T = gratingFlipPtb + (wanted - K), so the delivered
% time is the requested one rather than the nearest flip. estim is scheduled
% on the last flip before T and its onsetFunction busy-waits the residual,
% capped at rMax so the frame loop is never blocked long enough to drop a
% frame.
ivl = 1000/c.screen.frameRate;
off = args.pulseTime + args.onsetCorrection - args.K;   % ms after grating FLIP
if args.subFrame
    nf = floor(off/ivl);
    rW = off - nf*ivl;
    if rW > args.rMax
        nf = nf + 1; rW = 0;
    elseif rW > args.rMax - 2
        warning('exp1C:rMaxMargin', ...
            'busy-wait %.2f ms is close to rMax %.1f; check the trigger width.', ...
            rW, args.rMax);
    end
else
    nf = round(off/(1000/args.nominalFrameRate));  rW = 0;
end
if ~isempty(args.pulseFrame), nf = args.pulseFrame; rW = 0; end
assert(nf >= 0, ...
    'exp1C: pulse would be %d frames BEFORE the grating. Raise pulseTime.', nf);

e.on   = nf*ivl;              % g.on = 0, so the grating is frame 1
gOffMs = nf*ivl + rW;         % fire this long after the grating's flip
gOn  = 0;                     % grating is frame 1 in this protocol
gPtb = NaN; gTrial = -1;      % shared with the nested functions below

fprintf(['exp1C timing: flip %.4f ms, K %.1f ms\n' ...
         '  pulse on frame %d + %.2f ms wait, commanded estim.on = %.3f ms\n' ...
         '  lands %.2f ms after VISIBLE onset (asked %g)\n'], ...
        ivl, args.K, nf, rW, e.on, gOffMs + args.K, args.pulseTime);

e.duration    = args.trigDuration;
e.intanChn    = args.intanChn;
e.returnChn      = args.returnChn;       % '' = monopolar
e.returnSep      = args.returnSep;       % [um], checked against the 400 floor
e.returnPolarity = args.returnPolarity;  % [] = inverted automatically
e.testMode    = args.noIntan;

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
    error('exp1C: DaqDeviceIndex found no DAQ, so no trigger can be sent.');
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

    % The photodiode patch must sit OUTSIDE every aperture. With a
    % full-screen grating the 'nw' corner is inside the stimulus, so the
    % patch is drawn on top of it and still marks onset correctly -- but
    % check on the bench that it is not being modulated by the grating
    % underneath it before trusting the burst edges in exp1C.
    addDiodeFlasher(c.diode,'enabled',false);
    addDiodeFlasher(g,'size',0.075,'location','nw');

c.trialDuration = '@grating.stopTime+50';
c.iti = args.iti*1e3;

%% design -- 8 conditions in one factor, fully interleaved
d = neurostim.design('myDesign');
d.fac1.grating.sigma    = sig;        % cell array: [inner outer] per level
d.fac1.grating.mask     = msk;        % per-condition mask; GAUSS3 for the
                                      % disc when edge = 'gauss'
d.fac1.grating.contrast = con;
d.fac1.estim.fpa        = fpa;
d.fac1.estim.spa        = fpa;

blk = block('myBlock',d);
blk.nrRepeats = args.nrRepeats;

c.subject = subject;
c.paradigm = 'exp1C';
c.addPropsToInform('grating.sigma','grating.mask','grating.contrast','estim.fpa');

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
