function c = expVisualType(subject,varargin)
% EXPVISUALTYPE  Visual-only comparison of stimulus TYPES: bright spot, dark
% spot and Gabor (several orientations), crossed with size and contrast. No
% electrical stimulation, no Intan, no trigger: only the screen and Open Ephys.
% The recorded spike statistics (Fano factor, count distribution, PSTH shape)
% are compared with those of the electrical runs afterwards, to see which
% visual stimulus looks most like the electrical response.
%
% STATUS: written but NOT yet run in MATLAB or on the rig. Test with
% 'debug',true first.
%
% Built on exp1A.m (rig setup, diode, timing) and the visual part of the
% expTrain* files (same gabor object, same onset time, same duration).
%
%   c = expVisualType('CJ286');
%   c = expVisualType('CJ286','diameter',[0.5 1 6],'rfX',2,'rfY',-1, ...
%                     'contrast',[0.25 0.5 1],'nReps',30);
%
% The first input is the SUBJECT (animal ID, a label only, as in exp1A).
% Only the SUBJECT must be given.
%
% STANDALONE: this one file holds the experiment and the condition builder.
% It needs only what the stimulation PC already has (neurostim, marmolab,
% Psychtoolbox).
%
% =========================================================================
% INPUT LABELS  (name-value pairs, e.g. 'contrast',[0.25 1]). Every label has
% a default.
%
%   label          meaning                                   default
%   -------------  ----------------------------------------  -------------
%  STIMULUS TYPES
%   'visType'      'spot' | 'gabor' | 'both'                 'both'
%   'polarity'     SPOT only: 'bright' | 'dark' | 'both'    'both'
%   'ori'          GABOR orientations, deg (a list)          [0 45 90 135]
%                  The spot has no orientation.
%   'sf'           gabor spatial frequency, cycles/deg       1
%   'gaborPhase'   gabor phase, deg (0 = sine, 90 = bright centre)   0
%  SIZE, CONTRAST, POSITION
%   'diameter'     patch diameter LIST, deg                  6
%                  (every size is run; each extra size doubles the
%                   number of conditions. Smaller sizes (0.25-0.5 deg) are
%                   suggested by the ICMS literature: set them from RF mapping)
%   'contrast'     contrast LEVELS, 0..1                     [0.06 0.12 0.25 0.5 1]
%   'visDur'       visual duration, ms (rounded to frames)   250
%   'rfX','rfY'    deg, centre of the patch. PLACEHOLDER 0, 0:   0, 0
%                  set from RF mapping (a warning prints if both are 0)
%  TIMING
%   'pulseTime'    ms from trial start to the patch onset    50
%                  (the same value as the electrical files, so the visual
%                   onset has the same timing in every experiment)
%   'trialDur'     trial length, ms                          500
%   'iti'          inter-trial interval, s (inert, exp1A)    1.5
%   'nReps'        repeats of every condition                20
%  BENCH SWITCHES (none change the science)
%   'debug'        rig debug mode                            false
%   'skipSync'     Psychtoolbox SkipSyncTests                false
%   'confirm'      print a summary and ask before running    true
%
% =========================================================================
% CONDITIONS (built by the local function buildVisConds below)
% =========================================================================
% One condition for every combination of
%   stimulus variant x diameter x contrast
% where the variants are
%   spot bright      (if visType is 'spot' or 'both' and polarity includes bright)
%   spot dark        (if polarity includes dark)
%   gabor, one per orientation in 'ori'  (if visType is 'gabor' or 'both')
% plus ONE blank condition (nothing shown; contrast 0).
% Defaults: (1 bright spot + 1 dark spot + 4 gabor orientations) = 6 variants
% x 1 diameter x 5 contrasts = 30 conditions, + 1 blank = 31 conditions,
% 20 repeats, about 21 min. Each extra diameter doubles the 30 (about +20 min);
% each extra contrast adds 6 (about +4 min); each extra orientation adds 5.
% THE SESSION LENGTH IS PRINTED BEFORE YOU CONFIRM.
%
% Spot = the same gabor object with spatial frequency 0 (a disc of uniform
% luminance); bright = phase 90, dark = phase 270. Gabor = grating in the SAME
% circular window, same position, size, onset and duration as the spot.
% PHASE (checked in the shader source, GLSLShaders/gabor.frag.txt): luminance =
% background x (1 + contrast x sin(2 pi f x + phase)). A spot has f = 0, so its
% pattern is sin(phase): bright = 90, dark = 270. Phase 0 or 180 would be
% INVISIBLE (sin = 0). The gabor keeps 0 / 180 (a grating and its inverse).
%
% The design writes contrast, phase, spatial frequency, orientation and the
% window size (sigma) every trial. The window size changing per trial is
% UNVERIFIED on the rig: if the patch size does not change between
% conditions, run one diameter per call instead.
%
% =========================================================================
% TRIAL STRUCTURE
% =========================================================================
%   0 ms        trial starts
%   pulseTime   patch appears (default 50 ms), for visDur (default 250 ms)
%   trialDur    trial ends (default 500 ms)
% Align rasters and PSTHs to the PHOTODIODE flash. There is no electrical
% stimulation and no Intan, so the Intan hole and the header pause do not
% apply here.
%
% =========================================================================
% PROVENANCE TAGS (as in exp1A)
%   [EXP1A] copied from exp1A.m, behaviour unchanged
%   [NEW]   written for this experiment
% =========================================================================

% [EXP1A] argument handling
if ~exist('subject','var')
    error('No subject name provided.');
end
validateattributes(subject,{'char'},{'nonempty'},'','subject',1);
p = inputParser();
p.KeepUnmatched = true;

% [NEW] stimulus types
p.addParameter('visType','both',@(x) ismember(lower(x),{'spot','gabor','both'}));
p.addParameter('polarity','both',@(x) ismember(lower(x),{'bright','dark','both'}));
p.addParameter('ori',[0 45 90 135],@(x) isnumeric(x) && ~isempty(x));
p.addParameter('sf',1,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('gaborPhase',0,@(x) isnumeric(x) && isscalar(x));
% [NEW] size, contrast, position
p.addParameter('diameter',6,@(x) isnumeric(x) && ~isempty(x) && all(x>0));
p.addParameter('contrast',[0.06 0.12 0.25 0.5 1],@(x) isnumeric(x) && ~isempty(x) && all(x>0) && all(x<=1));
p.addParameter('visDur',250,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('rfX',0,@isnumeric);
p.addParameter('rfY',0,@isnumeric);
% [NEW] timing
p.addParameter('pulseTime',50,@(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('trialDur',500,@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('iti',1.5,@isnumeric);
p.addParameter('nReps',20,@(x) isnumeric(x) && isscalar(x) && x>=1);
% [NEW] bench switches
p.addParameter('debug',false,@islogical);
p.addParameter('skipSync',false,@islogical);
p.addParameter('confirm',true,@islogical);

p.parse(varargin{:});
args = p.Results;

if args.pulseTime + args.visDur > args.trialDur
    error('expVisualType:timing','The visual patch (%g ms from %g ms) does not fit in trialDur %g ms.', ...
        args.visDur,args.pulseTime,args.trialDur);
end
if args.rfX == 0 && args.rfY == 0
    warning('expVisualType:rf','rfX and rfY are both 0: this is a PLACEHOLDER. Set them from RF mapping.');
end

% [NEW] conditions
vis = struct('type',lower(args.visType),'polarity',lower(args.polarity), ...
             'ori',args.ori(:)','sf',args.sf,'gaborPhase',args.gaborPhase, ...
             'diameter',args.diameter(:)','contrast',args.contrast(:)');
cond = buildVisConds(vis);
nCond = numel(cond);
nTrials = nCond * args.nReps;

fprintf('\nexpVisualType: %d conditions x %d repeats = %d trials (about %.0f min at ~2 s per trial)\n', ...
    nCond,args.nReps,nTrials,nTrials*2/60);
fprintf('  %d visual + %d blank conditions\n', ...
    sum(~strcmp({cond.kind},'blank')),sum(strcmp({cond.kind},'blank')));
fprintf('  types %s, polarity %s, gabor orientations [%s] deg, sf %g c/deg\n', ...
    vis.type,vis.polarity,num2str(vis.ori),vis.sf);
fprintf('  diameters [%s] deg, contrasts [%s], %g ms from %g ms, at (%g,%g)\n', ...
    num2str(vis.diameter),num2str(vis.contrast),args.visDur,args.pulseTime,args.rfX,args.rfY);
for k = 1:nCond
    fprintf('  %3d  %s\n',k,cond(k).label);
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
    warning('expVisualType:skipSync','SkipSyncTests is ON. Flip timestamps are unreliable.');
end

% [EXP1A] rigcfg builds cic, opens the screen, and creates the oephys plugin
% and the diode stimulus.
c = marmolab.rigcfg('debug',args.debug,'subject',subject,p.Unmatched);

% [EXP1A] oephys.fakeConnection MUST stay false.
c.oephys.fakeConnection = false;

% -------------------------------------------------------------------------
% [NEW] THE VISUAL PATCH: one gabor object for spot and gabor alike. Spot =
% spatial frequency 0. The window is a circle whose radius is sigma. The
% object is made big enough for the LARGEST diameter; the design shrinks the
% window per trial (sigma). Contrast 0 hides the patch (blank condition).
% -------------------------------------------------------------------------
dMax = max(vis.diameter);
ivl = 1000/c.screen.frameRate;
nf  = round(args.pulseTime/ivl);

g = stimuli.gabor(c,'flash');
g.X = args.rfX;
g.Y = args.rfY;
g.mask   = 'CIRCLE';
g.sigma  = dMax/2;               % CIRCLE: sigma is the radius; the design overrides it
g.width  = ceil(dMax) + 1;
g.height = ceil(dMax) + 1;
g.orientation = 0;               % the design overrides this per trial
g.phaseSpeed  = 0;
g.frequency   = 0;               % the design overrides this per trial
g.phase       = 0;               % the design overrides this per trial
g.contrast    = NaN;             % the design writes this every trial
g.on       = nf*ivl;             % frame-grid onset
g.duration = args.visDur;
fprintf('expVisualType timing: patch requested at %g ms, on frame %d, shown at %.2f ms after trial start\n', ...
    args.pulseTime,nf+1,g.on);

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

% [EXP1A] diode flasher: c.diode is disabled; the flasher is attached to the
% patch so it marks the visual onset.
addDiodeFlasher(c.diode,'enabled',false);
addDiodeFlasher(g,'size',0.075,'location','nw');

c.trialDuration = args.trialDur;      % ms
c.iti = args.iti*1e3;                 % s -> ms

% -------------------------------------------------------------------------
% [NEW] DESIGN: one factor, the condition. Each property has one entry per
% condition.
% -------------------------------------------------------------------------
d = neurostim.design('myDesign');
d.fac1.flash.contrast    = [cond.contrast];
d.fac1.flash.phase       = [cond.phase];
d.fac1.flash.frequency   = [cond.freq];
d.fac1.flash.orientation = [cond.ori];
d.fac1.flash.sigma       = [cond.diameter]/2;
% The patch, and so its photodiode flash, is on for visual trials only; the
% blank trial shows nothing (on = Inf means 'never on').
visOn = inf(1,nCond);
visOn(strcmp({cond.kind},'visual')) = g.on;
d.fac1.flash.on          = visOn;

blk = block('myBlock',d);
blk.nrRepeats = args.nReps;

c.subject  = subject;
c.paradigm = 'expVisualType';

% [NEW] keep the condition table with the data
stamp = datestr(now,'yyyymmdd_HHMMSS');
condFile = sprintf('%s_expVisualType_conds_%s.mat',subject,stamp);
save(condFile,'cond','args');
fprintf('Condition table saved to %s\n',fullfile(pwd,condFile));

% [EXP1A] plugin order: oephys first
c.setPluginOrder('oephys','flash');

try
    c.run(blk);
catch me
    try; sca; end                                  %#ok<TRYNC>
    try; ListenChar(0); end                        %#ok<TRYNC>
    rethrow(me);
end
try; sca; end                                      %#ok<TRYNC>
try; ListenChar(0); end                            %#ok<TRYNC>

end

% =========================================================================
% LOCAL FUNCTIONS (so this file is standalone)
% =========================================================================

function cond = buildVisConds(vis)
% BUILDVISCONDS  Visual-only conditions: every stimulus variant x diameter x
% contrast, plus one blank.
%
%   cond = buildVisConds(vis)
%
% vis   struct: type ('spot'|'gabor'|'both'), polarity ('bright'|'dark'|
%       'both', spot only), ori (list, gabor only), sf (cycles/deg),
%       gaborPhase (deg), diameter (list, deg), contrast (list)
%
% cond  struct array, one entry per condition, fields
%       kind ('visual'|'blank'), stim ('spot'|'gabor'|'blank'),
%       polarity ('bright'|'dark'|'-'), ori (deg; 0 for spot and blank),
%       diameter (deg; the largest one for blank, hidden), contrast, phase, freq (spatial
%       frequency; 0 for spot), label.
%       The blank has contrast 0, so the patch is hidden.

% stimulus variants: stim, polarity, ori, phase, freq
variants = struct('stim',{},'polarity',{},'ori',{},'phase',{},'freq',{});
if any(strcmp(vis.type,{'spot','both'}))
    if any(strcmp(vis.polarity,{'bright','both'}))
        variants(end+1) = struct('stim','spot','polarity','bright','ori',0,'phase',90,'freq',0);
    end
    if any(strcmp(vis.polarity,{'dark','both'}))
        variants(end+1) = struct('stim','spot','polarity','dark','ori',0,'phase',270,'freq',0);
    end
end
if any(strcmp(vis.type,{'gabor','both'}))
    for o = vis.ori
        variants(end+1) = struct('stim','gabor','polarity','-','ori',o, ...
            'phase',vis.gaborPhase,'freq',vis.sf);
    end
end

cond = struct('kind',{},'stim',{},'polarity',{},'ori',{},'diameter',{}, ...
              'contrast',{},'phase',{},'freq',{},'label',{});
for v = 1:numel(variants)
    for dm = vis.diameter
        for ct = vis.contrast
            k = numel(cond)+1;
            cond(k).kind     = 'visual';
            cond(k).stim     = variants(v).stim;
            cond(k).polarity = variants(v).polarity;
            cond(k).ori      = variants(v).ori;
            cond(k).diameter = dm;
            cond(k).contrast = ct;
            cond(k).phase    = variants(v).phase;
            cond(k).freq     = variants(v).freq;
            if strcmp(variants(v).stim,'spot')
                cond(k).label = sprintf('spot %s %gdeg c%g',variants(v).polarity,dm,ct);
            else
                cond(k).label = sprintf('gabor ori%g %gdeg c%g',variants(v).ori,dm,ct);
            end
        end
    end
end

% one blank, nothing shown
k = numel(cond)+1;
cond(k).kind = 'blank';  cond(k).stim = 'blank';  cond(k).polarity = '-';
cond(k).ori = 0;  cond(k).diameter = max(vis.diameter);  cond(k).contrast = 0;
cond(k).phase = 0;  cond(k).freq = 0;  cond(k).label = 'blank';
end
