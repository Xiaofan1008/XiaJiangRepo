function diagExp1A(stage,dbg)
% DIAGEXP1A  Stage-by-stage diagnostic. Never leaves the screen open.
%
%   diagExp1A(1)   rig config only, NO window opens at all  <- start here
%   diagExp1A(2)   gabor, fixed contrast, no design, 3 trials
%   diagExp1A(3)   gabor + design on contrast, 6 trials
%   diagExp1A(4)   + estimIntan, no TCP and no DAQ
%   diagExp1A(5)   + DAQ trigger, still no Intan
%
% Every stage is wrapped so that whatever happens, the screen is closed and
% the keyboard released BEFORE anything is printed. Any error is written to
% diagExp1A_error.txt in the current folder -- send me that file.
%
% Second argument is debug mode, default true (small window, less likely to
% take over the display). diagExp1A(2,false) for full screen.
%
% If MATLAB ever does lock up: Ctrl-C, then type  sca  and press enter even
% if you cannot see what you are typing.

if nargin < 1; stage = 1; end
if nargin < 2; dbg = true; end

err = [];
c = [];
try
    switch stage
        case 1
            c = marmolab.rigcfg('debug',dbg,'subject','diag');
            report(c);
            fprintf('\nStage 1 done. No window was opened.\n');
            return                                  % nothing to close

        case 2
            c = baseRig(dbg);
            g = gabor(c);
            g.contrast = 1;                         % fixed, no design
            c.trialDuration = 520; c.iti = 500;
            d = neurostim.design('d'); d.fac1.grating.on = 0;   % dummy factor
            run3(c,d,3);

        case 3
            c = baseRig(dbg);
            gabor(c);                               % contrast set by design
            c.trialDuration = 520; c.iti = 500;
            d = neurostim.design('d');
            d.fac1.grating.contrast = [0 0.25 1];
            run3(c,d,2);

        case 4
            c = baseRig(dbg);
            gabor(c);
            e = estimIntan(c,'estim');
            e.on = 250; e.duration = 100;
            e.testMode = true;                      % no TCP
            e.daq = [];                             % no DAQ
            c.trialDuration = 520; c.iti = 500;
            d = neurostim.design('d');
            d.fac1.grating.contrast = [0 0.25 1];
            d.fac2.estim.fpa = [-1 50];
            d.fac2.estim.spa = [-1 50];
            c.setPluginOrder('grating','estim');
            run3(c,d,2);

        case 5
            c = baseRig(dbg);
            gabor(c);
            e = estimIntan(c,'estim');
            e.on = 250; e.duration = 100;
            e.testMode = true;                      % still no TCP
            DAQ = DaqDeviceIndex;
            fprintf('DaqDeviceIndex = %s\n',mat2str(DAQ));
            e.daq = DAQ; e.aoutChannel = 0;
            e.onsetFunction  = @(o,t) DaqAOut(DAQ,0,1);
            e.offsetFunction = @(o,t) DaqAOut(DAQ,0,0);
            c.trialDuration = 520; c.iti = 500;
            d = neurostim.design('d');
            d.fac1.grating.contrast = [0 1];
            d.fac2.estim.fpa = [-1 50];
            d.fac2.estim.spa = [-1 50];
            c.setPluginOrder('grating','estim');
            run3(c,d,2);
    end
catch err %#ok<CTCH>
end

% ---- always clean up BEFORE printing anything --------------------------
try; sca; end                      %#ok<TRYNC>  close all PTB screens
try; ListenChar(0); end            %#ok<TRYNC>  give the keyboard back
try; ShowCursor; end               %#ok<TRYNC>
try; Priority(0); end              %#ok<TRYNC>

if isempty(err)
    fprintf('\nStage %d completed with no error.\n',stage);
    return
end

fid = fopen('diagExp1A_error.txt','w');
out = @(varargin) cellfun(@(f) fprintf(f,varargin{:}),{1,fid});
out('diagExp1A stage %d, debug=%d, %s\n\n',stage,dbg,datestr(now));
out('%s\n\n',err.message);
for k = 1:numel(err.stack)
    out('  %s  (line %d)\n',err.stack(k).name,err.stack(k).line);
end
fclose(fid);
fprintf('\nError written to diagExp1A_error.txt -- send me that file.\n');
end

% ------------------------------------------------------------------------
function c = baseRig(dbg)
c = marmolab.rigcfg('debug',dbg,'subject','diag');
report(c);
end

function g = gabor(c)
g = neurostim.stimuli.gabor(c,'grating');
g.mask = 'CIRCLE';
g.sigma = 3; g.width = 7; g.height = 7;
g.frequency = 0.8;
g.orientation = 45;
g.phaseSpeed = 360*4/c.screen.frameRate;
g.contrast = NaN;
g.on = 0; g.duration = 500;
end

function run3(c,d,reps)
b = neurostim.block('b',d);
b.nrRepeats = reps;
c.subject = 'diag';
c.paradigm = 'diag';
c.run(b);
end

function report(c)
fprintf('\n---- rig ----\n');
fprintf('screen.number : %s\n',mat2str(c.screen.number));
fprintf('pixels        : %d x %d\n',c.screen.xpixels,c.screen.ypixels);
fprintf('width/height  : %g x %g  (screen units)\n',c.screen.width,c.screen.height);
fprintf('frameRate     : %g Hz\n',c.screen.frameRate);
fprintf('colorMode     : %s\n',c.screen.colorMode);
if ~isempty(c.plugins)
    fprintf('plugins       : %s\n',strjoin({c.plugins.name},', '));
else
    fprintf('plugins       : (none)\n');
end
if ~isempty(c.stimuli)
    fprintf('stimuli       : %s\n',strjoin({c.stimuli.name},', '));
else
    fprintf('stimuli       : (none)\n');
end
fprintf('PTB screens   : %s\n',mat2str(Screen('Screens')));
fprintf('-------------\n\n');
end
