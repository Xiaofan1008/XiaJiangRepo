classdef estimIntan < neurostim.stimulus
% ESTIMINTAN  Electrical stimulus that talks to Intan with YOUR functions.
%
% No neurostim Intan plugin. This class calls, unchanged, the functions on
% your Intan path, in the same order and with the same retry structure as
% RunMultiElect_exp.m. It adds nothing to that exchange.
%
%   beforeExperiment   local_TCPIP_server()  ->  NOMCS  ->  RECORD  ->  pause
%   beforeTrial        postPulsePause
%                      trigger line low
%                      pause 5 s every 500 trials (Intan file indexing)
%                      parameter_update -> checkTCPOK -> SET -> verifyParameters
%                      retried until it succeeds, recoverTCP on failure
%   afterExperiment    pause  ->  STOP  ->  DALL  ->  fclose
%
% =========================================================================
% THE INTAN RECORDING HOLE -- what is known, 261015-261019
% =========================================================================
% Intan stops and restarts the board once per SET. The restart lasts ~650 ms
% and nothing is written while it runs, so the Intan file is a concatenation
% of per-trial segments and time.dat resets at every boundary. It happens on
% catch trials too (94/94 and 12/12 in two separate tests), so it follows the
% SET command, not the delivered current.
%
% The Neuropixels recording is continuous and unaffected, in every run.
%
% WHERE the hole lands, measured as ack -> hole with the ack on the trial's
% first frame:
%
%   config                                   ack -> hole      trial     result
%   no postPulsePause, estim first             357 ms          550      tail lost
%   postPulsePause 0.15, estim last       559/571/576/653      550      0-103 ms margin
%
% So the hole sits AFTER the trial with this configuration, by 9 to 103 ms.
% The margin is not large and it is not controlled -- the centre moved ~95 ms
% across four consecutive runs on 261019 with everything identical. Watch
% "bursts complete" on every session; if it starts failing, shorten
% stimDuration rather than changing anything structural.
%
% WHY THE PAUSE IS HERE AND NOT SOMEWHERE ELSE. It separates the previous
% trial's pulse from this trial's parameter_update, which is the PULSE
% recovery (500 ms) that designMultielectExp_O.m sets and
% RunMultiElect_exp.m enforces with a timer. Without it -- and with the
% oephys HTTP messages also disabled -- the exchange stalled from 513 ms to
% 1849 ms and the run dropped 5 of 360 triggers (261016).
%
% Things that were tried and did NOT help: raising c.iti, reordering the
% plugins, lengthening the trial, moving the exchange into afterFrame (that
% one never fired at all -- neurostim only dispatches a stimulus's frame
% hooks while the stimulus is on, and estim is on for 100 ms).
%
% setTime / setOK / pulseClock are all on cic.clockTime, so differences
% between them are valid. Do NOT difference them against cic.prms.firstFrame,
% which has a different origin.
%
% The -1 sentinel works as it always has: parameter_update sets enabled = 0
% whenever firstPhaseAmplitude or firstPhaseDuration is -1.

    properties
        tcp = [];              % socket from local_TCPIP_server()
        headerPause = 30;      % s after RECORD and before STOP
        indexPause = 5;        % s every 500 trials, as RunMultiElect_exp.m waits
        testMode = false;      % true = no TCP at all, visual + trigger line only

        % [NEW] Skip parameter_update/SET when the parameters are identical
        % to the last ones actually sent. No upload means no board restart
        % and no hole in the Intan file.
        %
        % FALSE in Experiment 1, where the current changes most trials: a
        % skip there would make some trials shorter than others and make
        % trial duration correlate with whether the current changed.
        %
        % TRUE in Experiment 2, where the parameters NEVER change -- the
        % current is fixed and the control trials simply do not fire the
        % trigger. One exchange on trial 1, none after, so the Intan file is
        % continuous and the cycle is limited only by the trial and the ITI.
        skipIfUnchanged = false;

        % [RUNMULTI] The PULSE recovery, in seconds. Top of beforeTrial, so it
        % sits between the previous trial's pulse and this trial's
        % parameter_update. See the header. 0.15 was measured working on
        % 261019; 0.2 is that with a little more room.
        postPulsePause = 0.2;
    end

    methods (Access = public)
        function s = estimIntan(c,name)
            s = s@neurostim.stimulus(c,name);

            % The Intan trigger goes out on the MCC ANALOG out, exactly as
            % RunMultiElect_exp.m does it: DaqAOut(DAQ,0,1) to fire and
            % DaqAOut(DAQ,0,0) to reset. Not a digital line, and not through
            % the neurostim mcc plugin -- mcc.analogOut is a no-op stub.
            s.addProperty('daq',[],'validate',@(x) isnumeric(x) || isempty(x));
            s.addProperty('aoutChannel',0,'validate',@isnumeric);

            % StimParams columns, same order and defaults as your .mat files
            s.addProperty('intanChn','A-007','validate',@ischar);        % 1

            % ---- [25 Oct] BIPOLAR -------------------------------------
            % returnChn  '' = monopolar, the behaviour of every run so far.
            %            Set it to a second Intan channel and the SAME pulse
            %            is uploaded to that channel with the polarity
            %            INVERTED, so current flows between the pair instead
            %            of from one site to a distant ground.
            %
            % The mechanism is RunMultiElect_exp.m lines 109-127: several
            % parameter_update calls, then ONE SET. Two rows, one SET, one
            % board restart -- so this costs one extra upload per trial and
            % does NOT add a second hole.
            %
            % returnSep  [um] centre-to-centre distance between the two
            %            sites. YOU supply it; it is not computed. Channel
            %            number is not position -- the NeuroNexus channel map
            %            is not monotonic (that is what E_Mapnumber is for in
            %            RunMultiElect_exp.m), so deriving a distance from
            %            'A-007' and 'A-015' would be a guess. Below 400 um
            %            the run is refused.
            %
            % returnPolarity  [] = the opposite of stimPolarity, which is what
            %            bipolar means. Only set it by hand to deliberately
            %            drive both sites the same way.
            s.addProperty('returnChn','','validate',@ischar);
            s.addProperty('returnSep',NaN,'validate',@isnumeric);
            s.addProperty('returnPolarity',[], ...
                'validate',@(x) isempty(x) || isnumeric(x));
            s.addProperty('enabled',1,'validate',@isnumeric);            % 2
            s.addProperty('triggerSource',0,'validate',@isnumeric);      % 3
            s.addProperty('triggerEdgeOrLevel',0,'validate',@isnumeric); % 4
            s.addProperty('triggerHighOrLow',0,'validate',@isnumeric);   % 5
            s.addProperty('postTriggerDelay',0,'validate',@isnumeric);   % 6
            s.addProperty('pot',0,'validate',@isnumeric);                % 7
            s.addProperty('nsp',1,'validate',@isnumeric);                % 8
            s.addProperty('pulseTrainPeriod',10000,'validate',@isnumeric); % 9
            s.addProperty('ptr',1000,'validate',@isnumeric);             % 10
            s.addProperty('stSH',1,'validate',@isnumeric);               % 11
            s.addProperty('stimPolarity',0,'validate',@isnumeric);       % 12
            s.addProperty('fpd',100,'validate',@isnumeric);              % 13
            s.addProperty('spd',100,'validate',@isnumeric);              % 14
            s.addProperty('ipi',60,'validate',@isnumeric);               % 15
            s.addProperty('fpa',-1,'validate',@isnumeric);               % 16  -1 = catch
            s.addProperty('spa',-1,'validate',@isnumeric);               % 17
            s.addProperty('enAS',1,'validate',@isnumeric);               % 18
            s.addProperty('prAS',200,'validate',@isnumeric);             % 19
            s.addProperty('poAS',160000,'validate',@isnumeric);          % 20
            s.addProperty('maAS',1,'validate',@isnumeric);               % 21
            s.addProperty('enCR',1,'validate',@isnumeric);               % 22
            s.addProperty('prCR',100000,'validate',@isnumeric);          % 23
            s.addProperty('poCR',150000,'validate',@isnumeric);          % 24

            % ---- DIAGNOSTIC -------------------------------------------
            % setTime    cic.clockTime when verifyParameters returned
            % setOK      1 if Intan confirmed on the first try, 0 if it had
            %            to be resent (then setTime marks the FAILED attempt,
            %            so filter on setOK == 1)
            % pulseClock cic.clockTime at the trigger, written by markPulse
            %
            % pulseClock - setTime = SET-ack to pulse. Add the trigger->hole
            % offset chk_hole reads out of the Intan file and you have
            % SET -> board restart, which is still the unmeasured quantity.
            s.addProperty('setTime',NaN,'validate',@isnumeric);
            s.addProperty('setOK',NaN,'validate',@isnumeric);
            s.addProperty('pulseClock',NaN,'validate',@isnumeric);

            % didSet  1 = an exchange ran this trial, 0 = skipped because the
            %         parameters were already loaded. Always 1 unless
            %         skipIfUnchanged is on.
            s.addProperty('didSet',NaN,'validate',@isnumeric);
        end

        function markPulse(s)
        % Called from the experiment file's onsetFunction, straight after
        % DaqAOut raises the line. Same clock as setTime.
            s.pulseClock = s.cic.clockTime;
        end

        function checkBipolar(s)
        % [25 Oct] Refuse to start a bipolar run that cannot be interpreted.
        % Called from beforeExperiment, so it fails before RECORD rather than
        % after an hour of pairing.
            if isempty(s.returnChn)
                fprintf('estim: MONOPOLAR, stim %s\n', s.intanChn);
                return;
            end
            if strcmpi(strtrim(s.returnChn), strtrim(s.intanChn))
                error('estimIntan:sameChannel', ...
                    ['returnChn (%s) is the stimulating channel. A bipolar ' ...
                     'pair needs two different sites.'], s.returnChn);
            end
            if ~isfinite(s.returnSep)
                error('estimIntan:noSep', ...
                    ['Bipolar run on %s/%s but returnSep was not set. Pass the ' ...
                     'centre-to-centre separation in um. It is not computed ' ...
                     'from the channel names: the NeuroNexus channel map is ' ...
                     'not monotonic, so a number derived from "%s" and "%s" ' ...
                     'would be a guess.'], ...
                     s.intanChn, s.returnChn, s.intanChn, s.returnChn);
            end
            if s.returnSep < 400
                error('estimIntan:tooClose', ...
                    ['Bipolar separation is %.0f um, under the 400 um floor. ' ...
                     'Closer than that the two sites share a population and ' ...
                     'the dipole shunts the volley it is meant to evoke. ' ...
                     'Pick a return site at least 400 um from %s.'], ...
                     s.returnSep, s.intanChn);
            end
            if isempty(s.returnPolarity)
                pol = 1 - (s.stimPolarity ~= 0);
            else
                pol = s.returnPolarity;
                if pol == s.stimPolarity
                    warning('estimIntan:samePolarity', ...
                        ['returnPolarity equals stimPolarity, so BOTH sites ' ...
                         'are driven the same way. That is parallel ' ...
                         'monopolar stimulation, not bipolar.']);
                end
            end
            fprintf(['estim: BIPOLAR  stim %s (polarity %d) / return %s ' ...
                     '(polarity %d), %.0f um apart\n'], ...
                     s.intanChn, s.stimPolarity, s.returnChn, pol, s.returnSep);
            fprintf(['       2 uploads + 1 SET per trial. One board restart, ' ...
                     'so still ONE Intan hole.\n']);
        end

        function beforeExperiment(s)
            checkBipolar(s);
            if s.testMode; return; end

            % The experiment file normally opens the socket before c.run, so
            % the connect happens at the command line rather than behind the
            % stimulus window. Fall back to opening it here if it did not.
            if isempty(s.tcp)
                s.tcp = local_TCPIP_server();
                fprintf(s.tcp,'NOMCS');
            end
            % Set the Intan recording directory and file format by hand in the
            % GUI first. Nothing is sent about the save path.
            fprintf(s.tcp,'RECORD');
            pause(s.headerPause);     % let Intan generate a header
        end

        function beforeTrial(s)
            % [RUNMULTI] PULSE recovery. Must stay BEFORE parameter_update.
            pause(s.postPulsePause);

            % Reset the line for the next pulse, exactly as
            % RunMultiElect_exp.m does at the top of each trial loop.
            if ~isempty(s.daq)
                DaqAOut(s.daq,s.aoutChannel,0);
            end

            if s.testMode; return; end

            if rem(s.cic.trial,500) == 0
                fprintf('Pausing %d seconds for indexing...\n',s.indexPause);
                pause(s.indexPause);
            end

            id = s.cic.trial;
            rows = buildParamSet(s);     % 1 row monopolar, 2 bipolar

            % [NEW] See the property comment. Only ever true in Experiment 2.
            if s.skipIfUnchanged && ~isempty(s.lastSent) && isequal(rows,s.lastSent)
                s.didSet = 0;
                return;
            end

            % [25 Oct] One upload per channel, then ONE SET.
            % RunMultiElect_exp.m wraps `while(~OK)` AROUND the electrode
            % loop and tests OK after each upload, so a failure on channel 1
            % followed by a success on channel 2 leaves OK = 1 and channel 1
            % is never resent. Retry per channel instead.
            uploadRows(s,rows,id);

            fprintf(s.tcp,'SET');
            S = verifyParameters(s.tcp,id);
            s.setTime = s.cic.clockTime;
            s.setOK   = S;
            while ~S
                fprintf('WARNING: COMMUNICATIONS INTERRUPTED. RESENDING PARAMETERS\n');
                uploadRows(s,rows,id);     % ALL rows, not just the first
                fprintf(s.tcp,'SET');
                S = verifyParameters(s.tcp,id);
            end

            s.lastSent = rows;
            s.didSet   = 1;
        end

        function uploadRows(s,rows,id)
        % Upload every channel's row, retrying each one until Intan
        % acknowledges IT, before any SET is sent.
            for r = 1:numel(rows)
                OK = 0;
                while ~OK
                    parameter_update(s.tcp,rows{r},id);
                    OK = checkTCPOK(s.tcp);
                    if OK ~= 1
                        disp('WARNING');
                        recoverTCP(s.tcp);
                        pause(0.01);
                    end
                end
            end
        end

        function afterExperiment(s)
            if ~isempty(s.daq)
                DaqAOut(s.daq,s.aoutChannel,0);
            end
            if s.testMode || isempty(s.tcp); return; end
            pause(s.headerPause);     % closing buffer
            fprintf(s.tcp,'STOP');
            fprintf(s.tcp,'DALL');
            fclose(s.tcp);
            s.tcp = [];
        end
    end

    properties (Access = private)
        lastSent = {};       % last params cell actually uploaded
    end

    methods (Access = protected)
        function rows = buildParamSet(s)
        % [25 Oct] The rows to upload this trial, one per channel.
        % Monopolar gives {stimRow}; bipolar gives {stimRow, returnRow}.
        % They go up back to back and then ONE SET covers both, which is the
        % RunMultiElect_exp.m multi-electrode pattern (lines 109-127).
            rows = {buildParams(s, s.intanChn, s.stimPolarity)};
            if ~isempty(s.returnChn)
                if isempty(s.returnPolarity)
                    pol = 1 - (s.stimPolarity ~= 0);   % 0 <-> 1
                else
                    pol = s.returnPolarity;
                end
                % Everything else is identical, fpa/spa included. Equal
                % amplitude with opposite polarity is what makes the pair a
                % dipole, and it keeps the -1 catch sentinel working on both
                % channels at once: parameter_update turns fpa = -1 into
                % enabled = 0, so a catch trial disables BOTH sites.
                rows{2} = buildParams(s, s.returnChn, pol);
            end
        end

        function params = buildParams(s, chn, pol)
            % One StimParams row. parameter_update reads 15 of these and
            % derives the rest, including enabled, postTriggerDelay and
            % refractoryPeriod. Nothing is validated here on purpose.
            if nargin < 2 || isempty(chn), chn = s.intanChn;      end
            if nargin < 3 || isempty(pol), pol = s.stimPolarity;  end
            params     = cell(1,25);
            params{1}  = chn;
            params{2}  = s.enabled;
            params{3}  = s.triggerSource;
            params{4}  = s.triggerEdgeOrLevel;
            params{5}  = s.triggerHighOrLow;
            params{6}  = s.postTriggerDelay;
            params{7}  = s.pot;
            params{8}  = s.nsp;
            params{9}  = s.pulseTrainPeriod;
            params{10} = s.ptr;
            params{11} = s.stSH;
            params{12} = pol;
            params{13} = s.fpd;
            params{14} = s.spd;
            params{15} = s.ipi;
            params{16} = s.fpa;
            params{17} = s.spa;
            params{18} = s.enAS;
            params{19} = s.prAS;
            params{20} = s.poAS;
            params{21} = s.maAS;
            params{22} = s.enCR;
            params{23} = s.prCR;
            params{24} = s.poCR;
            params{25} = [];   % PULSE (post-pulse recovery) in your .mat.
                               % RunMultiElect_exp.m waits on it with a timer;
                               % here postPulsePause plus the oephys messages
                               % in beforeTrial sit between the pulse and
                               % parameter_update and serve the same purpose.
        end
    end
end
