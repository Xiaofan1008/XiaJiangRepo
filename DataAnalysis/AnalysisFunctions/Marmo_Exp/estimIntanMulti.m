classdef estimIntanMulti < neurostim.stimulus
% ESTIMINTANMULTI  Multi-electrode version of estimIntan.
%
% Same job as estimIntan.m (talks to Intan with YOUR functions, in the same
% order and with the same retry structure as RunMultiElect_exp.m), but one
% trial can now send SEVERAL StimParams rows, one per electrode, each with
% its own amplitude, post-trigger delay (PTD) and pulse count, followed by
% ONE SET and ONE trigger. This is what sequential / simultaneous stimulation,
% ISI sweeps, A-only / B-only recovery conditions and multi-electrode trains
% need.
%
% STATUS: written but NOT yet run in MATLAB or on the rig. Test with testMode
% first, then with Intan.
%
% =========================================================================
% HOW IT IS USED
% =========================================================================
%   e = estimIntanMulti(c,'estim');
%   e.chnList = {'A-007','A-012'};          % Intan names of ALL channels this
%                                           % run uses. Every one of them is
%                                           % sent on EVERY trial.
%   e.setConditions(cond);                  % see CONDITION TABLE below
%   d.fac1.estim.stimCond = 1:numel(cond);  % the design picks a condition
%
% This class only delivers stimulation. Charge-limit checks and saving the
% condition table belong in the experiment file (or a separate helper).
%
% CONDITION TABLE  (struct array, one element per condition)
%   cond(k).amp   1 x nCh   uA per channel. <= 0 means that channel is OFF
%                           (sent with amplitude -1, so parameter_update sets
%                           enabled = 0). 0 uA is treated as OFF on purpose:
%                           exp1A found -1 better than a 0 uA pulse, which
%                           still fires the stimulator.
%   cond(k).ptd   1 x nCh   post-trigger delay in MICROSECONDS per channel,
%                           written straight into StimParams column 6, as
%                           design_exp_v7 does. (Needs the 1_Shank_flexible
%                           parameter_update on the path -- see below.)
%   cond(k).nsp   1 x nCh   pulses per channel. Optional, default all 1.
%                           nsp = 0 turns that channel OFF (v7 event-level
%                           recovery levels use this).
%   cond(k).label           optional text, e.g. 'A->B 5ms'.
%   Pulse trains use the shared trainPeriod_us (1e6/frequency).
%
% WHY EVERY CHANNEL IS SENT EVERY TRIAL
%   1. Intan may keep an earlier trial's channel enabled if it is not
%      mentioned again. Sending all channels, with unused ones disabled,
%      removes that risk. (A safe default; not verified inside Intan.)
%   2. The exchange then takes the same time for every condition, so trial
%      timing does not depend on the condition (the estimIntan header explains
%      why that matters for the Intan recording hole).
%
% =========================================================================
% WHAT IS UNCHANGED FROM estimIntan
% =========================================================================
%   TCP socket, NOMCS, RECORD / STOP / DALL, postPulsePause, trigger line low
%   at the top of each trial, the 5 s index pause every 500 trials, the
%   setTime / setOK / pulseClock diagnostics, markPulse, testMode, and the
%   trigger itself (fired from the experiment file's onsetFunction via
%   DaqAOut; this class does not fire it).
%
% WHAT IS NEW
%   * chnList + setConditions: one row per channel, built from a condition.
%   * beforeTrial loops parameter_update / checkTCPOK over ALL rows, then ONE
%     SET + verifyParameters. If anything fails, ALL rows are resent. (The
%     resend path in RunMultiElect_exp.m resends only row n+1.)
%   * exchangeMs: how long the whole parameter exchange took, per trial, in the
%     same spirit as setTime. It can only be measured here, inside beforeTrial.
%
% DEPENDENCIES (on the stimulation PC path)
%   parameter_update, checkTCPOK, recoverTCP, verifyParameters,
%   local_TCPIP_server (or the Linux equivalent), DaqAOut.
%   parameter_update MUST be the 1_Shank_flexible version, which uses the
%   per-row postTriggerDelay (params{6}). The CodeSab version ignores column 6
%   and forces 210 us for every row, which would turn every sequential
%   condition into a simultaneous one.

    properties
        tcp = [];              % socket from local_TCPIP_server()
        headerPause = 30;      % s after RECORD and before STOP
        indexPause = 5;        % s every 500 trials
        testMode = false;      % true = no TCP at all
        skipIfUnchanged = false;   % skip the exchange if all rows are identical
                                   % to the last ones sent. Leave false unless
                                   % conditions rarely change (see estimIntan).
        postPulsePause = 0.2;  % s, PULSE recovery, top of beforeTrial

        % ---- channels and conditions ---------------------------------
        chnList = {};          % cellstr of Intan channel names, e.g. {'A-007','A-012'}
        conds = struct('amp',{},'ptd',{},'nsp',{},'label',{});

        % ---- shared waveform and Intan settings ----------------------
        % Each value below is one StimParams column (column number in
        % brackets). Notes marked [code] are stated in your own files;
        % notes marked [inferred] come from the parameter names and have not
        % been checked inside Intan. Times are in MICROSECONDS (us).
        %
        % PULSE WIDTH IS NOT DECIDED YET. 100 us per phase is a placeholder.
        fpd = 100;             % [13] first phase duration, us. This is the cathodic phase when stimPolarity = 0.      
        spd = 100;             % [14] second phase duration, us. Keep = fpd and equal amplitudes, so the pulse is charge    
                               %      balanced.
        ipi = 60;              % [15] interphase delay: gap between the two  phases, us. Only used because stSH = 1.   
        stSH = 1;              % [11] pulse shape. 1 = biphasic WITH an interphase delay [code: parameter_update
                               %      forces stimShape=1, whatever is set here]. 0 = biphasic without delay, 2 = triphasic
                               %      [inferred, not used].
        stimPolarity = 0;      % [12] which phase comes first.  0 = cathodic (negative) first [code: exp1A]
                               %      1 = anodic (positive) first [inferred]
        trainPeriod_us = 16667;% [9] time between pulses of a train, us  (= 1e6 / frequency; 16667 us = 60 Hz,
                               %      3333 us = 300 Hz). Ignored for single pulses.
        ptr = 1000;            % [10] post-stimulation refractory period, us [inferred from the name]. The 1_Shank_
                               %      flexible parameter_update sends it as  postStimRefractory.
        enabledFlag = 1;       % [2] NOT USED: parameter_update decides enabled itself (0 if the amplitude is
                               %      -1, else 1). Kept only so the columns line up with your .mat files.
        triggerSource = 0;     % [3] NOT SENT by parameter_update. Alignment only.
        triggerEdgeOrLevel = 0;% [4] SENT as triggerEdgeOrLevel by parameter_update. 0 is the value your tested files use.
                               %      Its meaning (edge vs level trigger) has  not been checked.
        triggerHighOrLow = 0;  % [5] NOT SENT by parameter_update. Alignment only.
        ptdOffset = 0;         % [6] us added to every channel's post-triggerdelay (PTD). 0 matches design_exp_v7,
                               %      which writes PTD straight into column 6.

        % Amp settle: Intan shorts the amplifier inputs around the pulse so the recording amplifiers recover quickly.
        enAS = 1;              % [18] 1 = amp settle on, 0 = off [code: exp1A]
        prAS = 200;            % [19] amp settle starts this long BEFORE the pulse, us    
        poAS = 160000;         % [20] amp settle stays on this long AFTER the  pulse, us (160000 = 160 ms). Kills the    
                               %      Intan LFP in that window; MUA is not affected [code: exp1A]   
        maAS = 1;              % [21] 1 = keep amp settle on through a whole pulse train [inferred from the name]
      
        % Charge recovery: Intan discharges the electrode after the pulse.
        enCR = 1;              % [22] 1 = charge recovery on, 0 = off [code: exp1A]
        prCR = 100000;         % [23] recovery switches ON this long after the pulse, us (100000 = 100 ms)     
        poCR = 150000;         % [24] recovery switches OFF this long after the  pulse, us (150000 = 150 ms)    
    end

    properties (Access = private)
        lastSent = {};         % last rows actually uploaded (nCh x 25 cell)
    end

    methods (Access = public)
        function s = estimIntanMulti(c,name)
            s = s@neurostim.stimulus(c,name);

            % Trigger on the MCC ANALOG out, as in estimIntan / RunMultiElect_exp
            s.addProperty('daq',[],'validate',@(x) isnumeric(x) || isempty(x));
            s.addProperty('aoutChannel',0,'validate',@isnumeric);

            % The design sets this each trial: index into s.conds
            s.addProperty('stimCond',1,'validate',@isnumeric);

            % ---- DIAGNOSTIC (same meaning as estimIntan) -------------
            s.addProperty('setTime',NaN,'validate',@isnumeric);
            s.addProperty('setOK',NaN,'validate',@isnumeric);
            s.addProperty('pulseClock',NaN,'validate',@isnumeric);
            s.addProperty('didSet',NaN,'validate',@isnumeric);
            % [NEW] ms taken by the whole parameter exchange (all rows +
            % SET + verify), measured with tic/toc on MATLAB's clock.
            s.addProperty('exchangeMs',NaN,'validate',@isnumeric);
        end

        function setConditions(s,cond)
            % Validate and store the condition table.
            assert(~isempty(s.chnList) && iscellstr(s.chnList), ...
                'estimIntanMulti: set chnList (cellstr) before setConditions.'); %#ok<ISCLSTR>
            nCh = numel(s.chnList);
            assert(isstruct(cond) && ~isempty(cond), ...
                'estimIntanMulti: cond must be a non-empty struct array.');
            out = struct('amp',{},'ptd',{},'nsp',{},'label',{});
            for k = 1:numel(cond)
                assert(isfield(cond,'amp') && isfield(cond,'ptd'), ...
                    'estimIntanMulti: each condition needs fields amp and ptd.');
                a = cond(k).amp(:)';  p = cond(k).ptd(:)';
                if isfield(cond,'nsp') && ~isempty(cond(k).nsp)
                    n = cond(k).nsp(:)';
                else
                    n = ones(1,nCh);
                end
                assert(numel(a)==nCh && numel(p)==nCh && numel(n)==nCh, ...
                    'estimIntanMulti: condition %d: amp, ptd and nsp must each have %d entries (one per channel).',k,nCh);
                if isfield(cond,'label') && ~isempty(cond(k).label)
                    lab = cond(k).label;
                else
                    lab = sprintf('cond %d',k);
                end
                out(k).amp = a;  out(k).ptd = p;  out(k).nsp = n;  out(k).label = lab; %#ok<AGROW>
            end
            s.conds = out;
            s.lastSent = {};
        end

        function markPulse(s)
            % Called from the experiment file's onsetFunction, straight after
            % DaqAOut raises the line. Same clock as setTime.
            s.pulseClock = s.cic.clockTime;
        end

        function beforeExperiment(s)
            if s.testMode; return; end

            if isempty(s.tcp)
                s.tcp = local_TCPIP_server();
                fprintf(s.tcp,'NOMCS');
            end
            fprintf(s.tcp,'RECORD');
            pause(s.headerPause);     % let Intan generate a header
        end

        function beforeTrial(s)
            % PULSE recovery. Must stay BEFORE parameter_update.
            pause(s.postPulsePause);

            if ~isempty(s.daq)
                DaqAOut(s.daq,s.aoutChannel,0);   % reset the trigger line
            end

            k = s.stimCond;
            assert(k >= 1 && k <= numel(s.conds) && k == round(k), ...
                'estimIntanMulti: stimCond %g is not a valid condition index.',k);
            rows = buildRows(s,k);

            if s.testMode; return; end

            if rem(s.cic.trial,500) == 0
                fprintf('Pausing %d seconds for indexing...\n',s.indexPause);
                pause(s.indexPause);
            end

            id = s.cic.trial;

            if s.skipIfUnchanged && ~isempty(s.lastSent) && isequal(rows,s.lastSent)
                s.didSet = 0;
                return;
            end

            tExchange = tic;

            % Send ALL rows, then ONE SET. If any row fails, start again from
            % row 1 so no electrode is left with stale parameters.
            sendAllRows(s,rows,id);
            fprintf(s.tcp,'SET');
            S = verifyParameters(s.tcp,id);
            s.setTime = s.cic.clockTime;
            s.setOK   = S;
            while ~S
                fprintf('WARNING: COMMUNICATIONS INTERRUPTED. RESENDING PARAMETERS\n');
                sendAllRows(s,rows,id);
                fprintf(s.tcp,'SET');
                S = verifyParameters(s.tcp,id);
            end

            s.exchangeMs = toc(tExchange)*1000;
            s.lastSent = rows;
            s.didSet   = 1;
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

    methods (Access = protected)
        function sendAllRows(s,rows,id)
            % parameter_update + checkTCPOK for every row. A failure on any
            % row triggers recoverTCP and a restart from the first row.
            allOK = false;
            while ~allOK
                allOK = true;
                for i = 1:size(rows,1)
                    parameter_update(s.tcp,rows(i,:),id);
                    OK = checkTCPOK(s.tcp);
                    if OK ~= 1
                        disp('WARNING');
                        recoverTCP(s.tcp);
                        pause(0.01);
                        allOK = false;
                        break;
                    end
                end
            end
        end

        function rows = buildRows(s,k)
            % nCh x 25 cell: one StimParams row per channel, same column
            % order as your .mat files / estimIntan.
            c   = s.conds(k);
            nCh = numel(s.chnList);
            rows = cell(nCh,25);
            for i = 1:nCh
                off = ~(c.amp(i) > 0) || ~(c.nsp(i) >= 1);
                if off
                    amp = -1;  nsp = 1;  pot = 0;  ptd = 0;   % -1 => enabled = 0
                else
                    amp = c.amp(i);
                    nsp = c.nsp(i);
                    pot = double(nsp > 1);                    % 0 pulse, 1 train
                    ptd = c.ptd(i) + s.ptdOffset;
                end
                rows{i,1}  = s.chnList{i};
                rows{i,2}  = s.enabledFlag;
                rows{i,3}  = s.triggerSource;
                rows{i,4}  = s.triggerEdgeOrLevel;
                rows{i,5}  = s.triggerHighOrLow;
                rows{i,6}  = ptd;               % post-trigger delay, us
                rows{i,7}  = pot;               % 0 = single pulse, 1 = pulse train
                rows{i,8}  = nsp;               % pulses in the train (1 for a single pulse)
                rows{i,9}  = s.trainPeriod_us;
                rows{i,10} = s.ptr;
                rows{i,11} = s.stSH;
                rows{i,12} = s.stimPolarity;
                rows{i,13} = s.fpd;
                rows{i,14} = s.spd;
                rows{i,15} = s.ipi;
                rows{i,16} = amp;               % first phase amplitude, uA; -1 = channel off
                rows{i,17} = amp;               % second phase amplitude, uA (= first: balanced)
                rows{i,18} = s.enAS;
                rows{i,19} = s.prAS;
                rows{i,20} = s.poAS;
                rows{i,21} = s.maAS;
                rows{i,22} = s.enCR;
                rows{i,23} = s.prCR;
                rows{i,24} = s.poCR;
                rows{i,25} = [];
            end
        end
    end
end