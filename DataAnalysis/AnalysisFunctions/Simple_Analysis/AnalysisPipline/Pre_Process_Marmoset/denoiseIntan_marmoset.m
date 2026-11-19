function denoiseIntan_marmoset(filepath,dName,T,par,varargin)
% DENOISEINTAN_MARMOSET  Copy of denoiseIntan_sab (2026-11-03). Only change:
% it calls simpleBlank_PTD_PulseTrain_Marmoset instead of ..._Xia_v1, so the
% marmoset recordings (long trial spacing) do not crash at a chunk end.
% Output file names are unchanged (<dName>_dn_sab.dat).
% Chunk overlap: instead of 1 s, 'carry' samples (>= shifttime + 3000) are
% kept between chunks, and abs0 / isLast / carry are passed to the blanking
% function so every trial is handled exactly once (see simpleBlank_PTD_PulseTrain_Marmoset).
dbstop if error
%% This function denoises an Intan datafile
%% Initial Setup
if size(varargin,2)==0
    Startpoint_analyse=0;
    Overall_time_to_analyse=0;
else
    Startpoint_analyse=cell2mat(varargin(1));
    Overall_time_to_analyse=cell2mat(varargin(2));
end
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;
if nChn>32
    E_Mapnumber=1;
else
    E_Mapnumber=0;
end
shortbytes=2;
%% Deal with the digital lines
if isempty(dir('*.trig.dat'))
cleanTrig_sabquick;
end
%% Blank stimulation artefact
trig = loadTrig(0);
theseTrig = trig;
% MARMOSET: chunk overlap in samples. shifttime is the same quantity the
% blanking function uses (min trial spacing - 0.101 s); add 3000 samples margin.
okTrig = sort(theseTrig(theseTrig > 0));
if numel(okTrig) > 1
    shiftSamp = min(diff(okTrig)) - 0.101*FS;
else
    shiftSamp = 0;
end
carry = max(FS, ceil(shiftSamp) + 3000);
carry = min(carry, floor(FS*T/2));      % never more than half a chunk
abs0 = 0;                               % absolute sample index of v(:,1) (0-based)
if isempty(dir([dName '_dn_sab.dat']))
    if ~isempty(theseTrig)
        info = dir([filepath filesep dName '.dat']);
        info = info.bytes/2;
        nL = (ceil(info / (nChn*FS*double(T)))+1);
        %lvFID = fopen([filepath '\' dName '.dat'],'r');
        vFID = fopen([filepath filesep dName '.dat'],'r');
        vdnFID = fopen([filepath filesep dName '_dn_sab.dat'],'W');
        %lbFID = fopen([filepath '\' dName '_lfpraw.dat'],'w');
        N = 1;
        CHK = split(filepath,'_');
        blank = 1;
        for i = 1:size(CHK,1)
            if strcmp(CHK{i},'CSD') || strcmp(dName,'analogin') || isempty(dir('*exp_datafile*.mat'))
                blank = 0;
            end
        end
        dispstat('','init');
        if (blank)
            dispstat(sprintf('Blanking artefact . . .'),'keepthis','n');
        else
            dispstat(sprintf('Skipping blanking . . .'),'keepthis','n');
        end
        BREAK = 1;
        while (N) && (BREAK)
            if strcmp(dName,'analogin')
                v = fread(vFID,[nChn, (FS * T)],'uint16');
                v = (v - 32768) .* 0.0003125;
            elseif strcmp(dName,'amplifier')
                if Startpoint_analyse~=0
                    offset=Startpoint_analyse*FS*nChn*shortbytes;%offset from beginning of file
                    ftell(vFID)
                    fseek(vFID,offset,'bof');
                    ftell(vFID)
                end
                v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195;
            end
            if ~size(v,2)
                BREAK = 0;
            end
            dispstat(sprintf('Progress %03.2f%%',100*((N-1)/nL)),'timestamp');
            % Append saved data to the start of the next loop
            if (N ~= 1)
                v = [hv v]; %#ok<*AGROW>
            end
            if (blank)
                if (par)
                    parfor iChn = 1:nChn
                        % v(iChn,:) = simpleBlank(v(iChn,:),N,T,theseTrig,1,FS);
                        % ---- 2 Electrode Seq stim ----
                        % v(iChn,:) = simpleBlank_PostTriggerDelay_Xia(v(iChn,:),N,T,theseTrig,1,FS); 
                
                        % ---- Multi Electrode Seq stim ----
                        % v(iChn,:) = simpleBlank_PTD_MultiElec_Xia_v1(v(iChn,:),N,T,theseTrig,1,FS);
                        % v(iChn,:) = simpleBlank_PTD_MultiElec_Xia_v2(v(iChn,:),N,T,theseTrig,1,FS);

                        % ---- Pulse Train -----
                        v(iChn,:) = simpleBlank_PTD_PulseTrain_Marmoset(v(iChn,:),N,T,theseTrig,1,FS,abs0,~BREAK,carry);
                    end
                else
                    d = Depth(E_Mapnumber);
                    for iChn = 1:nChn
                        if d(iChn) == 13
                            pause(0.01);
                        end
                        % v(iChn,:) = simpleBlank(v(iChn,:),N,T,theseTrig,1,FS);
                        % ---- 2 electrode stim ----
                        % v(iChn,:) = simpleBlank_PostTriggerDelay_Xia(v(iChn,:),N,T,theseTrig,1,FS);
                        % ---- Multi Electrodes Seq Stim ---- 
                        % v(iChn,:) = simpleBlank_PTD_MultiElec_Xia_v1(v(iChn,:),N,T,theseTrig,1,FS);
                        % v(iChn,:) = simpleBlank_PTD_MultiElec_Xia_v2(v(iChn,:),N,T,theseTrig,1,FS);
                    
                        % ---- Pulse Train ----
                        v(iChn,:) = simpleBlank_PTD_PulseTrain_Marmoset(v(iChn,:),N,T,theseTrig,1,FS,abs0,~BREAK,carry);
                    end
                end
            end
            % Save a little data at the start of each loop
            nKeep = min(carry, size(v,2));
            hv = v(:,end-nKeep+1:end);
            % Don't save the last 'carry' samples, we'll deal with them next time
            if (BREAK)
                v = v(:,1:end-nKeep);
            end
            abs0 = abs0 + size(v,2);          % where the next array starts
            if Overall_time_to_analyse~=0
                BREAK=0;
            end
            fwrite(vdnFID,v./0.195,'int16');
            dispstat(sprintf('Progress %03.2f%%',100*((N)/nL)),'timestamp');
            N = N + 1;
        end
    end
    %% Clean up
    fclose(vFID);
    fclose(vdnFID);
end
end