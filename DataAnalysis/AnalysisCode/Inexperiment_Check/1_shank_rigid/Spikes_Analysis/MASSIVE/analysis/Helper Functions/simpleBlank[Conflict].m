function data = simpleBlank(data,N,T,trig,mode)
TrialParams = loadTrialParams; TrialParams = cell2mat(TrialParams(:,2))';
nChn = size(data,1); FS = 30000;
if mode == 1 % for HPF
    BIN_b = -30; % Amount to blank in samples
    missed_trials = TrialParams(trig > ((N-1)*T*FS)-200 & trig < ((N-1)*T*FS)+200);
    missed_trig = trig(trig > ((N-1)*T*FS)-200 & trig < ((N-1)*T*FS)+200);
    trials = TrialParams(trig >= ((N-1)*T*FS)+200 & trig <= (N*T*FS)-200);
    trig = trig(trig >= ((N-1)*T*FS)+200 & trig <= (N*T*FS)-200);
    if ~isempty(missed_trig)
        trials = [missed_trials trials];
        trig = [missed_trig trig];        
    end
    if (N == 1)
        trig = trig - ((N-1)*T*FS);
    else
        trig = trig - ((N-1)*T*FS) + FS;
    end    
    for t = 1:length(trig)
        % Don't blank zero trials
        if trials(t) == 1
            continue;
        end
        thisTrig = trig(t);
        for c = 1:nChn
            BIN_a = 2*abs(BIN_b)+30;
            % Calculate the targeted voltage
            target_V = mean(data(c,thisTrig+210:thisTrig+270));
            target_D = diff([data(c,thisTrig+BIN_b), target_V]);
            blank_tolerance = 50;
            % Create temporary data structure
            tmp_data = data(c,thisTrig+(2*BIN_b):thisTrig+900);
            % Set up the while loop
            tmp_data(abs(BIN_b):BIN_a) = interpolate(tmp_data(abs(BIN_b):BIN_a),1);
            while abs(range(tmp_data(abs(BIN_b):240+2*(abs(BIN_b))))) > abs(target_D)+blank_tolerance
                BIN_a = BIN_a + 1;
                if tmp_data(BIN_a) >= tmp_data(abs(BIN_b))
                    BIN_a = find(tmp_data(BIN_a:end) <= target_V+blank_tolerance,1) + BIN_a;
                else
                    BIN_a = find(tmp_data(BIN_a:end) >= target_V-blank_tolerance,1) + BIN_a;
                end
                if BIN_a > (2*abs(BIN_b) + 150)
                    break;
                end
                tmp_data(abs(BIN_b):BIN_a) = interpolate(tmp_data(abs(BIN_b):BIN_a),1);
            end
            tmp_data(abs(BIN_b):BIN_a+20) = interpolate(tmp_data(abs(BIN_b):BIN_a+20),1);
            data(c,thisTrig+(2*BIN_b):thisTrig+900) = tmp_data;
%             % Check the size of the deviation
%             chk1 = 1000;
%             while (abs(chk1) > blankMax) % Magic number - maximum acceptable voltage deviation
%                 BIN_a = BIN_a + 1;
%                 %data(c,trig(t) + BIN_b:trig(t) + BIN_a) = interpolate(data(c,trig(t) + BIN_b:trig(t) + BIN_a),1);
%                 chk2 = zeros(1,BIN_a+30);
%                 for cc = BIN_a:BIN_a+30
%                     chk2(cc-BIN_a+1) = diff([data(c,thisTrig + BIN_a),data(c,thisTrig + cc)]);
%                 end
%                 chk1 = max(abs(chk2));
%                 if BIN_a >= 150
%                     break;
%                 end
%                 chk3 = chk1;
%                 while (abs(chk3) > 100)
%                     BIN_a = BIN_a + 1;
%                     chk3 = diff([data(c,thisTrig + BIN_b),data(c,thisTrig + BIN_a)]);
%                     if BIN_a >= 150
%                         break;
%                     end
%                 end
%             end
%             % Interpolate the data once.
%             data(c,thisTrig + BIN_b:thisTrig + BIN_a) = interpolate(data(c,thisTrig + BIN_b:thisTrig + BIN_a),1);
        end
    end
elseif mode == 2 % for LPF
    bBIN = [5, 250*30]; % Amount to blank in samples
    % Check for missed trigger lines
    missed_trig = trig(trig > ((N-1)*T*1000)-550 & trig < ((N-1)*T*1000)+550);    
    trig = trig(trig >= ((N-1)*T*1000)+550 & trig <= (N*T*1000)-550);    
    if ~isempty(missed_trig)
        trig = [missed_trig trig];
    end
    if (N == 1)
        trig = trig - ((N-1)*T*1000);
    else
        trig = trig - ((N-1)*T*1000) + 1000;
    end
    trig = cast(trig,'int64');
    for t = 1:length(trig)
        for c = 1:size(data,1)
            %distant_channel = 13 + c; % Get as far away as possible
            %if distant_channel > 26
            %    distant_channel = distant_channel - 26;
            %end
            %if distant_channel == 2 || distant_channel == 9
            %    distant_channel = distant_channel + 1; % Deal with dead channels
            %end            
            trend1 = interpolate(data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)),1);
            wave = (data(c,trig(t)*30 + bBIN(2):trig(t)*30 + bBIN(2) + diff(bBIN)));
            trend2 = interpolate(wave,1);
            wave = wave - trend2 + trend1;
            data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)) = wave;
            %data(c,trig(t)*30 - 5:trig(t)*30 + 5) = interpolate(data(c,trig(t)*30 - 5:trig(t)*30 + 5),1);
        end
    end
end
end