ID = [6,10];
X = 0:5:20;
for i = ID(1):ID(2)
    selectedTrials(1,:) = find(allTrials == i);
    trials(i,:) = time_stamps(1,selectedTrials(1,:));
end
for i = ID(1):ID(2)
    for r = 1:n_REP
        OFFSET = cast(nChn*2*(FS/1e3)*(trials(i,r)+BIN(1)),'int64');
        fseek(v_fid,OFFSET,'bof'); % Find the relevant timestamp
        v = fread(v_fid,[nChn,(FS/1e3)*6000],'int16') * 0.195;
        fseek(v_fid,0,'bof'); % Return to the start of the file
        %% Extract spikes from v. OUTPUT will contain n_REP cell arrays
        % Blank. Note that BlankArte always takes a window of [-500, 500].
        blank_v = v;
        blank_v(:,(501*30):(1500*30)) = BlankArte(v(:,(501*30):(1500*30)),nChn);
        % Downsample blank-v to a smaller window
        blank_v(:,999*30:1004*30) = 0;
        blank_v = blank_v(:,(901*30):(1200*30));
        % Filter
        tmp = conv(blank_v(16,:),Mufilt);
        filt_v = tmp(1,MuNf/2:nData+MuNf/2-1);
        % Generate a threshold
        sd = median(abs(filt_v(200*30:299*30)))./0.6745;
        thresh = threshfac*sd;
        %thresh = 35;
        % Plot for debugging
        %X = 1/30:1/30:length(filt_v)/30;
        %figure;
        %plot(X,filt_v);
        %line([0 X(end)],[thresh thresh],'Color','red');
        % Extract spikes
        [~,Sp_tmp] = spikeextract(filt_v,thresh,FS);
        Sp_tmp(Sp_tmp < 100) = [];
        Sp_tmp(Sp_tmp > 120) = [];
        MEAN(i,r) = length(Sp_tmp);
    end
end
for i = ID(1):ID(2)
    MEAN(i,1) = mean(MEAN(i,:));
    MEAN(i,1) = MEAN(i,1) * 50;
end
figure
hold on
scatter(X,MEAN(ID(1):ID(2),1),50,[0,0,0],'filled');
xlabel('Current Level (uA)');
ylabel('Response (Spikes/sec)');
title('Stimulation/Response - Channel 16 | Recording 1');
