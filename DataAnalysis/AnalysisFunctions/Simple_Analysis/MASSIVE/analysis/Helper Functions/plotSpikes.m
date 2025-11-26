H = figure; hold on
for i = 1:n_REP
    plot(SPIKES{i}');
end
AVERAGE = zeros(1,49);
COUNT = 0;
nSpikes = 0;
COUNT_all = 0;
nSpikes_all = 0;
for i = 1:n_REP
    if ~isempty(mean(SPIKES{i},1))
        if ~isnan(mean(SPIKES{i},1))
            AVERAGE = AVERAGE + mean(SPIKES{i},1);
            COUNT = COUNT + 1;
            nSpikes = nSpikes + size(SPIKES{i},1);
        end
    end
    if ~isempty(mean(OUTPUT_all{i},1))
        if ~isnan(mean(OUTPUT_all{i},1))
            COUNT_all = COUNT_all + 1;
            nSpikes_all = nSpikes_all + size(OUTPUT_all{i},1);
        end
    end
end
text([30 30], [-200 -200],['n = ' num2str(nSpikes)]);
AVERAGE = AVERAGE ./ COUNT;
plot(AVERAGE','Color','k','LineWidth',4);
ylim([-250 150]);
xlabel('Time (spikeextract samples)');
ylabel('Voltage (uV)');
title([num2str(n_REP) ' Trials | Stimulus ' num2str(thisAmp) ' uA | ' num2str(thisDur) ' us']);
disp(['There are: ' num2str(nSpikes) ' spikes in these trials.']);
disp(['The spiking rate is: ' num2str(nSpikes/(32*10/1000)) ' Spikes/sec.']);
disp(['There are: ' num2str(nSpikes_all) ' spikes in this dataset.']);
disp(['The overall spiking rate is: ' num2str(nSpikes_all/(32*(diff(BIN))/1000)) ' Spikes/sec.']);
savepath = [filepath '\' num2str(CHANNEL) '\' num2str(STIMCHN) '\'];
if ~exist(savepath,'dir')
    mkdir(savepath);
end
savepath = [savepath  num2str(thisAmp) '-' num2str(thisDur) '.png'];
saveas(H,savepath);
pause(1.000);
close(H);