function generate_Raster(chn,ID,nT,PHASE)
dbstop if error
%% Load in data from current directory
d = Depth; chn = d(chn);
sp = loadSpikes(chn);
sp = denoiseSpikes(sp,chn);
trig = loadTrig(0);
TP = loadTrialParams;
tID = cell2mat(TP(cell2mat(TP(:,2)) == ID,1));
theseTrig = trig(tID)./30;
trig = trig(tID);
%lfp = loadLFP(chn);
%% Set up the raster data structure
BIN = [-200 200]; offset = [-251,-200,-251,-200,-100,-25];
SPACING = 200; SMOOTHING = 2; MAX = 400;
xdata = [];
ydata = [];
% Grab Data for the Inset Figure
FS = 30000; nChn = 32; file = pwd;
mDIR = dir([file '\*mu_sab.dat']);
mNAME = mDIR.name;
mFID = fopen([file '\' mNAME],'r');
X = BIN(1):1/30:BIN(2) - 1/30;
mu = zeros(nT,(FS/1e3)*diff(BIN));
for tr = 1:nT
    theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
    for i = 1:length(theseSp)
        xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
        ydata = [ydata, tr*(MAX/nT)];
    end
    % Also, extract MUA waveforms
    OFFSET = cast(nChn*2*(trig(tr)+(BIN(1)*FS/1e3)),'int64');
    fseek(mFID,OFFSET,'bof');
    v = fread(mFID,[nChn, (FS/1e3)*diff(BIN)],'short') ./ 10;
    mu(tr,:) = v(chn,:);
    % Also, calculate phase for each trial
%     phase = generatePhaseVector(lfp,[4,5,6,10,20,80]);
%     for b = 1:6
%         thisPhase(tr,b) = phase{b}(cast(theseTrig(tr)/30 + offset(b),'int64'));
%     end
end
%% Plot the raster
if PHASE
    C = define_Colormap;
    % Need to sort the raster by phase
    [~,index] = sort(thisPhase(:,2));
    sYdata = unique(ydata);
    sYdata = sYdata(index');
    newYdata = [];
    for n = 1:nT
        count = size(ydata(ydata == sYdata(n)),2);
        newYdata = [newYdata, repmat(sYdata(n),[1,count])];
    end
    ydata = newYdata;
    figure; hold on; axM = gca;
    yScale = MAX./nT;
    for i = 1:size(xdata,2)
        line([xdata(i) xdata(i)],[ydata(i)-yScale/2.3 ydata(i)+yScale/2.3],'Color','k','LineWidth',2);
    end
    line([200 200],[0 MAX],'Color','b');
    % Convert coordinates for phase plot
    radius = yScale/3; b = 2; 
    len = 0;
    xScale = MAX/150;
    for p = 1:6
        index = (sort(thisPhase(:,2)) >= -180 + 60*(p-1) & sort(thisPhase(:,2)) <= -180 + 60*(p));
        index(index == 0) = [];
        newlen = length(index);
        line(axM,[abs(BIN(1))-5 abs(BIN(1))-5+(radius*xScale*cosd(-180 + (p-1)*60))],[(len+newlen/2)*yScale + radius (len+newlen/2)*yScale + radius + (radius*yScale*sind(-180 + (p-1)*60))],'Color','r','LineWidth',2);
        line(axM,[abs(BIN(1))-5 abs(BIN(1))-5+(radius*cosd(-180 + (p)*60))],[(len+newlen/2)*yScale + radius (len+newlen/2)*yScale + radius + (radius*yScale*sind(-180 + (p)*60))],'Color','r','LineWidth',2);
        len = len + length(index);
        line([150 300],[len*(MAX/nT)-yScale/2.3 len*(MAX/nT)-yScale/2.3],'Color',C{p});
    end
    % Add the convolved spikerate
    Z = hist(xdata,0:400); %#ok<HIST>
    window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
    rate = (1000/nT)*conv(Z,window);
    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
    plot(axM,rate,'b','LineWidth',2);
    ylim(axM,[0 MAX]); xlim(axM,[150 300]);
    % Tidy Up
    xlabel(axM,'Time (ms)');
    ylabel(axM,'Firing Rate (Sp/s)');
    xticks(axM,150:50:300);
    xticklabels(axM,-50:50:100);
    beautifyPlot(30,axM);
else
    figure; hold on; axM = gca;
    yScale = MAX./nT;
    for i = 1:size(xdata,2)
        line([xdata(i) xdata(i)],[ydata(i)-yScale/2.3 ydata(i)+yScale/2.3],'Color','k','LineWidth',2);
    end
    line([200 200],[0 MAX],'Color','b');
    % Add an inset to show it's not artefact
    axI = axes('Position',[.175 .7 .2 .2]); hold on;
    xlim(axI,[-10 20]); %xlabel('Time (ms)');    
    for tr = 1:nT
        plot(axI,X,mu(tr,:)+(tr-1)*SPACING,'Color','k');
        %line(axI,X,thresh*ones(1,length(X))+(tr-1)*SPACING,'Color','r');
    end
    ylim(axI,[3500 5000]); yticks(''); xticklabels('');
    line(axI,[0 0],[0 nT*SPACING],'Color','b');
    % Add the convolved spikerate
    Z = hist(xdata,0:400); %#ok<HIST>
    window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
    rate = (1000/nT)*conv(Z,window);
    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
    plot(axM,rate,'b','LineWidth',2);
    ylim(axM,[0 MAX]); xlim(axM,[150 300]);
    % Tidy Up
    xlabel(axM,'Time (ms)');
    ylabel(axM,'Firing Rate (Sp/s)');
    xticks(axM,150:50:300);
    xticklabels(axM,-50:50:100);
    beautifyPlot(30,axI);
    beautifyPlot(30,axM);
end
end