% This function extracts relevant information from a .nev file
function NEVresponse(MODE)
if ~nargin
    MODE = 1;
end
% MODE controls for NEV or NS6 approach.
% Variables
digLinesStillBroken = 1;
nChn = 96;
dead_channels = 1001;
INDEX = 1:100; INDEX(100) = []; INDEX(91) = []; INDEX(10) = []; INDEX(1) = [];
XLABEL = 92:1:99;
YLABEL = 11:20:91;
% Initial Logic
[file, path, ~] = ...
    uigetfile('*.nev', 'Select a .nev data file', 'MultiSelect', 'off');
here = pwd;
cd(path);
% Load the data into an output structure
init_struct = openNEV([path '\' file],'nosave');
cd(here);
if (MODE == 1)    
    % Extract time_stamps
    t_stamps = init_struct.Data.SerialDigitalIO.TimeStamp;
    if (digLinesStillBroken)
        t_stamps = t_stamps(cast([1 diff(t_stamps) > 1000],'logical'));
    end
    trig = t_stamps;
    trig = trig(1:2:end);
    % Extract spiking
    sp = cell(nChn,1);
    for t = 1:length(init_struct.Data.Spikes.TimeStamp)
        ind = init_struct.Data.Spikes.Electrode(t);
        if (ind <= nChn)
            sp{ind} = [sp{ind}; init_struct.Data.Spikes.TimeStamp(t)];
        end
    end
elseif (MODE == 2)
    sp = loadSpikesBR();
    trig = loadTrigBR(0);
end
% Extract channel mapping
elec_id = zeros(nChn,1);
for i = 1:nChn
    fDig = str2double(init_struct.ElectrodesInfo(i).ElectrodeLabel(5));
    sDig = str2double(init_struct.ElectrodesInfo(i).ElectrodeLabel(6));
    if isnan(sDig)
        elec_id(i,:) = fDig;
    else
        elec_id(i,:) = str2double(init_struct.ElectrodesInfo(i).ElectrodeLabel([5,6]));
    end
end
% Run a quick response check
n_REP = [];
baselineBIN = [-575, -21];
stimulationBIN = [1, 41];
responseFig = figure('units','normalized','outerposition',[0 0 1 1]);
s_axes = tight_subplot(10,10,[.035 .025],[.075 .080],[.05 .01]);
delete(s_axes([1,10,91,100]));
%% Load stimulus variables
[~, path, ~] = ...
    uigetfile('*.mat', 'Select a .mat experimental data file', 'MultiSelect', 'off');
cd(path);
AMP = loadAMP; loadCHN;
TrialParams = loadTrialParams; TrialParams = cell2mat(TrialParams);
loadNREP;
for tChn = 1:nChn
    nTr = length(AMP);
    nT = zeros(n_REP,nTr*length(CHN));
    for i = 1:nTr*length(CHN)
        tmp = TrialParams(TrialParams(:,2) == i,1);
        nT(:,i) = tmp(1:n_REP);
    end
    preStim = zeros(n_REP,nTr*length(nChn)); postStim = zeros(n_REP,nTr*length(nChn));
    TICKLABELS = cell(1,length(AMP)*length(CHN));
    for iC = 1:length(CHN)
        for i = 1:length(AMP)
            TICKLABELS{i+(iC-1)*length(AMP)} = AMP(i);
        end
    end
    TICKLABELS{1} = 'N.S.';
    % Select the right channel
    spChn = elec_id(tChn);
    spt = sp{spChn}(:,1);
    for r = 1:n_REP %#ok<*BDSCI>
        for t = 1:nTr*length(CHN)
            try
                thisTrig = trig(nT(r,t));
                preStim(r,t) = sum(spt >= (thisTrig + baselineBIN(1)) & spt <= (thisTrig + baselineBIN(2)));
                postStim(r,t) = sum(spt >= (thisTrig + stimulationBIN(1)) & spt <= (thisTrig + stimulationBIN(2)));
            catch
                preStim(r,t) = NaN;
                postStim(r,t) = NaN;
            end
        end
    end
    preMean = nanmean(preStim,1) ./ (diff(baselineBIN)/1e3);
    postMean = nanmean(postStim,1) ./ (diff(stimulationBIN)/1e3);
    preSte = (nanstd(preStim,1) ./ (diff(baselineBIN)/1e3)) ./ sqrt(n_REP);
    postSte = (nanstd(postStim,1) ./ (diff(stimulationBIN)/1e3)) ./ sqrt(n_REP);
    figure(responseFig); hold on;
    ax = s_axes(INDEX(tChn));
    hold(ax,'on');
    dotSize = 40;
    xdata = zeros(1,length(AMP)*length(CHN));
    for iC = 1:length(CHN)
        for iA = 1:length(AMP)
            xdata(1,iA+(iC-1)*length(AMP)) = AMP(iA) + (iC-1)*(AMP(end)+6);
        end
    end
    if any(dead_channels == tChn)
        xlim(ax,[0 1]);
        xticks(ax,'');
        ylim(ax,[0 1]);
        yticks(ax,'');
        line(ax,[0 1],[0 1],'Color','k');
        line(ax,[0 1],[1 0],'Color','k');
        beautifyPlot(30,ax);
        title(ax,num2str(tChn));
        drawnow;
    else
        scatter(ax,xdata,preMean,dotSize,[0,0,0]);
        errorbar(ax,xdata,preMean,preSte,'Color',[0,0,0],'LineStyle','none','LineWidth',1,'HandleVisibility','off');
        scatter(ax,xdata,postMean,dotSize,[0,0,0],'filled','r');
        errorbar(ax,xdata,postMean,postSte,'Color',[0,0,0],'LineStyle','none','LineWidth',1,'HandleVisibility','off');
        if find(XLABEL == INDEX(tChn),1)
            xlabel(ax,'Current (uA)');
            xlim(ax,[-2 xdata(end)+1]);
            xticks(ax,xdata);
            xticklabels(ax,TICKLABELS);
        else
            xlim(ax,[-2 xdata(end)+1]);
            xticks(ax,xdata);
            xticklabels(ax,'');
        end
        if find(YLABEL == INDEX(tChn),1)
            ylabel(ax,'Response (Sp/s)');
        end
        if max(postMean) == 0
            ylim(ax,[min(postMean)-3*max(postSte) 20]);
            yticks(ax,0:10:20);
        else
            ylim(ax,[min(postMean)-3*max(postSte) max(postMean)+3*max(postSte)]);
            y = ceil(max(postMean)+3*max(postSte)) / 3;
            if (y < 100)
                if (y < 10)
                    y = round(y);
                else
                    y = round(y,-1);
                end
            else
                y = round(y,-2);
            end
            if y == 0
                y = 1;
            end
            yticks(ax,round(min(postMean)-3*max(postSte)):y:round(max(postMean)+3*max(postSte)));
        end
        %ylim(ax,[0 100]);
        %yticks(ax,0:25:100);
        title(ax,num2str(tChn));
        beautifyPlot(14,ax);
        drawnow;
    end
end
end