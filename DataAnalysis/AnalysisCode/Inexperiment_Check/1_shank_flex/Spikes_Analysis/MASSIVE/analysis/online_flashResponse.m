function online_flashResponse(TYPE)
dbstop if error
if nargin<1
    TYPE = 1;
end
if ~(any(TYPE == [1,2]))
    TYPE = 1;
end
%% Basic Variables
BIN = [-100 2100];
DEPTH = Depth;
%% Setup Logic
% Navigate to the right place
exp = dir('amplifier.dat');
while isempty(exp)
    % No amplifier data here
    [~, path, ~] = ...
        uigetfile('*.dat','Select a .dat amplifier data file','MultiSelect','off');
    cd(path);
    exp = dir('amplifier.dat');
end
% Check whether the datafile has been processed
lfp = dir('*.lfp.dat');
if isempty(lfp)
    fprintf('Has this directory been processed?\n');
    fprintf('We will run processIntan(1,1) in this directory.\n');
    processIntan(1,1);
end
if (TYPE == 1)
    lfp = loadLFP;
    nChn = size(lfp,1);
elseif (TYPE == 2)
    sp = loadSpikes;
    nChn = size(sp,2);
    SMOOTHING = 3;
end
% Check trigger lines
trig = dir('*.trig.dat');
if isempty(trig)
    fprintf('No trigger lines detected.\n');
    fprintf('Please enter trigger values now.\n');
    WF = input('How many white frames per flash?    ');
    BF = input('How many black frames per flash?    ');
    nF = input('How many flashes in total?    ');
    FR = input('What is the screen framerate?    ');
    % Construct the trigger lines
    msPerFlash = 1000*(WF+BF)/FR;
    trig = abs(BIN(1))+1:msPerFlash:(abs(BIN(1)))+nF*msPerFlash;
else
    trig = loadTrig(1);
end
nTrig = length(trig);
%% Construction Logic
if (TYPE == 1)
    trig = cast(trig,'int64');
    LFP = zeros(nChn,diff(BIN)+1);
    for t = 1:nTrig
        for c = 1:nChn
            LFP(c,:) = LFP(c,:) + lfp(DEPTH(c),trig(t)+BIN(1):trig(t)+BIN(2));
        end
    end
    LFP = LFP ./ nTrig;
    %% Plotting Logic
    RANGE = zeros(1,nChn);
    for c = 1:32
        RANGE(c) = range(LFP(c,:));
    end
    SPACING = min(RANGE);
    figure; hold on;
    for c = 1:nChn
        plot(LFP(c,:)+(c-1)*SPACING);
    end
elseif (TYPE == 2)
    xdata = cell(nChn,1);   
    for t = 1:nTrig
        for c = 1:nChn
            spt = sp{c}(:,1);
            theseSp = (spt(spt > trig(t)+BIN(1) & spt < trig(t)+BIN(2)) - trig(t));
            for i = 1:length(theseSp)
                xdata{c} = [xdata{c}, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
            end
        end
    end
    %% Plotting Logic
    figure;
    s_axes = tight_subplot(6,6,[.025 .025],[.075 .080],[.05 .01]);
    delete(s_axes([1,6,31,36]));
    for c = 1:nChn
        A = zeros(6,6);
        for i = 1:32
            A(i+3) = i;
        end
        for i = 1:28
            A(i+2) = i;
        end
        A(1,6) = 0;
        for i = 1:4
            A(i+1) = i;
        end
        A(6,1) = 0;
        INDEX = find(A' == c);
        ax = s_axes(INDEX); %#ok<FNDSB>
        hold(ax,'on');
        Z = hist(xdata{c},0:diff(BIN)); %#ok<HIST>
        window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
        rate = (1000/nTrig)*conv(Z,window);
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING) ./ nTrig;
        plot(ax,rate+(c-1)*40,'k');
        xlim(ax,[abs(BIN(1))+200 abs(BIN(2))-200]);
        ylim([-10 40]);
    end
end