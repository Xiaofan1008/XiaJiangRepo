function stimConcatenate(tID,Chosen_trig,chan,type,Binstart,Binend,varargin)
dbstop if error
%% This function concatenates twenty-one milliseconds of high-pass filtered data from all stimulation events
if ~nargin
    file = pwd;
    chan = [1 2 3 4];
    tID = [1];
    type='MU';
    Binstart=-250;
    Binend=600;
end
if isempty(varargin)
    file = pwd;
end
Blackrock = 0;
%% Variables
if (Blackrock)
    nChn = 96;
else
    filepath = pwd;
    fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
    fileinfo = dir([filepath,filesep, 'info.rhs']);
    if (datetime(fileinfo.date) < fourShank_cutoff)
        nChn=32;
        E_Mapnumber=0;
    else
        E_Mapnumber=loadMapNum;
        if E_Mapnumber>0
            nChn=64;
        else
            nChn=32;
        end
    end
    E_MAP = Depth(E_Mapnumber);
end
FS = 30000;
BIN = [Binstart Binend]; % ms
SPACING = 200;

%% Loading in data
%SP = dir('*.sp.mat');
if (Blackrock)
    d = 1:96; %#ok<*UNRCH>
else
    d = Depth(E_Mapnumber);
end
stimChn = 1; loadStimChn; %thresh = load(SP.name,'thresh'); thresh = thresh.thresh; 
loadNREP; % number of repeats
figure; hold on;
for tt = tID    
    if (Blackrock)
        trig = loadTrigBR(0);
        mDIR = dir([file filesep '*muBR.dat']);
    else
        trig = loadTrig(0);
%         trig(550/2:588/2)=0;
%         trig(7244/2:7246/2)=0;
%         b = [0]; 
%         k = 550/2; %row position, can be 0,1,2 or 3 in this case
%         trig = [trig(1,1:k) b trig(1,k+1:end)];
%         k = 7244/2; %row position, can be 0,1,2 or 3 in this case
%         trig = [trig(1,1:k) b trig(1,k+1:end)];
        if strcmp(type,'MU')
            mDIR = dir([file filesep '*mu_sab.dat']);
        elseif strcmp(type,'MUA')
            mDIR = dir('MUAdata.dat');
            SPACING = 50;
        else
            mDIR = dir('amplifier.dat');
            SPACING = 500;
        end
    end
    mNAME = mDIR.name;
    mFID = fopen([file filesep mNAME],'r');
    TrialParams = loadTrialParams;
    if ~isempty(TrialParams)
        TrialParams = find(cell2mat(TrialParams(:,2)) == tt);
    else
        TrialParams = 1:length(trig);
    end
    TrialParams = TrialParams(2:2:n_REP*2);
    trig = trig(TrialParams./2);
    %output = trig_helper([6 -251],7);
    %ind = output(output(:,2)==tt,:);
    %ind = ind(ind(:,5)==1,1);
    %trig = trig(ind);
    nTrig = length(trig);    
    setupGraphsBasic;
    %% Initialise data matrices
    X = BIN(1):1/30:BIN(2) - 1/30;
    MEAN = zeros(1,diff(BIN)*30);
    for t=Chosen_trig(1):Chosen_trig(end)
        OFFSET = cast(nChn*2*(trig(t)+(BIN(1)*FS/1e3)),'int64');
        fseek(mFID,OFFSET,'bof'); % used to point to specific point in the data indicated by offset from begginning of file
        v = fread(mFID,[nChn, (FS/1e3)*diff(BIN)],'short') ./ 10;
        a = 0;
        for c = chan
            a = a + 1;
            if c == stimChn
                if (length(Chosen_trig)==1)&&(length(chan)==1)
                    plot(X,v(d(c),:) + (a-1)*SPACING,'Color','r');
                    ylabel('uV')
                else 
                    plot(X,v(d(c),:) + (a-1)*SPACING + (t-1)*SPACING,'Color','r');
                end
            else
                if (length(Chosen_trig)==1)&&(length(chan)==1)
                   plot(X,v(d(c),:) + (a-1)*SPACING,'Color','k')
                    ylabel('uV')
                else
                    plot(X,v(d(c),:) + (a-1)*SPACING + (t-1)*SPACING,'Color','k');
                end
                %line([X(1) X(end)],[thresh{d(c)} + (a-1)*SPACING + (t-1)*SPACING thresh{d(c)} + (a-1)*SPACING + (t-1)*SPACING],'Color','r');
            end
            %text(0,(a-1)*SPACING+150,num2str(c));
            MEAN = MEAN + v(d(c),:);
        end
    end
    MEAN = MEAN ./ nTrig;
    if length(chan) == 1 && length(tID) == 1
        %plot(X,MEAN + (a-1)*SPACING + (tt-1)*SPACING,'Color','r','LineWidth',2);
    end
    disp(['A total of ' num2str(t) ' trials were included, of ' num2str(nTrig)]);
end
%line([3.2 3.2],[-1000 1000],'Color','r');
%% Try something new
% figure; hold on;
% nT = t;
% for t = 1:1
%     OFFSET = cast(nChn*2*(FS/1e3)*(trig(t)+BIN(1)),'int64');
%     fseek(mFID,OFFSET,'bof');
%     v = fread(mFID,[nChn, (FS/1e3)*diff(BIN)],'short') ./ 10;
%     if any(range(v,2) > ARTMAX)
%         continue;
%     end
%     for c = chan
%         if (t == nTrig)
%             text(-5,5,num2str(c),'FontSize',12,'Color','r');
%         end        
%         plot(X,v(d(c),:)-MEAN,'Color','k');   
%     end
% end
%% Close gracefully
fclose(mFID);
end