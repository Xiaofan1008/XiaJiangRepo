function stimConcatenate(file)
dbstop if error
%% This function concatenates twenty-one milliseconds of high-pass filtered data from all stimulation events
if ~nargin
    file = pwd;
end
Blackrock = 0;
%% Variables
if (Blackrock)
    nChn = 96;
else
    nChn = 32;
end
FS = 30000;
BIN = [-250 300]; % ms
chan = [3];
SPACING = 200;
tID = [5];
%% Loading in data
%SP = dir('*.sp.mat');
if (Blackrock)
    d = 1:96; %#ok<*UNRCH>
else
    d = Depth;
end
stimChn = 1; loadStimChn; %thresh = load(SP.name,'thresh'); thresh = thresh.thresh; 
loadNREP;
figure; hold on;
for tt = tID    
    if (Blackrock)
        trig = loadTrigBR(0);
        mDIR = dir([file '\*muBR.dat']);
    else
        trig = loadTrig(0);
        mDIR = dir([file '\*mu.dat']);
    end
    mNAME = mDIR.name;
    mFID = fopen([file '\' mNAME],'r');
    TrialParams = loadTrialParams;
    if ~isempty(TrialParams)
        TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == tt));
    else
        TrialParams = 1:length(trig);
    end
    TrialParams = TrialParams(1:n_REP);
    trig = trig(TrialParams);
    %output = trig_helper([6 -251],7);
    %ind = output(output(:,2)==tt,:);
    %ind = ind(ind(:,5)==1,1);
    %trig = trig(ind);
    nTrig = length(trig);    
    setupGraphsBasic;
    %% Initialise data matrices
    X = BIN(1):1/30:BIN(2) - 1/30;
    MEAN = zeros(1,diff(BIN)*30);
    for t = 1:21
        OFFSET = cast(nChn*2*(trig(t)+(BIN(1)*FS/1e3)),'int64');
        fseek(mFID,OFFSET,'bof');
        v = fread(mFID,[nChn, (FS/1e3)*diff(BIN)],'short') ./ 10;
        a = 0;
        for c = chan
            a = a + 1;
            if c == stimChn
                plot(X,v(d(c),:) + (a-1)*SPACING + (t-1)*SPACING,'Color','k');
            else
                plot(X,v(d(c),:) + (a-1)*SPACING + (t-1)*SPACING,'Color','k');
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