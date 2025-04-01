function spikes = plotRawdata(chsp,tID,DataType,plotspktime)
%DataType= DT dn mu or amplifier
%plotspktime= 1 to plot spikes
%chsp is unsorted channel
%tID= trial ID
%spike is the spike time
nChn=64;
FS=30000;
timebefore=0.09;%in seconds 0.09
timeafter=0.09;%in seconds0.09
trig=loadTrig;
TrialParams=loadTrialParams; % 1: Trial number 2: Trial ID 3: Channel
numelect=find(diff(cell2mat(TrialParams(:,1)))~=0,1,'first'); % find the first non-zero value of the adjacent elements, to determine how many rows corresponsding to one trial. 
% make each trial match one single row 
TrialParamstID = find(cell2mat(TrialParams(1:numelect:end,2)) == tID); %identifies trial row matching trial ID
savespktimes=[];
trigtID = trig(TrialParamstID);
trigtID(trigtID==-500)=[];
spikes = cell(size(length(trigtID),2));

filepath = pwd;
[filepathm,name,ext] = fileparts(filepath);
name = name(1:end-14); % delect last 14 bits (date of data)
if plotspktime==1
    sp=load([name '.sp.mat'],'sp');
    sp=sp.sp{chsp};
    thresh=load([name '.sp.mat'],'thresh');
    thresh=thresh.thresh{chsp};
end
if strcmp(DataType,'amplifier') % Raw voltage recording
    fileID=fopen('amplifier.dat','r');
elseif strcmp(DataType,'dn')
    fileID=fopen('amplifier_dn_sab.dat','r'); % Downsampled data
elseif strcmp(DataType,'mu')
    fileID=fopen([name '.mu_sab.dat'],'r'); % Multi-Unit Activity data
elseif strcmp(DataType,'DT')
    fileID=fopen([name '_DT.mu.dat'],'r'); % Sorted spikes (Single-unit activity)
end

try
    ftell(fileID);
catch
    return
end

shortbytes=2;
for indT=1:length(trigtID)
    offset=trigtID(indT)*nChn*shortbytes-timebefore*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
    fseek(fileID,offset,'bof');
    %fseek(fileID2,offset,'bof');
    ftell(fileID);
    if strcmp(DataType,'amplifier') || strcmp(DataType,'dn')
        vblankmu = fread(fileID,[nChn, (timebefore*FS+timeafter*FS)],'int16') .* 0.195;
    else 
        vblankmu = fread(fileID,[nChn, (timebefore*FS+timeafter*FS)],'short')./10; %plots one second from trigger and 250ms brefore
    end
    % *0.195 or /10 used to convert ADC value to the real voltage value
    % based on the hardware amplifier gain
    

    % if any(vblankmu(chsp,:)<-50)
    %         figure
    %         plot (-0.02*1000:1000/FS:0.02*1000-1000/FS,vbladnkmu2(chsp,:))
    %         title(['Channel ' num2str(chsp)])
    %         xlabel('Time (ms)')
    %         ylabel('Voltage (uV)')

    % ------- plot the MUA activity (filtered raw voltage) ------% 
    % figure(1)
    % hold on
    % plot (-timebefore*1000:1000/FS:timeafter*1000-1000/FS,vblankmu(chsp,:)) % *1000 convert second to millisecond 
    % title(['Channel ' num2str(chsp) ', Trig ' num2str(indT)])
    % xlabel('Time (ms)')
    % ylabel('Voltage (uV)')


    if plotspktime==1
        yline(thresh);
        offsetspk=trigtID(indT)*1000./FS; % spike starting time in ms
        spktimes=sp(((sp(:,1)>(offsetspk-timebefore*1000))&(sp(:,1)<(offsetspk+timeafter*1000))),:);
        store = [];
        % mark individual spike occurance
        for numspk=1:size(spktimes,1)
            scatter(spktimes(numspk,1)-offsetspk,min(spktimes(numspk,2:end)),'r');
            store(end+1) = spktimes(numspk,1)-offsetspk;
        end
        spikes{indT} = store;
        % plot spike waveform
        if ~isempty(spktimes)
            savespktimes=[savespktimes; (spktimes(:,1)-offsetspk)];
            % figure(100)
            % hold on
            % plot(1*1000/FS:1000/FS:49*1000/FS,spktimes(:,2:end)')
            % title(['Channel ' num2str(chsp) ', Trig ' num2str(indT)])
            % xlabel('Time (ms)')
            % ylabel('Spike amplitude (uV)')
        else
            close
        end
        
    end

    %end
end
fclose(fileID);
end