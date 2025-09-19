function [Sp,Spktimes] = spikeextract(Mu,thresh,SAMPLING,artefact_low,artefact_high)
%  SPIKEEXTRACT 
%  
%   [Sp,Spktimes] = spikeextract(Mu,thresh)
%
%  Inputs:      MU  =   Array.  Multiunit activity
%               THRESH  =   Scalar.  Threshold for extraction.
%                               If two elements - extract positive and
%                               negative spikes
%               SAMPLING in Hz
%
%               artefact = Scalar. Ignores spikes below this - for artefact
%               removal
%
%  Outputs:     SP                =  Array.  Spike waveforms.
%                                               1ms upsampled 4x.
%                   SPKTIMES    =  Vector.  Spike times at 1ms.
%   

if(nargin < 3)
    SAMPLING = 20;
else
    SAMPLING=SAMPLING/1e3;
end
DT = SAMPLING*1.6;

flag = 0;
if thresh > 0
    flag = 1;
end
if length(thresh)==2
    thresh = sort(thresh,'ascend');
    flag = 2; 
end
if nargin >= 4
    flag = 3; % used to ignore artefact
end

% % extract waveforms            
tmu = zeros(1,size(Mu,2));
switch flag
    case 1
        %disp('Positive threshold')
        tmu(Mu > thresh)=1;   %  Positive threshold
    case 0
        %disp('Negative threshold')
        tmu(Mu < thresh)=1;   %  Negative threshold
    case 2
        %disp('Positive and negative threshold')
        tmu(find(Mu < thresh(1) | Mu > thresh(2)))=1;
    case 3
        tmu1=(find(Mu < artefact_low));
        tmu1=tmu1(diff(tmu1)~=1);
        if ~isempty(tmu1)
            if (tmu1(1)-50)<=0
                for i=1:length(tmu1)
                    if (tmu1(i)-50)<=0
                        Mu(1:tmu1(i)+20)=0;
                    else
                        Mu(tmu1(i)-50:tmu1(i)+20)=0;
                    end
                 end
            elseif (tmu1(end)+20)>length(Mu)
                Mu(tmu1(end)-50:end)=0;
                for i=1:length(tmu1)-1
                    Mu(tmu1(i)-50:tmu1(i)+20)=0;
                end
            else
                for i=1:length(tmu1)
                    Mu(tmu1(i)-50:tmu1(i)+20)=0;
                end
            end
        end
        tmu(find(Mu < thresh(1)))=1;
        flag=0;
end

Spktimes = find(diff(tmu)>0)';
Sp = [];

n = find(Spktimes > DT & Spktimes < size(Mu,2)-DT);
Spktimes = Spktimes(n);
nsp = length(n);
%disp([num2str(nsp) ' spikes']);
if nsp>0
    base = ones(nsp,1)*[-DT:DT];
    align = Spktimes(:,ones(1,size(base,2)));
    intspike = resample(reshape(Mu(base+align),size(align))',4,1);  % Upsampling
    switch flag
        case 1
        [a,ind] = max(intspike(3*DT:5*DT,:));    %  Positive threshold
        ind = ind + 3*DT-1;
        case 0
        [a,ind] = min(intspike(3*DT:5*DT,:));    %  Negative threshold
        ind = ind+3*DT-1; % correct the offset
        case 2
            neg = find(intspike(4*DT,:)<0);
            [a,negind] = min(intspike(3*DT:5*DT,neg));    %  Negative threshold
            pos = find(intspike(4*DT,:)>0);
            [a,posind] = max(intspike(3*DT:5*DT,pos));    %  Positive threshold
            [dum,interleaveind]=sort([neg pos],'ascend');
            negposind = [negind posind];
            negposind = negposind(interleaveind);
            ind = negposind + 3*DT-1;
    end
    ind = ind';
    % Sp = intspike;
    %  The following is to avoid indexing spikes at edge of spike extraction
    %  window
    n = find(ind > DT+1 & ind < size(intspike,1)-2*DT);  %  What are the right numbers here for 8x oversampling?
    %disp(['  ' num2str(length(n)) ' spikes']);
    if n > 0
        base = ones(length(n),1)*[-DT:3*DT];
        x = 1:length(n);
        xx = x(ones(1,size(base,2)),:)*size(intspike,1);
        xx = xx'-size(intspike,1);
        align = ind(n,ones(1,size(base,2)));
        index = base+align+xx;
        index = index';
        intspike = intspike(:,n);
        Sp = intspike(index(1:4:end,:));            %  Downsampling
        % Adjustment to how this works
        % We know the minimum of each wave - this should be centered
        % And spktime should reflect the minima of each wave
        % This will allow for removal of duplicate threshold crossings        
        % spktimes are taken from index DT*4+1 of intspike
        % ind has the intspike minimum
        Spktimes = Spktimes + (ind - (DT*4+1))./4; % splike time: compensate the upsampling offfset 
        Spktimes = Spktimes(n')./SAMPLING; % convert index to seconds
        if nargin ==5 % remove potential artifacts (if artifact threshold_high is provided)
            [r,c]=find(Sp>artefact_high);
            Sp(:,c)=[];
            Spktimes(c)=[];
        end
    else
        Sp = zeros(1,DT+1);
        Spktimes = zeros(1,1);
    end
else
    %Sp = zeros(1,DT+1);
    Sp = [];
    Spktimes = zeros(1,1);
end

Sp =Sp';
