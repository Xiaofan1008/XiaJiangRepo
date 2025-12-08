function flag=checkTrials(trig,timeafter,timebefore,nChn,FS,chsp, spikesafter,varargin)
flag=0;
if nargin==1
    baselinspikes=varargin{1};
else
    baselinspikes=0;
end
if isempty(spikesafter)
    spikesafter=0; 
end
if any(abs(diff(spikesafter(:,1),5))<4)|| any(abs(diff(baselinspikes(:,1),5))<4) % way of comparing num of spks per second
    fileID=fopen('amplifier.dat','r');
    shortbytes=2;
    offset=trig*nChn*shortbytes-timebefore*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
    fseek(fileID,offset,'bof');
    vblankmu = fread(fileID,[nChn, (timeafter*FS+timebefore*FS)],'short')./10; %plots one second from trigger and 250ms brefore
    if any(vblankmu(chsp,1:round(timebefore-5/1000)*FS)>1000) || any(vblankmu(chsp,1:round(timebefore-5/1000)*FS)<-1000) || any(vblankmu(chsp,round((timebefore+5/1000)*FS):end)>1000) || any(vblankmu(chsp,round((timebefore+5/1000)*FS):end)<-1000)
        flag=1;
    end
    fclose(fileID);
end
end