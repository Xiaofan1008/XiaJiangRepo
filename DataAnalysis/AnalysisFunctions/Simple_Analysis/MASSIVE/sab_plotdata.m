close all
clc

tic
nChn=32;
fs=30000;
nosamp=30000*100;
skipdigi=1;
if skipdigi==1
    fileID=fopen('digitalin.dat','r');
    data = fread(fileID,[fs*1000,1], 'uint16');
    % C = textscan(fileID,'%s %s %f32 %d8 %u %f %f %s %f');
    fclose(fileID);
    whos C
    [r,c]=find(data>0.5);
    changestate=diff(data);
    [r1,c1]=find(changestate~=0);
    rt=r1/fs;
    nospikes=length(r1)/2;
    time=length(data)/(fs);
    ratio=nospikes/time;
    samptriallength=diff(r1);
    timetriallength=samptriallength/fs;
    timespike=timetriallength(2:2:length(timetriallength));
    timebetweenspike=timetriallength(1:2:length(timetriallength)-1);
    avspiketime=sum(timespike)/length(timespike);
    avsbetspiketime=sum(timebetweenspike)/length(timebetweenspike);
    figure (2)
    hold on
    subplot(2,1,1)
    plot(timespike)
    title('Time high')
    ylabel('Time(s)')
    xlim([0 length(timespike)])
    hold off
    hold on
    subplot(2,1,2)
    plot(timebetweenspike)
    xlim([0 length(timebetweenspike)])
    title('Time low')
    xlabel('Spike number')
    ylabel('Time(s)')
    hold off
    figure
    px=1486371;
    py=1486450;
    figure (1)
    plot(px:py,data(px:py))
    toc
end

T=length(data)/fs;
tic
fileID=fopen('estim_pen2_002.lfp.dat','r');
ampdata = fread(fileID,[32, T*fs], 'int16');
fclose(fileID);
toc
%ampdatasep=reshape(ampdata,32,[]);

px=722125;
py=722200;
x=px/(fs):1/(fs):py/(fs);
figure (4)
hold on
xlabel('Time(s)')
% X=[rt(88) rt(88) rt(89) rt(89)];
% Y=[0.5*10^4 -0.5*10^4 -0.5*10^4 0.5*10^4];
% basevalue = -0.5;
% area(X,Y)
plot(x,ampdata(1, px:py),'b')
hold off



px=722125;
py=722200;
x=px/(fs):1/(fs):py/(fs);
figure 
hold on

for i=1:3:30
    plot(x,ampdata(i, px:py))
end
%legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32')
legend('1','5','9','13','17','21','25','29','30')
 
plot(x,ampdata(8, px:py),'k')
% plot(x,ampdata(3, px:py),'g')
 plot(x,ampdata(18, px:py),'b')
xlabel('Time(s)')
ylabel('Amplitude')


%6 to 9mins is okay data

%%
% figure (3)
% for i=1:32
% subplot(4,8,i)
% plot(x,ampdata(i,px:py))
% xlabel('Time(s)')
% ylim([-1*10^4 1*10^4])
% end
% 
% 
% px=401.6*fs;
% py=402.6*fs;
% x=px/(fs):1/(fs):py/(fs);
% figure (5)
% hold on
% xlabel('Time(s)')
% plot(x,ampdata(24, px:py),'b')
% ylim([-0.4*10^4 0.2*10^4])
% hold off



%%
[MUAfilt, LFPfilt] = generate_Filters;
MuNf = length(MUAfilt);
LfpNf = length(LFPfilt);
Ndata = size(ampdata,2);
MuOut = cell(1,32);
for iChn = 1:32
    flip_data = fliplr(ampdata(iChn,:));
    tmp = conv(flip_data,MUAfilt);
    mu = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
    MuOut{iChn} = mu;
end
n=3;%filterorder are 2n so n=3 6th order filter
lowcutfreqlfp=60/(fs/2);%0.5hz
highcutfreqlfp=500/(fs/2);%500hz

lowcutfreqspike=100/(fs/2);
highcutfreqspike=1000/(fs/2);


[bspike,aspike] = butter(n,[lowcutfreqspike highcutfreqspike],'bandpass');

filtereddataspike=filter(bspike,aspike,ampdata(22, px:py));

[blfp,alfp] = butter(n,[lowcutfreqlfp highcutfreqlfp],'bandpass');

filtereddatalfp=filter(blfp,alfp,ampdata(22, px:py));

figure
hold on
plot(ampdata(32,px:py),'b')
plot(mu(px:py),'r')
%ylim([-3000 1500])
hold off

figure
hold on
plot(filtereddataspike,'b')
plot(ampdata(22, px:py),'r')
%ylim([-3000 1500])
hold off

