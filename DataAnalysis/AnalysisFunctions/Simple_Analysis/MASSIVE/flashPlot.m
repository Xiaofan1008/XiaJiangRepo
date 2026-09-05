chk = 1; N = 0; Ninput = 1e6; time = 0; SCALEFACTOR = 10; T=256;
nChn=64;
E_Mapnumber=1;
par=0;
E_MAP = Depth(E_Mapnumber);
FS=30000;
filepath = pwd;
[~,name,~] = fileparts(filepath);
name = name(1:end-14);
[~,Lfpfilt] = generate_Filters;
LfpNf = length(Lfpfilt);
dName='amplifier';
fileinfo = dir([filepath filesep dName '.dat']);
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);


%if isempty(dir('meanlfp.mat'))
fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);
stimDig = flip(find(digital_in == 0)); % Fix for finding the falling line instead of the rising line

visDig = flip(find(digital_in == 2)); % Fix for finding the falling line instead of the rising line
dt = diff(stimDig);
kill = dt == -1;
stimDig(kill) = [];
dt = diff(visDig);
kill = dt == -1;
visDig(kill) = [];
nStamps = max([length(stimDig); length(visDig)]);
time_stamps = nan(1,nStamps);
time_stamps(1,1:length(stimDig)) = flip(stimDig);
time_stamps(1)=[];
trig=time_stamps;


%load([name '.lfp_sab.mat'])
LFPstruct=[];
n_REP=nStamps-1;

        nTrig = n_REP;
        fprintf('Percentage trig processed: \n')
        for indT=1:nTrig
                fileID=fopen([name '.lfp_sab.dat'],'r');
                shortbytes=2;
                offset=trig(indT)*nChn*shortbytes-0.25*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
                %ftell(fileID)
                fseek(fileID,offset,'bof');
                %ftell(fileID)
                vlfp = fread(fileID,[nChn, (0.25*FS+1.5*FS)],'short')./10; %plots 500ms from trigger and 250ms brefore
                fclose(fileID);
                vlfp = vlfp(:,1:30:end);
                fprintf('%0.8f%% \n', (indT*100/nTrig))
           for chsp=1:1:nChn
                check=['T', num2str(chsp)];
                if isfield(LFPstruct,check)
                    LFPstruct.(check) = [LFPstruct.(check); vlfp(chsp,:)];
                else
                    LFPstruct.(check)= vlfp(chsp,:);
                end
            end
        end
        meanlfpstruct=zeros(nChn,size(vlfp,2));
        for chsp=1:1:nChn
            check=['T', num2str(chsp)];
            meanlfpstruct(chsp,:)=mean(LFPstruct.(check),1);
        end




%%
n=800;
figure
        for chsp=1:1:nChn
            hold on
            check=['T', num2str(chsp)];
            plot(-250:1500-1,LFPstruct.(check)(n,:));
        end
xline(0,'r')
xline(1111.11,'r')
xlim([-250 1500])
title(['Trial ' num2str(n) ' filtered'])
xlabel('Time (ms)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING
%%
E_Mapnumber=1;
E_MAP = Depth(E_Mapnumber);
load('meanlfp.mat','meanlfpstruct')
%load('meanresp.mat','meanrespstruct')
%%
avgnospT=meanlfpstruct(E_MAP,:);
figure
subplot(1,4,1)
hold on
title('Shank 1')
xlabel('Time (ms)')
%ylabel('uV')
sepdist=100;
for i=16:-1:1
    hold on
plot(avgnospT(i,:)+i*sepdist,'k')
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xline(1111.11,'r')
xlim([-250 1500])

subplot(1,4,4)
hold on
title('Shank 4')
xlabel('Time (ms)')
%ylabel('uV')
for i=32:-1:17
plot(-250:1500-1,avgnospT(i,:)+i*sepdist,'k')
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xline(1111.11,'r')
xlim([-250 1500])

subplot(1,4,2)
hold on
title('Shank 2')
for i=48:-1:33
plot(-250:1500-1,avgnospT(i,:)+i*sepdist,'k')
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xlabel('Time (ms)')
%ylabel('uV')
xline(1111.11,'r')
xlim([-250 1500])

subplot(1,4,3)
hold on
title('Shank 3')
xlabel('Time (ms)')
%ylabel('uV')
for i=64:-1:49
plot(-250:1500-1,avgnospT(i,:)+i*sepdist,'k')
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xline(1111.11,'r')
xlim([-250 1500])



% figure
% for chsp=1:1:nChn
%             check=['T', num2str(chsp)];
%             
%             hold on
%             plot(-250:1000/30000:500-1000/30000,LFPstruct.(check)(1:200:end,:),'k')
% end
% figure
% for chsp=1:1:nChn
%             check=['T', num2str(chsp)];
%             
%             hold on
%             plot(-250:1000/30000:500-1000/30000,meanlfpstruct.(check),'k')
% end
% xlabel('Time (ms)')
% ylabel('uV')
% title('Average LFP waveform of each electrode timelocked to flash stimuli n_r_e_p=900')
% xline(0,'r')